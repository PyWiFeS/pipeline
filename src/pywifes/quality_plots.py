import matplotlib.pyplot as plt
import pickle
from photutils.aperture import RectangularAperture
from astropy.io import fits
import numpy as np

def slitlet_cutout(image, aperture):
    mask = aperture.to_mask(method="center")
    cutout = mask.cutout(image)
    return cutout

def slitlet_aperture(boundaries, bin_x=1, bin_y=2):
    xmin, xmax, ymin, ymax = boundaries
    xmin = np.round(xmin // bin_x)
    xmax = np.round(xmax // bin_x)
    ymin = np.round(ymin // bin_y)
    ymax = np.round(ymax // bin_y)
    center_x = (xmin + xmax) / 2
    center_y = (ymin + ymax) / 2
    width = xmax - xmin
    height = ymax - ymin
    aperture = RectangularAperture((center_x, center_y), width, height)
    return aperture

def read_pkl(path_pkl):
    with open(path_pkl, 'rb') as file:
        data = pickle.load(file)
    return data

def plot_slitlet(ax, path_slitlet):
    slitlets = read_pkl(path_slitlet)
    # Toy plot for labeling
    ax.plot(0,0, color ='red', lw=1, ls='--',label='Slitlet boundary',zorder=-1)

    for slit_number in slitlets:
        boundaries = slitlets[slit_number]
        aperture = slitlet_aperture(boundaries, bin_x=1, bin_y=2)
        aperture.plot(ax=ax, color='red', lw=1, ls='--')

def plot_fits(ax, image_path, min=5, max=95, cmap="nipy_spectral"):
    data = fits.getdata(image_path)
    vmin, vmax = np.percentile(data, (min, max))
    img = ax.imshow(data, vmin=vmin, vmax=vmax, cmap=cmap, origin="lower",aspect='auto',interpolation='None')
    return img

def plot_collapsed_slitlets(ax, path_slitlet, image_path):
    image = fits.getdata(image_path)
    slitlets = read_pkl(path_slitlet)
    for slit_number in slitlets:
        boundaries = slitlets[slit_number]
        aperture = slitlet_aperture(boundaries, bin_x=1, bin_y=2)
        cutout = slitlet_cutout(image, aperture)
        median_slitlet = np.median(cutout, axis=0)
        ax.plot(median_slitlet)

def slitlet_yticks(path_slitlet,bin_x=1, bin_y=2):
    slitlets = read_pkl(path_slitlet)
    y_centers = []
    slit_numbers = []
    for slit_number in slitlets:
        boundaries = slitlets[slit_number]
        xmin, xmax, ymin, ymax = boundaries
        ymin = np.round(ymin // bin_y)
        ymax = np.round(ymax // bin_y)
        center_y = (ymin + ymax) / 2

        y_centers.append(center_y)
        slit_numbers.append(slit_number)
    return y_centers, slit_numbers



def flatfiled_plot(flat_image_path, slitlet_path, title, output_plot):
    # Predefined list of colors (20 colors from 'viridis' colormap)
    colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        '#393b79', '#637939', '#8c6d31', '#843c39', '#7b4173', 
        '#5254a3', '#6b6ecf', '#9c9ede', '#8ca252', '#b5cf6b',
        '#cedb9c', '#bd9e39', '#e7ba52', '#e7969c', '#de9ed6'
    ]
    # Create figure and subplots with shared x-axis and different heights
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 10), sharex=True, 
                                   gridspec_kw={'height_ratios': [1.5, 1]})
    
    # Plot Title
    plt.suptitle(title,size=20)

    # Plot slitlets over 2D flatfield on the top
    img = plot_fits(ax1, flat_image_path, cmap='viridis')
    x_limits = ax1.get_xlim()
    plot_slitlet(ax1, slitlet_path)
    ax1.label_outer()  # Hide x-tick labels for the top subplot

    # Add color bar to the top of ax1
    cbar = fig.colorbar(img, ax=ax1, orientation='horizontal', pad=0.01, location='top', fraction=0.05, aspect=50)
    cbar.set_label('Count [e$^{-1}$]', size=15)

    # Set yticks and labels for ax1
    y_centers, slit_numbers = slitlet_yticks(slitlet_path, bin_x=1, bin_y=2)
    ax1.set_yticks(y_centers)
    ax1.set_yticklabels(slit_numbers)
    ax1.set_ylabel('Slitlet number', size=15)

    # Plot collapsed slitlet on the bottom
    slitlets = read_pkl(slitlet_path)
    for i, slit_number in enumerate(slitlets):
        boundaries = slitlets[slit_number]
        aperture = slitlet_aperture(boundaries, bin_x=1, bin_y=2)
        cutout = slitlet_cutout(fits.getdata(flat_image_path), aperture)
        median_slitlet = np.median(cutout, axis=0)
        color = colors[i % len(colors)]  # Use color from the list, cycle if more slitlets than colors
        ax2.plot(median_slitlet, color=color)

        # Set the corresponding tick label color in ax1
        tick_labels = ax1.get_yticklabels()
        tick_labels[i].set_color(color)
        ax1.set_yticklabels(tick_labels,size=12)

    ax2.set_ylabel('Count [e$^{-1}$]', size=15)

    y2_limits = ax2.get_ylim()

    # Add a light gray shaded area in ax2
    ax2.fill_between(x_limits, 0, 20000, lw=0,facecolor='none', hatch='////', edgecolor='lightgray', label='Optimal range')

    # Ensure the x-limits of the 2D image are preserved
    ax1.set_xlim(x_limits)
    ax2.set_xlim(x_limits)
    ax2.set_ylim(y2_limits)

    
    ax2.set_xlabel('X-axis [pixel]', size=15)
    ax2.legend(loc='upper right', framealpha=1.0)
    ax1.legend(loc='upper right', framealpha=1.0)

    # Adjust layout and show the plot
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    plt.close()

# Example usage
flat_image_path = 'data_products/master_calib/wifes_red_super_domeflat_raw.fits'
slitlet_path = 'data_products/master_calib/wifes_red_slitlet_defs.pkl'
output_plot = 'prueba.png'

flatfiled_plot(flat_image_path, slitlet_path, 'Titulo', output_plot)

