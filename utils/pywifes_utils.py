"""
Utilities functions to assist with reductions.

Note: When running without ssh -Y, get RuntimeError: Invalid DISPLAY variable
"""
from __future__ import print_function, division
import glob
import numpy as np
from astropy.io import fits
import sys
import os
import process_stellar as ps
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import pickle
from collections import Counter
from astropy.table import Table
import wifes_calib


def extract_stellar_spectra_ascii(root, night, steps = ["08", "09", "10"]):
    """Generate ascii spectra for each of the data reduction steps in steps.
    Saves ascii spectra in a new folder in the night directory called 'ascii'.

    Parameters:
    -----------
    root: string
        The base path to the reduced data (i.e. where the nightly folders are
        stored.)
    night: string
        The night to extract spectra from (in the form 201XYYZZ).
    steps: list of strings
        The PyWiFeS data reduction steps to extract and convert to 1D spectra.
    """
    # Sort out directory structures
    data_dir = os.path.join(root, night)
    out_dir = os.path.join(data_dir, 'ascii')

    if not os.path.isdir(out_dir) and not os.path.exists(out_dir):
        os.mkdir(out_dir)

    print('Converting to ascii:', data_dir)

    # Extract all files
    for path, subdirs, files in os.walk(data_dir):
        for name in files:
            fl=os.path.join(path, name)
            step = fl.split(".")[-2][1:]

            # Only run on specified data reduction outputs
            if step in steps and fl.endswith('%s.fits' % step):
                print(fl)
                f = fits.open(fl)
                header = f[0].header
                objectid = header['OBJNAME']
                run = header['RUN']
                f.close()

                # Extract spectrum
                flux, wave = ps.read_and_find_star_p08(fl)
                spectrum, sig = ps.weighted_extract_spectrum(flux)

                # Determine output format depending on spectral arm
                if 'T2m3wr' in name:
                    filename = '%s_%s_%s_r.dat'%(night, objectid, step)
                elif 'T2m3wb' in name:
                    filename = '%s_%s_%s_b.dat'%(night, objectid, step)
                
                filename = filename.replace(' ', '_')
                
                # Save output
                fln = os.path.join(out_dir, filename)
                print(fln)
                np.savetxt(fln, np.transpose([wave, spectrum]))


def make_extracted_stellar_cube(root, night, steps = ["08", "09", "10"]):
    """Make a fits data cube of the reduced and extracted stellar spectra, with
    a different HDU for each data reduction step. Store this new cube in a 
    folder called 'extracted_1D_cubes' in the night folder. By default:
     - 8: Final (non-fluxed or telluric corrected)
     - 9: Final (fluxed)
     -10: Final (fluxed and telluric corrected)

    Parameters:
    -----------
    root: string
        The base path to the reduced data (i.e. where the nightly folders are
        stored.)
    night: string
        The night to extract spectra from (in the form 201XYYZZ).
    steps: list of strings
        The PyWiFeS data reduction steps to extract and convert to 1D spectra.
    """
    # Sort out directory structures
    data_dir = os.path.join(root, night)
    out_dir = os.path.join(data_dir, 'extracted_cubes_1d')

    if not os.path.isdir(out_dir) and not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Get a list of the of the reduced files. Assume files have the naming
    # convention T2m3w[r/b]-20190828.104313-0031.p08.fits' where [r/b] will be
    # a single character indicating the arm of WiFeS used.
    unique_obs = glob.glob(os.path.join(data_dir, "*", "*%s.fits" % steps[0]))

    # Generalise by removing the suffixes
    unique_obs = [ob.replace(".p%s.fits" % steps[0], "") for ob in unique_obs]

    # Now go through one observation at a time, and create a data cube 
    # containing the original header information, and a HDU for each step of
    # the data reduction listed in steps.
    for ob in unique_obs:
        # Get the list of fits files to extract 1D spectra from
        fits_files = [ob + ".p%s.fits" % step for step in steps]

        # Get the header information from the first
        header = fits.getheader(fits_files[0])

        obj_id = header['OBJNAME'].replace(" ","")

        # Construct a new fits file
        hdus = []
        hdus.append(fits.PrimaryHDU(header=header))

        print("Making cube for %s%s" % (ob, ".pXX.fits"))

        for fits_file, step in zip(fits_files, steps):
            # Extract the 1D spectrum
            print("\tExtracting 1D spectrum for step %s" % step)
            flux, wave = ps.read_and_find_star_p08(fits_file)
            spectrum, sig = ps.weighted_extract_spectrum(flux)

            # Make fits table from numpy record array
            data = np.array([wave, spectrum, sig]).T.tolist()
            rec = np.rec.array(data, names=["wave", "spectrum", "sigma"])

            hdus.append(fits.BinTableHDU.from_columns(rec))
            hdus[-1].header["REDSTEP"] = step
        
        # Determine output format depending on spectral arm
        if 'T2m3wr' in ob:
            output_filename = '%s_%s_r.fits' % (night, obj_id)
        elif 'T2m3wb' in ob:
            output_filename = '%s_%s_b.fits' % (night, obj_id)

        output_filepath = os.path.join(out_dir, output_filename)

        # Write the fits file
        print("Writing cube to %s \n" % output_filename)
        hdu_list = fits.HDUList(hdus)
        hdu_list.writeto(output_filepath, overwrite=True)

def exclude_bad_frames(all_files, excluderun_filename):
    """
    Take a list of files (any kind: science, arc, cal...)
    that are bad. Exclude these files from the list.
    The list should include the dates and run numbers, e.g.
    20190824 345
    20190218 1
    etc.
    
    Returns:
    'all_files' where bad frames had been filtered out.
    """

    # File with dates and run numbers
    excluderun = np.loadtxt(excluderun_filename, comments='#', dtype=int)
    
    # Filter bad files
    all_files2=[]
    for x in all_files:
        # Data for frames of this night
        date = x.split('/')[-2]
        run = int(x.split('/')[-1].split('-')[-1].replace('.fits', ''))
        
        # Go through the list of bad files
        for xx in excluderun:
            d=xx[0]
            r=xx[1]
            #~ print('/%d/'%d, '%04d.fits'%r)
            if '/%d/'%d in x and '%04d.fits'%r in x:
                print('*** exclude me! ', x, date, run)
            else:
                if x not in all_files2:
                    all_files2.append(x)
    
    n = len(all_files)-len(all_files2)
    print('%d bad files excluded.'%n)
    all_files = all_files2

    return all_files

def update_object_header(root, night, print_only=False):
    """Quick function to sort out instances where OBJNAME header keyword is
    correct, but OBJECT keyword header is wrong (but has priority).

    Currently this looks for "exp" (i.e. exposure time) in OBJECT, due to a 
    TAROS copy-paste error

    Parameters
    ----------
    root: string
        Root directory where all the night folders are stored.

    night: string
        Name of the folder for a night's data, e.g. 20190828

    print_only: boolean
        Whether to actually change the headers, or just pring what they are
    """
    # Get a list of the files
    fits_files = glob.glob(os.path.join(root, night, "*.fits"))
    fits_files.sort()

    # For all files, open and check the state of the headers
    for fits_file in fits_files:
        with fits.open(fits_file, "update") as ff:
            if ("OBJNAME" in ff[0].header 
                and ff[0].header["IMAGETYP"]=="object"):
                if ("exp" in ff[0].header["OBJNAME"] 
                    and ff[0].header["OBJECT"] != ""):
                    print("%s: \t OBJECT: %s\t OBJNAME: %s" 
                           % (fits_file, ff[0].header["OBJECT"], 
                              ff[0].header["OBJNAME"]))
                    if not print_only:
                        ff[0].header["OBJNAME"] = ff[0].header["OBJECT"]

def list_of_all_observed_files():
    root = '/data/mash/marusa/2m3data/wifes/'
    young_stars_input_catalog_filename = 'young_stars_input_catalog.dat'

    keywords = ['FILENAME', 'EXPTIME', 'INSTRUME', 'DATE-OBS', 'PROPID', 'OBJNAME']

    imagetypes=[]

    filenames = []
    exptimes = []
    instrumes = []
    dateobss = []
    propids = []
    objnames = []
    dates = []
    detectors = []
    airmass = []
    ras = []
    decs = []

    # wifes
    gratingb=[]
    gratingr=[]
    beamsplt=[]
    ccdsum=[]
    ccdsec=[]
    
    runs=[]

    for r, d, f in os.walk(root):
        for filename in f:
            if filename.endswith('.fits') and 'T2m3ag' not in filename:
                filename=os.path.join(r, filename)
                try:
                    imagetyp = fits.getheader(filename, 0)['IMAGETYP']
                    #~ exptime = fits.getheader(filename, 0)['EXPTIME']

                    RA = fits.getheader(filename, 0)['RA']
                    #~ if imagetyp=='object' and exptime>0:
                        
                    v1 = fits.getheader(filename, 0)['FILENAME']
                    date = int(v1.split('-')[1].split('.')[0])
                    v2 = fits.getheader(filename, 0)['INSTRUME']
                    v3 = fits.getheader(filename, 0)['DATE-OBS']
                    v4 = fits.getheader(filename, 0)['PROPID']
                    v5 = fits.getheader(filename, 0)['OBJNAME']
                    v6 = fits.getheader(filename, 0)['EXPTIME']
                    v7 = fits.getheader(filename, 0)['DETECTOR']
                    v8 = fits.getheader(filename, 0)['AIRMASS']
                    v9 = fits.getheader(filename, 0)['DEC']
                    
                    v11 = fits.getheader(filename, 0)['GRATINGB']
                    v22 = fits.getheader(filename, 0)['GRATINGR']
                    v33 = fits.getheader(filename, 0)['BEAMSPLT']
                    v44 = fits.getheader(filename, 0)['CCDSUM']
                    v55 = fits.getheader(filename, 0)['CCDSEC']
                    
                    run = int(v1.split('-')[-1].replace('.fits', ''))

                    imagetypes.append(imagetyp)
                    ras.append(RA)
                    dates.append(date)
                    filenames.append(v1)
                    instrumes.append(v2)
                    dateobss.append(v3)
                    propids.append(int(v4))
                    objnames.append(v5)
                    exptimes.append(v6)
                    detectors.append(v7)
                    airmass.append(v8)
                    decs.append(v9)
                    runs.append(run)
                    

                    gratingb.append(v11)
                    gratingr.append(v22)
                    beamsplt.append(v33)
                    ccdsum.append(v44)
                    ccdsec.append(v55)


                except:
                    continue


    tab = Table([filenames], names=['filename'])
    tab['OBJNAME'] = objnames
    tab['RA'] = ras
    tab['DEC'] = decs
    tab['DATE'] = dates
    tab['RUN'] = runs
    tab['IMAGETYP'] = imagetypes    
    tab['INSTRUME'] = instrumes
    tab['DETECTOR'] = detectors
    tab['DATE-OBS'] = dateobss
    tab['AIRMASS'] = airmass
    tab['EXPTIME'] = exptimes
    tab['PROPID'] = propids

    # WiFeS
    tab['GRATINGB'] = gratingb
    tab['GRATINGR'] = gratingr
    tab['BEAMSPLT'] = beamsplt
    tab['CCDSUM'] = ccdsum
    tab['CCDSEC'] = ccdsec
    
    print(tab)
    filename = os.path.join(root, 'files.fits')
    tab.write(filename, format='fits', overwrite=True)
    print(filename, 'written.')


    ### PROGRAM ##########
    ys = np.loadtxt(young_stars_input_catalog_filename, comments='#', dtype=str)
    gaia_id = [int(x) for x in ys[:,0]]
    tmass_id = ys[:,1]
    mask = np.logical_or(np.in1d(tab['OBJNAME'], gaia_id), np.in1d(tab['OBJNAME'], tmass_id))







    #~ print(t)
    #~ print(len(t))

    objnames = [x[1] for x in t]
    print(set(objnames))
    print(len(set(objnames)))

def list_of_all_observed_young_star_files():
    root = '/data/mash/marusa/2m3data/wifes/'
    young_stars_input_catalog_filename = 'young_stars_input_catalog.dat'

    filename = os.path.join(root, 'files.fits')
    tab=Table.read(filename)


    ### PROGRAM ##########
    ys = np.loadtxt(young_stars_input_catalog_filename, comments='#', dtype=str)
    gaia_id = [int(x) for x in ys[:,0]]
    tmass_id = ys[:,1]
    mask = np.logical_or(np.in1d(tab['OBJNAME'], gaia_id), np.in1d(tab['OBJNAME'], tmass_id))

    t=tab[mask]




    #~ tab.write(filename, format='fits', overwrite=True)
    #~ print(filename, 'written.')

    print(t)

def bias_stats():
    """
    MZ: Write the doc!!
    
    Goal: Identify bad biases.
    """
    
    # Folder with fits files
    #~ root = sys.argv[1]
    root='/data/mash/marusa/2m3data/wifes/'
    print('root', root)
    
    # FILENAMES
    try:
        filenames = np.loadtxt('biases_raw.dat', dtype=str)
        print('Read')
    except:
        filenames=[]
    
        for path, subdirs, files in os.walk(root):
            for name in files:
                # Read header
                fl=os.path.join(path, name)
                try:
                    header = fits.getheader(fl, 0)
                    imagetyp = header['IMAGETYP']
                except:
                    continue
                if imagetyp.lower()!='bias' and imagetyp.lower()!='zero':
                    continue

                filenames.append(fl)
        print(filenames)
        np.savetxt('biases_raw.dat', filenames, fmt='%s')

    filenames = sorted(filenames)
    # Compute
    resultb=[]
    resultr=[]
    for fl in filenames:
        if 'ag' in fl:
            continue
        
        # We don't have WiFeS data in 2018
        if '2m3data/wifes/2018' in fl:
            continue
        
        header = fits.getheader(fl, 0)
        ccdsec = header['CCDSEC']
        ccdsum = header['CCDSUM']
        
        gratingb = header['GRATINGB']
        gratingr = header['GRATINGR']
        beamsplt = header['BEAMSPLT']

        image_data = fits.getdata(fl, ext=0)
        
        # Take one line
        line = image_data#[50,:]
        
        med = int(np.median(image_data))
        #~ med = int(np.median(line[10:-10]))
        
        if 'T2m3wb' in fl:
            #~ resultb.append(med)
            resultb.append(line)
            
            if med>3300:
                #~ print(fl, ccdsec, ccdsum, med, gratingb, gratingr, beamsplt)
                print(fl)
            
        elif 'T2m3wr' in fl:
            #~ resultr.append(med)
            resultr.append(line)
            
            if med>4900:
                #~ print(fl, ccdsec, ccdsum, med, gratingb, gratingr, beamsplt)
                print(fl)
    
    resultb=np.array(resultb)
    resultr=np.array(resultr)
    print(resultb)
    x=range(resultb.shape[1])
    fig=plt.figure()
    ax=fig.add_subplot(211)
    for y in resultb:
        ax.plot(x, y, c='k', alpha=0.5)
    ax=fig.add_subplot(212)
    for y in resultr:
        ax.plot(x, y, c='k', alpha=0.5)
    
    
    plt.show()

def dark_stats():
    """
    MZ: Write the doc!!
    
    Goal: Identify bad biases.
    """
    
    # Folder with fits files
    #~ root = sys.argv[1]
    root='/data/mash/marusa/2m3data/wifes/'
    print('root', root)
    
    # FILENAMES
    try:
        filenames = np.loadtxt('darks_raw.dat', dtype=str)
        print('Read')
    except:
        filenames=[]
    
        for path, subdirs, files in os.walk(root):
            for name in files:
                # Read header
                fl=os.path.join(path, name)
                try:
                    header = fits.getheader(fl, 0)
                    imagetyp = header['IMAGETYP']
                except:
                    continue
                
                if imagetyp.lower()!='dark':
                    continue

                filenames.append(fl)
        print(filenames)
        np.savetxt('darks_raw.dat', filenames, fmt='%s')
                    


    # Compute
    resultb=[]
    resultr=[]
    for fl in filenames:
        header = fits.getheader(fl, 0)

        image_data = fits.getdata(fl, ext=0)
        exptime = header['EXPTIME']
        image_data = image_data / exptime
        
        # Take one line
        line = image_data[50,:]
        
        #~ med = np.median(image_data)
        med = np.median(line)
        print(fl, med)
        
        if 'T2m3wb' in fl:
            #~ resultb.append(med)
            resultb.append(line)
        elif 'T2m3wr' in fl:
            #~ resultr.append(med)
            resultr.append(line)

    
    resultb=np.array(resultb)
    resultr=np.array(resultr)
    print(resultb)
    x=range(resultb.shape[1])
    fig=plt.figure()
    ax=fig.add_subplot(211)
    for y in resultb:
        ax.plot(x, y, c='k', alpha=0.5)
    ax=fig.add_subplot(212)
    for y in resultr:
        ax.plot(x, y, c='k', alpha=0.5)
    
    
    plt.show()

def flats_stats(camera='Blue'):
    """
    Check how stable flats are over time. Plot p02 files.  
    """
    
    # Folder with fits files
    #~ root = sys.argv[1]
    #~ root='/data/mash/marusa/2m3reduced/wifes/'
    root='/data/mash/marusa/wifes_old/'
    #~ root='/data/mash/marusa/2m3data/wifes/'
    print('root', root)
    
    # Find filenames
    try:
        filenames = np.loadtxt('filenames_flats_p02_red??.dat', dtype=str)
        #~ doesnt_work=True
    except:
        print('READING ALL FITS FILES')
        filenames=[]
        for path, subdirs, files in os.walk(root):
            for name in files:
                fl=os.path.join(path, name)
                
                if camera=='Blue':
                    if 'ag' in fl or '/reduced_b/T2m3wb' not in fl or '.fits' not in fl or '.pkl' in fl or 'p02' not in fl:
                        continue
                elif camera=='Red':
                    if 'ag' in fl or '/reduced_r/T2m3wr' not in fl or '.fits' not in fl or '.pkl' in fl or 'p02' not in fl:
                        continue
                
                #~ print(fl)
                
                #~ if 'wifesB_super_domeflat.fits' in fl:
                    #~ print(fl)
                    
                # Read data
                header = fits.getheader(fl, 0)
                
                try:
                    imagetyp = header['IMAGETYP']
                except:
                    continue
                
                if imagetyp != 'flat':
                    continue
                #~ print(imagetyp)
                
                ccdsec = header['CCDSEC']
                ccdsum = header['CCDSUM']
                beamsplt = header['BEAMSPLT']
                gratingb = header['GRATINGB']
                exptime = int(header['EXPTIME'])
                lamp = header['M1ARCLMP']
                #~ print('EXPTIME', exptime)
                #~ if exptime!=5:
                    #~ continue
                
                #~ if beamsplt!='RT480' or gratingb!='B3000':
                if beamsplt!='RT480':
                    continue                
                
                if lamp!='QI-1':
                    continue

                filenames.append(fl)
        
        
        # Sort filenames to get time sequence
        filenames=sorted(filenames)
        #~ for x in filenames:
            #~ print(x)
        
        np.savetxt('filenames_flats.dat', filenames, fmt='%s')
    
    # Exclude bad filenames
    bad = ['/data/mash/marusa/2m3reduced/wifes/20191013/reduced_b/T2m3wb-20190805.201127-0690.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20191013/reduced_b/T2m3wb-20190805.201202-0691.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20190805/reduced_r/T2m3wr-20190805.201127-0690.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20190805/reduced_r/T2m3wr-20190805.201202-0691.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20191013/reduced_r/T2m3wr-20190805.201127-0690.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20191013/reduced_r/T2m3wr-20190805.201202-0691.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20191014/reduced_r/T2m3wr-20191014.190602-0112.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20191014/reduced_r/T2m3wr-20191014.190742-0113.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20190731/reduced_r/T2m3wr-20190731.071750-0126.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20190731/reduced_r/T2m3wr-20190731.071825-0127.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20190731/reduced_r/T2m3wr-20190731.071900-0128.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20190731/reduced_r/T2m3wr-20190731.071934-0129.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20190805/reduced_r/T2m3wr-20190805.201714-0696.p02.fits', '/data/mash/marusa/2m3reduced/wifes/20191013/reduced_r/T2m3wr-20190805.201714-0696.p02.fits'] # mirror was not in, or overexposed
    #~ filenames = [x for x in filenames if x not in bad]
    

    # Colors
    start = 0.0
    stop = 1.0
    number_of_lines = len(filenames)
    cm_subsection = np.linspace(start, stop, number_of_lines) 
    colors = [cm.jet(x) for x in cm_subsection]
    
    alpha=1.0

    # both
    fig=plt.figure()
    ax=fig.add_subplot(111)
    
    counts=np.zeros(4)
    medians=[]
    labels=[]

    data=[]

    for i, fl in enumerate(filenames):
        # Read data
        #~ print(fl)
        header = fits.getheader(fl, 0)
        image_data = fits.getdata(fl, ext=0)
        ccdsec = header['CCDSEC']
        ccdsum = header['CCDSUM']
        
        # Extract one line
        # This depends on whether it is full/stellar frame and the binning!!
        if ccdsec == '[1:4202,2057:4112]' and ccdsum == '1 1': # stellar and ybin 1; OK mostly
            line = image_data[446:520,:]
            c='g'
            label='stellar 1'
            counts[0]=counts[0]+1
        elif ccdsec == '[1:4202,2057:4112]' and ccdsum == '1 2': # stellar and ybin 2
            line = image_data[225:258,:]
            c='k'
            label='stellar 2'
            counts[1]=counts[1]+1
        elif ccdsec == '[1:4202,1:4112]' and ccdsum == '1 1': # full frame and ybin 1; OK
            line = image_data[2500:2575,:]
            c='r'
            label='full 1'
            counts[2]=counts[2]+1
            #~ continue
        elif ccdsec == '[1:4202,1:4112]' and ccdsum == '1 2': # full frame and ybin 2; OK
            line = image_data[1251:1287,:]
            c='b'
            label='full 2'
            counts[3]=counts[3]+1
        
        # Combine all rows of slitlet into one
        line = np.nanmedian(line, axis=0)
        
        # Normalize
        #~ m = np.max(line)
        m=np.percentile(line, 90)
        medians.append(m)
        line = line/m

        # FILTERS: RED
        if np.median(line[3000:3500])<0.5 or line[2000]>0.9:
        #~ if line[0]<0.23:
            print(fl, m, label)
            continue
        if line[0]<0.235:
            print(fl, m, label, '******')
            #~ continue

        date=fl.split('/')[-3]
        
        #~ print(fl, label, len(line))
        
        x=range(len(line)) # xbin is always 1!
        ax.plot(x, line, c='k', alpha=0.1)
        
        data.append(line)

    print(counts)

    # An average flat    
    data=np.array(data)
    overall_median = np.median(data, axis=0)
    np.savetxt('average_flat.dat', overall_median)
    
    ax.plot(x, overall_median, c='r')

    plt.show()

def settings_stats():
    """
    How many images were taken with a give set of settings
    """
    keywords=['GRATINGB', 'GRATINGR', 'BEAMSPLT', 'CCDSUM', 'CCDSEC']

    root='/data/mash/marusa/2m3data/wifes/'
    
    result=dict()
    
    # Find filenames
    print('READING ALL FITS FILES')
    filenames=[]
    for path, subdirs, files in os.walk(root):
        for name in files:
            fl=os.path.join(path, name)

            if 'ag' in fl or '.fits' not in fl:
                continue
            
            # Read data
            f=fits.open(fl)
            header = f[0].header
            f.close()
            
            try:
                imagetyp = header['IMAGETYP']
            except:
                continue
            
            try:
                tmp=[]
                for k in keywords:
                    v = header[k]
                    if v=='[1:4202,2057:4112]':
                        v='stel'
                    elif v=='[1:4202,1:4112]':
                        v='full'
                    elif v=='1 1':
                        v=1
                    elif v=='1 2':
                        v=2
                    tmp.append(v)
                
                tmp=tuple(tmp)
                
                try:
                    result[tmp].append(imagetyp)
                except:
                    result[tmp]=[imagetyp]
            except:
                continue


    r=[]
    for k, v in result.iteritems():
        # Count imagetypes
        values = Counter(v).keys()
        numbers = Counter(v).values()
        
        r.append(['%04d'%len(v), k, ['(%d, %s)'%(x, y) for x, y in zip(numbers, values)]])
        
        #~ print(k, len(v), ['(%d, %s)'%(x, y) for x, y in zip(numbers, values)])
    
    r=sorted(r, key=lambda x: x[0])

    for x in r:
        print(x)

def compare_sens_func_best_fit_and_ratio():
    """
    
    """

    best_fit = [7.58898732e-33, -2.89682376e-28, 4.13047385e-24, -2.17367092e-20, -5.97792537e-17, 1.05761370e-12, -3.65871852e-10, -4.60135001e-05, 2.84268254e-01, -7.26833856e+02, 7.13182151e+05]

    func=np.poly1d(best_fit)

    with open('/data/mash/marusa/2m3reduced/wifes/20199999/reduced_r/wifesR_calib.pkl', 'rb') as h:
        calib = pickle.load(h)
    
    plt.plot(calib['wave'], calib['cal'], c='k')
    plt.plot(calib['wave'], func(calib['wave']), c='r')
    
    # New best fit
    best_fit = np.polyfit(calib['wave'], calib['cal'], 40)
    func=np.poly1d(best_fit)
    plt.plot(calib['wave'], func(calib['wave']), c='blue')

    plt.show()

def fits2png():
    """
    Diagnostics: Convert all fits files in the folder into png/jpg and save them into separate folders based on their imagetype. This is to have a quick look if all imagetypes are correct.
    """
    # Folder with fits files
    root = sys.argv[1]
    print('root', root)

    for path, subdirs, files in os.walk(root):
        for name in files:
            fl=os.path.join(path, name)
            if fl.endswith('.fits') or fl.endswith('.fit') or fl.endswith('.FITS') or fl.endswith('.FIT'):
                print(fl)
                
                # Read data
                f=fits.open(fl)
                header = f[0].header
                #~ objectid=header['OBJNAME']
                #~ run=header['RUN']
                imagetype=header['IMAGETYP']
                f.close()
                image_data = fits.getdata(fl, ext=0)

                # Set colorscale
                minimum = np.percentile(image_data, 30)
                maximum = np.percentile(image_data, 95)
                
                print(minimum, maximum)
                

                # Make an output folder
                out_dir = os.path.join(root, 'png')
                if not os.path.isdir(out_dir) and not os.path.exists(out_dir):
                    os.mkdir(out_dir)
                out_dir = os.path.join(root, 'png', imagetype.lower())
                print('out_dir', out_dir)
                if not os.path.isdir(out_dir) and not os.path.exists(out_dir):
                    os.mkdir(out_dir)

                # Plot and save figure
                plt.figure()
                #~ plt.imshow(image_data, cmap='gray', vmin=minimum, vmax=maximum)
                plt.imshow(image_data, cmap='gray', norm=LogNorm(vmin=minimum, vmax=maximum))
                #~ plt.colorbar()
                fl_out = os.path.join(out_dir, name.replace('.fits', '.png'))
                print('fl_out', fl_out)
                plt.savefig(fl_out)

def delete_ag_images():
    """
    Some fits files are from the Imager (?) from daytime testing??

    Delete them in all subfolders.
    """

    import os
    import sys


    path = sys.argv[1]

    for r, d, f in os.walk(path):
        for file in f:
            if 'T2m3ag' in file:
                filename = os.path.join(r, file)
                print(filename)
                os.remove(filename)

def find_all_standard_star_frames(root='/data/mash/marusa/2m3data/wifes/', exclude_bad_files=True, exclude_bad_filename='/data/mash/marusa/2m3data/wifes/list_of_bad_exposures_that_we_shouldn_use.dat'):
    """
    Find all flux and telluric standards in the folder and its subfolders.
    """

    # List of standard stars
    stdstar_list = wifes_calib.ref_fname_lookup.keys()
    stdstar_is_flux_cal = wifes_calib.ref_flux_lookup
    stdstar_is_telluric = wifes_calib.ref_telluric_lookup

    keywords=['GRATINGB', 'GRATINGR', 'BEAMSPLT', 'CCDSUM', 'CCDSEC']

    print('READING ALL FITS FILES')
    all_files=[]
    for path, subdirs, files in os.walk(root):
        for name in files:
            if name.endswith('.fits') and 'T2m' in name:
                fl=os.path.join(path, name)
                all_files.append(fl)
    
    if exclude_bad_files:
        all_files = exclude_bad_frames(all_files, exclude_bad_filename)
    
    
    filenames_flux=dict()
    filenames_tell=dict()
    for fl in all_files:
        try:
            header = fits.getheader(fl, 0)
            imagetyp = header['IMAGETYP']
        except:
            continue
        
        #~ if imagetyp.lower() != 'science' and imagetyp.lower() != 'object':
            #~ continue

        try: 
            objname = header['OBJNAME'].replace(' ', '')
        except:
            continue
        
        # Flux standards
        try:
            flux = stdstar_is_flux_cal[objname]

            if flux:
                mode = tuple([header[key] for key in keywords])

                try:
                    filenames_flux[mode].append(fl)
                except:
                    filenames_flux[mode]=[fl]
        except:
            continue
        
        # Telluric standards
        try:
            tell = stdstar_is_telluric[objname]
            
            if tell:
                mode = tuple([header[key] for key in keywords])
                
                try:
                    filenames_tell[mode].append(fl)
                except:
                    filenames_tell[mode]=[fl]
        except:
            continue
 

    print('Flux')
    for k, v in filenames_flux.iteritems():
        print(k, len(v))

    print()
    print()
    print()
    print('Tell')
    for k, v in filenames_tell.iteritems():
        print(k, len(v))

    print()  
    for k, v in filenames_flux.iteritems():
        #~ if k==('B3000', 'R7000', 'RT480', '1 1', '[1:4202,1:4112]'):
        if k==('B3000', 'R7000', 'RT480', '1 1', '[1:4202,2057:4112]'):
            print(k)
            for x in v:
                print(x)
  
