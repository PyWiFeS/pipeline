#! /usr/bin/env python

##########################################
import pygtk
pygtk.require('2.0')
import gtk

import numpy
import os
import re
from astropy.io import fits as pyfits

import pylab
import matplotlib.cm
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtk import NavigationToolbar2GTK as NavigationToolbar

from mjc_snid import compare_all_snid, retrieve_snid_template_data, clean_sn_spectrum

pylab.rc('figure.subplot',
         top=0.95,
         bottom=0.06,
         left=0.06,
         right=0.95)

############################################
class GalaxyGrabber:
    def __init__(self):
        #-----------------------------------------
        # internal variables to be initialized
        #-----------------------------------------
        # filename info
        self.blue_cube_fn = None
        self.red_cube_fn = None
        self.object_name = None

        # spaxel indices and save selectors
        self.blue_gal_spax = None
        self.blue_sky_spax = None
        self.blue_save_choice = 'Subbed Spaxels'
        self.red_gal_spax = None
        self.red_sky_spax = None
        self.red_save_choice = 'Subbed Spaxels'

        # cube data storage
        self.blue_cube_data = None
        self.blue_var_cube_data = None
        self.blue_wavelength_array = None
        self.red_cube_data = None
        self.red_var_cube_data = None
        self.red_wavelength_array = None

        # display limit variables
        self.blue_disp_wmin = None
        self.blue_disp_wmax = None
        self.red_disp_wmin = None
        self.red_disp_wmax = None

        # correlator spectrum storage
        self.corr_blue_wave = None
        self.corr_blue_flux = None
        self.corr_red_wave = None
        self.corr_red_flux = None
        self.corr_wave = None
        self.corr_flux = None

        # header info
        self.blue_cube_header = None
        self.blue_header_menu = None
        self.red_cube_header = None
        self.red_header_menu = None

        # dictionary of color schemes
        self.chosen_color_scheme = 'Jet'
        #self.chosen_color_scheme = 'Hot'
        #self.chosen_color_scheme = 'Plant'
        #self.chosen_color_scheme = 'Blue'
        
        self.color_schemes = {}
        self.color_schemes['Jet'] = {'cmap':'jet',
                                     'gal_color':'w',
                                     'sky_color':'k',
                                     'dead_color':'r',
                                     'gal_spec_color':'b',
                                     'sky_spec_color':'r',
                                     'sub_spec_color':'g'}
        self.color_schemes['Hot'] = {'cmap':'hot',
                                     'gal_color':'b',
                                     'sky_color':'m',
                                     'dead_color':'w',
                                     'gal_spec_color':'b',
                                     'sky_spec_color':'m',
                                     'sub_spec_color':'g'}
        self.color_schemes['Plant'] = {'cmap':'RdYlGn_r',
                                       'gal_color':'b',
                                       'sky_color':'k',
                                       'dead_color':'r',
                                       'gal_spec_color':'b',
                                       'sky_spec_color':'m',
                                       'sub_spec_color':'g'}
        self.color_schemes['Blue'] = {'cmap':'Blues_r',
                                      'gal_color':'g',
                                      'sky_color':'k',
                                      'dead_color':'r',
                                      'gal_spec_color':'g',
                                      'sky_spec_color':'k',
                                      'sub_spec_color':'r'}
        self.color_schemes['Green'] = {'cmap':'Greens_r',
                                       'gal_color':'b',
                                       'sky_color':'k',
                                       'dead_color':'r',
                                       'gal_spec_color':'b',
                                       'sky_spec_color':'k',
                                       'sub_spec_color':'r'}
        self.color_schemes['Purple'] = {'cmap':'Purples_r',
                                        'gal_color':'b',
                                        'sky_color':'k',
                                        'dead_color':'r',
                                        'gal_spec_color':'b',
                                        'sky_spec_color':'k',
                                        'sub_spec_color':'r'}
        
        #------------------------------------------
        # window setup and kill signal definitions
        #------------------------------------------
        # create a window
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_title('WiFeS Data Cube Explorer')

        # what to do when 'delete event' signal comes
        # defined in self function below
        self.window.connect('delete_event', self.delete_event)
        # what to do with a destroy event
        self.window.connect('destroy', self.destroy)

        # set border width of window
        self.window.set_border_width(10)

        #-----------------------------------------
        # details of what's in the window
        #-----------------------------------------
        # start with Vbox to hold everything
        self.vbox = gtk.VBox(False, 0)
        self.window.add(self.vbox)
        self.vbox.show()
        # hold everything in a Notebook!
        self.notebook = gtk.Notebook()
        self.notebook.set_tab_pos(gtk.POS_TOP)
        self.vbox.pack_start(self.notebook, False, False, 0)
        self.notebook.show()
        
        #------------------------------------------------------------
        #------------------------------------------------------------
        #   PAGE 1 - RED CUBE DATA
        #------------------------------------------------------------
        #------------------------------------------------------------
        self.red_vbox = gtk.VBox(False, 10)
        self.red_vbox.show()
        red_vbox_label = gtk.Label('Red Data Cube')
        self.notebook.append_page(self.red_vbox, red_vbox_label)
        
        #------------------------------------------------
        # Section 1: data selection
        red_data_select_frame = gtk.Frame('Red Data Selection')
        self.red_vbox.pack_start(red_data_select_frame, False, False, 0)
        red_data_select_frame.set_border_width(10)
        red_data_select_frame.show()
        
        # it will all be one vbox
        vbox1 = gtk.VBox(False, 10)
        red_data_select_frame.add(vbox1)
        vbox1.show()
                
        # make a hbox1c for:
        # -red file names
        hbox1c = gtk.HBox(False, 0)
        vbox1.pack_start(hbox1c, True, False, 0)
        hbox1c.show()

        # 1c1 - red file name
        red_fn_label = gtk.Label('Red Cube Filename ')
        hbox1c.pack_start(red_fn_label, False, False, 5)
        red_fn_label.show()
        
        red_fn_box = gtk.Entry()
        self.red_fn_entry = red_fn_box
        red_fn_box.set_max_length(500)
        red_fn_box.connect('activate', self.update_red_fn)
        hbox1c.pack_start(red_fn_box, True, True, 5)
        red_fn_box.show()

        # button to choose red filename
        red_filename_button = gtk.Button('Choose...')
        red_filename_button.connect('clicked', self.select_red_filename)
        hbox1c.pack_start(red_filename_button, False, False, 5)
        red_filename_button.show()
        
        # hbox for header info
        header_hbox = gtk.HBox(False, 0)
        vbox1.pack_start(header_hbox, True, False, 0)
        header_hbox.show()

        red_header_label = gtk.Label('Red Header Field:')
        header_hbox.pack_start(red_header_label, False, False, 10)
        red_header_label.show()

        red_header_selector = gtk.Entry()
        red_header_selector.set_max_length(8)
        red_header_selector.set_width_chars(9)
        red_header_selector.connect('activate', self.update_red_header_info)
        header_hbox.pack_start(red_header_selector, False, False, 10)
        red_header_selector.show()

        red_header_label2 = gtk.Label('Contents:')
        header_hbox.pack_start(red_header_label2, False, False, 10)
        red_header_label2.show()

        self.red_header_entry = gtk.Entry()
        self.red_header_entry.set_editable(False)
        header_hbox.pack_start(self.red_header_entry, True, True, 10)
        self.red_header_entry.show()

        #------------------------------------------------
        # section2: red cube image and spectra
        red_cube_frame = gtk.Frame('Red Cube')
        self.red_vbox.pack_start(red_cube_frame, False, False, 0)
        red_cube_frame.show()

        # the frame will be managed in an hbox
        self.red_hbox = gtk.HBox(False, 10)
        red_cube_frame.add(self.red_hbox)
        self.red_hbox.show()

        #-------------------
        #-------------------
        # add a place for cube image
        red_cube_vbox = gtk.VBox(False, 0)
        self.red_hbox.pack_start(red_cube_vbox, False, False, 0)
        red_cube_vbox.show()
        
        red_cube_fig = Figure(figsize=(3,5), dpi=50)
        self.red_cube_plot = red_cube_fig.add_subplot(1,1,1)
        self.red_cube_plot.hold(True)
        self.red_cube_image = FigureCanvas(red_cube_fig)
        self.red_cube_image.mpl_connect('button_press_event',
                                        self.click_red_cube_image)
        self.red_cube_image.set_size_request(300,500)
        red_cube_vbox.pack_start(self.red_cube_image, False, False, 0)
        self.red_cube_image.show()

        # colorbar scale adjustor!
        self.red_im_adj = gtk.Adjustment(1.0, 0.0, 1.0, 0.01, 0.01)
        self.red_im_adj.connect('value_changed', self.replot_red_cube_image)
        red_im_scaler = gtk.HScrollbar(self.red_im_adj)
        red_im_scaler.set_update_policy(gtk.UPDATE_CONTINUOUS)
        red_cube_vbox.pack_start(red_im_scaler, False, False, 0)
        red_im_scaler.show()

        # and the corresponding expand button
        red_cube_expand_button = gtk.Button('Expand')
        red_cube_expand_button.connect('clicked', self.expand_red_cube)
        red_cube_vbox.pack_start(red_cube_expand_button, False, False, 0)
        red_cube_expand_button.show()

        #--------------------------------
        #--------------------------------
        # overall vbox for spaxel selection and saving frame
        red_spax_big_vbox = gtk.VBox(False, 10)
        self.red_hbox.pack_start(red_spax_big_vbox, False, False, 0)
        red_spax_big_vbox.show()

        #-+-+-+-+-+-+-+-+-+
        # frame for setting display wavelength bounds for cube image
        red_disp_frame = gtk.Frame('Red Cube Display Options')
        red_spax_big_vbox.pack_start(red_disp_frame, False, False, 0)
        red_disp_frame.show()

        # which will be managed by a vbox
        red_disp_vbox1 = gtk.VBox(False, 0)
        red_disp_frame.add(red_disp_vbox1)
        red_disp_vbox1.show()

        #-+-+-+
        # hbox for wavelength min
        red_wmin_hbox = gtk.HBox(False, 0)
        red_disp_vbox1.pack_start(red_wmin_hbox, False, False, 0)
        red_wmin_hbox.show()

        # label
        red_wmin_label = gtk.Label('Lambda Min:')
        red_wmin_hbox.pack_start(red_wmin_label, False, False, 0)
        red_wmin_label.show()

        # entry
        self.red_wmin_entry = gtk.Entry()
        red_wmin_hbox.pack_start(self.red_wmin_entry, False, False, 0)
        self.red_wmin_entry.set_width_chars(8)
        self.red_wmin_entry.show()

        #-+-+-+        
        # hbox for wavelength max
        red_wmax_hbox = gtk.HBox(False, 0)
        red_disp_vbox1.pack_start(red_wmax_hbox, False, False, 0)
        red_wmax_hbox.show()

        # label
        red_wmax_label = gtk.Label('Lambda Max:')
        red_wmax_hbox.pack_start(red_wmax_label, False, False, 0)
        red_wmax_label.show()

        # entry
        self.red_wmax_entry = gtk.Entry()
        red_wmax_hbox.pack_start(self.red_wmax_entry, False, False, 0)
        self.red_wmax_entry.set_width_chars(8)
        self.red_wmax_entry.show()

        #-+-+-+
        # button to update plot
        red_disp_button = gtk.Button('Set')
        red_disp_button.connect('clicked', self.set_red_disp_wave)
        red_disp_vbox1.pack_start(red_disp_button, True, False, 0)
        red_disp_button.show()
        
        #-+-+-+-+-+-+-+-+-+
        # spaxel selection tools will go in a frame
        red_spax_frame = gtk.Frame('Red Spaxel Selection')
        red_spax_big_vbox.pack_start(red_spax_frame, False, False, 0)
        red_spax_frame.show()

        # which will be managed by a vbox
        red_spax_vbox1 = gtk.VBox(False, 0)
        red_spax_frame.add(red_spax_vbox1)
        red_spax_vbox1.show()

        #-+-+-+-+-+
        # first make a vbox for red gal spax
        red_gal_spax_vbox = gtk.VBox(False, 0)
        red_spax_vbox1.pack_start(red_gal_spax_vbox, False, False, 10)
        red_gal_spax_vbox.show()

        #-+-+
        # label in an hbox
        red_gal_spax_hbox1 = gtk.HBox(False, 0)
        red_gal_spax_vbox.pack_start(red_gal_spax_hbox1, False, False, 0)
        red_gal_spax_hbox1.show()
        
        red_gal_spax_label = gtk.Label('Red Galaxy Spaxel(s)')
        red_gal_spax_hbox1.pack_start(red_gal_spax_label, True, False, 0)
        red_gal_spax_label.show()

        #-+-+
        # then entry box
        red_gal_spax_hbox2 = gtk.HBox(False, 0)
        red_gal_spax_vbox.pack_start(red_gal_spax_hbox2, False, False, 0)
        red_gal_spax_hbox2.show()
        
        self.red_gal_spax_entry = gtk.Entry()
        #self.red_gal_spax_entry.set_max_length(3)
        self.red_gal_spax_entry.connect('activate', self.update_red_gal_spax)
        red_gal_spax_hbox2.pack_start(self.red_gal_spax_entry, False, True, 0)
        self.red_gal_spax_entry.show()

        #-+-+
        # buttons in an hbox
        red_gal_spax_hbox3 = gtk.HBox(False, 0)
        red_gal_spax_vbox.pack_start(red_gal_spax_hbox3, False, False, 0)
        red_gal_spax_hbox3.show()
        
        red_gal_spax_button = gtk.Button('Set')
        red_gal_spax_button.connect('clicked', self.update_red_gal_spax)
        red_gal_spax_hbox3.pack_start(red_gal_spax_button, True, False, 0)
        red_gal_spax_button.show()
        
        red_gal_spax_button = gtk.Button('Clear')
        red_gal_spax_button.connect('clicked', self.clear_red_gal_spax)
        red_gal_spax_hbox3.pack_start(red_gal_spax_button, True, False, 0)
        red_gal_spax_button.show()

        #-+-+
        # copy button in another hbox
        red_gal_spax_hbox4 = gtk.HBox(False, 0)
        red_gal_spax_vbox.pack_start(red_gal_spax_hbox4, False, False, 0)
        red_gal_spax_hbox4.show()
        
        red_gal_copy_button = gtk.Button('Copy Blue Gal Spax')
        red_gal_copy_button.connect('clicked', self.copy_blue_to_red_gal_spax)
        red_gal_spax_hbox4.pack_start(red_gal_copy_button, True, False, 0)
        red_gal_copy_button.show()

        #-+-+-+-+-+
        # then make a vbox for red sky spax
        red_sky_spax_vbox = gtk.VBox(False, 0)
        red_spax_vbox1.pack_start(red_sky_spax_vbox, False, False, 5)
        red_sky_spax_vbox.show()

        #-+-+
        # label in an hbox
        red_sky_spax_hbox1 = gtk.HBox(False, 0)
        red_sky_spax_vbox.pack_start(red_sky_spax_hbox1, False, False, 0)
        red_sky_spax_hbox1.show()
        
        red_sky_spax_label = gtk.Label('Red Sky Spaxel(s)')
        red_sky_spax_hbox1.pack_start(red_sky_spax_label, True, False, 0)
        red_sky_spax_label.show()
        
        #-+-+
        # then entry in an hbox2
        red_sky_spax_hbox2 = gtk.HBox(False, 0)
        red_sky_spax_vbox.pack_start(red_sky_spax_hbox2, False, False, 0)
        red_sky_spax_hbox2.show()
        
        self.red_sky_spax_entry = gtk.Entry()
        #self.red_sky_spax_entry.set_max_length(3)
        self.red_sky_spax_entry.connect('activate', self.update_red_sky_spax)
        red_sky_spax_hbox2.pack_start(self.red_sky_spax_entry, False, True, 0)
        self.red_sky_spax_entry.show()

        #-+-+
        # buttons in an hbox
        red_sky_spax_hbox3 = gtk.HBox(False, 0)
        red_sky_spax_vbox.pack_start(red_sky_spax_hbox3, False, False, 0)
        red_sky_spax_hbox3.show()
        
        red_sky_spax_button = gtk.Button('Set')
        red_sky_spax_button.connect('clicked', self.update_red_sky_spax)
        red_sky_spax_hbox3.pack_start(red_sky_spax_button, True, False, 0)
        red_sky_spax_button.show()
        
        red_sky_spax_button = gtk.Button('Clear')
        red_sky_spax_button.connect('clicked', self.clear_red_sky_spax)
        red_sky_spax_hbox3.pack_start(red_sky_spax_button, True, False, 0)
        red_sky_spax_button.show()
        
        #red_sky_spax_button = gtk.Button('Autoset')
        #red_sky_spax_button.connect('clicked', self.autoset_red_sky_spax)
        #red_sky_spax_hbox3.pack_start(red_sky_spax_button, True, False, 0)
        #red_sky_spax_button.show()

        #-+-+
        # copy button in another hbox
        red_sky_spax_hbox4 = gtk.HBox(False, 0)
        red_sky_spax_vbox.pack_start(red_sky_spax_hbox4, False, False, 0)
        red_sky_spax_hbox4.show()
        
        red_sky_copy_button = gtk.Button('Copy Blue Sky Spax')
        red_sky_copy_button.connect('clicked', self.copy_blue_to_red_sky_spax)
        red_sky_spax_hbox4.pack_start(red_sky_copy_button, True, False, 0)
        red_sky_copy_button.show()

        #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        # put instructions in their own frame
        #red_spax_instructions = gtk.Frame('Instructions')
        #red_spax_big_vbox.pack_start(red_spax_instructions, False, False, 0)
        #red_spax_instructions.show()
        #red_instructions_string = 'Left click a spaxel to select\n'
        #red_instructions_string += '(or deselect) a galaxy spaxel,\n'
        #red_instructions_string += 'right click a spaxel to select\n'
        #red_instructions_string += '(or deselect) a sky spaxel.'
        #red_instructions_label = gtk.Label(red_instructions_string)
        #red_spax_instructions.add(red_instructions_label)
        #red_instructions_label.show()

        #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        # saving will have its own frame
        red_save_frame = gtk.Frame('Red Data Saving')
        red_spax_big_vbox.pack_start(red_save_frame)
        red_save_frame.show()

        # 'save to file' managed in a vbox
        red_spec_save_vbox = gtk.VBox(False, 5)
        red_save_frame.add(red_spec_save_vbox)
        red_spec_save_vbox.show()
        
        #-+-+
        # hbox1 will have 'save (spax name) to'
        red_spec_save_hbox1 = gtk.HBox(False, 0)
        red_spec_save_vbox.pack_start(red_spec_save_hbox1, False, False, 0)
        red_spec_save_hbox1.show()

        rs_label1 = gtk.Label('Save ')
        red_spec_save_hbox1.pack_start(rs_label1, True, False, 0)
        rs_label1.show()

        # red save menu
        self.red_save_menu = gtk.Menu()

        red_save_gal = gtk.MenuItem('Galaxy Spaxel')
        self.red_save_menu.append(red_save_gal)
        red_save_gal.connect('activate',
                             self.set_red_save_selector, 'Galaxy Spaxel')
        red_save_gal.show()

        red_save_sky = gtk.MenuItem('Sky Spaxel')
        self.red_save_menu.append(red_save_sky)
        red_save_sky.connect('activate',
                             self.set_red_save_selector, 'Sky Spaxel')
        red_save_sky.show()

        red_save_sub = gtk.MenuItem('Subbed Spaxels')
        self.red_save_menu.append(red_save_sub)
        red_save_sub.connect('activate',
                             self.set_red_save_selector, 'Subbed Spaxels')
        red_save_sub.show()

        self.red_save_selector = gtk.Button('Subbed Spaxels')
        self.red_save_selector.connect('pressed', self.popup_red_save_menu)
        red_spec_save_hbox1.pack_start(self.red_save_selector, True, True, 0)
        self.red_save_selector.show()

        rs_label2 = gtk.Label(' to:')
        red_spec_save_hbox1.pack_start(rs_label2, True, False, 0)
        rs_label2.show()
        
        #-+-+
        # then save entry in an hbox2
        red_spec_save_hbox2 = gtk.HBox(False, 0)
        red_spec_save_vbox.pack_start(red_spec_save_hbox2, False, False, 0)
        red_spec_save_hbox2.show()
        
        self.red_spec_save_entry = gtk.Entry()
        self.red_spec_save_entry.set_max_length(50)
        red_spec_save_hbox2.pack_start(self.red_spec_save_entry, False, True, 0)
        self.red_spec_save_entry.show()
        
        #-+-+
        # finally a save button in hbox4
        red_spec_save_hbox4 = gtk.HBox(False, 0)
        red_spec_save_vbox.pack_start(red_spec_save_hbox4, False, False, 0)
        red_spec_save_hbox4.show()
        
        red_spec_save_button = gtk.Button('Save')
        red_spec_save_button.connect('clicked', self.save_red_spec)
        red_spec_save_hbox4.pack_start(red_spec_save_button, True, True, 0)
        red_spec_save_button.show()
        
        #-+-+
        # export to correlator button in hbox5
        red_spec_save_hbox5 = gtk.HBox(False, 0)
        red_spec_save_vbox.pack_start(red_spec_save_hbox5, False, False, 0)
        red_spec_save_hbox5.show()
        
        red_spec_corr_button = gtk.Button('Export to Correlator')
        red_spec_corr_button.connect('clicked', self.export_red_to_prep)
        red_spec_save_hbox5.pack_start(red_spec_corr_button, True, True, 0)
        red_spec_corr_button.show()
        
        #--------------------------------------------------
        #--------------------------------------------------
        # and finally spectra images
        red_spec_vbox = gtk.VBox(False, 0)
        self.red_hbox.pack_start(red_spec_vbox, False, False, 0)
        red_spec_vbox.show()

        # spectra figure
        red_spec_fig = Figure(figsize=(5,5), dpi=50)
        self.red_gal_spec_plot = red_spec_fig.add_subplot(3,1,1)
        self.red_gal_spec_plot.hold(False)
        self.red_sky_spec_plot = red_spec_fig.add_subplot(3,1,2)
        self.red_sky_spec_plot.hold(False)
        self.red_sub_spec_plot = red_spec_fig.add_subplot(3,1,3)
        self.red_sub_spec_plot.hold(False)
        
        self.red_spec_image = FigureCanvas(red_spec_fig)
        self.red_spec_image.set_size_request(500,500)
        red_spec_vbox.pack_start(self.red_spec_image, False, False, 0)
        self.red_spec_image.show()

        # and the corresponding expand button
        red_spec_expand_button = gtk.Button('Expand')
        red_spec_expand_button.connect('clicked', self.expand_red_specs)
        red_spec_vbox.pack_start(red_spec_expand_button, False, False, 0)
        red_spec_expand_button.show()

        #------------------------------------------------------------
        #------------------------------------------------------------
        #   PAGE 2 - BLUE CUBE DATA
        #------------------------------------------------------------
        #------------------------------------------------------------
        self.blue_vbox = gtk.VBox(False, 10)
        self.blue_vbox.show()
        blue_vbox_label = gtk.Label('Blue Data Cube')
        self.notebook.append_page(self.blue_vbox, blue_vbox_label)
        
        #------------------------------------------------
        # Section 1: data selection
        blue_data_select_frame = gtk.Frame('Blue Data Selection')
        self.blue_vbox.pack_start(blue_data_select_frame, False, False, 0)
        blue_data_select_frame.set_border_width(10)
        blue_data_select_frame.show()
        
        # it will all be one vbox
        vbox1 = gtk.VBox(False, 10)
        blue_data_select_frame.add(vbox1)
        vbox1.show()
                
        # make a hbox1c for:
        # -blue file names
        hbox1c = gtk.HBox(False, 0)
        vbox1.pack_start(hbox1c, True, False, 0)
        hbox1c.show()

        # 1c1 - blue file name
        blue_fn_label = gtk.Label('Blue Cube Filename ')
        hbox1c.pack_start(blue_fn_label, False, False, 5)
        blue_fn_label.show()
        
        blue_fn_box = gtk.Entry()
        self.blue_fn_entry = blue_fn_box
        blue_fn_box.set_max_length(500)
        blue_fn_box.connect('activate', self.update_blue_fn)
        hbox1c.pack_start(blue_fn_box, True, True, 5)
        blue_fn_box.show()

        # button to choose blue filename
        blue_filename_button = gtk.Button('Choose...')
        blue_filename_button.connect('clicked', self.select_blue_filename)
        hbox1c.pack_start(blue_filename_button, False, False, 5)
        blue_filename_button.show()

        # hbox for header info
        header_hbox = gtk.HBox(False, 0)
        vbox1.pack_start(header_hbox, True, False, 0)
        header_hbox.show()

        blue_header_label = gtk.Label('Blue Header Field:')
        header_hbox.pack_start(blue_header_label, False, False, 10)
        blue_header_label.show()

        blue_header_selector = gtk.Entry()
        blue_header_selector.set_max_length(8)
        blue_header_selector.set_width_chars(9)
        blue_header_selector.connect('activate', self.update_blue_header_info)
        header_hbox.pack_start(blue_header_selector, False, False, 10)
        blue_header_selector.show()

        blue_header_label2 = gtk.Label('Contents:')
        header_hbox.pack_start(blue_header_label2, False, False, 10)
        blue_header_label2.show()

        self.blue_header_entry = gtk.Entry()
        self.blue_header_entry.set_editable(False)
        header_hbox.pack_start(self.blue_header_entry, True, True, 10)
        self.blue_header_entry.show()

        #------------------------------------------------
        # section2: blue cube image and spectra
        blue_cube_frame = gtk.Frame('Blue Cube')
        self.blue_vbox.pack_start(blue_cube_frame, False, False, 0)
        blue_cube_frame.show()

        # the frame will be managed in an hbox
        self.blue_hbox = gtk.HBox(False, 10)
        blue_cube_frame.add(self.blue_hbox)
        self.blue_hbox.show()

        #-------------------
        #-------------------
        # add a place for cube image
        blue_cube_vbox = gtk.VBox(False, 0)
        self.blue_hbox.pack_start(blue_cube_vbox, False, False, 0)
        blue_cube_vbox.show()
        
        blue_cube_fig = Figure(figsize=(3,5), dpi=50)
        self.blue_cube_plot = blue_cube_fig.add_subplot(1,1,1)
        self.blue_cube_plot.hold(True)        
        self.blue_cube_image = FigureCanvas(blue_cube_fig)
        self.blue_cube_image.mpl_connect('button_press_event',
                                         self.click_blue_cube_image)
        self.blue_cube_image.set_size_request(300,500)
        blue_cube_vbox.pack_start(self.blue_cube_image, False, False, 0)
        self.blue_cube_image.show()

        # colorbar scale adjustor!
        self.blue_im_adj = gtk.Adjustment(1.0, 0.0, 1.0, 0.01, 0.01)
        self.blue_im_adj.connect('value_changed', self.replot_blue_cube_image)
        blue_im_scaler = gtk.HScrollbar(self.blue_im_adj)
        blue_im_scaler.set_update_policy(gtk.UPDATE_CONTINUOUS)
        blue_cube_vbox.pack_start(blue_im_scaler, False, False, 0)
        blue_im_scaler.show()

        # and the corresponding expand button
        blue_cube_expand_button = gtk.Button('Expand')
        blue_cube_expand_button.connect('clicked', self.expand_blue_cube)
        blue_cube_vbox.pack_start(blue_cube_expand_button, False, False, 0)
        blue_cube_expand_button.show()

        #--------------------------------
        #--------------------------------
        # overall vbox for spaxel selection and saving frame
        blue_spax_big_vbox = gtk.VBox(False, 10)
        self.blue_hbox.pack_start(blue_spax_big_vbox, False, False, 0)
        blue_spax_big_vbox.show()

        #-+-+-+-+-+-+-+-+-+
        # frame for setting display wavelength bounds for cube image
        blue_disp_frame = gtk.Frame('Blue Cube Display Options')
        blue_spax_big_vbox.pack_start(blue_disp_frame, False, False, 0)
        blue_disp_frame.show()

        # which will be managed by a vbox
        blue_disp_vbox1 = gtk.VBox(False, 0)
        blue_disp_frame.add(blue_disp_vbox1)
        blue_disp_vbox1.show()

        #-+-+-+
        # hbox for wavelength min
        blue_wmin_hbox = gtk.HBox(False, 0)
        blue_disp_vbox1.pack_start(blue_wmin_hbox, False, False, 0)
        blue_wmin_hbox.show()

        # label
        blue_wmin_label = gtk.Label('Lambda Min:')
        blue_wmin_hbox.pack_start(blue_wmin_label, False, False, 0)
        blue_wmin_label.show()

        # entry
        self.blue_wmin_entry = gtk.Entry()
        blue_wmin_hbox.pack_start(self.blue_wmin_entry, False, False, 0)
        self.blue_wmin_entry.set_width_chars(8)
        self.blue_wmin_entry.show()

        #-+-+-+        
        # hbox for wavelength max
        blue_wmax_hbox = gtk.HBox(False, 0)
        blue_disp_vbox1.pack_start(blue_wmax_hbox, False, False, 0)
        blue_wmax_hbox.show()

        # label
        blue_wmax_label = gtk.Label('Lambda Max:')
        blue_wmax_hbox.pack_start(blue_wmax_label, False, False, 0)
        blue_wmax_label.show()

        # entry
        self.blue_wmax_entry = gtk.Entry()
        blue_wmax_hbox.pack_start(self.blue_wmax_entry, False, False, 0)
        self.blue_wmax_entry.set_width_chars(8)
        self.blue_wmax_entry.show()

        #-+-+-+
        # button to update plot
        blue_disp_button = gtk.Button('Set')
        blue_disp_button.connect('clicked', self.set_blue_disp_wave)
        blue_disp_vbox1.pack_start(blue_disp_button, True, False, 0)
        blue_disp_button.show()
        
        #-+-+-+-+-+-+-+-+-+
        # spaxel selection tools will go in a frame
        blue_spax_frame = gtk.Frame('Blue Spaxel Selection')
        blue_spax_big_vbox.pack_start(blue_spax_frame, False, False, 0)
        blue_spax_frame.show()

        # which will be managed by a vbox
        blue_spax_vbox1 = gtk.VBox(False, 0)
        blue_spax_frame.add(blue_spax_vbox1)
        blue_spax_vbox1.show()

        #-+-+-+-+-+
        # first make a vbox for blue gal spax
        blue_gal_spax_vbox = gtk.VBox(False, 0)
        blue_spax_vbox1.pack_start(blue_gal_spax_vbox, False, False, 10)
        blue_gal_spax_vbox.show()

        #-+-+
        # label in an hbox
        blue_gal_spax_hbox1 = gtk.HBox(False, 0)
        blue_gal_spax_vbox.pack_start(blue_gal_spax_hbox1, False, False, 0)
        blue_gal_spax_hbox1.show()
        
        blue_gal_spax_label = gtk.Label('Blue Galaxy Spaxel(s)')
        blue_gal_spax_hbox1.pack_start(blue_gal_spax_label, True, False, 0)
        blue_gal_spax_label.show()

        #-+-+
        # then entry box
        blue_gal_spax_hbox2 = gtk.HBox(False, 0)
        blue_gal_spax_vbox.pack_start(blue_gal_spax_hbox2, False, False, 0)
        blue_gal_spax_hbox2.show()
        
        self.blue_gal_spax_entry = gtk.Entry()
        #self.blue_gal_spax_entry.set_max_length(3)
        self.blue_gal_spax_entry.connect('activate', self.update_blue_gal_spax)
        blue_gal_spax_hbox2.pack_start(self.blue_gal_spax_entry, False, True, 0)
        self.blue_gal_spax_entry.show()

        #-+-+
        # buttons in an hbox
        blue_gal_spax_hbox3 = gtk.HBox(False, 0)
        blue_gal_spax_vbox.pack_start(blue_gal_spax_hbox3, False, False, 0)
        blue_gal_spax_hbox3.show()
        
        blue_gal_spax_button = gtk.Button('Set')
        blue_gal_spax_button.connect('clicked', self.update_blue_gal_spax)
        blue_gal_spax_hbox3.pack_start(blue_gal_spax_button, True, False, 0)
        blue_gal_spax_button.show()
        
        blue_gal_spax_button = gtk.Button('Clear')
        blue_gal_spax_button.connect('clicked', self.clear_blue_gal_spax)
        blue_gal_spax_hbox3.pack_start(blue_gal_spax_button, True, False, 0)
        blue_gal_spax_button.show()

        #-+-+
        # copy button in another hbox
        blue_gal_spax_hbox4 = gtk.HBox(False, 0)
        blue_gal_spax_vbox.pack_start(blue_gal_spax_hbox4, False, False, 0)
        blue_gal_spax_hbox4.show()
        
        blue_gal_copy_button = gtk.Button('Copy Red Gal Spax')
        blue_gal_copy_button.connect('clicked', self.copy_red_to_blue_gal_spax)
        blue_gal_spax_hbox4.pack_start(blue_gal_copy_button, True, False, 0)
        blue_gal_copy_button.show()

        #-+-+-+-+-+
        # then make a vbox for blue sky spax
        blue_sky_spax_vbox = gtk.VBox(False, 0)
        blue_spax_vbox1.pack_start(blue_sky_spax_vbox, False, False, 5)
        blue_sky_spax_vbox.show()

        #-+-+
        # label in an hbox
        blue_sky_spax_hbox1 = gtk.HBox(False, 0)
        blue_sky_spax_vbox.pack_start(blue_sky_spax_hbox1, False, False, 0)
        blue_sky_spax_hbox1.show()
        
        blue_sky_spax_label = gtk.Label('Blue Sky Spaxel(s)')
        blue_sky_spax_hbox1.pack_start(blue_sky_spax_label, True, False, 0)
        blue_sky_spax_label.show()
        
        #-+-+
        # then entry in an hbox2
        blue_sky_spax_hbox2 = gtk.HBox(False, 0)
        blue_sky_spax_vbox.pack_start(blue_sky_spax_hbox2, False, False, 0)
        blue_sky_spax_hbox2.show()
        
        self.blue_sky_spax_entry = gtk.Entry()
        #self.blue_sky_spax_entry.set_max_length(3)
        self.blue_sky_spax_entry.connect('activate', self.update_blue_sky_spax)
        blue_sky_spax_hbox2.pack_start(self.blue_sky_spax_entry, False, True, 0)
        self.blue_sky_spax_entry.show()

        #-+-+
        # buttons in an hbox
        blue_sky_spax_hbox3 = gtk.HBox(False, 0)
        blue_sky_spax_vbox.pack_start(blue_sky_spax_hbox3, False, False, 0)
        blue_sky_spax_hbox3.show()
        
        blue_sky_spax_button = gtk.Button('Set')
        blue_sky_spax_button.connect('clicked', self.update_blue_sky_spax)
        blue_sky_spax_hbox3.pack_start(blue_sky_spax_button, True, False, 0)
        blue_sky_spax_button.show()
        
        blue_sky_spax_button = gtk.Button('Clear')
        blue_sky_spax_button.connect('clicked', self.clear_blue_sky_spax)
        blue_sky_spax_hbox3.pack_start(blue_sky_spax_button, True, False, 0)
        blue_sky_spax_button.show()
        
        #blue_sky_spax_button = gtk.Button('Autoset')
        #blue_sky_spax_button.connect('clicked', self.autoset_blue_sky_spax)
        #blue_sky_spax_hbox3.pack_start(blue_sky_spax_button, True, False, 0)
        #blue_sky_spax_button.show()

        #-+-+
        # copy button in another hbox
        blue_sky_spax_hbox4 = gtk.HBox(False, 0)
        blue_sky_spax_vbox.pack_start(blue_sky_spax_hbox4, False, False, 0)
        blue_sky_spax_hbox4.show()
        
        blue_sky_copy_button = gtk.Button('Copy Red Sky Spax')
        blue_sky_copy_button.connect('clicked', self.copy_red_to_blue_sky_spax)
        blue_sky_spax_hbox4.pack_start(blue_sky_copy_button, True, False, 0)
        blue_sky_copy_button.show()

        #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        # put instructions in their own frame
        #blue_spax_instructions = gtk.Frame('Instructions')
        #blue_spax_big_vbox.pack_start(blue_spax_instructions, False, False, 0)
        #blue_spax_instructions.show()
        #blue_instructions_string = 'Left click a spaxel to select\n'
        #blue_instructions_string += '(or deselect) a galaxy spaxel,\n'
        #blue_instructions_string += 'right click a spaxel to select\n'
        #blue_instructions_string += '(or deselect) a sky spaxel.'
        #blue_instructions_label = gtk.Label(blue_instructions_string)
        #blue_spax_instructions.add(blue_instructions_label)
        #blue_instructions_label.show()

        #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        # saving will have its own frame
        blue_save_frame = gtk.Frame('Blue Data Saving')
        blue_spax_big_vbox.pack_start(blue_save_frame)
        blue_save_frame.show()

        # 'save to file' managed in a vbox
        blue_spec_save_vbox = gtk.VBox(False, 5)
        blue_save_frame.add(blue_spec_save_vbox)
        blue_spec_save_vbox.show()

        #-+-+
        # hbox1 will have 'save (spax name) to'
        blue_spec_save_hbox1 = gtk.HBox(False, 0)
        blue_spec_save_vbox.pack_start(blue_spec_save_hbox1, False, False, 0)
        blue_spec_save_hbox1.show()

        rs_label1 = gtk.Label('Save ')
        blue_spec_save_hbox1.pack_start(rs_label1, True, False, 0)
        rs_label1.show()

        # blue save menu
        self.blue_save_menu = gtk.Menu()

        blue_save_gal = gtk.MenuItem('Galaxy Spaxel')
        self.blue_save_menu.append(blue_save_gal)
        blue_save_gal.connect('activate',
                              self.set_blue_save_selector, 'Galaxy Spaxel')
        blue_save_gal.show()

        blue_save_sky = gtk.MenuItem('Sky Spaxel')
        self.blue_save_menu.append(blue_save_sky)
        blue_save_sky.connect('activate',
                              self.set_blue_save_selector, 'Sky Spaxel')
        blue_save_sky.show()

        blue_save_sub = gtk.MenuItem('Subbed Spaxels')
        self.blue_save_menu.append(blue_save_sub)
        blue_save_sub.connect('activate',
                              self.set_blue_save_selector, 'Subbed Spaxels')
        blue_save_sub.show()

        self.blue_save_selector = gtk.Button('Subbed Spaxels')
        self.blue_save_selector.connect('pressed', self.popup_blue_save_menu)
        blue_spec_save_hbox1.pack_start(self.blue_save_selector, True, True, 0)
        self.blue_save_selector.show()

        rs_label2 = gtk.Label(' to:')
        blue_spec_save_hbox1.pack_start(rs_label2, True, False, 0)
        rs_label2.show()
        
        #-+-+
        # then save entry in an hbox2
        blue_spec_save_hbox2 = gtk.HBox(False, 0)
        blue_spec_save_vbox.pack_start(blue_spec_save_hbox2, False, False, 0)
        blue_spec_save_hbox2.show()
        
        self.blue_spec_save_entry = gtk.Entry()
        self.blue_spec_save_entry.set_max_length(50)
        blue_spec_save_hbox2.pack_start(
            self.blue_spec_save_entry, False, True, 0)
        self.blue_spec_save_entry.show()
        
        #-+-+
        # finally a save button in hbox4
        blue_spec_save_hbox4 = gtk.HBox(False, 0)
        blue_spec_save_vbox.pack_start(blue_spec_save_hbox4, False, False, 0)
        blue_spec_save_hbox4.show()
        
        blue_spec_save_button = gtk.Button('Save')
        blue_spec_save_button.connect('clicked', self.save_blue_spec)
        blue_spec_save_hbox4.pack_start(blue_spec_save_button, True, True, 0)
        blue_spec_save_button.show()
        
        #-+-+
        # export to correlator button in hbox5
        blue_spec_save_hbox5 = gtk.HBox(False, 0)
        blue_spec_save_vbox.pack_start(blue_spec_save_hbox5, False, False, 0)
        blue_spec_save_hbox5.show()
        
        blue_spec_corr_button = gtk.Button('Export to Correlator')
        blue_spec_corr_button.connect('clicked', self.export_blue_to_prep)
        blue_spec_save_hbox5.pack_start(blue_spec_corr_button, True, True, 0)
        blue_spec_corr_button.show()

        #--------------------------------------------------
        #--------------------------------------------------
        # and finally spectra images
        blue_spec_vbox = gtk.VBox(False, 0)
        self.blue_hbox.pack_start(blue_spec_vbox, False, False, 0)
        blue_spec_vbox.show()

        # spectra figure
        blue_spec_fig = Figure(figsize=(5,5), dpi=50)
        self.blue_gal_spec_plot = blue_spec_fig.add_subplot(3,1,1)
        self.blue_gal_spec_plot.hold(False)
        self.blue_sky_spec_plot = blue_spec_fig.add_subplot(3,1,2)
        self.blue_sky_spec_plot.hold(False)
        self.blue_sub_spec_plot = blue_spec_fig.add_subplot(3,1,3)
        self.blue_sub_spec_plot.hold(False)
        
        self.blue_spec_image = FigureCanvas(blue_spec_fig)
        self.blue_spec_image.set_size_request(500,500)
        blue_spec_vbox.pack_start(self.blue_spec_image, False, False, 0)
        self.blue_spec_image.show()

        # and the corresponding expand button
        blue_spec_expand_button = gtk.Button('Expand')
        blue_spec_expand_button.connect('clicked', self.expand_blue_specs)
        blue_spec_vbox.pack_start(blue_spec_expand_button, False, False, 0)
        blue_spec_expand_button.show()

        #------------------------------------------------------------
        #------------------------------------------------------------
        #   PAGE 3 - SPECTRUM PREPARATION
        #------------------------------------------------------------
        #------------------------------------------------------------
        self.prep_vbox = gtk.VBox(False, 10)
        self.prep_vbox.show()
        prep_vbox_label = gtk.Label('Spectrum Preparation')
        self.notebook.append_page(self.prep_vbox, prep_vbox_label)        

        #---------------------------
        # allow user to select blue and red spectra
        # BLUE
        blue_prep_hbox = gtk.HBox(False, 0)
        self.prep_vbox.pack_start(blue_prep_hbox, False, False, 0)
        blue_prep_hbox.show()
        # label
        blue_prep_label = gtk.Label('Blue Spectrum File')
        blue_prep_hbox.pack_start(blue_prep_label, False, False, 0)
        blue_prep_label.show()
        # entry
        self.blue_prep_entry = gtk.Entry()
        self.blue_prep_entry.set_max_length(500)
        self.blue_prep_entry.connect('activate', self.update_blue_prep)
        blue_prep_hbox.pack_start(self.blue_prep_entry, True, True, 10)
        self.blue_prep_entry.show()        
        # button
        blue_prep_button = gtk.Button('Choose...')
        blue_prep_button.connect('clicked', self.select_blue_prep)
        blue_prep_hbox.pack_start(blue_prep_button, False, False, 5)
        blue_prep_button.show()
        # RED
        red_prep_hbox = gtk.HBox(False, 0)
        self.prep_vbox.pack_start(red_prep_hbox, False, False, 0)
        red_prep_hbox.show()
        # label
        red_prep_label = gtk.Label('Red Spectrum File')
        red_prep_hbox.pack_start(red_prep_label, False, False, 0)
        red_prep_label.show()
        # entry
        self.red_prep_entry = gtk.Entry()
        self.red_prep_entry.set_max_length(500)
        self.red_prep_entry.connect('activate', self.update_red_prep)
        red_prep_hbox.pack_start(self.red_prep_entry, True, True, 10)
        self.red_prep_entry.show()        
        # button
        red_prep_button = gtk.Button('Choose...')
        red_prep_button.connect('clicked', self.select_red_prep)
        red_prep_hbox.pack_start(red_prep_button, False, False, 5)
        red_prep_button.show()

        #------------------------------------------------------
        # hbox for displaying data and scaling red/blue data
        prep_disp_hbox = gtk.HBox(False, 0)
        self.prep_vbox.pack_start(prep_disp_hbox, False, False, 0)
        prep_disp_hbox.show()

        # spectrum image display
        prep_spec_vbox = gtk.VBox(False, 0)
        prep_disp_hbox.pack_start(prep_spec_vbox, False, False, 0)
        prep_spec_vbox.show()
        
        prep_spec_fig = Figure(figsize=(3,6), dpi=50)
        self.prep_spec_plot = prep_spec_fig.add_subplot(1,1,1)
        self.prep_spec_plot.hold(True)        
        self.prep_spec_image = FigureCanvas(prep_spec_fig)
        self.prep_spec_image.set_size_request(800,600)
        prep_spec_vbox.pack_start(self.prep_spec_image, False, False, 0)
        self.prep_spec_image.show()

        # and the corresponding expand button
        prep_spec_expand_button = gtk.Button('Expand')
        prep_spec_expand_button.connect('clicked', self.expand_prep_spec)
        prep_spec_vbox.pack_start(prep_spec_expand_button, False, False, 0)
        prep_spec_expand_button.show()        
        
        #---------------------------
        # scale blue/red channels!!!
        prep_scale_vbox = gtk.VBox(False, 0)
        prep_disp_hbox.pack_start(prep_scale_vbox, False, False, 10)
        prep_scale_vbox.show()

        prep_scale_label = gtk.Label('Red / Blue Scaling:')
        prep_scale_vbox.pack_start(prep_scale_label, False, False, 0)
        prep_scale_label.show()

        prep_scale_hbox = gtk.HBox(False, 0)
        prep_scale_vbox.pack_start(prep_scale_hbox, False, False, 0)
        prep_scale_hbox.show()

        # entry
        self.prep_scale_entry = gtk.Entry()
        self.prep_scale_entry.set_width_chars(8)
        self.prep_scale_entry.set_text('1.0')
        self.prep_scale_entry.connect('activate', self.adjust_prep_scaling)
        prep_scale_hbox.pack_start(self.prep_scale_entry, False, False, 0)
        self.prep_scale_entry.show()

        # button
        prep_scale_button = gtk.Button('Set')
        prep_scale_button.show()
        prep_scale_button.connect('clicked', self.adjust_prep_scaling)
        prep_scale_hbox.pack_start(prep_scale_button, False, False, 0)

        #------------------------------------------------------------
        #------------------------------------------------------------
        #   PAGE 4 - SPECTRUM CORRELATION
        #------------------------------------------------------------
        #------------------------------------------------------------
        self.corr_vbox = gtk.VBox(False, 10)
        self.corr_vbox.show()
        corr_vbox_label = gtk.Label('Spectrum Correlation')
        self.notebook.append_page(self.corr_vbox, corr_vbox_label)        

        #------------------------------------
        # output showing continuum-subtracted, apodized spectrum
        # ALSO with current SNID template of interest!

        #++++ PLOT / ETC. HBOX
        corr_plot_hbox = gtk.HBox(False, 0)
        self.corr_vbox.pack_start(corr_plot_hbox, False, False, 5)
        corr_plot_hbox.show()
        #+--- stored in vbox
        corr_plot_vbox = gtk.VBox(False, 0)
        corr_plot_hbox.pack_start(corr_plot_vbox, False, False, 0)
        corr_plot_vbox.show()
        #+- Figure for plotting
        corr_plot_fig = Figure(figsize=(9,5), dpi=50)
        self.corr_plot = corr_plot_fig.add_subplot(1,1,1)
        self.corr_plot.hold(True)
        self.corr_image = FigureCanvas(corr_plot_fig)
        self.corr_image.set_size_request(700,450)
        corr_plot_vbox.pack_start(self.corr_image, False, False, 0)
        self.corr_image.show()
        #+- Navigation bar
        corr_toolbar = NavigationToolbar(self.corr_image, corr_plot_fig)
        corr_plot_vbox.pack_start(corr_toolbar, False, False, 0)
        corr_toolbar.show()
        #+- Expand button
        corr_expand_button = gtk.Button('Expand')
        corr_expand_button.connect('clicked', self.expand_corr_plot)
        corr_plot_vbox.pack_start(corr_expand_button, False, False, 0)
        corr_expand_button.show()
        #------------------------
        #+--- Vbox to hold it all!
        corr_logistics_vbox = gtk.VBox(False, 0)
        corr_plot_hbox.pack_start(corr_logistics_vbox, False, False, 0)
        corr_logistics_vbox.show()
        #+-- Frame for zmax and 'run snid' button
        corr_run_frame = gtk.Frame('SNID Input')
        corr_logistics_vbox.pack_start(corr_run_frame, False, False, 0)
        corr_run_frame.set_border_width(10)
        corr_run_frame.show()
        #+-- Vbox...
        corr_run_vbox = gtk.VBox(False, 5)
        corr_run_frame.add(corr_run_vbox)
        corr_run_vbox.show()
        #+- Hbox for setting zmax (label + entry)
        corr_zmax_hbox = gtk.HBox(False, 5)
        corr_run_vbox.pack_start(corr_zmax_hbox, False, False, 5)
        corr_zmax_hbox.show()
        corr_zmax_label = gtk.Label('z Max: ')
        corr_zmax_hbox.pack_start(corr_zmax_label, False, False, 5)
        corr_zmax_label.show()
        self.corr_zmax_entry = gtk.Entry()
        self.corr_zmax_entry.set_width_chars(10)
        corr_zmax_hbox.pack_start(self.corr_zmax_entry, False, False, 5)
        self.corr_zmax_entry.show()
        #+- 'Run SNID' Button
        run_snid_button = gtk.Button('Run SNID')
        run_snid_button.connect('clicked', self.run_snid)
        corr_run_vbox.pack_start(run_snid_button, False, False, 5)
        run_snid_button.show()
        #+- 'Clear SNID Results' Button
        clear_snid_button = gtk.Button('Clear SNID Results')
        clear_snid_button.connect('clicked', self.clear_snid_results)
        corr_run_vbox.pack_start(clear_snid_button, False, False, 5)
        clear_snid_button.show()
        #------------------------
        #+-- Frame SNID OUTPUT
        corr_out_frame = gtk.Frame('SNID Output')
        corr_logistics_vbox.pack_start(corr_out_frame, False, False, 0)
        corr_out_frame.set_border_width(10)
        corr_out_frame.show()
        #+-- Vbox...
        corr_out_vbox = gtk.VBox(False, 5)
        corr_out_frame.add(corr_out_vbox)
        corr_out_vbox.show()
        #+- display output best redshift and error
        corr_zbest_hbox = gtk.HBox(False, 5)
        corr_out_vbox.pack_start(corr_zbest_hbox, False, False, 5)
        corr_zbest_hbox.show()
        corr_zbest_label = gtk.Label('Best Redshift:')
        corr_zbest_hbox.pack_start(corr_zbest_label, False, False, 5)
        corr_zbest_label.show()
        self.corr_zbest_entry = gtk.Entry()
        self.corr_zbest_entry.set_width_chars(20)
        self.corr_zbest_entry.set_editable(False)
        corr_zbest_hbox.pack_start(self.corr_zbest_entry, False, False, 5)
        self.corr_zbest_entry.show()
        #+- display output best phase and error
        corr_tbest_hbox = gtk.HBox(False, 5)
        corr_out_vbox.pack_start(corr_tbest_hbox, False, False, 5)
        corr_tbest_hbox.show()
        corr_tbest_label = gtk.Label('Best Phase:')
        corr_tbest_hbox.pack_start(corr_tbest_label, False, False, 5)
        corr_tbest_label.show()
        self.corr_tbest_entry = gtk.Entry()
        self.corr_tbest_entry.set_width_chars(20)
        self.corr_tbest_entry.set_editable(False)
        corr_tbest_hbox.pack_start(self.corr_tbest_entry, False, False, 5)
        self.corr_tbest_entry.show()

        #------------------------------------
        #++++ RESULTS SCROLLBOX
        corr_results_sbox = gtk.ScrolledWindow()
        self.corr_vbox.pack_start(corr_results_sbox, False, False, 5)
        corr_results_sbox.set_policy(gtk.POLICY_NEVER, gtk.POLICY_ALWAYS)
        corr_results_sbox.set_size_request(500, 200)
        corr_results_sbox.show()
        # SNID results will use TreeView widget!!
        self.snid_treestore = gtk.ListStore(
            str, str, float, float, float, float)
        self.snid_treeview = gtk.TreeView(self.snid_treestore)

        # set up columns and allow ordering
        # - Col 0: SN Name
        self.snid_col0 = gtk.TreeViewColumn('SN Name')
        self.snid_treeview.append_column(self.snid_col0)
        self.snid_cell0 = gtk.CellRendererText()
        self.snid_col0.pack_start(self.snid_cell0, True)
        self.snid_col0.add_attribute(self.snid_cell0, 'text', 0)
        self.snid_col0.set_sort_column_id(0)
        # - Col 1: SN Type
        self.snid_col1 = gtk.TreeViewColumn('SN Type')
        self.snid_treeview.append_column(self.snid_col1)
        self.snid_cell1 = gtk.CellRendererText()
        self.snid_col1.pack_start(self.snid_cell1, True)
        self.snid_col1.add_attribute(self.snid_cell1, 'text', 1)
        self.snid_col1.set_sort_column_id(1)
        # - Col 2: Best Epoch
        self.snid_col2 = gtk.TreeViewColumn('Epoch')
        self.snid_treeview.append_column(self.snid_col2)
        self.snid_cell2 = gtk.CellRendererText()
        self.snid_col2.pack_start(self.snid_cell2, True)
        self.snid_col2.add_attribute(self.snid_cell2, 'text', 2)
        self.snid_col2.set_sort_column_id(2)
        # - Col 3: Best z
        self.snid_col3 = gtk.TreeViewColumn('Best z')
        self.snid_treeview.append_column(self.snid_col3)
        self.snid_cell3 = gtk.CellRendererText()
        self.snid_col3.pack_start(self.snid_cell3, True)
        self.snid_col3.add_attribute(self.snid_cell3, 'text', 3)
        self.snid_col3.set_sort_column_id(3)
        # - Col 4: z Error
        self.snid_col4 = gtk.TreeViewColumn('z Err')
        self.snid_treeview.append_column(self.snid_col4)
        self.snid_cell4 = gtk.CellRendererText()
        self.snid_col4.pack_start(self.snid_cell4, True)
        self.snid_col4.add_attribute(self.snid_cell4, 'text', 4)
        self.snid_col4.set_sort_column_id(4)
        # - Col 5: rlap
        self.snid_col5 = gtk.TreeViewColumn('rlap')
        self.snid_treeview.append_column(self.snid_col5)
        self.snid_cell5 = gtk.CellRendererText()
        self.snid_col5.pack_start(self.snid_cell5, True)
        self.snid_col5.add_attribute(self.snid_cell5, 'text', 5)
        self.snid_col5.set_sort_column_id(5)

        self.snid_treeview.set_search_column(5)
        self.snid_treeview.set_reorderable(True)
        
        # when user changes selection in TreeView, update the plot!
        self.snid_selection = self.snid_treeview.get_selection()
        self.snid_selection.connect('changed', self.update_snid_plot)

        # finally add to the hbox!
        corr_results_sbox.add_with_viewport(self.snid_treeview)
        self.snid_treeview.show()

        #------------------------------------------------------------
        #------------------------------------------------------------
        #   FOOTER: ERROR REPORTING
        #------------------------------------------------------------
        #------------------------------------------------------------
        error_report_frame = gtk.Frame('Error Reporting')
        self.vbox.pack_start(error_report_frame, False, False, 0)
        error_report_frame.show()

        # make an error report entry box
        error_msg_hbox = gtk.HBox(False, 0)
        error_report_frame.add(error_msg_hbox)
        error_msg_hbox.show()
        
        self.error_message_box = gtk.Entry()
        self.error_message_box.set_editable(False)
        error_msg_hbox.pack_start(self.error_message_box, True, True, 10)
        self.error_message_box.show()

        error_msg_clear = gtk.Button('Clear')
        error_msg_clear.connect('clicked', self.clear_error_box)
        error_msg_hbox.pack_start(error_msg_clear, False, False, 0)
        error_msg_clear.show()
        
        #---------------------
        # and finally show window
        self.window.show()
        
    #--------------------------------------------
    # death instructions
    #--------------------------------------------
    def delete_event(self, widget, event, data = None):
        #print 'take me to your leader!'
        # if False will also run destroy event
        return False
    
    def delete_expanded_plot(self, widget, event, data = None):
        #print 'take me to your leader!'
        # if False will also run destroy event
        widget.hide()
        return True

    def destroy(self, widget, data = None):
        gtk.main_quit()

    #-------------------------------------
    # real functions
    #-------------------------------------
    # filename functions
    def update_blue_fn(self, widget, input_blue_fn = None):
        if input_blue_fn == None:
            new_blue_fn = self.blue_fn_entry.get_text()
        else:
            new_blue_fn = input_blue_fn
        # now check that the fn is for real
        if os.path.isfile(new_blue_fn):
            self.blue_cube_fn = new_blue_fn
            self.blue_fn_entry.set_text(new_blue_fn)
            #self.update_object_name()
            self.update_blue_data(widget)
        else:
            error_msg = 'Invalid Blue Cube Filename ' + new_blue_fn
            self.error_announce(error_msg)

    def update_red_fn(self, widget, input_red_fn = None):
        if input_red_fn == None:
            new_red_fn = self.red_fn_entry.get_text()
        else:
            new_red_fn = input_red_fn
        # now check that the fn is for real
        if os.path.isfile(new_red_fn):
            self.red_cube_fn = new_red_fn
            self.red_fn_entry.set_text(new_red_fn)
            #self.update_object_name()
            self.update_red_data(widget)
        else:
            error_msg = 'Invalid Red Cube Filename  ' + new_red_fn
            self.error_announce(error_msg)

    def select_blue_filename(self, widget):
        # open file browser dialog
        dialog = gtk.FileChooserDialog(
            title='Select WiFeS Blue Data Cube',
            action=gtk.FILE_CHOOSER_ACTION_OPEN,
            buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                     gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_CANCEL)
        response = dialog.run()
        # update the filename if a filename was chosen
        if response == gtk.RESPONSE_OK:
            new_blue_fn = dialog.get_filename()
            self.update_blue_fn(widget, new_blue_fn)
        dialog.destroy()
        return

    def select_red_filename(self, widget):
        # open file browser dialog
        dialog = gtk.FileChooserDialog(
            title='Select WiFeS Red Data Cube',
            action=gtk.FILE_CHOOSER_ACTION_OPEN,
            buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                     gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_CANCEL)
        response = dialog.run()
        # update the filename if a filename was chosen
        if response == gtk.RESPONSE_OK:
            new_red_fn = dialog.get_filename()
            self.update_red_fn(widget, new_red_fn)
        dialog.destroy()
        return

    #-----------------------------
    def toggle_red_data_display(self, red_toggle_widget):
        if red_toggle_widget.get_active():
            self.red_hbox.hide()
            red_toggle_widget.set_label('Show Red Data?')
        else:
            self.red_hbox.show()
            red_toggle_widget.set_label(' Hide Red Data?')
            
    def toggle_blue_data_display(self, blue_toggle_widget):
        if blue_toggle_widget.get_active():
            self.blue_hbox.hide()
            blue_toggle_widget.set_label('Show Blue Data?')
        else:
            self.blue_hbox.show()
            blue_toggle_widget.set_label('Hide Blue Data?')

    #-----------------------------
    def update_red_header_info(self, entry_widget):
        requested_key = entry_widget.get_text()
        try:
            header_value = self.red_cube_header[requested_key]
            self.red_header_entry.set_text(str(header_value))
        except:
            error_msg = 'Red header key ' + requested_key + 'not known.'
            
    def update_blue_header_info(self, entry_widget):
        requested_key = entry_widget.get_text()
        try:
            header_value = self.blue_cube_header[requested_key]
            self.blue_header_entry.set_text(str(header_value))
        except:
            error_msg = 'Blue header key ' + requested_key + 'not known.'

    #----------------------------------------------
    # Color Scheme Selection Functions
    #----------------------------------------------
    #def setup_color_scheme_menu(self):
    #    keylist = self.color_schemes.keys()
    #    for scheme in keylist:
    #        color_scheme_item = gtk.MenuItem(scheme)
    #        self.color_scheme_menu.append(color_scheme_item)
    #        color_scheme_item.connect('activate', self.set_color_scheme_selector, scheme)
    #        color_scheme_item.show()
    #
    #def popup_color_scheme_menu(self, widget):
    #    self.color_scheme_menu.popup(None, None, None, 0, 0)
    #    
    #def set_color_scheme_selector(self, widget, choice_string):
    #    self.chosen_color_scheme = choice_string
    #    self.color_scheme_selector.set_label(self.chosen_color_scheme)
    #    # replot all the plots
    #    self.replot_red_cube_image()
    #    self.replot_red_specs()
    #    self.replot_blue_cube_image()
    #    self.replot_blue_specs()

    #------------------------------    
    # WIFES DATA GRABBING FUNCTIONS
    #------------------------------
    def extract_wifes_cube_data(self, cube_fn):
        f = pyfits.open(cube_fn)
        if len(f) == 76:
            ny,nlam = numpy.shape(f[1].data)
            nx=25
            # get wavelength array
            lam0 = f[1].header['CRVAL1']
            dlam = f[1].header['CDELT1']
            lam_array = lam0+dlam*numpy.arange(nlam,dtype='d')
            # get data and variance
            obj_cube_data = numpy.zeros([nlam,ny,nx],dtype='d')
            obj_cube_var  = numpy.zeros([nlam,ny,nx],dtype='d')
            for i in range(nx):
                curr_data = f[i+1].data
                curr_var = f[i+26].data
                obj_cube_data[:,:,nx-i-1] = curr_data.T
                obj_cube_var[:,:,nx-i-1] = curr_var.T
        else:
            nlam, ny, nx = numpy.shape(f[1].data)
            # get wavelength array
            lam0 = f[1].header['CRVAL3']
            dlam = f[1].header['CDELT3']
            lam_array = lam0+dlam*numpy.arange(nlam,dtype='d')
            # get data and variance
            obj_cube_data = numpy.zeros([nlam,ny,nx],dtype='d')
            obj_cube_var  = numpy.zeros([nlam,ny,nx],dtype='d')
            for i in range(nx):
                curr_data = f[1].data[:,:,i]
                curr_var = f[2].data[:,:,i]
                obj_cube_data[:,:,nx-i-1] = curr_data
                obj_cube_var[:,:,nx-i-1] = curr_var
            
        f.close()
        # return flux, variance, wavelength
        return obj_cube_data, obj_cube_var, lam_array

    #------------------------------    
    # Red Data Update Functions
    #------------------------------    
    def update_red_data(self, widget):
        # first resave internal data from the red filename
        f1 = pyfits.open(self.red_cube_fn)
        self.red_cube_header = dict(f1[1].header.items())
        f1.close()
        data, var, lam = self.extract_wifes_cube_data(self.red_cube_fn)
        self.red_wavelength_array = lam
        self.red_cube_data = data
        self.red_var_cube_data = var
        # put a default save name
        obj_name = self.red_cube_header['OBJECT']
        red_save_fn = obj_name + '_wifes_R.dat'
        self.red_spec_save_entry.set_text(red_save_fn)
        # then resave the spaxel numbers
        self.resave_red_gal_spax()
        self.resave_red_sky_spax()
        # finally replot the data
        self.replot_red_cube_image()
        self.replot_red_specs()

    def update_red_gal_spax(self, widget):
        self.resave_red_gal_spax()
        self.replot_red_cube_image()
        self.replot_red_gal_spec()
        self.replot_red_sub_spec()

    def update_red_sky_spax(self, widget):
        self.resave_red_sky_spax()
        self.replot_red_cube_image()
        self.replot_red_sky_spec()
        self.replot_red_sub_spec()

    #--------------------
    def clear_red_gal_spax(self, widget):
        self.red_gal_spax_entry.set_text('')
        self.update_red_gal_spax(widget)
        
    def clear_red_sky_spax(self, widget):
        self.red_sky_spax_entry.set_text('')
        self.update_red_sky_spax(widget)

    #--------------------
    def copy_blue_to_red_gal_spax(self, widget):
        blue_gal_spax = self.blue_gal_spax_entry.get_text()
        self.red_gal_spax_entry.set_text(blue_gal_spax)
        self.update_red_gal_spax(widget)
        
    def copy_blue_to_red_sky_spax(self, widget):
        blue_sky_spax = self.blue_sky_spax_entry.get_text()
        self.red_sky_spax_entry.set_text(blue_sky_spax)
        self.update_red_sky_spax(widget)

    #--------------------
    def autofind_red_sky_spax(self):
        cube_data = self.red_cube_data
        cube_var  = self.red_var_cube_data
        return self.autoselect_sky_spaxels(cube_data, cube_var)

    def autoset_red_sky_spax(self, widget):
        red_sky_spaxels = self.autofind_red_sky_spax()
        red_sky_spax_string = ';'.join([str(x) for x in
                                        red_sky_spaxels])
        self.red_sky_spax_entry.set_text(red_sky_spax_string)
        self.update_red_sky_spax(widget)

    #--------------------
    def resave_red_gal_spax(self):
        # okay now to get multiple spaxels i am gonna hope
        # the user doesn't input garbage
        red_gal_entry_contents = self.red_gal_spax_entry.get_text()
        if len(red_gal_entry_contents) == 0:
            self.red_gal_spax = None
            return
        red_gal_spax_list = re.split(';', red_gal_entry_contents)
        try:
            #self.red_gal_spax = [int(x) for x in red_gal_spax_list]
            self.red_gal_spax = red_gal_spax_list
        except:
            self.red_gal_spax = None

    def resave_red_sky_spax(self):
        # okay now to get multiple spaxels i am gonna hope
        # the user doesn't input garbage
        red_sky_entry_contents = self.red_sky_spax_entry.get_text()
        if len(red_sky_entry_contents) == 0:
            self.red_sky_spax = None
            return
        red_sky_spax_list = re.split(';', red_sky_entry_contents)
        try:
            #self.red_sky_spax = [int(x) for x in red_sky_spax_list]
            self.red_sky_spax = red_sky_spax_list
        except:
            self.red_sky_spax = None        

    #------------------------------    
    # Blue Data Update Functions
    #------------------------------
    def update_blue_data(self, widget):
        # first resave internal data from the blue filename
        f1 = pyfits.open(self.blue_cube_fn)
        self.blue_cube_header = dict(f1[1].header.items())
        f1.close()
        data, var, lam = self.extract_wifes_cube_data(self.blue_cube_fn)
        self.blue_wavelength_array = lam
        self.blue_cube_data = data
        self.blue_var_cube_data = var
        # put a default save name
        obj_name = self.blue_cube_header['OBJECT']
        blue_save_fn = obj_name + '_wifes_B.dat'
        self.blue_spec_save_entry.set_text(blue_save_fn)
        # then resave the spaxel numbers
        self.resave_blue_gal_spax()
        self.resave_blue_sky_spax()
        # finally replot the data
        self.replot_blue_cube_image()
        self.replot_blue_specs()

    def update_blue_gal_spax(self, widget):
        self.resave_blue_gal_spax()
        self.replot_blue_cube_image()
        self.replot_blue_gal_spec()
        self.replot_blue_sub_spec()

    def update_blue_sky_spax(self, widget):
        self.resave_blue_sky_spax()
        self.replot_blue_cube_image()
        self.replot_blue_sky_spec()
        self.replot_blue_sub_spec()

    #--------------------
    def clear_blue_gal_spax(self, widget):
        self.blue_gal_spax_entry.set_text('')
        self.update_blue_gal_spax(widget)
        
    def clear_blue_sky_spax(self, widget):
        self.blue_sky_spax_entry.set_text('')
        self.update_blue_sky_spax(widget)

    #--------------------
    def copy_red_to_blue_gal_spax(self, widget):
        red_gal_spax = self.red_gal_spax_entry.get_text()
        self.blue_gal_spax_entry.set_text(red_gal_spax)
        self.update_blue_gal_spax(widget)
        
    def copy_red_to_blue_sky_spax(self, widget):
        red_sky_spax = self.red_sky_spax_entry.get_text()
        self.blue_sky_spax_entry.set_text(red_sky_spax)
        self.update_blue_sky_spax(widget)

    #--------------------
    def autofind_blue_sky_spax(self):
        cube_data = self.blue_cube_data
        cube_var  = self.blue_var_cube_data
        return self.autoselect_sky_spaxels(cube_data, cube_var)

    def autoset_blue_sky_spax(self, widget):
        blue_sky_spaxels = self.autofind_blue_sky_spax()
        blue_sky_spax_string = ';'.join([str(x) for x in
                                         blue_sky_spaxels])
        self.blue_sky_spax_entry.set_text(blue_sky_spax_string)
        self.update_blue_sky_spax(widget)

    #--------------------
    def resave_blue_gal_spax(self):
        # okay now to get multiple spaxels i am gonna hope
        # the user doesn't input garbage
        blue_gal_entry_contents = self.blue_gal_spax_entry.get_text()
        if len(blue_gal_entry_contents) == 0:
            self.blue_gal_spax = None
            return
        blue_gal_spax_list = re.split(';', blue_gal_entry_contents)
        try:
            #self.blue_gal_spax = [int(x) for x in blue_gal_spax_list]
            self.blue_gal_spax = blue_gal_spax_list
        except:
            self.blue_gal_spax = None

    def resave_blue_sky_spax(self):
        # okay now to get multiple spaxels i am gonna hope
        # the user doesn't input garbage
        blue_sky_entry_contents = self.blue_sky_spax_entry.get_text()
        if len(blue_sky_entry_contents) == 0:
            self.blue_sky_spax = None
            return
        blue_sky_spax_list = re.split(';', blue_sky_entry_contents)
        try:
            #self.blue_sky_spax = [int(x) for x in blue_sky_spax_list]
            self.blue_sky_spax = blue_sky_spax_list
        except:
            self.blue_sky_spax = None

    #-----------------------------------------------------------
    # Common data functions
    #-----------------------------------------------------------
    def find_sky_indices(self, data, var):
        # first remove places where i set var=0.0 (dead spaxels)
        data = data[numpy.where(var!= 0.0)[0]]
        var  =  var[numpy.where(var!= 0.0)[0]]
        # for 1d data compute indices that fall within
        # 4sigma of each other
        ind0 = numpy.arange(len(data))
        ind = numpy.arange(len(data))
        prev_npts = 0
        niter_max=10
        niter = 0
        while(len(ind) != prev_npts and niter<niter_max):
            prev_npts = len(ind)
            ave_flux = numpy.sum(data/var)/numpy.sum(1/var)
            sigma = numpy.median(var**0.5)
            ind = numpy.where(abs(data-ave_flux)<4*sigma)[0]
            data = data[ind]
            var = var[ind]
            ind0 = ind0[ind]
            niter += 1
        return ind0

    def autoselect_sky_spaxels(self, cube_data, cube_var):
        # this first part finds the frequency with which
        # each spaxel is selected as sky in the above function
        cube_shape = numpy.shape(cube_data)
        npts = cube_shape[1]
        nspax = cube_shape[0]
        selection_sum = numpy.zeros(nspax, dtype = 'f')
        selector = numpy.zeros(nspax, dtype = 'f')
        for i in range(npts):
            selector *= 0
            sky_indices = self.find_sky_indices(cube_data[:,i],cube_var[:,i])
            selector[sky_indices] = 1.0
            selection_sum += selector
        sel_freq = selection_sum / npts
        # the next part finds spaxels that self-consistently fall
        # within 3sigma of 1 (always a sky spax)
        ind0 = numpy.arange(len(sel_freq))
        ind  = numpy.arange(len(sel_freq))
        temp_sf = sel_freq[:]
        prev_npts = 0
        niter_max=10
        niter = 0
        try:
            while(len(ind) != prev_npts and niter<niter_max):
                prev_npts = len(ind)
                sigma = numpy.median(1.0 - temp_sf)
                ind = numpy.where(abs(1.0 - temp_sf)<3.0*sigma)[0]
                temp_sf = temp_sf[ind]
                ind0 = ind0[ind]
                niter += 1
        except:
            sigma = 0.0
        # return the list of sky spaxels with 2 sigma
        nsig = 2.0
        return numpy.where(sel_freq >= (1.0 - nsig*sigma))[0]
    
    #--------------------------------------------
    # Red Cube image fucntions
    #--------------------------------------------
    def set_red_disp_wave(self, widget=None):
        # get wave wmin and wmax
        try:
            init_red_wmin = float(self.red_wmin_entry.get_text())
        except:
            init_red_wmin = 0.0
        red_wmin = max(init_red_wmin, numpy.min(self.red_wavelength_array))
        try:
            init_red_wmax = float(self.red_wmax_entry.get_text())
        except:
            init_red_wmax = 20000.0
        red_wmax = min(init_red_wmax, numpy.max(self.red_wavelength_array))
        # set internal variables
        self.red_disp_wmin = red_wmin
        self.red_disp_wmax = red_wmax
        # replot!
        self.replot_red_cube_image()
    
    def replot_red_cube_image(self, widget=None):
        # first clear the figure
        self.red_cube_plot.clear()
        self.plot_red_cube_image(self.red_cube_plot)
                
        # finally tell it to draw
        self.red_cube_image.draw()

    def expand_red_cube(self, widget):
        expanded_rcube = gtk.Window(gtk.WINDOW_TOPLEVEL)
        expanded_rcube.set_title('Red Cube Summed Spaxels')

        # what to do when 'delete event' signal comes
        # defined in self function below
        expanded_rcube.connect('delete_event', self.delete_expanded_plot)
        # what to do with a destroy event
        expanded_rcube.connect('destroy', self.delete_expanded_plot)

        # first make a the figure
        expanded_red_fig = Figure(figsize = (3,6), dpi = 50)
        ex_red_im = FigureCanvas(expanded_red_fig)
        expanded_rcube.add(ex_red_im)
        ex_red_im.set_size_request(600,900)
        ex_red_im.show()

        # now make a plot space and plot in that widget
        expanded_red_cube_plot = expanded_red_fig.add_subplot(1,1,1)
        self.plot_red_cube_image(expanded_red_cube_plot)
        
        # finally tell it to draw
        ex_red_im.draw()
        #pylab.show()
        expanded_rcube.show()

    def plot_red_cube_image(self, target_plot_widget):
        # if not data, skip out
        if self.red_cube_data == None:
            return
        # now set wavelength limits
        if self.red_disp_wmin == None:
            red_disp_imin = 0
        else:
            red_disp_imin = numpy.nonzero(
                self.red_wavelength_array >= self.red_disp_wmin)[0][0]
        if self.red_disp_wmax == None:
            red_disp_imax = len(self.red_wavelength_array)
        else:
            red_disp_imax = numpy.nonzero(
                self.red_wavelength_array <= self.red_disp_wmax)[0][-1]+1
        # assume the spaxel numbers have been correctly set
        # do the replotting of the image
        red_spax_im = numpy.sum(
            self.red_cube_data[red_disp_imin:red_disp_imax,:,:],axis=0)
        init_vmax = numpy.max(red_spax_im)
        set_vmax = init_vmax*self.red_im_adj.value
        ny, nx = numpy.shape(red_spax_im)
        edge_box = [0, nx, 0, ny]

        # finally plot the image
        chosen_cmap = self.color_schemes[self.chosen_color_scheme]['cmap']
        target_plot_widget.imshow(red_spax_im,
                                  cmap = matplotlib.cm.get_cmap(chosen_cmap),
                                  interpolation = 'nearest',
                                  extent=edge_box,
                                  vmax=set_vmax,
                                  vmin=0.0,
                                  origin='lower')

        # then put boxes for galaxy spaxels
        # again i am gonna assume the user input real stuff
        red_gal_spax = self.red_gal_spax
        if red_gal_spax != None:
            for gspax in red_gal_spax:
                # goes from spax number to [i,j] pair
                xstr, ystr = gspax[1:-1].split(',')
                gi = int(float(ystr))
                gj = int(float(xstr))
                #print gi, gj
                red_gal_box_x = [gj, gj,   gj+1, gj+1, gj]
                red_gal_box_y = [gi, gi+1, gi+1, gi,   gi]
                gal_color = self.color_schemes[self.chosen_color_scheme]['gal_color']
                target_plot_widget.plot(red_gal_box_x, red_gal_box_y,
                                        color = gal_color, linewidth = 2)

         # ... and sky spax
        red_sky_spax = self.red_sky_spax
        if red_sky_spax != None:    
            for sspax in red_sky_spax:
                # goes from spax number to [i,j] pair
                xstr, ystr = sspax[1:-1].split(',')
                si = int(float(ystr))
                sj = int(float(xstr))
                #print si, si
                red_sky_box_x = [sj, sj,   sj+1, sj+1, sj]
                red_sky_box_y = [si, si+1, si+1, si,   si]
                sky_color = self.color_schemes[self.chosen_color_scheme]['sky_color']
                target_plot_widget.plot(red_sky_box_x, red_sky_box_y,
                                        color = sky_color, linewidth = 2)

        # set the axis
        target_plot_widget.axis(edge_box)

    #--------------------------------------------
    # Blue Cube image fucntions
    #--------------------------------------------
    def set_blue_disp_wave(self, widget=None):
        # get wave wmin and wmax
        try:
            init_blue_wmin = float(self.blue_wmin_entry.get_text())
        except:
            init_blue_wmin = 0.0
        blue_wmin = max(init_blue_wmin, numpy.min(self.blue_wavelength_array))
        try:
            init_blue_wmax = float(self.blue_wmax_entry.get_text())
        except:
            init_blue_wmax = 20000.0
        blue_wmax = min(init_blue_wmax, numpy.max(self.blue_wavelength_array))
        # set internal variables
        self.blue_disp_wmin = blue_wmin
        self.blue_disp_wmax = blue_wmax
        # replot!
        self.replot_blue_cube_image()

    def replot_blue_cube_image(self, widget=None):
        # first clear the figure
        self.blue_cube_plot.clear()

        # plot in that space
        self.plot_blue_cube_image(self.blue_cube_plot)
                
        # finally tell it to draw
        self.blue_cube_image.draw()

    def expand_blue_cube(self, widget):
        expanded_rcube = gtk.Window(gtk.WINDOW_TOPLEVEL)
        expanded_rcube.set_title('Blue Cube Summed Spaxels')

        # what to do when 'delete event' signal comes
        # defined in self function below
        expanded_rcube.connect('delete_event', self.delete_expanded_plot)
        # what to do with a destroy event
        expanded_rcube.connect('destroy', self.delete_expanded_plot)
        # first make a the figure
        expanded_blue_fig = Figure(figsize = (3,6), dpi = 50)
        
        ex_blue_im = FigureCanvas(expanded_blue_fig)
        expanded_rcube.add(ex_blue_im)
        ex_blue_im.set_size_request(600,900)
        ex_blue_im.show()
        
        expanded_blue_cube_plot = expanded_blue_fig.add_subplot(1,1,1)
        self.plot_blue_cube_image(expanded_blue_cube_plot)
        
        # finally tell it to draw
        ex_blue_im.draw()
        #pylab.show()
        expanded_rcube.show()

    def plot_blue_cube_image(self, target_plot_widget):
        # if no data skip out
        if self.blue_cube_data == None:
            return
        # now set wavelength limits
        if self.blue_disp_wmin == None:
            blue_disp_imin = 0
        else:
            blue_disp_imin = numpy.nonzero(
                self.blue_wavelength_array >= self.blue_disp_wmin)[0][0]
        if self.blue_disp_wmax == None:
            blue_disp_imax = len(self.blue_wavelength_array)
        else:
            blue_disp_imax = numpy.nonzero(
                self.blue_wavelength_array <= self.blue_disp_wmax)[0][-1]+1
        # assume the spaxel numbers have been correctly set
        # do the replotting of the image
        blue_spax_im = numpy.sum(
            self.blue_cube_data[blue_disp_imin:blue_disp_imax,:,:],axis=0)
        init_vmax = numpy.max(blue_spax_im)
        set_vmax = init_vmax*self.blue_im_adj.value
        ny, nx = numpy.shape(blue_spax_im)
        edge_box = [0, nx, 0, ny]

        # plot it
        chosen_cmap = self.color_schemes[self.chosen_color_scheme]['cmap']
        target_plot_widget.imshow(blue_spax_im,
                                  cmap = matplotlib.cm.get_cmap(chosen_cmap),
                                  interpolation = 'nearest',
                                  extent=edge_box,
                                  vmin=0.0,
                                  vmax=set_vmax,
                                  origin='lower')
        
        # then put boxes for galaxy spaxels
        # again i am gonna assume the user input real stuff
        blue_gal_spax = self.blue_gal_spax
        if blue_gal_spax != None:    
            for gspax in blue_gal_spax:          
                # goes from spax number to [i,j] pair
                xstr, ystr = gspax[1:-1].split(',')
                gi = int(float(ystr))
                gj = int(float(xstr))
                blue_gal_box_x = [gj, gj,   gj+1, gj+1, gj]
                blue_gal_box_y = [gi, gi+1, gi+1, gi,   gi]
                gal_color = self.color_schemes[self.chosen_color_scheme]['gal_color']
                target_plot_widget.plot(blue_gal_box_x, blue_gal_box_y,
                                        color = gal_color, linewidth = 2)

        # ... and sky spax
        blue_sky_spax = self.blue_sky_spax
        if blue_sky_spax != None:          
            for sspax in blue_sky_spax:
                # goes from spax number to [i,j] pair
                xstr, ystr = sspax[1:-1].split(',')
                si = int(float(ystr))
                sj = int(float(xstr))
                blue_sky_box_x = [sj, sj,   sj+1, sj+1, sj]
                blue_sky_box_y = [si, si+1, si+1, si,   si]
                sky_color = self.color_schemes[self.chosen_color_scheme]['sky_color']
                target_plot_widget.plot(blue_sky_box_x, blue_sky_box_y,
                                        color = sky_color, linewidth = 2)

        # set the axis
        target_plot_widget.axis(edge_box)

    #----------------------------------------------------
    # Cube Image Mouse click events
    #----------------------------------------------------
    def click_red_cube_image(self, event):
        # if there is no data don't do anything
        if self.red_cube_data == None:
            return
        # check that it was clicked in the viewing area
        try:
            click_x = int(event.xdata)
            click_y = int(event.ydata)
        except:
            return
        # now figure out the spaxel number
        new_spax_str = '[%d,%d]' % (click_x, click_y)

        # figure out if it is galaxy or sky, and add to the appropriate list
        # or remove if already on that list
        if event.button == 1:
            # galaxy
            curr_gal_entry = self.red_gal_spax_entry.get_text()
            if len(curr_gal_entry) == 0:
                self.red_gal_spax_entry.set_text(new_spax_str)
            else:
                # check if it is already in the galaxy list:
                try:
                    spax_pos = self.red_gal_spax.index(new_spax_str)
                    # if it is here then it needs to be deleted
                    spax_string = ''
                    for spax in self.red_gal_spax:
                        if spax == new_spax_str:
                            continue
                        else:
                            spax_string += spax + ';'
                    self.red_gal_spax_entry.set_text(spax_string[:-1])
                # if not then add it to what is there
                except:
                    self.red_gal_spax_entry.set_text(curr_gal_entry+';'+new_spax_str)

            # finally check if it is in the sky list:
            try:
                spax_pos = self.red_sky_spax.index(new_spax_str)
                # if it is here then it needs to be deleted
                spax_string = ''
                for spax in self.red_sky_spax:
                    if spax == new_spax_str:
                        continue
                    else:
                        spax_string += spax + ';'
                self.red_sky_spax_entry.set_text(spax_string[:-1])
                self.update_red_sky_spax(self.red_sky_spax_entry)
            # if not then do nothing
            except:
                pass

            # this dumb this has to be sent a widget, grrr
            self.update_red_gal_spax(self.red_gal_spax_entry)

        elif event.button == 3:
            # sky
            curr_sky_entry = self.red_sky_spax_entry.get_text()
            if len(curr_sky_entry) == 0:
                self.red_sky_spax_entry.set_text(new_spax_str)
            else:
                # check if it is already in the sky list
                try:
                    spax_pos = self.red_sky_spax.index(new_spax_str)
                    # if it is here then it needs to be deleted
                    spax_string = ''
                    for spax in self.red_sky_spax:
                        if spax == new_spax_str:
                            continue
                        else:
                            spax_string += spax + ';'
                    self.red_sky_spax_entry.set_text(spax_string[:-1])
                # if not then add to what is already there
                except:
                    self.red_sky_spax_entry.set_text(curr_sky_entry+';'+new_spax_str)


            # finally need to if it is in the galaxy list:
            try:
                spax_pos = self.red_gal_spax.index(new_spax_str)
                # if it is here then it needs to be deleted
                spax_string = ''
                for spax in self.red_gal_spax:
                    if spax == new_spax_str:
                        continue
                    else:
                        spax_string += spax + ';'
                self.red_gal_spax_entry.set_text(spax_string[:-1])
                self.update_red_gal_spax(self.red_gal_spax_entry)
            # if not then do nothing
            except:
                pass

            # this dumb this has to be sent a widget, grrr
            self.update_red_sky_spax(self.red_sky_spax_entry)

        # otherwise it must be the other mouse button
        else:
            error_msg = 'Mouse Button 2 not bound to event'
            self.error_announce(error_msg)

    def click_blue_cube_image(self, event):
        # if there is no data don't do anything
        if self.blue_cube_data == None:
            return
        # check that it was clicked in the viewing area
        try:
            click_x = int(event.xdata)
            click_y = int(event.ydata)
        except:
            return
        # now figure out the spaxel number
        new_spax_str = '[%d,%d]' % (click_x, click_y)

        # figure out if it is galaxy or sky, and add to the appropriate list
        # or remove if already on that list
        if event.button == 1:
            # galaxy
            curr_gal_entry = self.blue_gal_spax_entry.get_text()
            if len(curr_gal_entry) == 0:
                self.blue_gal_spax_entry.set_text(new_spax_str)
            else:
                # check if it is already in the galaxy list:
                try:
                    spax_pos = self.blue_gal_spax.index(new_spax_str)
                    # if it is here then it needs to be deleted
                    spax_string = ''
                    for spax in self.blue_gal_spax:
                        if spax == new_spax_str:
                            continue
                        else:
                            spax_string += spax + ';'
                    self.blue_gal_spax_entry.set_text(spax_string[:-1])
                # if not then add it to what is there
                except:
                    self.blue_gal_spax_entry.set_text(curr_gal_entry+';'+new_spax_str)

            # finally check if it is in the sky list:
            try:
                spax_pos = self.blue_sky_spax.index(new_spax_str)
                # if it is here then it needs to be deleted
                spax_string = ''
                for spax in self.blue_sky_spax:
                    if spax == new_spax_str:
                        continue
                    else:
                        spax_string += spax + ';'
                self.blue_sky_spax_entry.set_text(spax_string[:-1])
                self.update_blue_sky_spax(self.blue_sky_spax_entry)
            # if not then do nothing
            except:
                pass

            # this dumb this has to be sent a widget, grrr
            self.update_blue_gal_spax(self.blue_gal_spax_entry)

        elif event.button == 3:
            # sky
            curr_sky_entry = self.blue_sky_spax_entry.get_text()
            if len(curr_sky_entry) == 0:
                self.blue_sky_spax_entry.set_text(new_spax_str)
            else:
                # check if it is already in the sky list
                try:
                    spax_pos = self.blue_sky_spax.index(new_spax_str)
                    # if it is here then it needs to be deleted
                    spax_string = ''
                    for spax in self.blue_sky_spax:
                        if spax == new_spax_str:
                            continue
                        else:
                            spax_string += spax + ';'
                    self.blue_sky_spax_entry.set_text(spax_string[:-1])
                # if not then add to what is already there
                except:
                    self.blue_sky_spax_entry.set_text(curr_sky_entry+';'+new_spax_str)


            # finally need to if it is in the galaxy list:
            try:
                spax_pos = self.blue_gal_spax.index(new_spax_str)
                # if it is here then it needs to be deleted
                spax_string = ''
                for spax in self.blue_gal_spax:
                    if spax == new_spax_num:
                        continue
                    else:
                        spax_string += spax + ';'
                self.blue_gal_spax_entry.set_text(spax_string[:-1])
                self.update_blue_gal_spax(self.blue_gal_spax_entry)
            # if not then do nothing
            except:
                pass

            # this dumb this has to be sent a widget, grrr
            self.update_blue_sky_spax(self.blue_sky_spax_entry)

        # otherwise it must be the other mouse button
        else:
            error_msg = 'Mouse Button 2 not bound to event'
            self.error_announce(error_msg)

    #------------------------------------
    # Red Spectrum data access functions
    #------------------------------------
    def retrieve_red_gal_flux(self):
        red_gal_spax = self.red_gal_spax
        # return None if no gal spax specified
        if red_gal_spax == None:
            return None
        # get red galaxy data
        temp_red_gal_wt_flux = numpy.zeros(len(self.red_wavelength_array), dtype = 'f')
        temp_red_gal_weights = numpy.zeros(len(self.red_wavelength_array), dtype = 'f')
        for gspax in red_gal_spax:
            # goes from spax number to [i,j] pair
            xstr, ystr = gspax[1:-1].split(',')
            gi = int(float(ystr))
            gj = int(float(xstr))
            curr_flux = self.red_cube_data[:,gi,gj]
            curr_var  = self.red_var_cube_data[:,gi,gj]
            temp_red_gal_wt_flux += curr_flux / curr_var
            temp_red_gal_weights += 1.0 / curr_var
        red_gal_flux = temp_red_gal_wt_flux / temp_red_gal_weights
        return red_gal_flux
    
    def retrieve_red_gal_var(self):
        red_gal_spax = self.red_gal_spax
        # return None if no gal spax specified
        if red_gal_spax == None:
            return None
        # get red galaxy data
        temp_red_gal_weights = numpy.zeros(len(self.red_wavelength_array), dtype = 'f')
        for gspax in red_gal_spax:
            xstr, ystr = gspax[1:-1].split(',')
            gi = int(float(ystr))
            gj = int(float(xstr))
            curr_var  = self.red_var_cube_data[:,gi,gj]
            temp_red_gal_weights += 1.0 / curr_var
        red_gal_var = 1.0 / temp_red_gal_weights
        return red_gal_var

    def retrieve_red_sky_flux(self):
        red_sky_spax = self.red_sky_spax
        # if no sky spaxels specified, return None
        if red_sky_spax == None:
            return None
        # get red sky data
        temp_red_sky_wt_flux = numpy.zeros(len(self.red_wavelength_array), dtype = 'f')
        temp_red_sky_weights = numpy.zeros(len(self.red_wavelength_array), dtype = 'f')
        for sspax in red_sky_spax:
            xstr, ystr = sspax[1:-1].split(',')
            si = int(float(ystr))
            sj = int(float(xstr))
            curr_flux = self.red_cube_data[:,si,sj]
            curr_var  = self.red_var_cube_data[:,si,sj]
            temp_red_sky_wt_flux += curr_flux / curr_var
            temp_red_sky_weights += 1.0 / curr_var
        red_sky_flux = temp_red_sky_wt_flux / temp_red_sky_weights
        return red_sky_flux

    def retrieve_red_sky_var(self):
        red_sky_spax = self.red_sky_spax
        # if no sky spaxels specified, return None
        if red_sky_spax == None:
            return None
        # get red sky data
        temp_red_sky_weights = numpy.zeros(len(self.red_wavelength_array), dtype = 'f')
        for sspax in red_sky_spax:
            xstr, ystr = sspax[1:-1].split(',')
            si = int(float(ystr))
            sj = int(float(xstr))
            curr_var  = self.red_var_cube_data[:,si,sj]
            temp_red_sky_weights += 1.0 / curr_var
        red_sky_var = 1.0 / temp_red_sky_weights
        return red_sky_var

    def retrieve_red_sub_flux(self):
        red_gal_flux = self.retrieve_red_gal_flux()
        red_sky_flux = self.retrieve_red_sky_flux()
        if red_gal_flux == None or red_sky_flux == None:
            return None
        else:
            return red_gal_flux - red_sky_flux

    def retrieve_red_sub_var(self):
        red_gal_var = self.retrieve_red_gal_var()
        red_sky_var = self.retrieve_red_sky_var()
        if red_gal_var == None or red_sky_var == None:
            return None
        else:
            return 1.0/(1.0/red_gal_var + 1.0/red_sky_var)
        
    #------------------------------------
    # Blue Spectrum data access functions
    #------------------------------------
    def retrieve_blue_gal_flux(self):
        blue_gal_spax = self.blue_gal_spax
        # return None if no gal spax specified
        if blue_gal_spax == None:
            return None
        # get blue galaxy data
        temp_blue_gal_wt_flux = numpy.zeros(len(self.blue_wavelength_array), dtype = 'f')
        temp_blue_gal_weights = numpy.zeros(len(self.blue_wavelength_array), dtype = 'f')
        for gspax in blue_gal_spax:
            xstr, ystr = gspax[1:-1].split(',')
            gi = int(float(ystr))
            gj = int(float(xstr))
            curr_flux = self.blue_cube_data[:,gi,gj]
            curr_var  = self.blue_var_cube_data[:,gi,gj]
            temp_blue_gal_wt_flux += curr_flux / curr_var
            temp_blue_gal_weights += 1.0 / curr_var
        blue_gal_flux = temp_blue_gal_wt_flux / temp_blue_gal_weights
        return blue_gal_flux
    
    def retrieve_blue_gal_var(self):
        blue_gal_spax = self.blue_gal_spax
        # return None if no gal spax specified
        if blue_gal_spax == None:
            return None
        # get blue galaxy data
        temp_blue_gal_weights = numpy.zeros(len(self.blue_wavelength_array), dtype = 'f')
        for gspax in blue_gal_spax:
            xstr, ystr = gspax[1:-1].split(',')
            gi = int(float(ystr))
            gj = int(float(xstr))
            curr_var  = self.blue_var_cube_data[:,gi,gj]
            temp_blue_gal_weights += 1.0 / curr_var
        blue_gal_var = 1.0 / temp_blue_gal_weights
        return blue_gal_var

    def retrieve_blue_sky_flux(self):
        blue_sky_spax = self.blue_sky_spax
        # if no sky spaxels specified, return None
        if blue_sky_spax == None:
            return None
        # get blue sky data
        temp_blue_sky_wt_flux = numpy.zeros(len(self.blue_wavelength_array), dtype = 'f')
        temp_blue_sky_weights = numpy.zeros(len(self.blue_wavelength_array), dtype = 'f')
        for sspax in blue_sky_spax:
            xstr, ystr = sspax[1:-1].split(',')
            si = int(float(ystr))
            sj = int(float(xstr))
            curr_flux = self.blue_cube_data[:,si,sj]
            curr_var  = self.blue_var_cube_data[:,si,sj]
            temp_blue_sky_wt_flux += curr_flux / curr_var
            temp_blue_sky_weights += 1.0 / curr_var
        blue_sky_flux = temp_blue_sky_wt_flux / temp_blue_sky_weights
        return blue_sky_flux

    def retrieve_blue_sky_var(self):
        blue_sky_spax = self.blue_sky_spax
        # if no sky spaxels specified, return None
        if blue_sky_spax == None:
            return None
        # get blue sky data
        temp_blue_sky_weights = numpy.zeros(len(self.blue_wavelength_array), dtype = 'f')
        for sspax in blue_sky_spax:
            xstr, ystr = sspax[1:-1].split(',')
            si = int(float(ystr))
            sj = int(float(xstr))
            curr_var  = self.blue_var_cube_data[:,si,sj]
            temp_blue_sky_weights += 1.0 / curr_var
        blue_sky_var = 1.0 / temp_blue_sky_weights
        return blue_sky_var

    def retrieve_blue_sub_flux(self):
        blue_gal_flux = self.retrieve_blue_gal_flux()
        blue_sky_flux = self.retrieve_blue_sky_flux()
        if blue_gal_flux == None or blue_sky_flux == None:
            return None
        else:
            return blue_gal_flux - blue_sky_flux

    def retrieve_blue_sub_var(self):
        blue_gal_var = self.retrieve_blue_gal_var()
        blue_sky_var = self.retrieve_blue_sky_var()
        if blue_gal_var == None or blue_sky_var == None:
            return None
        else:
            return 1.0/(1.0/blue_gal_var + 1.0/blue_sky_var)

    #------------------------------------
    # Red Spectrum plotting functions
    #------------------------------------
    def expand_red_specs(self, widget):
        expanded_rspec = gtk.Window(gtk.WINDOW_TOPLEVEL)
        expanded_rspec.set_title('Red Cube Spectra')

        # what to do when 'delete event' signal comes
        # defined in self function below
        expanded_rspec.connect('delete_event', self.delete_expanded_plot)
        # what to do with a destroy event
        expanded_rspec.connect('destroy', self.delete_expanded_plot)

        # let's make a vbox to put it all in
        exp_rspec_vbox = gtk.VBox(False, 0)
        expanded_rspec.add(exp_rspec_vbox)
        exp_rspec_vbox.show()
        
        # first make the figure
        expanded_red_specs = Figure(figsize = (3,6), dpi = 50)
        exp_red_spec_im = FigureCanvas(expanded_red_specs)
        exp_red_spec_im.set_size_request(600,900)
        exp_rspec_vbox.add(exp_red_spec_im)
        exp_red_spec_im.show()

        ## now what to put in the figure
        # first find out what spaxels are set
        red_gal_flux = self.retrieve_red_gal_flux()
        red_sky_flux = self.retrieve_red_sky_flux()
        red_sub_flux = self.retrieve_red_sub_flux()
        
        # if galaxy spaxel set, plot it
        gal_subplot = expanded_red_specs.add_subplot(3,1,1)
        if red_gal_flux != None:
            gal_spec_color = self.color_schemes[self.chosen_color_scheme]['gal_spec_color']
            gal_subplot.plot(self.red_wavelength_array, red_gal_flux,
                             color = gal_spec_color)
            gal_subplot.set_title('Galaxy Spaxel')
            gal_subplot.set_xlabel('Wavelength (AA)')
            gal_subplot.set_ylabel('Flux (ergs)')

        # if sky spaxel set, plot it
        sky_subplot = expanded_red_specs.add_subplot(3,1,2)
        if red_sky_flux != None:
            sky_spec_color = self.color_schemes[self.chosen_color_scheme]['sky_spec_color']
            sky_subplot.plot(self.red_wavelength_array, red_sky_flux,
                             color = sky_spec_color)
            sky_subplot.set_title('Sky Spaxel')
            sky_subplot.set_xlabel('Wavelength (AA)')
            sky_subplot.set_ylabel('Flux (ergs)')

        # if both spaxels_set, plot sub
        sub_subplot = expanded_red_specs.add_subplot(3,1,3)
        if red_sub_flux != None:
            sub_spec_color = self.color_schemes[self.chosen_color_scheme]['sub_spec_color']
            sub_subplot.plot(self.red_wavelength_array, red_sub_flux,
                             color = sub_spec_color)
            sub_subplot.set_title('Galaxy Spaxel (sky subtracted)')
            sub_subplot.set_xlabel('Wavelength (AA)')
            sub_subplot.set_ylabel('Flux (ergs)')
        
        # draw the figure
        exp_red_spec_im.draw()
        
        # add a bar for tools
        toolbar = NavigationToolbar(exp_red_spec_im, expanded_rspec)
        exp_rspec_vbox.pack_start(toolbar, False, False, 0)
        toolbar.show

        # finally show the window
        expanded_rspec.show()
        
    def replot_red_specs(self):
        # again assume spaxels have been resaved correctly
        self.replot_red_gal_spec()
        self.replot_red_sky_spec()
        self.replot_red_sub_spec()

    def replot_red_gal_spec(self):
        red_gal_flux = self.retrieve_red_gal_flux()
        # skip failure modes
        if red_gal_flux == None:
            self.red_gal_spec_plot.clear()
            self.red_spec_image.draw()
            return
        else:
            gal_spec_color = self.color_schemes[self.chosen_color_scheme]['gal_spec_color']
            self.red_gal_spec_plot.plot(self.red_wavelength_array, red_gal_flux,
                                        color = gal_spec_color)
            self.red_spec_image.draw()

    def replot_red_sky_spec(self):
        red_sky_flux = self.retrieve_red_sky_flux()
        # skip failure modes
        if red_sky_flux == None:
            self.red_sky_spec_plot.clear()
            self.red_spec_image.draw()
            return
        else:
            sky_spec_color = self.color_schemes[self.chosen_color_scheme]['sky_spec_color']
            self.red_sky_spec_plot.plot(self.red_wavelength_array, red_sky_flux,
                                        color = sky_spec_color)
            self.red_spec_image.draw()

    def replot_red_sub_spec(self):
        red_sub_flux = self.retrieve_red_sub_flux()
        # skip failure modes
        if red_sub_flux == None:
            self.red_sub_spec_plot.clear()
            self.red_spec_image.draw()
            return
        else:
            sub_spec_color = self.color_schemes[self.chosen_color_scheme]['sub_spec_color']
            self.red_sub_spec_plot.plot(self.red_wavelength_array, red_sub_flux,
                                        color = sub_spec_color)
            self.red_spec_image.draw()

    #------------------------------------
    # Blue Spectrum plotting functions
    #------------------------------------
    def expand_blue_specs(self, widget):
        expanded_bspec = gtk.Window(gtk.WINDOW_TOPLEVEL)
        expanded_bspec.set_title('Blue Cube Spectra')

        # what to do when 'delete event' signal comes
        # defined in self function below
        expanded_bspec.connect('delete_event', self.delete_expanded_plot)
        # what to do with a destroy event
        expanded_bspec.connect('destroy', self.delete_expanded_plot)

        # let's make a vbox to put it all in
        exp_bspec_vbox = gtk.VBox(False, 0)
        expanded_bspec.add(exp_bspec_vbox)
        exp_bspec_vbox.show()
        
        # first make a the figure
        expanded_blue_specs = Figure(figsize = (3,6), dpi = 50)
        exp_blue_spec_im = FigureCanvas(expanded_blue_specs)
        exp_blue_spec_im.set_size_request(600,900)
        exp_bspec_vbox.add(exp_blue_spec_im)
        exp_blue_spec_im.show()

        ## now what to put in the figure
        # first find out what spaxels are set
        blue_gal_flux = self.retrieve_blue_gal_flux()
        blue_sky_flux = self.retrieve_blue_sky_flux()
        blue_sub_flux = self.retrieve_blue_sub_flux()
        
        # if galaxy spaxel set, plot it
        gal_subplot = expanded_blue_specs.add_subplot(3,1,1)
        if blue_gal_flux != None:
            gal_spec_color = self.color_schemes[self.chosen_color_scheme]['gal_spec_color']
            gal_subplot.plot(self.blue_wavelength_array, blue_gal_flux,
                             color = gal_spec_color)
            gal_subplot.set_title('Galaxy Spaxel')
            gal_subplot.set_xlabel('Wavelength (AA)')
            gal_subplot.set_ylabel('Flux (ergs)')

        # if sky spaxel set, plot it
        sky_subplot = expanded_blue_specs.add_subplot(3,1,2)
        if blue_sky_flux != None:
            sky_spec_color = self.color_schemes[self.chosen_color_scheme]['sky_spec_color']
            sky_subplot.plot(self.blue_wavelength_array, blue_sky_flux,
                             color = sky_spec_color)
            sky_subplot.set_title('Sky Spaxel')
            sky_subplot.set_xlabel('Wavelength (AA)')
            sky_subplot.set_ylabel('Flux (ergs)')

        # if both spaxels_set, plot sub
        sub_subplot = expanded_blue_specs.add_subplot(3,1,3)
        if blue_sub_flux != None:
            sub_spec_color = self.color_schemes[self.chosen_color_scheme]['sub_spec_color']
            sub_subplot.plot(self.blue_wavelength_array, blue_sub_flux,
                             color = sub_spec_color)
            sub_subplot.set_title('Galaxy Spaxel (sky subtracted)')
            sub_subplot.set_xlabel('Wavelength (AA)')
            sub_subplot.set_ylabel('Flux (ergs)')
        
        # draw the figure
        exp_blue_spec_im.draw()
        
        # add a bar for tools
        toolbar = NavigationToolbar(exp_blue_spec_im, expanded_bspec)
        exp_bspec_vbox.pack_start(toolbar, False, False, 0)
        toolbar.show

        # finally show the window
        expanded_bspec.show()
        
    def replot_blue_specs(self):
        # again assume spaxels have been resaved correctly
        self.replot_blue_gal_spec()
        self.replot_blue_sky_spec()
        self.replot_blue_sub_spec()

    def replot_blue_gal_spec(self):
        blue_gal_flux = self.retrieve_blue_gal_flux()
        # skip failure modes
        if blue_gal_flux == None:
            self.blue_gal_spec_plot.clear()
            self.blue_spec_image.draw()
            return
        else:
            gal_spec_color = self.color_schemes[self.chosen_color_scheme]['gal_spec_color']
            self.blue_gal_spec_plot.plot(self.blue_wavelength_array, blue_gal_flux,
                                         color = gal_spec_color)
            self.blue_spec_image.draw()

    def replot_blue_sky_spec(self):
        blue_sky_flux = self.retrieve_blue_sky_flux()
        # skip failure modes
        if blue_sky_flux == None:
            self.blue_sky_spec_plot.clear()
            self.blue_spec_image.draw()
            return
        else:
            sky_spec_color = self.color_schemes[self.chosen_color_scheme]['sky_spec_color']
            self.blue_sky_spec_plot.plot(self.blue_wavelength_array, blue_sky_flux,
                                         color = sky_spec_color)
            self.blue_spec_image.draw()

    def replot_blue_sub_spec(self):
        blue_sub_flux = self.retrieve_blue_sub_flux()
        # skip failure modes
        if blue_sub_flux == None:
            self.blue_sub_spec_plot.clear()
            self.blue_spec_image.draw()
            return
        else:
            sub_spec_color = self.color_schemes[self.chosen_color_scheme]['sub_spec_color']
            self.blue_sub_spec_plot.plot(self.blue_wavelength_array, blue_sub_flux,
                                         color = sub_spec_color)
            self.blue_spec_image.draw()

    #----------------------------------------------
    # Red Saving Functions
    #----------------------------------------------
    def popup_red_save_menu(self, red_save_button_widget):
        self.red_save_menu.popup(None, None, None, 0, 0)
        
    def set_red_save_selector(self, widget, choice_string):
        self.red_save_choice = choice_string
        self.red_save_selector.set_label(self.red_save_choice)

    def retrieve_red_save_flux(self):
        if self.red_save_choice == 'Galaxy Spaxel':
            return self.retrieve_red_gal_flux()
        elif self.red_save_choice == 'Sky Spaxel':
            return self.retrieve_red_sky_flux()
        elif self.red_save_choice == 'Subbed Spaxels':
            return self.retrieve_red_sub_flux()
        else:
            error_msg = 'Somehow the spaxel select is broken.'
            self.error_announce(error_msg)

    def retrieve_red_save_var(self):
        if self.red_save_choice == 'Galaxy Spaxel':
            return self.retrieve_red_gal_var()
        elif self.red_save_choice == 'Sky Spaxel':
            return self.retrieve_red_sky_var()
        elif self.red_save_choice == 'Subbed Spaxels':
            return self.retrieve_red_sub_var()
        else:
            error_msg = 'Somehow the spaxel select is broken.'
            self.error_announce(error_msg)
        
    def save_red_spec(self, widget):
        save_flux = self.retrieve_red_save_flux()
        save_var  = self.retrieve_red_save_var()
        
        # failure modes with error messages should go here
        if save_flux == None:
            # appropriate error msg here
            error_msg = 'No flux to save!'
            self.error_announce(error_msg)
            return

        # save wavelength, flux, (and var if requested)
        red_wl = self.red_wavelength_array

        # save variance in the third column
        good_inds = numpy.nonzero(save_flux == save_flux)[0]
        save_lines = [str(red_wl[i]) + '\t' + str(save_flux[i])
                      + '\t' + str(save_var[i])
                      + '\n' for i in good_inds]

        # finally write out the file
        init_save_fn = self.red_spec_save_entry.get_text()
        # if save file not specified just make a temp one
        if len(init_save_fn) == 0:
            init_save_fn = 'red_gal_temp.dat'
        # if you are not me then don't write to my directory
        save_dir = os.getcwd() + '/'
        save_fn = save_dir + init_save_fn
        # now write it out
        f = open(save_fn, 'w')
        f.writelines(save_lines)
        f.close()
        # NEW! - same yaml file with gg info
        #self.save_red_yaml(save_fn)

    #def save_red_yaml(self, save_fn):
    #    # ok pyyaml doesn't seem to be working, so i
    #    # have to just save this text-wise
    #    out_yaml = ['---\n',
    #                'name: \'%s\'\n'      % self.object_name,
    #                'snifs_run: \'%s\'\n' % self.snifs_run,
    #                'channel: \'R\'\n',
    #                'cube_fn: \'%s\'\n'   % self.red_cube_fn,
    #                'gal_spax: %s\n'  % str(self.red_gal_spax),
    #                'sky_spax: %s\n'  % str(self.red_sky_spax),
    #                'save_fn: \'%s\''     % save_fn]
    #    # save it to appropriate file
    #    save_extn = save_fn.split('.')[-1]
    #    yaml_fn = save_fn[:-len(save_extn)] + 'yaml'
    #    f = open(yaml_fn, 'w')
    #    f.writelines(out_yaml)
    #    f.close()

    #----------------------------------------------
    # Blue Saving Functions
    #----------------------------------------------
    def popup_blue_save_menu(self, blue_save_button_widget):
        self.blue_save_menu.popup(None, None, None, 0, 0)
        
    def set_blue_save_selector(self, widget, choice_string):
        self.blue_save_choice = choice_string
        self.blue_save_selector.set_label(self.blue_save_choice)

    def retrieve_blue_save_flux(self):
        if self.blue_save_choice == 'Galaxy Spaxel':
            return self.retrieve_blue_gal_flux()
        elif self.blue_save_choice == 'Sky Spaxel':
            return self.retrieve_blue_sky_flux()
        elif self.blue_save_choice == 'Subbed Spaxels':
            return self.retrieve_blue_sub_flux()
        else:
            error_msg = 'Somehow the spaxel select is broken.'
            self.error_announce(error_msg)

    def retrieve_blue_save_var(self):
        if self.blue_save_choice == 'Galaxy Spaxel':
            return self.retrieve_blue_gal_var()
        elif self.blue_save_choice == 'Sky Spaxel':
            return self.retrieve_blue_sky_var()
        elif self.blue_save_choice == 'Subbed Spaxels':
            return self.retrieve_blue_sub_var()
        else:
            error_msg = 'Somehow the spaxel select is broken.'
            self.error_announce(error_msg)
        
    def save_blue_spec(self, widget):
        save_flux = self.retrieve_blue_save_flux()
        save_var  = self.retrieve_blue_save_var()
        
        # failure modes with error messages should go here
        if save_flux == None:
            # appropriate error msg here
            error_msg = 'No flux to save!'
            self.error_announce(error_msg)
            return

        # save wavelength, flux, (and var if requested)
        blue_wl = self.blue_wavelength_array

        # save variance in the third column
        save_lines = [str(blue_wl[i]) + '\t' + str(save_flux[i])
                      + '\t' + str(save_var[i])
                      + '\n' for i in range(len(save_flux))]
        
        # finally write out the file
        save_fn = self.blue_spec_save_entry.get_text()
        # if save file not specified just make a temp one
        if len(save_fn) == 0:
            save_fn = 'blue_gal_temp.dat'
        # if you are not me then don't write to my directory
        save_dir = os.getcwd() + '/'
        save_fn = save_dir + save_fn
        # now write it out
        f = open(save_fn, 'w')
        f.writelines(save_lines)
        f.close()
        # NEW! - same yaml file with gg info
        #self.save_blue_yaml(save_fn)

    #def save_blue_yaml(self, save_fn):
    #    # ok pyyaml doesn't seem to be working, so i
    #    # have to just save this text-wise
    #    out_yaml = ['---\n',
    #                'name: \'%s\'\n'      % self.object_name,
    #                'snifs_run: \'%s\'\n' % self.snifs_run,
    #                'channel: \'B\'\n',
    #                'cube_fn: \'%s\'\n'   % self.blue_cube_fn,
    #                'gal_spax: %s\n'  % str(self.blue_gal_spax),
    #                'sky_spax: %s\n'  % str(self.blue_sky_spax),
    #                'save_fn: \'%s\''     % save_fn]
    #    # save it to appropriate file
    #    save_extn = save_fn.split('.')[-1]
    #    yaml_fn = save_fn[:-len(save_extn)] + 'yaml'
    #    f = open(yaml_fn, 'w')
    #    f.writelines(out_yaml)
    #    f.close()


    #----------------------------------------------
    # Correlator Preparation Function
    #----------------------------------------------
    def update_blue_prep(self, widget=None, data=None):
        if data == None:
            blue_prep_fn = self.blue_prep_entry.get_text()
            data = numpy.loadtxt(blue_prep_fn)
        self.corr_blue_wave = data[:,0]
        self.corr_blue_flux = data[:,1]
        self.update_prepped_spec()
        self.replot_prep_spec()
        return

    def export_blue_to_prep(self, widget=None):
        # grab the data
        wave = self.blue_wavelength_array
        flux = self.retrieve_blue_sub_flux()
        data = numpy.zeros([len(flux), 2], dtype='d')
        data[:,0] = wave
        data[:,1] = flux
        # update it!
        self.update_blue_prep(data=data)
        return

    def select_blue_prep(self, widget=None):
        # filename selection
        dialog = gtk.FileChooserDialog(
            title='Select Blue Data For Correlation',
            action=gtk.FILE_CHOOSER_ACTION_OPEN,
            buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                     gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_CANCEL)
        response = dialog.run()
        # update the filename if a filename was chosen
        if response == gtk.RESPONSE_OK:
            new_blue_fn = dialog.get_filename()
            self.blue_prep_entry.set_text(new_blue_fn)
            self.update_blue_prep()
        dialog.destroy()
        return

    def update_red_prep(self, widget=None, data=None):
        if data == None:
            red_prep_fn = self.red_prep_entry.get_text()
            data = numpy.loadtxt(red_prep_fn)
        self.corr_red_wave = data[:,0]
        self.corr_red_flux = data[:,1]
        self.update_prepped_spec()
        self.replot_prep_spec()
        return

    def export_red_to_prep(self, widget=None):
        # grab the data
        wave = self.red_wavelength_array
        flux = self.retrieve_red_sub_flux()
        data = numpy.zeros([len(flux), 2], dtype='d')
        data[:,0] = wave
        data[:,1] = flux
        # update it!
        self.update_red_prep(data=data)
        return

    def select_red_prep(self, widget=None):
        # filename selection
        dialog = gtk.FileChooserDialog(
            title='Select Red Data For Correlation',
            action=gtk.FILE_CHOOSER_ACTION_OPEN,
            buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                     gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_CANCEL)
        response = dialog.run()
        # update the filename if a filename was chosen
        if response == gtk.RESPONSE_OK:
            new_red_fn = dialog.get_filename()
            self.red_prep_entry.set_text(new_red_fn)
            self.update_red_prep()
        dialog.destroy()
        return

    def adjust_prep_scaling(self, widget=None):
        self.update_prepped_spec()
        self.replot_prep_spec()
        return

    def update_prepped_spec(self):
        # all cases:
        # 1 - neither blue or red set
        if self.corr_blue_flux == None and self.corr_red_flux == None:
            self.corr_flux = None
            self.corr_wave = None
        # 2 - only 1 is set
        elif self.corr_blue_flux == None:
            self.corr_flux = self.corr_red_flux
            self.corr_wave = self.corr_red_wave
        elif self.corr_red_flux == None:
            self.corr_flux = self.corr_blue_flux
            self.corr_wave = self.corr_blue_wave
        # 3 - both are set, with a scaling
        else:
            # get blue/red scaling
            br_scale = float(self.prep_scale_entry.get_text())
            wave = numpy.concatenate((self.corr_blue_wave,
                                      self.corr_red_wave))
            flux = numpy.concatenate((br_scale*self.corr_blue_flux,
                                      self.corr_red_flux))
            order = wave.argsort()
            self.corr_flux = flux[order]
            self.corr_wave = wave[order]
        # exit
        return 

    def plot_prep_spec(self, target_plot_widget):
        target_plot_widget.clear()
        target_plot_widget.plot(self.corr_wave, self.corr_flux)
        return

    def replot_prep_spec(self):
        self.plot_prep_spec(self.prep_spec_plot)
        self.prep_spec_image.draw()
        return

    def expand_prep_spec(self, widget):
        expanded_prep = gtk.Window(gtk.WINDOW_TOPLEVEL)
        expanded_prep.set_title('Blue Cube Summed Spaxels')

        # what to do when 'delete event' signal comes
        # defined in self function below
        expanded_prep.connect('delete_event', self.delete_expanded_plot)
        # what to do with a destroy event
        expanded_prep.connect('destroy', self.delete_expanded_plot)

        # must hold all in a vbox
        exp_vbox = gtk.VBox(False, 0)
        expanded_prep.add(exp_vbox)
        exp_vbox.show()

        # first make a the figure
        expanded_prep_fig = Figure(figsize = (6,4), dpi = 50)
        
        ex_prep_im = FigureCanvas(expanded_prep_fig)
        exp_vbox.pack_start(ex_prep_im, False, False, 0)
        ex_prep_im.set_size_request(600,400)
        ex_prep_im.show()
        
        expanded_prep_cube_plot = expanded_prep_fig.add_subplot(1,1,1)
        self.plot_prep_spec(expanded_prep_cube_plot)
        # tell it to draw
        ex_prep_im.draw()

        # add a toolbar
        toolbar = NavigationToolbar(ex_prep_im, expanded_prep)
        exp_vbox.pack_start(toolbar, False, False, 0)
        toolbar.show()

        # finally show it
        expanded_prep.show()
        return

    #----------------------------------------------
    # Function for Correlation Tab
    #----------------------------------------------
    def run_snid(self, widget):
        try:
            zmax = float(self.corr_zmax_entry.get_text())
        except:
            zmax = 0.1
            self.corr_zmax_entry.set_text(str(zmax))
        # 1 - run snid on the data
        snid_results = compare_all_snid(
            self.corr_wave,
            self.corr_flux,
            zmax)
        # 2 - calculate and display best redshift and phase
        epochs = numpy.array([float(x) for x in snid_results[2]])
        zbests = snid_results[3]
        zerrs  = snid_results[4]
        rlaps  = snid_results[5]
        # set up weights
        weights = numpy.zeros(len(snid_results[0]),dtype='d')
        weights[numpy.nonzero(rlaps >= 4.0)[0]] = 1.0
        weights[numpy.nonzero(rlaps >= 5.0)[0]] = 3.0
        weights[numpy.nonzero(rlaps >= 6.0)[0]] = 5.0
        wt_sum = numpy.sum(weights)
        # final weighted values
        zbest      = numpy.sum(weights*zbests)/wt_sum
        zbest_err  = (numpy.sum(weights/zerrs**2)/wt_sum)**-0.5
        tbest      = numpy.sum(weights*epochs)/wt_sum
        #tbest_err  = (numpy.sum(weights*((tbest-epochs)**2))/wt_sum)**0.5
        #tbest_err  = numpy.sum(weights*numpy.abs(tbest-epochs))/wt_sum
        tbest_err  = (numpy.sum(weights*((tbest-epochs)**-2))/wt_sum)**-0.5
        # make strings and display
        zbest_str = '%.06f +- %.06f' % (zbest, zbest_err)
        tbest_str = '%.02f +- %.02f' % (tbest, tbest_err)
        self.corr_zbest_entry.set_text(zbest_str)
        self.corr_tbest_entry.set_text(tbest_str)
        # 3 - populate the results tree
        self.snid_treestore.clear()
        n_sne = len(snid_results[0])
        # sort by rlap
        sort_order = snid_results[5].argsort()[::-1]
        for q in range(n_sne):
            i = sort_order[q]
            self.snid_treestore.append([x[i] for x in snid_results])
        self.window.show_all()
        return

    def clear_snid_results(self, widget=None):
        self.snid_treestore.clear()
        return

    def get_current_snid_spec(self):
        try:
            curr_tree, curr_iter = self.snid_selection.get_selected()
            curr_sn = curr_tree.get_value(curr_iter, 0)
            curr_epoch = float(curr_tree.get_value(curr_iter, 2))
            curr_z = float(curr_tree.get_value(curr_iter, 3))
            # get the data for that SN and plot it!
            init_temp_wave, temp_flux = retrieve_snid_template_data(
                curr_sn, curr_epoch)
            temp_wave = init_temp_wave*(1.0+curr_z)
        except:
            temp_wave = []
            temp_flux = []
        return temp_wave, temp_flux

    def update_snid_plot(self, widget=None):
        temp_wave, temp_flux = self.get_current_snid_spec()
        # reformat the correlation input data
        data_wave = self.corr_wave
        data_flux = clean_sn_spectrum(self.corr_wave,
                                      self.corr_flux,
                                      fft_filter = False)
        # plot the data!
        self.corr_plot.clear()
        self.corr_plot.plot(data_wave, data_flux, color='k')
        self.corr_plot.plot(temp_wave, temp_flux, color='r', lw=2)
        self.corr_plot.set_xlim([data_wave.min(), data_wave.max()])
        self.corr_image.draw()
        return

    def expand_corr_plot(self, widget):
        temp_wave, temp_flux = self.get_current_snid_spec()
        # reformat the correlation input data
        data_wave = self.corr_wave
        data_flux = clean_sn_spectrum(self.corr_wave,
                                      self.corr_flux,
                                      fft_filter = False)
        # make a new window!
        expanded_snid = gtk.Window(gtk.WINDOW_TOPLEVEL)
        expanded_snid.set_title('Red Cube Spectra')
        # what to do when 'delete event' signal comes
        # defined in self function below
        expanded_snid.connect('delete_event', self.delete_expanded_plot)
        # what to do with a destroy event
        expanded_snid.connect('destroy', self.delete_expanded_plot)
        # let's make a vbox to put it all in
        exp_snid_vbox = gtk.VBox(False, 0)
        expanded_snid.add(exp_snid_vbox)
        exp_snid_vbox.show()        
        # make the figure
        expanded_red_specs = Figure(figsize = (8,6), dpi = 50)
        exp_snid_im = FigureCanvas(expanded_red_specs)
        exp_snid_im.set_size_request(800,600)
        exp_snid_vbox.add(exp_snid_im)
        exp_snid_im.show()
        # plot the data!
        snid_subplot = expanded_red_specs.add_subplot(1,1,1)
        snid_subplot.plot(data_wave, data_flux, color='k')
        snid_subplot.plot(temp_wave, temp_flux, color='r', lw=2)
        snid_subplot.set_xlim([data_wave.min(), data_wave.max()])
        # draw the figure
        exp_snid_im.draw()        
        # add a bar for tools
        toolbar = NavigationToolbar(exp_snid_im, expanded_snid)
        exp_snid_vbox.pack_start(toolbar, False, False, 0)
        toolbar.show()
        # finally show the window
        expanded_snid.show()
        return        

    #-----------------------------------------------------
    # Error reporting
    def error_announce(self, error_string):
        self.error_message_box.set_text(error_string)

    def clear_error_box(self, widget):
        self.error_announce('')
        
    def be_annoying(self, widget):
        try:
            self.n_clicks += 1
        except:
            self.n_clicks = 1
            
        if self.n_clicks < 5:
            error_msg = 'Stop being impatient!'
        elif self.n_clicks < 10:
            error_msg = 'Seriously just chill out homes'
        elif self.n_clicks < 15:
            error_msg = 'You must really like to click buttons.'
        else:
            error_msg = 'Silly ' + os.getlogin() + ', you must have other work you should be doing!'

        self.error_announce(error_msg)

    #------------------------------------
    def main(self):
        gtk.main()

################################
# run it
if __name__ == '__main__':
    ggrab = GalaxyGrabber()
    ggrab.main()
else:
    print __name__
