#! /usr/bin/env python

import sys
import os
from astropy.io import fits as pyfits

from . import wifes_calib






def classifier(obs,data_dir):
    stdstar_list = wifes_calib.ref_fname_lookup.keys()

    # classify each obs
    bias = []
    domeflat = []
    twiflat = []
    dark = []
    arc = []
    wire = []
    stdstar = {}
    science = {}

    for ob in obs:
        fn = data_dir + ob + '.fits'
        f = pyfits.open(fn)
        imagetype = f[0].header['IMAGETYP'].upper()
        obj_name = f[0].header['OBJECT']
        f.close()
        #---------------------------
        # check if it is within a close distance to a standard star
        # if so, fix the object name to be the good one from the list!
        try:
            near_std, std_dist = wifes_calib.find_nearest_stdstar(fn)
            if std_dist < 100.0:
                obj_name = near_std
        except:
            pass
        #---------------------------
        # 1 - bias frames
        if imagetype == 'BIAS':
            bias.append(ob)
        # 2 - quartz flats
        if imagetype == 'FLAT':
            domeflat.append(ob)
        # 3 - twilight flats
        if imagetype == 'SKYFLAT':
            twiflat.append(ob)
        # 4 - dark frames
        if imagetype == 'DARK':
            dark.append(ob)
        # 5 - arc frames
        if imagetype == 'ARC':
            arc.append(ob)
        # 6 - wire frames
        if imagetype == 'WIRE':
            wire.append(ob)
        # 7 - standard star
        if imagetype == 'STANDARD':
            # group standard obs together!
            if obj_name in stdstar.keys():
                stdstar[obj_name].append(ob)
            else:
                stdstar[obj_name] = [ob]
        
        # all else are science targets (also consider standar star in imagety = OBJECT)
        if imagetype == 'OBJECT':
            if obj_name in stdstar_list:
                # group standard obs together!
                if obj_name in stdstar.keys():
                    stdstar[obj_name].append(ob)
                else:
                    stdstar[obj_name] = [ob]
            else:
                # group science obs together!
                if obj_name in science.keys():
                    science[obj_name].append(ob)
                else:
                    science[obj_name] = [ob]

    # #------------------
    # science dictionay

    sci_obs = []

    for obj_name in science.keys():
        obs_list = science[obj_name]
        sci_obs.append({'sci':obs_list, 'sky':[]})


    #------------------
    # stdstars dictionary
    std_obs = []

    for obj_name in stdstar.keys():
        obs_list = stdstar[obj_name]
        std_obs.append({'sci':obs_list, 'name':obj_name,'type':['flux', 'telluric']})


    obs_metadata = {
        'bias' : bias,
        'domeflat' : domeflat,
        'twiflat' : twiflat,
        'dark' : dark,
        'wire' : wire,
        'arc'  : arc,
        'sci'  : sci_obs,
        'std'  : std_obs}
    
    return obs_metadata 






def classify(data_dir):

    # What is this option for? TODO change to a default variable in the function
    try:
        naxis2_use = int(sys.argv[2])
    except:
        naxis2_use = 0

    # Get list of all fits files in directory
    all_files = os.listdir(data_dir)

    # Filtering the data as per blue and red arm
    blue_obs = []
    red_obs = []

    obs_date = None
    for fn in all_files:
        obs = fn.replace('.fits', '')

        # Date of the observations: Is really needed?
        if obs_date == None:
            try:
                f = pyfits.open(data_dir+fn)
                obs_date = f[0].header['DATE-OBS'].split('T')[0].replace('-', '')
                f.close()
            except:
                continue
        # ------------------------------------------------
        try:
            f = pyfits.open(data_dir+fn)
            camera = f[0].header['CAMERA']
            naxis2 = f[0].header['NAXIS2']
            f.close()
        except:
            continue
        if naxis2_use!=0:
            if naxis2 != naxis2_use:
                continue
        if camera == 'WiFeSBlue':
            if obs in blue_obs:
                continue
            else:
                blue_obs.append(obs)
        if camera == 'WiFeSRed':
            if obs in red_obs:
                continue
            else:
                red_obs.append(obs)


    blue_obs_metadata = classifier(blue_obs,data_dir)
    red_obs_metadata = classifier(red_obs,data_dir)


    return {"blue":blue_obs_metadata,"red": red_obs_metadata}
