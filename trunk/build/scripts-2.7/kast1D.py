#!/dark/anaconda/anaconda27/envs/halenv/bin/python
#
#  python 3 script for deimos reduction
#
###############################################################
description = "> reduce kast data, run the script in a directory with kast data "
usage = "%prog   [--iraf (reduce with iraf )\n --directory (reduce data in a different directory)\n  --interactive (reduce in interactive mode)"
import kast
from kast import kastutil
from optparse import OptionParser, OptionGroup
import pyds9
import os
import pickle
import time
import re
import sys
import numpy as np
from matplotlib import pylab as plt
from astropy.io import fits
import pyds9
import glob
                
pyversion = sys.version_info[0]

# check that ds9 is open 
plt.ion()
ds9 = pyds9.DS9(str(time.time()))
ds9.set('frame 1')
ds9.set('scale zscale');

if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description, version="%prog 1.0")
    parser.add_option("-d", "--directory", dest="directory", default=None, type="str",
                      help='reduce data in this directory \t [%default]')
    parser.add_option("--stage", dest="stage", default=None, type="str",
                      help='reduce data a single stage \t [%default, trim,trim+,sky,sky+,trace,trace+,extract,extract+,wave,wave+,flux,flux+]')   
    parser.add_option("-i", "--interactive", action="store_true",
                      dest='interactive', default=False, help='Interactive \t\t\ [%default]')
    parser.add_option("-F", "--force", dest="force", action="store_true",default=False)
    
    option, args = parser.parse_args()
    _interactive = option.interactive
    _directory = option.directory
    _force = option.force
    stage = option.stage

    if _interactive:
        _verbose= True
        _interiraf= 'yes'
    else:
        _verbose= False
        _interiraf= 'no'
    #
    #  initialize the dictionary with all the infos
    #
    dictionary, setup_object, setup_flat, setup_arc, setup_bias, setup_standard= kast.kastutil.checkalldata(directory=_directory,verbose=_verbose, all=True)
    #
    #  check that all data are in the directory and
    # 


    readaxi = {'kastr':'column', 'kastb': 'line'}
    trimsec = {'kastb':'[1:2048,20:320]', 'kastr': '[135:435,1:2700]'}
    specredaxis = {'kastb': 1, 'kastr': 2}
    

##################################################################################
            
    list = dictionary.keys()
    header= '                   '
    for key in kast.kastutil.listhd:
        if key in ['DATE-OBS']:
            header = header + '%25s' % (key)
        else:
            header = header + '%15s' % (key)

    print('#'*10 + '  STANDARDS ' +'#'*10 + '\n')
    print(header)
    for img in setup_standard['kastr']:
        value = img + ' '
        for key in kast.kastutil.listhd:
            if key in ['DATE-OBS']:
                value = value + '%25s' % (str(dictionary[img][key]))
            else:
                value = value + '%15s' % (str(dictionary[img][key]))
        print(value)
        
    print('#'*10 + '  OBJECT ' +'#'*10 + '\n')
    print(header)
    for img in setup_object['kastr']:
        value = img + ' '
        for key in kast.kastutil.listhd:
            if key in ['DATE-OBS']:
                value = value + '%25s' % (str(dictionary[img][key]))
            else:
                value = value + '%15s' % (str(dictionary[img][key]))
        print(value)

    print('#'*10 + '  ARCS ' +'#'*10 + '\n')
    print(header)
    for img in setup_arc['kastr']:
        value = img + ' '
        for key in kast.kastutil.listhd:
            if key in ['DATE-OBS']:
                value = value + '%25s' % (str(dictionary[img][key]))
            else:
                value = value + '%15s' % (str(dictionary[img][key]))
        print(value)

    objectlist={'obj':{},'std':{},'arc':{},'arcnoslit':{}}
    ###########
    for arm in setup_object.keys():
        for img in setup_object[arm]:
            _slit = dictionary[img]['SLIT_N']
            if arm == 'kastr':
                _disp = dictionary[img]['GRATNG_N']
            elif arm == 'kastb':
                _disp = dictionary[img]['GRISM_N']
            _bsplit = dictionary[img]['BSPLIT_N']
            if (arm,_disp,_bsplit,_slit) not in objectlist['obj']:
                objectlist['obj'][arm,_disp,_bsplit,_slit] = []
            objectlist['obj'][arm,_disp,_bsplit,_slit].append(img)
    ###########
    for arm in setup_standard.keys():
        for img in setup_standard[arm]:
            _slit = dictionary[img]['SLIT_N']
            if arm == 'kastr':
                _disp = dictionary[img]['GRATNG_N']
            elif arm == 'kastb':
                _disp = dictionary[img]['GRISM_N']
            _bsplit = dictionary[img]['BSPLIT_N']
            if (arm,_disp,_bsplit,_slit) not in objectlist['std']:
                objectlist['std'][arm,_disp,_bsplit,_slit] = []
            objectlist['std'][arm,_disp,_bsplit,_slit].append(img)
    ###########
    for arm in setup_arc.keys():
        for img in setup_arc[arm]:
            _slit = dictionary[img]['SLIT_N']
            if arm == 'kastr':
                _disp = dictionary[img]['GRATNG_N']
            elif arm == 'kastb':
                _disp = dictionary[img]['GRISM_N']
            _bsplit = dictionary[img]['BSPLIT_N']
            if (arm,_disp,_bsplit,_slit) not in objectlist['arc']:
                objectlist['arc'][arm,_disp,_bsplit,_slit] = []
            objectlist['arc'][arm,_disp,_bsplit,_slit].append(img)
            
            if (arm,_disp,_bsplit) not in objectlist['arcnoslit']:
                objectlist['arcnoslit'][arm,_disp,_bsplit] = []
            objectlist['arcnoslit'][arm,_disp,_bsplit].append(img)            

    ###########################################################
    #for key in objectlist['obj'].keys():
    #    skip=False
    #    print('########## Setup: '+ '-'.join(key))
    #    if key not in objectlist['std'].keys():
    #        print('WARNING not standard with this setup ')
    #        skip = True
    #    else:
    #        print('STD found')
    #    if key not in objectlist['arc'].keys():
    #        print('WARNING not ARC with this setup ')
#   #         skip = True
    #    else:
    #        print('ARC found')
    #        
    #    if skip is False:
    print('# extract spectra  and arcs ....')
    ######### OBJECT ###############
    for key in objectlist['obj'].keys():
                print(key)
                _arm = key[0]
                _disp = key[1]
                _slit = key[3]

                if key in objectlist['arc'].keys():
                    observedarc = objectlist['arc'][key]                    
                elif key[0:3] in objectlist['arcnoslit'].keys():
                    print('arc with different slit')
                    observedarc = objectlist['arcnoslit'][key[0:3]]                    
                else:
                    observedarc = []
                    
                for img in objectlist['obj'][key]:
                    imgex = os.path.splitext(img)[0] + '_ex.fits'
                    _reference = ''
                    _trace = 'yes'
                    _find = 'yes'
                    _fittrac = 'yes'
                    _recenter = 'yes'
                    _edit = 'yes'
                    _gain =1
                    _rdnoise = 1
                    dv = kast.kastutil.dvex()
                    _interactive = 'yes'
                    _review = 'yes'
                    _mode = 'q'
                    _type = 'obj'
                    kast.kastutil.extractspectrum(img,imgex,_reference,_trace,
                                                  _fittrac,_find,_recenter,_edit,
                                                  _gain,_rdnoise,dv,_interactive,
                                                  _review,_type,_mode, _arm, key, _force)
                    
                    if observedarc:
                        # need to find a better way to select the arc (close in JD)
                        arcfile = observedarc[0]
                        print(img,imgex,arcfile,_arm)
                        arcfileex =  kast.kastutil.arcextraction(arcfile, img, imgex, _arm,dv,True)
                    else:
                        print('no arc with this setup')
                        

                    
                    #raw_input('stop')
#    ######### STD ###############
    print('# extract spectra and arcs ....')
    for key in objectlist['std'].keys():
                _arm = key[0]
                _disp = key[1]
                _slit = key[3]

                if key in objectlist['arc'].keys():
                    observedarc = objectlist['arc'][key]                    
                elif key[0:3] in objectlist['arcnoslit'].keys():
                    print('arc with different slit')
                    observedarc = objectlist['arcnoslit'][key[0:3]]                    
                else:
                    observedarc = []
                
                for img in objectlist['std'][key]:
                    imgex = os.path.splitext(img)[0] + '_ex.fits'
                    _reference = ''
                    _trace = 'yes'
                    _find = 'yes'
                    _fittrac = 'yes'
                    _recenter = 'yes'
                    _edit = 'yes'
                    _gain =1
                    _rdnoise = 1
                    dv = kast.kastutil.dvex()
                    _interactive = 'yes'
                    _review = 'yes'
                    _mode = 'q'
                    _type = 'std'
                    kast.kastutil.extractspectrum(img,imgex,_reference,_trace,
                                                  _fittrac,_find,_recenter,_edit,
                                                  _gain,_rdnoise,dv,_interactive,
                                                  _review,_type,_mode, _arm, key, _force)

                    if observedarc:
                        # need to find a better way to select the arc (close in JD)
                        arcfile = observedarc[0]
                        print(img,imgex,arcfile,_arm)
                        arcfileex =  kast.kastutil.arcextraction(arcfile, img, imgex, _arm,dv,_force)
                    else:
                        print('no arc with this setup')
################################################################################
    #####   wavelengh calibraton objects
    for key in objectlist['obj'].keys():
                print(key)
                _arm = key[0]
                _disp = key[1]
                _dicroic = key[2]
                _slit = key[3]
                for img in objectlist['obj'][key]:
                    imgex = os.path.splitext(img)[0] + '_ex.fits'
                    arcfilex = 'arc_' + imgex
                    directory = kast.__path__[0] + '/archive/' + str(_arm) + '/arc/' + _disp + '/' + _dicroic 
                    print(directory)
                    listarc = glob.glob(directory + '/*fits')                    
                    print(listarc)
                    if not listarc:
                        kast.kastutil.identify(arcfilex, img, _arm, dv, arcref = False, force=_force)
                    else:
                        _arcref = listarc[0]
                        print('#######',_arcref)
                        kast.kastutil.identify(arcfilex, img, _arm, dv, arcref = _arcref, force=_force)

################################################################################
    #####   wavelengh calibraton objects
    for key in objectlist['std'].keys():
                print(key)
                _arm = key[0]
                _disp = key[1]
                _dicroic = key[2]
                _slit = key[3]
                for img in objectlist['std'][key]:
                    imgex = os.path.splitext(img)[0] + '_ex.fits'
                    arcfilex = 'arc_' + imgex
                    directory = kast.__path__[0] + '/archive/' + str(_arm) + '/arc/' + _disp + '/' + _dicroic 
                    print(directory)
                    listarc = glob.glob(directory + '/*fits')                    
                    print(listarc)
                    if not listarc:
                        kast.kastutil.identify(arcfilex, img, _arm, dv, arcref = False, force=_force)
                    else:
                        _arcref = listarc[0]
                        print('#######',_arcref)
                        kast.kastutil.identify(arcfilex, img, _arm, dv, arcref = _arcref, force=_force)
                        
#######################################################################################################
     #########   sens function
    for key in objectlist['std'].keys():
        print(key)
#        _output = '_'.join(key[:-1])
        _output = re.sub('/','','_'.join(key[:-1]))
        for img in objectlist['std'][key]:
            imgl = os.path.splitext(img)[0] + '_l.fits'
            if os.path.isfile(imgl):
                
                imgclean,atmofile = kast.kastutil.make_atmo(imgl)
                kast.kastutil.sensfunc(imgclean, _output ,_key=key,_split= True,  _function='spline3', _order=8, interactive='yes')
                
#######################################################################################################

    for key in objectlist['obj'].keys():
        print(key)
        _output = re.sub('/','','_'.join(key[:-1]))
        senslist = glob.glob('sens*'+_output+'*fits')
        print(senslist)
        if len(senslist):
            sensfile= senslist[0]
            for img in objectlist['obj'][key]:
                imgl = os.path.splitext(img)[0] + '_l.fits'
                imgf = os.path.splitext(img)[0] + '_f.fits'
                kast.kastutil.calibrate(imgl, imgf, sensfile ,  force=True, interactive='yes')
