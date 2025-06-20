#!/usr/bin/env python
#
#  python 3 script for deimos reduction
#
###############################################################
description = "> reduce kast data, run the script in a directory with kast data "
usage = "%prog   [--iraf (reduce with iraf )\n --directory (reduce data in a different directory)\n  --interactive (reduce in interactive mode)"
import kast
from kast import kastutil
from kast import cosmics
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
from dateutil.parser import parse
import shutil

pyversion = sys.version_info[0]

# check that ds9 is open 
#plt.ion()
#ds9 = pyds9.DS9(str(time.time()))
#ds9.set('frame 1')
#ds9.set('scale zscale');

if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description, version="%prog 1.0")
    parser.add_option("-d", "--directory", dest="directory", default=None, type="str",
                      help='reduce data in this directory \t [%default]')
    parser.add_option("--stage", dest="stage", default=None, type="str",
                      help='reduce data a single stage \t [%default, cosmic, cosmic+, extract, extract+,\
                      wave, wave+, sens, sens+, flux, flux+, combine, combine+, atmo]')   
    parser.add_option("-i", "--interactive", action="store_true",
                      dest='interactive', default=False, help='Interactive \t\t\ [%default]')
    parser.add_option("-F", "--force", dest="force", action="store_true",default=False)
    parser.add_option("--showspectra", dest="showspectra", action="store_true",default=False)
    parser.add_option("--showsens", dest="showsens", action="store_true",default=False)
    parser.add_option("--senslist", dest="senslist", default=None, type="str",
                      help='use external sensitivity list \t [%default]')

    
#    parser.add_option("--merge", dest="merge", action="store_true",default=False)
#    parser.add_option("--combine", dest="combine", action="store_true",default=False)
#    parser.add_option("--redocal", dest="redocal", action="store_true",default=False)
#    parser.add_option("--docosmic", dest="docosmic", action="store_true",default=False)

    
    option, args = parser.parse_args()
    _interactive = option.interactive
    _directory = option.directory
    _senslist = option.senslist
    _force = option.force
    _showsens = option.showsens
    _showspectra = option.showspectra
    stage = option.stage
    if stage not in [None,'all','cosmic','cosmic+', 'extract','extract+','wave','wave+',\
                     'sens','sens+','flux','flux+','combine','combine+','atmo']:
        sys.argv.append('--help')
        option, args = parser.parse_args()
        
    #  0000010  2   cosmic
    #  0000100  4   extract
    #  0001000  8   wave
    #  0001000  16  sens
    #  0010000  32  flux
    #  0100000  64  combine
    #  0100000  128  atmo

    _run=0
    if stage == 'cosmic':      _run = _run + 2
    elif stage == 'cosmic+':   _run = _run + 2 + 4 + 8 + 16 + 32 + 64 + 128
    elif stage == 'extract':   _run = _run + 4
    elif stage == 'extract+':  _run = _run + 4  + 8 + 16 + 32 + 64 + 128
    elif stage == 'wave':      _run = _run + 8
    elif stage == 'wave+':     _run = _run + 8 + 16 + 32 + 64 + 128 
    elif stage == 'sens':      _run = _run + 16
    elif stage == 'sens+':     _run = _run + 16 + 32 + 64 + 128
    elif stage == 'flux':      _run = _run + 32 
    elif stage == 'flux+':     _run = _run + 32 + 64 + 128
    elif stage == 'combine':   _run = _run + 64
    elif stage == 'combine+':  _run = _run + 64 + 128
    elif stage == 'atmo':      _run = _run + 128
    elif stage == 'all':       _run = _run + 2 + 4 + 8 + 16 + 32 + 64 + 128

    print(_run)    

    if _interactive:
        _verbose= True
        _interiraf= 'yes'
    else:
        _verbose= False
        _interiraf= 'no'

    print(_interiraf)
    #
    #  initialize the dictionary with all the infos
    #
    dictionary, setup_object, setup_flat, setup_arc, setup_bias, setup_standard= kast.kastutil.checkalldata(directory=_directory,verbose=_verbose, all=True)
    #
    #  check that all data are in the directory and
    # 

    readaxi = {'kastr':'column', 'kastb': 'line'}
    trimsec = {'kastb':'[1:1900,20:320]', 'kastr': '[135:435,60:2200]'}
    specredaxis = {'kastb': 1, 'kastr': 2}
    dv = kast.kastutil.dvex()
    _gain = 1
    _rdnoise = 1

##################################################################################
            
    list = dictionary.keys()
    header= '                   '
    for key in kast.kastutil.listhd:
        if key in ['DATE-OBS','EXPTIME','OBJECT']:
            if key in ['DATE-OBS']:
                header = header + '%25s' % (key)
            else:
                header = header + '%15s' % (key)
    header = header + '%25s  %25s  %25s\n' % ('extracted','wave','flux')


    for arm in setup_object.keys():
        print('#'*10 + '  STANDARDS ' +'#'*10 + '\n')
        print(header)
        for img in setup_standard[arm]:
            imgex = re.sub('.fits','_ex.fits',img)
            imgl = re.sub('.fits','_l.fits',img)
            imgf = re.sub('.fits','_f.fits',img)
            value = img + ' '
            for key in kast.kastutil.listhd:
                if key in ['DATE-OBS','EXPTIME','OBJECT']:
                    if key in ['DATE-OBS']:
                        value = value + '%25s' % (str(dictionary[img][key]))
                    else:
                        value = value + '%16s' % (str(dictionary[img][key]))
                        
            if os.path.isfile(imgex):   value = value + '%25s' % (imgex)
            else: value = value + '%25s' % ('X')
            if os.path.isfile(imgl):   value = value + '%25s' % (imgl)
            else: value = value + '%25s' % ('X')
            if os.path.isfile(imgf):   value = value + '%25s' % (imgf)
            else: value = value + '%25s' % ('X')
            print(value)
            
        print('#'*10 + '  OBJECT ' +'#'*10 + '\n')
        print(header)
        for img in setup_object[arm]:
            imgex = re.sub('.fits','_ex.fits',img)
            imgl = re.sub('.fits','_l.fits',img)
            imgf = re.sub('.fits','_f.fits',img)
            value = img + ' '
            for key in kast.kastutil.listhd:
                if key in ['DATE-OBS','EXPTIME','OBJECT']:
                    if key in ['DATE-OBS']:
                        value = value + '%25s' % (str(dictionary[img][key]))
                    else:
                        value = value + '%15s' % (str(dictionary[img][key]))
            if os.path.isfile(imgex):   value = value + '%25s' % (imgex)
            else: value = value + '%25s' % ('X')
            if os.path.isfile(imgl):   value = value + '%25s' % (imgl)
            else: value = value + '%25s' % ('X')
            if os.path.isfile(imgf):   value = value + '%25s' % (imgf)
            else: value = value + '%25s' % ('X')
            print(value)
        
        print('#'*10 + '  ARCS ' +'#'*10 + '\n')
        print(header)
        for img in setup_arc[arm]:
            imgex = re.sub('.fits','_ex.fits',img)
            value = img + ' '
            for key in kast.kastutil.listhd:
                if key in ['DATE-OBS','EXPTIME','OBJECT']:
                    if key in ['DATE-OBS']:
                        value = value + '%25s' % (str(dictionary[img][key]))
                    else:
                        value = value + '%15s' % (str(dictionary[img][key]))
                        
            if os.path.isfile(imgex):   value = value + '%25s' % (imgex)
            else: value = value + '%25s' % ('X')
            print(value)

   ##############################################

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
    ############################################################
    if _run & 2 ==2:
        for key in objectlist['obj'].keys():
            for img in objectlist['obj'][key]:
                print(img)
                imgdata, imghdr = fits.getdata(img, header=True)
                if imghdr.get('LACOSMIC'):
                    print('already done')
                else:
                    if os.path.isfile('c' + img):
                        os.remove('c' + img)
                    kast.cosmics.lacos(img, output='c' + img, gain=_gain, readn=_rdnoise, xorder=9, yorder=9,
                                     sigclip=3, sigfrac=0.5, objlim=1, verbose=True, interactive=True)
                    kast.kastutil.updateheader('c' + img, 0, {'LACOSMIC': [True, 'Laplacian cosmic ray rejection']})
                    os.remove(img)
                    os.rename('c'+img,img)
                    print('\n### cosmic rays rejections ........ done ')
                
    print('# extract spectra  and arcs ....')
    ######### OBJECT ###############
    if _run & 4 == 4:
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
                    _review = 'yes'
                    _mode = 'q'
                    _type = 'obj'
                    kast.kastutil.extractspectrum(img,imgex,_reference,_trace,
                                                  _fittrac,_find,_recenter,_edit,
                                                  _gain,_rdnoise,dv,_interiraf,
                                                  _review,_type,_mode, _arm, key, _force)
                    
                    if observedarc:
                        # need to find a better way to select the arc (close in JD)
                        arcfile = observedarc[0]
                        print(img,imgex,arcfile,_arm)
                        arcfileex =  kast.kastutil.arcextraction(arcfile, img, imgex, _arm,dv,_force)
                    else:
                        print('no arc with this setup')
                        

                    
                    #raw_input('stop')
#    ######### STD ###############
    if _run & 4 == 4:
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
                    #_interactive = 'yes'
                    _review = 'yes'
                    _mode = 'q'
                    _type = 'std'
                    kast.kastutil.extractspectrum(img,imgex,_reference,_trace,
                                                  _fittrac,_find,_recenter,_edit,
                                                  _gain,_rdnoise,dv,_interiraf,
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
    if _run & 8 == 8:
        for key in objectlist['obj'].keys():
            print(key)
            _arm = key[0]
            _disp = key[1]
            _dicroic = key[2]
            _slit = key[3]
            for img in objectlist['obj'][key]:
                imgex = os.path.splitext(img)[0] + '_ex.fits'
                imgl = os.path.splitext(img)[0] + '_l.fits'
                arcfilex = 'arc_' + imgex
                run = True
                if os.path.isfile(imgex) is False:
                    run = False
                    print(img,'  run first extraction')
                else:
                    if os.path.isfile(imgl):
                        if _force:
                            os.remove(imgl)
                        else:
                            print('already wavelength calibrated')
                            run = False
                if os.path.isfile(arcfilex) is False:
                    print('arcfile not found')
                    directory = kast.__path__[0] + '/archive/' + str(_arm) + '/arc/' + _disp + '/' + _dicroic 
                    listarc = glob.glob(directory + '/*fits')                    
                    if not listarc:
                        print('no arc with this setup in the archive, stop reduction ')
                        run = False
                    else:
                        shutil.copy(listarc[0], os.path.basename(listarc[0]))
                        if not os.path.isdir('database'): os.mkdir('database')
                        shutil.copy(directory +'/database/id' + re.sub('.fits','',os.path.basename(listarc[0])),\
                                    'database/id'+ re.sub('.fits','',os.path.basename(listarc[0])))
                        arcfilex = os.path.basename(listarc[0])
                        
                if run is True:
                        directory = kast.__path__[0] + '/archive/' + str(_arm) + '/arc/' + _disp + '/' + _dicroic 
                        listarc = glob.glob(directory + '/*fits')                    
                        if not listarc:
                            imgl = kast.kastutil.identify(arcfilex, img, _arm, dv, arcref = False, force=_force, interactive=_interiraf)
                        else:
                            _arcref = listarc[0]
                            print('#######',_arcref)
                            imgl= kast.kastutil.identify(arcfilex, img, _arm, dv, arcref = _arcref, force=_force, interactive=_interiraf)
                        
                        if _arm == 'kastr':
                            _skyfile = kast.__path__[0]+'/standard/ident/sky_red.fits'
                        else:
                            _skyfile = kast.__path__[0]+'/standard/ident/sky_blu.fits'
                        kast.kastutil.checkwavelength_obj(imgl, _skyfile, _interiraf,True, arm = _arm)
                        
################################################################################
    #####   wavelengh calibration standard
    if _run & 8 == 8:
        for key in objectlist['std'].keys():
            print(key)
            _arm = key[0]
            _disp = key[1]
            _dicroic = key[2]
            _slit = key[3]
            for img in objectlist['std'][key]:
                imgex = os.path.splitext(img)[0] + '_ex.fits'
                imgl = os.path.splitext(img)[0] + '_l.fits'
                arcfilex = 'arc_' + imgex
                run = True
                if os.path.isfile(imgex) is False:
                    print(img,'  run first extraction')
                    run = False
                if os.path.isfile(arcfilex) is False:
                    print(img,' arc not found')
                    run = False

                if run is True:                                
                    if os.path.isfile(imgl):
                        if _force:
                            os.remove(imgl)
                        else:
                            print('already wavelength calibrated')
                            run = False
                            
                    if run is True:            
                            directory = kast.__path__[0] + '/archive/' + str(_arm) + '/arc/' + _disp + '/' + _dicroic 
                            listarc = glob.glob(directory + '/*fits')                    
                            if not listarc:
                                imgl = kast.kastutil.identify(arcfilex, img, _arm, dv, arcref = False, force=_force, interactive=_interiraf)
                            else:
                                _arcref = listarc[0]
                                print('#######',_arcref)
                                imgl = kast.kastutil.identify(arcfilex, img, _arm, dv, arcref = _arcref, force=_force, interactive=_interiraf)
            
                            if _arm == 'kastr':
                                _skyfile = kast.__path__[0] + '/standard/ident/sky_new_0.fits'
                                kast.kastutil.checkwavestd(imgl, _skyfile, _interiraf, True, arm = _arm)
                            else:
                                _skyfile = kast.__path__[0]+'/standard/ident/sky_blu.fits'
                                kast.kastutil.checkwavelength_obj(imgl, _skyfile, _interiraf, True, arm = _arm)
#######################################################################################################
    #########   sens function
    if _run & 16 == 16:
        if _senslist is not None:
            print('use external sensitivity function, skip sens stage')
        else:
            for key in objectlist['std'].keys():
                print(key)
                if key[1] in ['600/7500','600/4310']:
                    split = False
                else:
                    split = True
                _output = re.sub('/','','_'.join(key[:-1]))
                for img in objectlist['std'][key]:
                    imgl = os.path.splitext(img)[0] + '_l.fits'
                    if os.path.isfile(imgl):
                        imgclean,atmofile = kast.kastutil.make_atmo(imgl)
                        kast.kastutil.sensfunc(imgclean, _output ,_key=key,_split= split,  _function='spline3',\
                                               _order=8, interactive=_interiraf,force=_force)
                    else:
                        print('run first wavelengh solution for standard')
                    
#######################################################################################################
    ######## calib spectra
    if _run & 32 == 32:
        for key in objectlist['obj'].keys():
            _arm = key[0]
            _disp = key[1]
            _slit = key[3]
            _dicroic = key[2]
            run = True
            if _senslist is not None:
                senslist = _senslist
            else:
                _output = re.sub('/','','_'.join(key[:-1]))
                senslist = glob.glob('sens*'+_output+'*fits')

            if len(senslist)==0:
                print('No senitivity function founded')
                print('search in the archive')
                directory = kast.__path__[0] + '/archive/' + str(_arm) + '/sens/' + _disp + '/' + _dicroic 
                senslist = glob.glob(directory+'/sens*fits')
            if len(senslist)==0:                
                run = False
            if run:
                if _interactive:
                    for i,j in enumerate(senslist):
                        print(i,j)
                    answ = kast.kastutil.ask('which sens function do you want to use ? [0] ')
                    if not answ:
                        answ = 0
                    sensfile = senslist[int(answ)]
                else:
                    sensfile= senslist[0]
                    
                for img in objectlist['obj'][key]:
                    imgl = os.path.splitext(img)[0] + '_l.fits'
                    imgf = os.path.splitext(img)[0] + '_f.fits'
                    kast.kastutil.calibrate(imgl, imgf, sensfile ,  force=_force, interactive=_interiraf)
            
######################################################################################
################################################
    if _run & 64 == 64:
        objectlist1 = {}
        for key in objectlist['obj'].keys():
            for img in objectlist['obj'][key]:
                imgdata, imghdr = fits.getdata(img, header=True)
                _object = imghdr.get('OBJECT')
                #_date = parse(imghdr.get('DATE-OBS')).strftime("%Y%m%d")
                if _object not in objectlist1:
                    objectlist1[_object]=[]
                objectlist1[_object].append(re.sub('.fits','_f.fits',img))
        
        objectlist2 = {}
        for _object in objectlist1.keys():
            if _object not in objectlist2:
                objectlist2[_object]=[]
                
            imgdata, imghdr = fits.getdata(objectlist1[_object][0], header=True)
            _object = imghdr.get('OBJECT')
            _date = parse(imghdr.get('DATE-OBS')).strftime("%Y%m%d")
        
            _output = 'kast_'+ _object +'_'+_date+'_blue.fits'
            if _output not in objectlist2:
                objectlist2[_object].append(_output)
                
            _output = 'kast_'+ _object +'_'+_date+'_red.fits'
            if _output not in objectlist2:
                objectlist2[_object].append(_output)
            
################################################
        for key in objectlist1:
            listcomb  = [i for i in objectlist1[key] if '_b' in i]
            _output  = [i for i in objectlist2[key] if '_blue' in i][0]
            kast.kastutil.combine_same_arm(listcomb, _output, _combine='average',
                                           _w1= 'INDEF',_w2= 5650,_scale = True, _sample= '4000:5000')
            
            listcomb  = [i for i in objectlist1[key] if '_r' in i]
            _output  = [i for i in objectlist2[key] if '_red' in i][0]
            kast.kastutil.combine_same_arm(listcomb, _output, _combine='average',_w1= 5450 ,
                                           _w2= 'INDEF',_scale = True, _sample= '6000:7000')

########################################################
        if _interactive:
            plt.figure()
            plt.ion()
        for key in objectlist2:
            if _interactive:
                plt.clf()
                ww0,ff0= kast.kastutil.readspectrum(objectlist2[key][0])
                plt.plot(ww0,ff0,'-b',label='blue')
                ww1,ff1= kast.kastutil.readspectrum(objectlist2[key][1])
                plt.plot(ww1,ff1,'-r',label='blue')
            try:
                _output = re.sub('blue','merge',objectlist2[key][0])
                output = kast.kastutil.combine_same_arm(objectlist2[key], _output, _combine='average',_w1= 'INDEF',\
                                                        _w2= 'INDEF',_scale = True, _sample= '5500:5600')
                if _interactive:
                    ww3,ff3= kast.kastutil.readspectrum(_output)
                    plt.plot(ww3,ff3,'-g',label='merge')
            except:
                print('error merging')
            if _interactive:
                plt.xlim(5200,6000)
                plt.ylim(np.percentile(ff1,1),np.percentile(ff1,99))
                plt.legend(ncol=1)
                answ = kast.kastutil.ask('stop here')



######################################################################################
    ######## atmo correction
    if _run & 128 == 128:
        run= True
        imglist = glob.glob('atmo*_r*fits')
        if len(imglist)==0:
            print('atmofile not found')
            run = False
        else:
            atmo = imglist[1]
            
        imglist = glob.glob('*merge.fits')
        if len(imglist)==0:
            run = False
            print('merge file not found')

        if run:
            for simg in imglist:
                kast.kastutil.correct_for_atmo(simg,atmo, _so2=None,_sh2o= None, _interactive = _interactive)
                
########################################################        
    ######## plot spectra
    if _showspectra:
        imglist = glob.glob('*merge_e.fits')
        if len(imglist)==0:
            imglist = glob.glob('*merge.fits')
        if len(imglist)!=0:
            output = 'kast_spectra.pdf'
            kast.kastutil.plotspectra(imglist,output, minmax=False)
        else:
            print('warning, final merge files not available yet')
########################################################        
    ######## plot sens
    if _showsens:
        imglist = glob.glob('*sens_*.fits')
        if len(imglist)!=0:
            output = 'kast_sensitivity.pdf'
            kast.kastutil.plotspectra(imglist,output, minmax=True)
        else:
            print('warning, final merge files not available yet')

            ####################
#        _color = ['b','r','g','m','c','k','orange','yellow','brown','b',(.3,.4,.5)]
#        plt.figure()
#        for key in objectlist1:
#            n = 0
#            plt.clf()
#            plt.title(key)
#            for img in objectlist1[key]:
#                n = n + 1
#                ww,ff= kast.kastutil.readspectrum(img)
#                plt.plot(ww,ff,'-',color='r',label=img)
#            plt.ylim(np.percentile(ff,1),np.percentile(ff,99))
#            plt.legend(ncol=1)
#            xx = kast.kastutil.ask('stop here')
####################################################################################


         
