#!/usr/bin/env python
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
import shutil
import glob
import numpy as np
from matplotlib import pylab as plt
from astropy.io import fits
import pyds9
                
pyversion = sys.version_info[0]

# check that ds9 is open 
plt.ion()
ds9 = pyds9.DS9(str(time.time()))
ds9.set('frame 1')
ds9.set('scale zscale');

def trytrim(img,show=False):
    global ax1, ax2
    hdu = fits.open(img)
    data = hdu[0].data
    if hdu[0].header['version'] == 'kastr':
            arm='red'
            axis = 0
            order=410
    else:
            arm = 'blu'
            axis =1
            order =110
         
    mean = data.mean(axis)
    x = np.arange(len(mean))
    meanmean= np.mean(mean)
    mean2 = mean[np.argsort(mean)] 
    if 'Flat' in hdu[0].header['object']:
        xmin =x[mean>meanmean][0]+10
        xmax =x[mean>meanmean][-1]-10
    elif 'Bias' in hdu[0].header['object'] or 'Arc' in hdu[0].header['object']:
        print('I can not use arc of bias to find the trim')
        xmin,xmax = None, None
    else:
        xmin =x[mean>mean2[order]][0]+10
        xmax =x[mean>mean2[order]][-1]-10

    if xmin and show:
        ax1.clear()
        ax1.plot(x,mean2,'-c')
        print(xmin)
        print(xmax)

        ax2.plot(x,mean,'-r', label = 'mean')
        ax2.plot([xmin,xmin],[np.min(mean),np.max(mean)],'k-')
        ax2.plot([xmax,xmax],[np.min(mean),np.max(mean)],'k-')
        raw_input(' stop '+ img )
        ax1.clear()
        ax2.clear()            
    return xmin,xmax


if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description, version="%prog 1.0")
    parser.add_option("-d", "--directory", dest="directory", default=None, type="str",
                      help='reduce data in this directory \t [%default]')
    parser.add_option("--stage", dest="stage", default=None, type="str",
                      help='reduce data a single stage \t [%default, trim,trim+,sky,sky+,trace,trace+,extract,extract+,wave,wave+,flux,flux+]')   
    parser.add_option("-i", "--interactive", action="store_true",
                      dest='interactive', default=False, help='Interactive \t\t\ [%default]')
    
    option, args = parser.parse_args()
    _interactive = option.interactive
    _directory = option.directory
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
    dictionary, setup_object, setup_flat, setup_arc, setup_bias, setup_standard= kast.kastutil.checkalldata(directory=_directory,verbose=_verbose)
    #
    #  check that all data are in the directory and
    #
    # add standard to objects, we do not need standard at this time
    for key in setup_standard:
        for img in setup_standard[key]:
            setup_object[key].append(img)
    ############################
    _sizeobject = {}
    for key in setup_object:
        _sizeobject[key]=[]
        for img in setup_object[key]:
            xsize = dictionary[img]['DATASEC'][1:-1]
            if xsize not in _sizeobject[key]:
                _sizeobject[key].append(xsize)
    ############################
    ############################
    _sizeflat = {}
    for key in setup_flat:
        _sizeflat[key]=[]
        for img in setup_flat[key]:
            xsize = dictionary[img]['DATASEC'][1:-1]
            if xsize not in _sizeflat[key]:
                _sizeflat[key].append(xsize)
    ############################
    ############################
    _sizebias = {}
    for key in setup_bias:
        _sizebias[key]=[]
        for img in setup_bias[key]:
            xsize = dictionary[img]['DATASEC'][1:-1]
            if xsize not in _sizebias[key]:
                _sizebias[key].append(xsize)
    ############################
    ############################
    _sizearc = {}
    for key in setup_arc:
        _sizearc[key]=[]
        for img in setup_arc[key]:
            xsize = dictionary[img]['DATASEC'][1:-1]
            if xsize not in _sizearc[key]:
                _sizearc[key].append(xsize)
    ############################
    print('object ', _sizeobject)
    print('flat ', _sizeflat)
    print('bias', _sizebias)
    print('arc  ', _sizearc)
    ############################
    for key in _sizeobject:
        for setup in _sizeobject[key]:
            if setup not in _sizeflat[key]:
                print('warning: not flat with the same size ' + str(setup))
            if setup not in _sizearc[key]:
                print('warning: not arc with the same size ' + str(setup))
            if setup not in _sizebias[key]:
                print('warning: not bias with the same size ' + str(setup))
    ############################
    if _verbose:
        for key in setup_bias:
            if len(setup_bias[key])>0:
                for img in setup_bias[key]:
                    print('bias', key, img)
            else:
                print('Warning: no bias found for arm: ', key)
    ############################
        for key in setup_flat:
            if len(setup_flat[key])>0:
                for img in setup_flat[key]:
                    print('flat ', key, img)
            else:
                print('Warning: no flat found for arm: ', key)
    ############################
    proceed = {}
    for key in _sizeobject:
            if len(_sizeobject[key]) ==1:
                proceed[key]= [True, 'only one size']
            else:
                proceed[key] = [False,'images with different size, please split them in different directories']


    readaxi = {'kastr':'column', 'kastb': 'line'}
    trimsec = {'kastb':'[1:1900,20:320]', 'kastr': '[150:420,60:2200]'}
    specredaxis = {'kastb': 1, 'kastr': 2}
    
    for arm in _sizeobject:
        if  proceed[arm][0] is True:
            
            ##################
            masterbias = 'masterbias_' + arm + '.fits'
            if len(setup_bias[arm])>0:
                if _verbose:
                    for img in setup_bias[arm]:
                        ds9.set_np2arr(dictionary[img]['fits'][0].data)
                        answer = kast.kastutil.ask('good [y/n] [y]?')
                        if not answer:
                            answer = 'y'
                        if answer in ['n','N','no']:
                            setup_bias[arm].remove(img)
            
                _rdnoise = 1
                _gain = 1
                print('combine bias for ' + arm)
                kast.kastutil.combinebias(setup_bias[arm], masterbias,_rdnoise,_gain, comb = 'median',rej = 'ccdclip')
            else:
                print('warning bias not there, using bias from archive')
                masterbias = 'masterbias_' + arm + '.fits'
                mast = kast.kastutil.searchbias(arm)
                if mast is not None:
                    shutil.copyfile(mast,masterbias)
                else:
                    print('skip bias')
            
            #########################

            masterflat = 'masterflat_' + arm + '.fits'
            if len(setup_flat[arm])>0:            
                if _verbose:
                    for img in setup_flat[arm]:
                        ds9.set_np2arr(dictionary[img]['fits'][0].data)
                        answer = kast.kastutil.ask('good [y/n] [y]?')
                        if not answer:
                            answer = 'y'
                        if answer in ['n','N','no']:
                            setup_flat[arm].remove(img)
        
                _rdnoise = 1
                _gain = 1
                _order = 80
                if len(setup_flat[arm]):
                    print('combine flat for '+arm)
                    kast.kastutil.combineflat(setup_flat[arm], masterflat,_rdnoise,_gain, comb = 'median',rej = 'ccdclip')
            else:
                print('Warning no flats for ' + arm)
                
            ###################################
            ###### define trima and trimb 
            if os.path.exists(masterflat):
                hdr=fits.open(masterflat)
                data = hdr[0].data
                if arm =='kastr':
                    y =data.mean(0)
                    x = np.arange(len(y))
                    trima, trimb = x[y>np.average(y)][0] + 10 ,x[y>np.average(y)][-1] -10
                else:
                    y =data.mean(1)
                    x = np.arange(len(y))
                    trima, trimb = x[y>np.average(y)][0] + 10 , x[y>np.average(y)][-1] - 10
            else:
                print('no flat, try define trim using science image')

                if _verbose:
                    fig, (ax1,ax2) = plt.subplots(1,2)
                xminvec = []
                xmaxvec = []
                for img in setup_object[arm]:
                    xmin,xmax = trytrim(img,show=_verbose)
                    if xmin:
                        xminvec.append(xmin)
                        xmaxvec.append(xmax)
                        
                print(xminvec)
                print(xmaxvec)
                xmin = int(np.mean(xminvec))
                xmax = int(np.mean(xmaxvec))
                print(xmin,xmax)
                if xmax -xmin <  100:
                    print('did not work, use default')
                
                    ######### define Trim ########
                    if arm =='kastb':
                        range = _sizeobject[arm][0].split(',')[1].split(':')
                    else:
                        range = _sizeobject[arm][0].split(',')[0].split(':')
                    trima = int(range[1])/2 -150
                    trimb = int(range[1])/2+150
                    if trima<=0:
                        trima = 0
                    if trimb >= float(range[1]):
                        trimb = int(range[1])
                else:
                        trima, trimb = xmin, xmax
                        
            print(trima,trimb)
            #############################
            
            if  dictionary[setup_object[arm][0]]['GRISM_N'] in ['600/7500','600/4310']:
                if arm =='kastb':
                    trimsec[arm] = '[1:2048,'+str(trima)+':'+str(trimb)+']'
                else:
                    trimsec[arm] =  '['+str(trima)+':'+str(trimb)+',24:2200]'
            elif  dictionary[setup_object[arm][0]]['GRISM_N'] in ['452/3306']:
                if arm =='kastb': # it should be, but I leave the kastr option just in case
                    # this is not perfect. if there are file with mirror and file with dicroit will do something bad
                    if dictionary[setup_object[arm][0]]['BSPLIT_N'] in ['mirror']:
                        # if there is no dicroit I have flux in the red
                        print('warning: first file does not have dicroit, I assume you want all the flux in the red part of the blue channel')
                        trimsec[arm] = '[1:2000,'+str(trima)+':'+str(trimb)+']'
                    else:
                        trimsec[arm] = '[1:1900,'+str(trima)+':'+str(trimb)+']'
                else:
                    trimsec[arm] =  '['+str(trima)+':'+str(trimb)+',60:2200]'
            else:    
                if arm =='kastb':
                    trimsec[arm] = '[1:1900,'+str(trima)+':'+str(trimb)+']'
                else:
                    trimsec[arm] =  '['+str(trima)+':'+str(trimb)+',60:2200]'
                    
            ############################
            if os.path.exists(masterbias):
                print('trim  bias '+ arm)
                _rdnoise = 1
                _gain = 1
                kast.kastutil.ccdprocimage(masterbias,'t'+masterbias,_trimcor='yes',_overscan='no',_zerocor='no',_flatcor='no', 
                                           _zero = '', _biassec='', _trimsec = trimsec[arm], _flat = '',
                                           _readaxi= readaxi[arm],direction = specredaxis[arm])
            else:
                print('no masterbias')
                
            ######################
            if os.path.exists(masterflat):
                print('trim flat and normalize')
                _rdnoise = 1
                _gain = 1
                _order = 80
                kast.kastutil.ccdprocimage(masterflat,'t'+masterflat,_trimcor='yes',_overscan='no',_zerocor='yes',_flatcor='no', 
                                           _zero ='tmasterbias_' + arm + '.fits', _biassec='', _trimsec = trimsec[arm], _flat = '',
                                           _readaxi= readaxi[arm],direction = specredaxis[arm])
            
                kast.kastutil.responseflat('t'+masterflat , 't'+masterflat,  'n'+masterflat,
                                           _order, function= 'spline3',direction = specredaxis[arm],_arm=arm, _interactive = _interiraf)
            else:
                print('no masterflat')

            ######################
            print('pre-reduce objects')
            if len(setup_object[arm])>0:            
                if _verbose:
                    for img in setup_object[arm]:
                        print(dictionary[img]['OBJECT'],dictionary[img]['EXPTIME'])
                        ds9.set_np2arr(dictionary[img]['fits'][0].data)
                        answer = kast.kastutil.ask('is this a science target [y/n] [y]?')
                        if not answer:
                            answer = 'y'
                        if answer in ['n','N','no']:
                            setup_object[arm].remove(img)
                for img in setup_object[arm]:
                    nameobj = dictionary[img]['OBJECT'] + '_' + img
                    if _verbose: print(nameobj)

                    if os.path.exists('nmasterflat_' + arm + '.fits'):
                        kast.kastutil.ccdprocimage(img,nameobj,_trimcor='yes',_overscan='no',_zerocor='yes',_flatcor='yes', 
                                                   _zero ='tmasterbias_' + arm + '.fits', _biassec='', _trimsec = trimsec[arm],
                                                   _flat = 'nmasterflat_' + arm + '.fits',
                                                   _readaxi= readaxi[arm],direction = specredaxis[arm])
                    else:
                        kast.kastutil.ccdprocimage(img,nameobj,_trimcor='yes',_overscan='no',_zerocor='yes',_flatcor='no', 
                                                   _zero ='tmasterbias_' + arm + '.fits', _biassec='', _trimsec = trimsec[arm],
                                                   _flat = '', _readaxi= readaxi[arm],direction = specredaxis[arm])
            else:
                print('warning no objects with this arm')
                
            print('pre-reduce arc')
            if len(setup_arc[arm])>0:            
                if _verbose:
                    for img in setup_arc[arm]:
                        ds9.set_np2arr(dictionary[img]['fits'][0].data)
                        answer = kast.kastutil.ask('is this an arc  [y/n] [y]?')
                        if not answer:
                            answer = 'y'
                        if answer in ['n','N','no']:
                            setup_arc[arm].remove(img)
                for img in setup_arc[arm]:
                    nameobj = dictionary[img]['OBJECT'] + '_' + img
                    if _verbose: print(nameobj)

                    if os.path.exists('nmasterflat_' + arm + '.fits'):
                        kast.kastutil.ccdprocimage(img,nameobj,_trimcor='yes',_overscan='no',_zerocor='yes',_flatcor='yes', 
                                                   _zero ='tmasterbias_' + arm + '.fits', _biassec='', _trimsec = trimsec[arm],
                                                   _flat = 'nmasterflat_' + arm + '.fits',
                                                   _readaxi= readaxi[arm],direction = specredaxis[arm])
                    else:
                        kast.kastutil.ccdprocimage(img,nameobj,_trimcor='yes',_overscan='no',_zerocor='yes',_flatcor='no', 
                                                   _zero ='tmasterbias_' + arm + '.fits', _biassec='', _trimsec = trimsec[arm],
                                                   _flat = '', _readaxi= readaxi[arm],direction = specredaxis[arm])
            else:
                print('warning no arc file with this  arm')
        else:
            print(proceed[arm])
##################################################################################
