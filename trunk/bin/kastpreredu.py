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
            if  dictionary[setup_object[arm][0]]['GRISM_N'] in ['600/7500','600/4310']:
                if arm =='kastb':
                    trimsec[arm] = '[1:2048,'+str(trima)+':'+str(trimb)+']'
                else:
                    trimsec[arm] =  '['+str(trima)+':'+str(trimb)+',24:2200]'
            else:    
                if arm =='kastb':
                    trimsec[arm] = '[1:1900,'+str(trima)+':'+str(trimb)+']'
                else:
                    trimsec[arm] =  '['+str(trima)+':'+str(trimb)+',60:2200]'
                    
            print('combine bias for ' + arm)
            if len(setup_bias[arm])>0:
                if _verbose:
                    for img in setup_bias[arm]:
                        ds9.set_np2arr(dictionary[img]['fits'][0].data)
                        answer = kast.kastutil.ask('good [y/n] [y]?')
                        if not answer:
                            answer = 'y'
                        if answer in ['n','N','no']:
                            setup_bias[arm].remove(img)
            
                masterbias = 'masterbias_' + arm + '.fits'
                _rdnoise = 1
                _gain = 1
                kast.kastutil.combinebias(setup_bias[arm], masterbias,_rdnoise,_gain, comb = 'median',rej = 'ccdclip')
                print('trim  bias '+ arm)
                kast.kastutil.ccdprocimage(masterbias,'t'+masterbias,_trimcor='yes',_overscan='no',_zerocor='no',_flatcor='no', 
                                           _zero = '', _biassec='', _trimsec = trimsec[arm], _flat = '',
                                           _readaxi= readaxi[arm],direction = specredaxis[arm])
            else:
                print('warning bias not there, using bias from archive')
                masterbias = 'masterbias_' + arm + '.fits'
                mast = kast.kastutil.searchbias(arm)
                if mast is not None:
                    shutil.copyfile(mast,masterbias)
                    print('trim  bias '+ arm)
                    kast.kastutil.ccdprocimage(masterbias,'t'+masterbias,_trimcor='yes',_overscan='no',_zerocor='no',_flatcor='no', 
                                               _zero = '', _biassec='', _trimsec = trimsec[arm], _flat = '',
                                               _readaxi= readaxi[arm],direction = specredaxis[arm])
                else:
                    print('skip bias')

            print('combine flat for '+arm)
            if len(setup_flat[arm])>0:            
                if _verbose:
                    for img in setup_flat[arm]:
                        ds9.set_np2arr(dictionary[img]['fits'][0].data)
                        answer = kast.kastutil.ask('good [y/n] [y]?')
                        if not answer:
                            answer = 'y'
                        if answer in ['n','N','no']:
                            setup_flat[arm].remove(img)
        
                masterflat = 'masterflat_' + arm + '.fits'
                _rdnoise = 1
                _gain = 1
                _order = 80
                if len(setup_flat[arm]):
                    kast.kastutil.combineflat(setup_flat[arm], masterflat,_rdnoise,_gain, comb = 'median',rej = 'ccdclip')
        
                    kast.kastutil.ccdprocimage(masterflat,'t'+masterflat,_trimcor='yes',_overscan='no',_zerocor='yes',_flatcor='no', 
                                               _zero ='tmasterbias_' + arm + '.fits', _biassec='', _trimsec = trimsec[arm], _flat = '',
                                               _readaxi= readaxi[arm],direction = specredaxis[arm])
            
                    kast.kastutil.responseflat('t'+masterflat , 't'+masterflat,  'n'+masterflat,
                                               _order, function= 'spline3',direction = specredaxis[arm],_arm=arm, _interactive = _interiraf)
            else:
                print('Warning no flats for ' + arm)

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
                    kast.kastutil.ccdprocimage(img,nameobj,_trimcor='yes',_overscan='no',_zerocor='yes',_flatcor='yes', 
                                               _zero ='tmasterbias_' + arm + '.fits', _biassec='', _trimsec = trimsec[arm],
                                               _flat = 'nmasterflat_' + arm + '.fits',
                                               _readaxi= readaxi[arm],direction = specredaxis[arm])
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
                    kast.kastutil.ccdprocimage(img,nameobj,_trimcor='yes',_overscan='no',_zerocor='yes',_flatcor='yes', 
                                               _zero ='tmasterbias_' + arm + '.fits', _biassec='', _trimsec = trimsec[arm],
                                               _flat = 'nmasterflat_' + arm + '.fits',
                                               _readaxi= readaxi[arm],direction = specredaxis[arm])
            else:
                print('warning no objects with this arm')
        else:
            print(proceed[arm])
##################################################################################
#####################   
#    print('pre-reduce standard')
#    for arm in setup_standard.keys():
#        if _verbose:
#            for img in setup_standard[arm]:
#                ds9.set_np2arr(dictionary[img]['fits'][0].data)
#                answer = kast.kastutil.ask('is this a science target [y/n] [y]?')
#                if not answer:
#                    answer = 'y'
#                if answer in ['n','N','no']:
#                    setup_standard[arm].remove(img)
#        for img in setup_standard[arm]:
#            nameobj = dictionary[img]['OBJECT'] + '_' + img
#            if _verbose: print(nameobj)
#            kast.kastutil.ccdprocimage(img,nameobj,_trimcor='yes',_overscan='no',_zerocor='yes',_flatcor='yes', 
#                                       _zero ='tmasterbias_' + arm + '.fits', _biassec='', _trimsec = trimsec[arm],
#                                       _flat = 'nmasterflat_' + arm + '.fits',
#                                       _readaxi= readaxi[arm],direction = specredaxis[arm])                    
            
#    list = dictionary.keys()
#    header= ''
#    #for key in dictionary[list[0]].keys():
#    for key in kast.kastutil.listhd:
#        if key in ['DATE-OBS']:
#            header = header + '%25s' % (key)
#        else:
#            header = header + '%15s' % (key)
#    print(header)
#    for img in setup_standard['kastr']:
#        value = ''
#        for key in kast.kastutil.listhd:
#            if key in ['DATE-OBS']:
#                value = value + '%25s' % (str(dictionary[img][key]))
#            else:
#                value = value + '%15s' % (str(dictionary[img][key]))
#        print(value)
        
    
   
