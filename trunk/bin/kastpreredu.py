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


    readaxi = {'kastr':'column', 'kastb': 'line'}
    trimsec = {'kastb':'[1:1900,20:320]', 'kastr': '[150:420,60:2200]'}
    specredaxis = {'kastb': 1, 'kastr': 2}
    
    ############################
    print('combine bias')
    for arm in setup_bias.keys():
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

    ############################        
    print('combine flat')    
    for arm in setup_flat.keys():
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
    ############################

    print('pre-reduce objects')
    for arm in setup_object.keys():
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
    
    print('pre-reduce standard')
    for arm in setup_standard.keys():
        if _verbose:
            for img in setup_standard[arm]:
                ds9.set_np2arr(dictionary[img]['fits'][0].data)
                answer = kast.kastutil.ask('is this a science target [y/n] [y]?')
                if not answer:
                    answer = 'y'
                if answer in ['n','N','no']:
                    setup_standard[arm].remove(img)
        for img in setup_standard[arm]:
            nameobj = dictionary[img]['OBJECT'] + '_' + img
            if _verbose: print(nameobj)
            kast.kastutil.ccdprocimage(img,nameobj,_trimcor='yes',_overscan='no',_zerocor='yes',_flatcor='yes', 
                                       _zero ='tmasterbias_' + arm + '.fits', _biassec='', _trimsec = trimsec[arm],
                                       _flat = 'nmasterflat_' + arm + '.fits',
                                       _readaxi= readaxi[arm],direction = specredaxis[arm])                    

    print('pre-reduce arc')
    for arm in setup_arc.keys():
        if _verbose:
            for img in setup_arc[arm]:
                ds9.set_np2arr(dictionary[img]['fits'][0].data)
                answer = kast.kastutil.ask('is this a science target [y/n] [y]?')
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
            
##################################################################################
            
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
        
    
   
