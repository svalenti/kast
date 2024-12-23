import astropy
from astropy.io import fits
from astropy import time
import glob
import os
import re
from matplotlib import pylab as plt
import numpy as np
import math
import sys
#import pyds9
from astropy.convolution import convolve, Box1DKernel
from astropy.stats import sigma_clip
from astropy.coordinates import SkyCoord
import astropy.units as u
import kast
from dateutil.parser import parse

pyversion = sys.version_info[0]


def ask(question):
    if pyversion>=3:
        answ = input(question)
    else:
        answ = raw_input(question)
    return answ

    
listhd = ['MJD','EXPTIME','AIRMASS','OBJECT',\
          'VERSION','RA','DEC','DATE-OBS','LAMPSTAJ','LAMPSTAD',\
          'BSPLIT_N','SLIT_N','GRATNG_N','GRISM_N']

def readstandard():
    path = kast.__path__[0]
    std = np.loadtxt(path+'/standard/standard_kast.txt',str)
    rastd = []
    decstd = []
    namestd= []
    for line in std:
        c = SkyCoord(line[1],line[2],frame='fk5',unit=(u.hourangle,u.deg))
        rastd.append(c.ra.value)
        decstd.append(c.dec.value)
        namestd.append(line[0])
        
    rastd = np.array(rastd)
    decstd = np.array(decstd)
    return rastd, decstd, namestd
    
def checkalldata(directory=False,verbose=False, all=False):
    import kast
    from astropy.io import ascii
    if all:
        if directory:
            imglist = glob.glob(directory + '*fits') 
        else:
            imglist = glob.glob('*fits')
        imglist = [i for i in imglist if '_ex.fits' not in i]
        imglist = [i for i in imglist if '_l.fits' not in i]
        imglist = [i for i in imglist if '_f.fits' not in i]
        imglist = [i for i in imglist if 'atmo_' not in i]
        imglist = [i for i in imglist if 'std_' not in i]
        imglist = [i for i in imglist if 'sens_' not in i]
        imglist = [i for i in imglist if '_clean.fits' not in i]
        imglist = [i for i in imglist if '_std' not in i]
        imglist = [i for i in imglist if '_sens' not in i]
        imglist = [i for i in imglist if '_obj' not in i]
        imglist = [i for i in imglist if '_sub' not in i]
        imglist = [i for i in imglist if 'kast_' not in i]
        imglist = [i for i in imglist if '_atlas' not in i]
        
    else:
        if directory:
            imglist = glob.glob(directory + 'r*fits') + glob.glob(directory + 'b*fits')
        else:
            imglist = glob.glob('r*fits') + glob.glob('b*fits') 


#              'SLSELNAM','TARGNAME','LAMPS','FLAMPS','OBSTYPE','IMTYPE','RA','DEC']
    dictionary={}
    for img in imglist:
        dictionary[img]={}
        hdu = fits.open(img)
        dictionary[img]['fits']= hdu
        for _header in listhd:
            if hdu[0].header.get(_header) is not None:
                dictionary[img][_header]= hdu[0].header.get(_header)
            else:
                dictionary[img][_header]= None
                
        if dictionary[img]['MJD'] is None:
            if  dictionary[img]['DATE-OBS'] is not None:
                _mjd = time.Time(dictionary[img]['DATE-OBS']).mjd
                dictionary[img]['MJD'] = _mjd

                
    setup_object={'kastb':[],
                  'kastr':[]}
    setup_bias={'kastb':[],
                  'kastr':[]}
    setup_arc={'kastb':[],
                  'kastr':[]}
    setup_flat={'kastb':[],
                'kastr':[]}
    skip={'kastb':[],
          'kastr':[]}

    rastd, decstd, namestd = readstandard()
    for img in dictionary:
        
        _ra = dictionary[img]['RA']
        _dec = dictionary[img]['DEC']
        c = SkyCoord(_ra,_dec,frame='fk5',unit=(u.hourangle,u.deg))
        _ra = c.ra.value
        _dec = c.dec.value
        distance = 3600 * ((rastd - _ra)**2 + (decstd - _dec)**2)**.5
        
        if dictionary[img]['EXPTIME'] == 0:
            setup_bias[dictionary[img]['VERSION']].append(img)
        elif dictionary[img]['GRATNG_N'].lower() in ['gmirror']:
            skip[dictionary[img]['VERSION']].append(img)
        elif dictionary[img]['OBJECT'].lower() in ['flat','flats']:
            setup_flat[dictionary[img]['VERSION']].append(img)
        elif dictionary[img]['LAMPSTAD'] in ['on']:
            setup_arc[dictionary[img]['VERSION']].append(img)            
        elif dictionary[img]['OBJECT'].lower() in ['arcs','hgcdarnehe_arc']:
            setup_arc[dictionary[img]['VERSION']].append(img)
        elif dictionary[img]['LAMPSTAJ'].lower() in ['on']:# or dictionary[img]['LAMPSTAG'].lower() in ['on']:
            setup_arc[dictionary[img]['VERSION']].append(img)
        elif min(distance)> 10:
            setup_object[dictionary[img]['VERSION']].append(img)
        else:                
            skip[dictionary[img]['VERSION']].append(img)

    return dictionary,    setup_object, setup_flat, setup_arc, setup_bias, skip

######################################################################

def combinebias(biaslist, masterbias,_rdnoise,_gain, comb = 'median',rej = 'ccdclip'):
    from pyraf import iraf
    iraf.noao(_doprint=0, Stdout=0)
    iraf.imred(_doprint=0, Stdout=0)
    iraf.ccdred(_doprint=0, Stdout=0)
    with open("biaslist", "w") as outfile:
        outfile.write("\n".join(biaslist))
        
    if os.path.isfile(masterbias):
        os.remove(masterbias)
    iraf.ccdred.zerocombine('@biaslist', output=masterbias, combine=comb,
                            reject= rej, ccdtype=' ', rdnoise=_rdnoise, gain=_gain,
                            process='no', Stdout=1)
    return masterbias


def combineflat(flatlist, masterflat,_rdnoise,_gain, comb = 'average',rej = 'none'):
    from pyraf import iraf
    iraf.noao(_doprint=0, Stdout=0)
    iraf.imred(_doprint=0, Stdout=0)
    iraf.ccdred(_doprint=0, Stdout=0)
    
    with open("flatlist", "w") as outfile:
        outfile.write("\n".join(flatlist))
    if os.path.isfile(masterflat):
        os.remove(masterflat)    
    iraf.ccdred.flatcombine('"@flatlist"', output=masterflat, combine=comb,
                            reject=rej, ccdtype=' ', rdnoise=_rdnoise, gain=_gain,
                            process='no', Stdout=1)
    return masterflat


def ccdprocimage(image, _output,_trimcor,_overscan,_zerocor,_flatcor, _readaxi, _zero,_trimsec,
                 _biassec, _flat, direction=2):
    from pyraf import iraf
    iraf.noao(_doprint=0, Stdout=0)
    iraf.imred(_doprint=0, Stdout=0)
    iraf.ccdred(_doprint=0, Stdout=0)
    iraf.specred(_doprint=0, Stdout=0)
    iraf.specred.dispaxis = direction
        
    if os.path.isfile(_output):
        os.remove(_output)
        
    iraf.ccdproc(image, output=_output,darkcor='no',fixpix='no', readaxi=_readaxi, ccdtype='',
                 overscan=_overscan,    biassec=_biassec,
                 trim=_trimcor,         trimsec=_trimsec,
                 zerocor=_zerocor,      zero=_zero,
                 flatcor=_flatcor,      flat=_flat,
                 Stdout=1)

    
def responseflat(_calibration , _normalization, _output, _order= 80, function= 'spline3',direction =2, _arm = 'kastb', _interactive = 'yes'):
    from pyraf import iraf
    iraf.noao(_doprint=0, Stdout=0)
    iraf.imred(_doprint=0, Stdout=0)
    iraf.ccdred(_doprint=0, Stdout=0)
    iraf.specred(_doprint=0, Stdout=0)
    iraf.specred.dispaxis = direction
    
    if os.path.isfile(_output):
        os.remove(_output)
        
    iraf.specred.response(calibration=_calibration , normalization = _normalization,  response=_output,
                          order= _order, function= 'spline3', interactive= _interactive)
    if _arm == 'kastb':
        iraf.imreplace(_output+'[1750:1900,*]',value =1, lower= 0.99, upper= 'INDEF')
        iraf.imreplace(_output+'[1750:1900,*]',value =1, lower= 'INDEF', upper= 1.01)
        iraf.imreplace(_output+'[1:250,*]',value =1, lower= 0.99, upper= 'INDEF')
        iraf.imreplace(_output+'[1:250,*]',value =1, lower= 'INDEF', upper= 1.01)
    elif _arm == 'kastr':
        iraf.imreplace(_output+'[*,1960:2140]',value =1, lower= 0.99, upper= 'INDEF')
        iraf.imreplace(_output+'[*,1960:2140]',value =1, lower= 'INDEF', upper= 1.01)
        iraf.imreplace(_output+'[*,1:170]',value =1, lower= 0.99, upper= 'INDEF')
        iraf.imreplace(_output+'[*,1:170]',value =1, lower= 'INDEF', upper= 1.01)
        iraf.imreplace(_output+'[*,1570:1670]',value =1, lower= 0.99, upper= 'INDEF')
        iraf.imreplace(_output+'[*,1570:1670]',value =1, lower= 'INDEF', upper= 1.01)
    else:
        print('arm not defined')
        
#############################################################################
#Setup: kastb-452/3306-d57-2.0 arcsec
#Setup: kastr-300/7500-d57-2.0 arcsec

def dvex():
    dv = {}
    dv['line'] = {'452/3306': 1000, '300/7500': 1200, }
    dv['std'] = {'_t_order': 6, '_t_niter': 50, '_t_sample': '*', '_t_nlost': 20, '_width': 10, '_radius': 10,
                 '_weights': 'variance',
                 '_nsum': 30, '_t_step': 10, '_t_nsum': 10, '_lower': -10, '_upper': 10, '_b_sample': '-40:-20,20:40',
                 '_resize': 'no'}
    dv['obj'] = {'_t_order': 3, '_t_niter': 50, '_t_sample': '*', '_t_nlost': 20, '_width': 10, '_radius': 10,
                 '_weights': 'variance',
                 '_nsum': -400, '_t_step': 10, '_t_nsum': 10, '_lower': -5, '_upper': 5, '_b_sample': '-25:-15,15:25',
                 '_resize': 'yes'}
    dv['dispaxis'] = {'kastb':1,'kastr':2}
    dv['ident']= {'cradius':10, 'fwhm':7, 'function':'legendre','order':5}
    return dv

def extractspectrum(img,imgex,_reference,_trace,_fittrac,_find,_recenter,_edit,
                    _gain,_rdnoise,dv,_interactive,_review,_type,_mode,arm, setup ,_force=False):

    run = True
    if os.path.isfile(imgex):
        if _force:
            os.remove(imgex)
        elif _interactive=='yes':
            hdr0 = fits.getheader(imgex)
            _ob = hdr0.get('OBJECT')
            print(_ob)
            answ = kast.kastutil.ask('Already extracted. do you want to extract again? [Y/[N]] ')
            if answ.lower() in ['no','n','']:
                run = False
        else:
            print('already extracted')
            run = False

    if run:
        from pyraf import iraf
        iraf.noao(_doprint=0, Stdout=0)
        iraf.imred(_doprint=0, Stdout=0)
        iraf.ccdred(_doprint=0, Stdout=0)
        iraf.specred(_doprint=0, Stdout=0)
        iraf.specred.dispaxis = dv['dispaxis'][arm]
        dist = dv['line'][setup[1]]
        iraf.specred.apall(img, output=imgex, referen=_reference, trace=_trace, fittrac=_fittrac, find=_find,
                           recenter=_recenter, edit=_edit,
                           nfind=1, extract='yes', backgro='fit', gain=_gain, readnoi=_rdnoise, lsigma=4, usigma=4,
                           format='multispec',
                           b_function='legendre', b_sample=dv[_type]['_b_sample'], clean='yes', pfit='fit1d',
                           lower=dv[_type]['_lower'], upper=dv[_type]['_upper'], t_niter=dv[_type]['_t_niter'],
                           width=dv[_type]['_width'],radius=dv[_type]['_radius'], line=dist,
                           nsum=dv[_type]['_nsum'], t_step=dv[_type]['_t_step'],
                           t_nsum=dv[_type]['_t_nsum'], t_nlost=dv[_type]['_t_nlost'], t_sample=dv[_type]['_t_sample'],
                           resize=dv[_type]['_resize'], t_order=dv[_type]['_t_order'],
                           weights=dv[_type]['_weights'], interactive=_interactive, review=_review, mode=_mode)


def arcextraction(arcfile, img, imgex, arm,dv,_force =False):
    run = True
    arcfilex = 'arc_' + imgex
    if os.path.isfile(arcfilex):
        if _force:
            os.remove(arcfilex)
        else:
            print('already extracted')
            run = False

    if run:
        from pyraf import iraf
        iraf.noao(_doprint=0, Stdout=0)
        iraf.imred(_doprint=0, Stdout=0)
        iraf.ccdred(_doprint=0, Stdout=0)
        iraf.specred(_doprint=0, Stdout=0)
        iraf.specred.dispaxis = dv['dispaxis'][arm]
        print(dv['dispaxis'][arm],arm)
        iraf.specred.apsum(arcfile, output='arc_' + imgex, referen=img, interac='no', find='no',
                                       recente='no', resize='no', edit='no', trace='no', fittrac='no',
                                       extract='yes', extras='no', review='no', backgro='none')

        return arcfilex
    



def identify(arcfilex, img, arm, dv, arcref=False, force =False, interactive = 'no'):
    imgl = os.path.splitext(img)[0] + '_l.fits'
    imgex = os.path.splitext(img)[0] + '_ex.fits'
    run = True

    if run is True:
        from pyraf import iraf
        iraf.noao(_doprint=0, Stdout=0)
        iraf.imred(_doprint=0, Stdout=0)
        iraf.ccdred(_doprint=0, Stdout=0)
        iraf.specred(_doprint=0, Stdout=0)
        iraf.specred.dispaxis = dv['dispaxis'][arm]
        iraf.set(direc=kast.__path__[0] + '/')

        print(dv['dispaxis'][arm],arm)
        if arcref is False:
            identific = iraf.specred.identify(images=arcfilex, section='middle line',
                                              coordli='direc$standard/ident/licklinelist.dat',
                                              nsum=10, fwidth=dv['ident']['fwhm'],
                                              cradius=dv['ident']['cradius'],
                                              function=dv['ident']['function'], order=dv['ident']['order'],
                                              mode='h', Stdout=1)
        else:
            _shift = 0
            os.system('cp '+ arcref + ' ./')
            databasename = os.path.dirname(arcref) +'/database/id' +re.sub('.fits','',os.path.basename(arcref))
            if not os.path.exists('database'):  os.makedirs('database/')
            os.system('cp '+ databasename + ' database/' )
            arcref0 = os.path.basename(arcref)            
            identific = iraf.specred.reidentify(referenc=arcref0, images=arcfilex,
                                                interac= interactive, section='middle line',
                                                shift=_shift, coordli='direc$standard/ident/licklinelist.dat',
                                                overrid='yes', cradius=dv['ident']['cradius'],
                                                step=0,
                                                newaps='no', nsum=5, nlost=2, mode='h',
                                                verbose='yes', Stdout=1)

        hedvec = {'REFSPEC1': [re.sub('.fits', '', arcfilex), ' reference arc']}
        updateheader(imgex, 0, hedvec)
        iraf.specred.dispcor(imgex, output=imgl, flux='yes')        
    return imgl
    
#def searcharc(img, listarc):
#    if not lisarc:
#        imglist = glob.glob(kast.__path__)
#    else:
#        pass
#    return imglist


def updateheader(filename, dimension, headerdict):
    from astropy.io import fits
    tupledict = {key: tuple(value) for key, value in headerdict.items()}
    try:
        hdulist = fits.open(filename, mode='update')
        header = hdulist[dimension].header
        header.update(tupledict)
        hdulist.close()
    except Exception as e:
        print 'header of', filename, 'not updated:'
        print e

#################################################################
        
def readspectrum(img):
    from numpy import array
    from astropy.io import fits
    import string

    fl = ''
    lam = ''
    graf = 1
    spec = fits.open(img)
    head = spec[0].header
    try:
        if spec[0].data.ndim == 1:
            fl = spec[0].data
        elif spec[0].data.ndim == 2:
            fl = spec[0].data[:, 0]
        elif spec[0].data.ndim == 3:
            fl = spec[0].data[0, 0, :]
    except:
        if spec[0].data.rank == 1:
            fl = spec[0].data
        elif spec[0].data.rank == 2:
            fl = spec[0].data[:, 0]
        elif spec[0].data.rank == 3:
            fl = spec[0].data[0, 0, :]
    naxis1 = head['naxis1']
    try:
        crpix1 = head['crpix1']
        crval1 = head['crval1']
        try:
            cdelt1 = head['cdelt1']
        except:
            cdelt1 = head['cd1_1']
        pix = array(range(1, naxis1 + 1, 1))
        pix = array(range(1, len(fl) + 1, 1))
        lam = (pix - crpix1) * cdelt1 + crval1
    except:
        try:
            WAT = head['WAT2_001']
            pix = array(range(1, naxis1 + 1, 1))
            crpix1 = string.split(string.split(WAT, '"')[1])[0]
            crval1 = string.split(string.split(WAT, '"')[1])[3]
            cdelt1 = string.split(string.split(WAT, '"')[1])[4]
            lam = (pix - float(crpix1)) * float(cdelt1) + float(crval1)
        except:
            graf = 0
    return lam, fl

##############################################################

def make_atmo(stdfile):
   xx,yy = kast.kastutil.readspectrum(stdfile)
   zz = (( 6820 > xx ) | ( xx> 7090)) &\
       (( 7140 > xx ) | ( xx> 7420)) &\
       (( 7570 > xx ) | ( xx> 7750)) &\
       (( 7820 > xx ) | ( xx> 8450)) &\
       (( 8910 > xx ) | ( xx> 9225)) &\
       (( 9267 > xx ) | ( xx> 9890))

   yynew = np.interp(xx,xx[zz],yy[zz])

   hdu = fits.open(stdfile)
   hdu[0].data[0] = yy/yynew
   atmofile = 'atmo_'+stdfile
   hdu.writeto(atmofile,overwrite=True)
   stdclean = re.sub('.fits','_clean.fits',stdfile)
   hdu[0].data[0] = yynew
   hdu.writeto(stdclean,overwrite=True)
   return stdclean,atmofile

#######################################################################


def sensfunc(standardfile, _output = None, _key=('kastb','x'), _split = False, _function= 'spline3',
             _order = 8, interactive = 'yes', force= True):
        from pyraf import iraf
        iraf.noao(_doprint=0, Stdout=0)
        iraf.imred(_doprint=0, Stdout=0)
        iraf.ccdred(_doprint=0, Stdout=0)
        iraf.specred(_doprint=0, Stdout=0)        
        iraf.set(direc=kast.__path__[0] + '/')
        _caldir = 'direc$standard/MAB/'
        _extinctdir = 'direc$standard/extinction/'
        
        _extinction = 'lickext.dat'
        _observatory = 'lick'
        hdu = fits.open(standardfile)
        _airmass = hdu[0].header['airmass']
        _exptime = hdu[0].header['exptime']
        _date = parse(hdu[0].header.get('DATE-OBS')).strftime("%Y%m%d%H%M")
        if _output is None:
            _outputstd= 'std_'+ str(_date) + '_' + standardfile
            _outputsens= 'sens_'+ str(_date) + '_' +  standardfile
        else:
            print(_output)
            _outputstd= 'std_'+ str(_date) + '_'+ _output
            _outputsens= 'sens_'+ str(_date) + '_'+ _output + '.fits'


        _run = True
        if os.path.isfile(_outputsens):
            if force:
                os.remove(_outputsens)
            else:
                _run = False
                print('sensitivity function already done')

        rastd, decstd, namestd = readstandard()        
        _ra = hdu[0].header['RA']
        _dec = hdu[0].header['DEC']
        c = SkyCoord(_ra,_dec,frame='fk5',unit=(u.hourangle,u.deg))
        _ra = c.ra.value
        _dec = c.dec.value
        distance = 3600 * ((rastd - _ra)**2 + (decstd - _dec)**2)**.5
        if np.min(distance)> 10:
            print('object not found in the list standard_kast.txt')
            refstar = 'INDEF'
            _run = False
        else:
            refstar = 'm'+namestd[np.argmin(distance)]
            dire = kast.__path__[0] + '/standard/MAB/'+ refstar + '.dat'
            print(dire)
            if not os.path.isfile(dire):
                print('standard table not found in the standard/MAB directory')
                _run = False
            
        if _run:                
            if os.path.isfile(_outputstd):
                    os.remove(_outputstd)
                    
            if _split:
                if _key[0]=='kastb':
                    _w01 = 3000
                    _w02 = 5300
                    _w11 = 5150
                    _w12 = 6000
                    sens0 = '_sensb0.fits'
                    sens1 = '_sensb1.fits'
                    _order0 = 10
                    _order1 = 25
                else:
                    _w01 = 5200
                    _w02 = 6000 
                    _w11 = 5850
                    _w12 = 12000
                    sens0 = '_sensr0.fits'
                    sens1 = '_sensr1.fits'
                    _order0 = 25
                    _order1 = 10
                    
                obj0 = '_obj0.fits'
                obj1 = '_obj1.fits'
                std1 = '_std1'
                std0 = '_std0'
                std1 = '_std1.fits'
        
                for _file in [obj0, obj1, std0, std1, sens0, sens1]:
                    if os.path.isfile(_file):
                        os.remove(_file)
        
                    
                iraf.scopy(standardfile, obj0, w1=_w01, w2=_w02)
                iraf.scopy(standardfile, obj1, w1=_w11, w2=_w12)
        
                ######
                iraf.specred.standard(input=obj0, output=std0, extinct=_extinctdir + _extinction,
                                      caldir=_caldir, observa=_observatory, star_nam=refstar, airmass=_airmass,
                                      exptime=_exptime, interac=interactive)
            
                iraf.specred.sensfunc(standard=std0, sensitiv=sens0, extinct=_extinctdir + _extinction,
                                      ignorea='yes', observa=_observatory, graphs='sri', functio=_function, order=_order0,
                                      interac=interactive)
                ######
                iraf.specred.standard(input=obj1, output=std1, extinct=_extinctdir + _extinction,
                                      caldir=_caldir, observa=_observatory, star_nam=refstar, airmass=_airmass,
                                      exptime=_exptime, interac=interactive)
            
                iraf.specred.sensfunc(standard=std1, sensitiv=sens1, extinct=_extinctdir + _extinction,
                                      ignorea='yes', observa=_observatory, graphs='sri', functio=_function, order=_order1 ,
                                      interac=interactive)
        
                combstring = sens0+','+sens1
                sss = iraf.specred.scombine(combstring, output=_outputsens, combine='average', reject='none',
                                    scale='none', weight='none', Stdout=1)
        
                updateheader(_outputsens, 0, {'standard': (refstar, '')})
                
                print('split')
            else:
                iraf.specred.standard(input=standardfile, output=_outputstd, extinct=_extinctdir + _extinction,
                                      caldir=_caldir, observa=_observatory, star_nam=refstar, airmass=_airmass,
                                      exptime=_exptime, interac=interactive)
                   
                iraf.specred.sensfunc(standard=_outputstd, sensitiv=_outputsens, extinct=_extinctdir + _extinction,
                                      ignorea='yes', observa=_observatory, graphs='sri', functio=_function, order=_order,
                                      interac=interactive)
                updateheader(_outputsens, 0, {'standard': (refstar, '')})

##############################################################            

def calibrate(imgl, imgf, sensfile, force = True, interactive = 'yes'):
    from pyraf import iraf
    iraf.noao(_doprint=0, Stdout=0)
    iraf.imred(_doprint=0, Stdout=0)
    iraf.ccdred(_doprint=0, Stdout=0)
    iraf.specred(_doprint=0, Stdout=0)

    run = True
    if os.path.isfile(imgf):
        if force:
            os.remove(imgf)
        else:
            print('already flux calibrated')
            run = False
            
    if run:   
        iraf.set(direc=kast.__path__[0] + '/')
        _caldir = 'direc$standard/MAB/'
        _extinctdir = 'direc$standard/extinction/'
        _extinction = 'lickext.dat'
        _observatory = 'lick'
        hdu = fits.open(imgl)
        _airmass = hdu[0].header['airmass']
        _exptime = hdu[0].header['exptime']
        print(imgl,imgf,sensfile)
        qqq = iraf.specred.calibrate(input=imgl, output=imgf, sensiti=sensfile, extinct='yes', flux='yes',
                                     extinction=_extinctdir + _extinction, observatory=_observatory,
                                     airmass=_airmass, ignorea='yes', exptime=_exptime, fnu='no')
        updateheader(imgf, 0, {'sensfun': (sensfile, '')})
                              
######################################################################################

def checkwavelength_arc(xx1, yy1, xx2, yy2, xmin, xmax, _interactive='yes'):
    from numpy import array, trapz, compress
    from numpy import interp as ninterp

    minimo = max(min(xx1), min(xx2)) + 60
    massimo = min(max(xx1), max(xx2)) - 60
    yy1 = [0 if e < 0 else e for e in array(yy1)]
    yy2 = [0 if e < 0 else e for e in array(yy2)]
    _shift, integral = [], []
    for shift in range(-600, 600, 1):
        xxnew = xx1 + shift / 10.
        yy2interp = ninterp(xxnew, xx2, yy2)
        yy2timesyy = yy2interp * yy1
        xxcut = compress((array(xxnew) >= minimo) & (array(xxnew) <= massimo), array(xxnew))
        yycut = compress((array(xxnew) >= minimo) & (array(xxnew) <= massimo), array(yy2timesyy))
        integrale = trapz(yycut, xxcut)
        integral.append(integrale)
        _shift.append(shift / 10.)
    result = _shift[integral.index(max(integral))]
    if _interactive in ['YES', 'y', 'Y', 'yes', 'Yes', True]:
        #   import matplotlib as mpl  
        #   mpl.use("TKAgg")  
        from pylab import plot, show, ion, clf, legend, xlim, ylim

        ion()
        clf()
        ratio = trapz(yy1, xx1) / trapz(yy2, xx2)
        yy3 = array(yy2) * float(ratio)
        xx4 = xx1 + result
        plot(xx1, yy1, label='spectrum')
        plot(xx2, yy3, label='reference sky')
        plot(xx4, yy1, label='shifted spectrum')
        legend(numpoints=1, markerscale=1.5)
        if xmin != '' and xmax != '':
            xlim(xmin, xmax)
    return result

###################################################################

def continumsub(imagefile, _order1=6, _order2=1,_low = 4, _high=3, _interactive = 'no' ):
    from pyraf import iraf
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    if os.path.isfile('tsky.fits'):
        os.remove('tsky.fits')
    toforget = ['specred.continuum']
    for t in toforget: iraf.unlearn(t)

    output = '_subtracted.fits'
    iraf.specred.continuum(imagefile, output='tsky.fits', type='difference', interact= _interactive, function='legendre',
                           niterat=300, low_rej=_low, high_re=_high, sample='*', order=_order1, ask='YES')
    iraf.specred.continuum('tsky.fits', output=output, type='difference', interact=_interactive, function='spline1',
                   overrid='yes', niterat=10, low_rej=4, high_re=2, sample='*', order=_order2, ask='YES')

    if os.path.isfile('tsky.fits'):
        os.remove('tsky.fits')
    return output

#######################################################

def checkwavelength_obj(fitsfile, skyfile, _interactive='yes', usethirdlayer=True, arm = 'kastr'):
    from astropy.io import fits
    import numpy as np
    from pyraf import iraf

    if arm=='kastr':
        minw = 5500
        maxw = 8000
    else:
        minw = 4000
        maxw = 6000
    
    if _interactive.lower() in ['yes', 'y']:
        do_shift = raw_input('### Do you want to check the wavelength calibration with telluric lines? [[y]/n] ')
    else:
        print '### Checking wavelength calibration with telluric lines'
        do_shift = ''
    if do_shift != 'n':
        if arm=='kastr':
            order1 = 8
            order2 = 1
            low = 4
            high = 3
        else:
            order1 = 30
            order2 = 1
            low = 5
            high = 3
        if usethirdlayer:
            iraf.scopy(fitsfile + '[*,1,3]', 'skylayer.fits',w1='INDEF',w2='INDEF') # iraf.continuum doesn't allow slices
            subtracted = continumsub('skylayer.fits', _order1=order1, _order2=order2, _low=low,_high=high, _interactive = 'no')
            os.remove('skylayer.fits')
        else:
            subtracted = continumsub(fitsfile, _order1=order1, _order2=order2, _low=low,_high=high, _interactive='no')
        sky_spec = fits.open(subtracted)[0]
        y1 = sky_spec.data
        crval1 = sky_spec.header['CRVAL1']
        x1 = crval1 + np.arange(len(y1)) * sky_spec.header['CD1_1']
        sky_arch = fits.open(skyfile)[0]
        y2 = sky_arch.data
        x2 = sky_arch.header['CRVAL1'] + np.arange(len(y2)) * sky_arch.header['CD1_1']
        shift = checkwavelength_arc(x1, y1, x2, y2, minw, maxw, _interactive)
        if _interactive.lower() in ['yes', 'y']:
            answ = raw_input('By how much do you want to shift the wavelength calibration? [{}] '.format(shift))
            if answ:
                shift = float(answ)
                
        updateheader(fitsfile, 0, {'CRVAL1': (crval1 + shift, ''), 'SHIFT': (shift, '')})
        if os.path.isfile(subtracted):
            os.remove(subtracted)
    else:
        shift = 0
    return shift

######################################################33
def checkwavestd(imgl, skyfile, _interactive='yes', _type=1, arm = 'kastr'):
    from astropy.io import fits
    from numpy import arange, array

    if arm=='kastr':
        minw = 5500
        maxw = 8000
    else:
        minw = 4000
        maxw = 6000
        
    print '\n### Warning: check in wavelength with sky lines not performed\n'
    if _interactive.upper() in ['YES', 'Y']:
        answ = raw_input('\n### Do you want to check the wavelength calibration with telluric lines [[y]/n]? ')
        if not answ: answ = 'y'
    else:
        answ = 'y'
    if answ in ['y', 'yes']:
        print('\n### check wavelength calibration with telluric lines \n')
        skydata, skyhdr = fits.getdata(skyfile, header=True)
        skyff = 1 - skydata
        crval1 = skyhdr['CRVAL1']
        cd1 = skyhdr['CD1_1']
        skyxx = np.arange(len(skyff))
        skyaa = crval1 + skyxx * cd1
        stdclean, atmofile = make_atmo(imgl)
        atmodata, atmohdr = fits.getdata(atmofile, header=True)
        if _type == 1:
            atmoff = 1 - atmodata[0][0]
        else:
            atmoff = 1 - atmodata
        crval1 = atmohdr['CRVAL1']
        cd1 = atmohdr['CD1_1']
        atmoxx = np.arange(len(atmoff))
        atmoaa = crval1 + atmoxx * cd1
        shift = checkwavelength_arc(atmoaa, atmoff, skyaa, skyff, minw, maxw,  _interactive)
        if _interactive.lower() in ['yes', 'y']:
            answ = raw_input('By how much do you want to shift the wavelength calibration? [{}] '.format(shift))
            if answ:
                shift = float(answ)
        if shift!=0:
            updateheader(imgl, 0, {'CRVAL1': (crval1 + shift, ''), 'SHIFT': (shift, '')})
#        if os.path.isfile(subtracted):
#            os.remove(subtracted)
#        floyds.util.delete('atmo2_' + _tel + '_' + imgex)
    else:
        shift = 0
    return shift

##############################################################
def delete(listfile):
    import os, string, re, glob

    if listfile[0] == '@':
        ff = open(listfile[1:])
        files = ff.readlines()
        imglist = []
        for ff in files:
            ff = re.sub(' ', '', ff)
            if not ff == '\n' and ff[0] != '#':
                ff = re.sub('\n', '', ff)
                imglist.append(ff)
    elif ',' in listfile:
        imglist = string.split(listfile, sep=',')
    else:
        imglist = [listfile]
    lista = []
    for _file in imglist:   lista = lista + glob.glob(_file)
    if lista:
        for _file in lista:
            if os.path.isfile(_file):
                try:
                    os.system('rm ' + _file)
                except:
                    pass


############################################################
def combine_same_arm(lista, _output, _combine='average',_w1= 'INDEF',_w2= 'INDEF',_scale = False, _sample= '4000:5000'):
    import re
    from pyraf import iraf
    iraf.noao(_doprint=0, Stdout=0)
    iraf.imred(_doprint=0, Stdout=0)
    iraf.ccdred(_doprint=0, Stdout=0)
    iraf.specred(_doprint=0, Stdout=0)
    iraf.unlearn(iraf.scombine)
    import time 
    hdr0 = fits.getheader(lista[0])
    if _scale: scomb_scale = 'median'
    else:     scomb_scale = 'none'
    datavecs = []
    hdrvec = []
    for n in range(0,4):
        list1= [i+'[*,1,'+str(n+1)+']' for i in lista]
        list2 = ','.join(list1)
        print(list2)
        output1 = '_'+str(n)+_output
        if os.path.isfile(output1):
            os.remove(output1)
        print(list2,output1,_combine,scomb_scale,_sample)
        iraf.specred.scombine(list2, w1='INDEF', w2='INDEF', output=output1, combine=_combine,\
                      scale=scomb_scale, sample=_sample)
        hdr = fits.open(output1)
        datavec, head = fits.getdata(output1, header=True) # these have shape (4440,)
        datavecs.append(datavec)
        hdrvec.append(head)
        os.remove(output1)        
        ####### this is just tkae the input file and replace the data
        #hdr0[0].data[n] = hdr[0].data
        #######
        time.sleep(1)

    for key in ['NAXIS1', 'CRVAL1', 'CD1_1', 'CRPIX1']:
        hdr0[key] = hdrvec[0][key]
    data3d = np.rollaxis(np.dstack(datavecs), 2)
    if os.path.isfile(_output):
        os.remove(_output)
    fits.writeto(_output, data3d, hdr0)
        
    #hdr0.writeto(_output,overwrite=True)
    if _w1!='INDEF' or _w2!='INDEF':
        kast.kastutil.delete(re.sub('.fits','_cut.fits',_output))
        iraf.scopy(_output,re.sub('.fits','_cut.fits',_output),w1=str(_w1),w2=str(_w2))
        kast.kastutil.delete(_output)
        os.rename(re.sub('.fits','_cut.fits',_output),_output)

    return _output
#############################################################
def plotspectra(imglist, output,minmax=False):
    fig,axs = plt.subplots(len(imglist),sharex=True)
    for j,img in enumerate(imglist):
        xx,yy = kast.kastutil.readspectrum(img)
        axs[j].plot(xx,yy,'-r',label=img.split('_')[1])
        minimo = np.min(yy[(xx>4000)&(xx<10000)])
        massimo = np.max(yy[(xx>4000)&(xx<10000)])
        if minmax is False:
            axs[j].set_ylim(np.percentile(yy[xx>4000],1),np.percentile(yy[xx>4000],99))
        else:
            axs[j].set_ylim(np.min(yy),np.max(yy))
#        axs[j].set_ylim(minimo,massimo)
        axs[j].legend(ncol=1)    

    fig.set_size_inches(6,11)
    plt.show()
    plt.savefig(output)

################################################################

def correct_for_atmo(simg,atmo, _so2=None,_sh2o= None, _interactive = False):
    import shutil
    from pyraf import iraf
    iraf.noao(_doprint=0, Stdout=0)
    iraf.imred(_doprint=0, Stdout=0)
    iraf.ccdred(_doprint=0, Stdout=0)
    iraf.specred(_doprint=0, Stdout=0)        
    iraf.set(direc=kast.__path__[0] + '/')
    print(simg,atmo)
    
    llatmo, ffatmo = kast.kastutil.readspectrum(atmo)
    llsci, ffsci = kast.kastutil.readspectrum(simg)
    
    llo2  = llatmo[(llatmo>=7550)&(llatmo<=7550)]
    ffo2  = ffatmo[(llatmo>=7550)&(llatmo<=7550)]
    llh2o = llatmo[(llatmo>=7100)&(llatmo<=7500)]
    ffh2o = ffatmo[(llatmo>=7100)&(llatmo<=7500)]
    _path = kast.__path__[0]

    _skyfileh2o = 'direc$standard/ident/ATLAS_H2O.fits'
    _skyfileo2 = 'direc$standard/ident/ATLAS_O2.fits'
    atlas_smooto2 = '_atlas_smoot_o2.fits'
    atlas_smooth2o = '_atlas_smoot_h2o.fits'
    _sigma = 200
    kast.kastutil.delete(atlas_smooto2)
    kast.kastutil.delete(atlas_smooth2o)
    iraf.imfilter.gauss(_skyfileh2o, output=atlas_smooth2o, sigma=_sigma)
    iraf.imfilter.gauss(_skyfileo2, output=atlas_smooto2, sigma=_sigma)
    llskyh2o, ffskyh2o = kast.kastutil.readspectrum(atlas_smooth2o)
    llskyo2, ffskyo2 = kast.kastutil.readspectrum(atlas_smooto2)
    
    ffskyo2cut = np.interp(llo2, llskyo2, ffskyo2)
    ffskyh2ocut = np.interp(llh2o, llskyh2o, ffskyh2o)
    _scaleh2o = []
    integral_h2o = []
    for i in range(1, 21):
            j = 0.6 + i * 0.04
            _ffskyh2ocut = (ffskyh2ocut * j) + 1 - j
            diff_h2o = abs(_ffskyh2ocut - ffh2o)
            integraleh2o = np.trapz(diff_h2o, llh2o)
            integral_h2o.append(integraleh2o)
            _scaleh2o.append(j)
    _scaleo2 = []
    integral_o2 = []
    for i in range(1, 21):
            j = 0.6 + i * 0.04
            _ffskyo2cut = (ffskyo2cut * j) + 1 - j
            diff_o2 = abs(_ffskyo2cut - ffo2)
            integraleo2 = np.trapz(diff_o2, llo2)
            integral_o2.append(integraleo2)
            _scaleo2.append(j)
    sh2o = _scaleh2o[np.argmin(integral_h2o)]
    so2 = _scaleo2[np.argmin(integral_o2)]
    if _so2 is not None and _sh2o is not None:
        telluric_features = (ffskyh2o * _sh2o) + 1 - _sh2o + (ffskyo2 * _so2) - _so2
    else:
        telluric_features = (ffskyh2o * sh2o) + 1 - sh2o + (ffskyo2 * so2) - so2
    telluric_features = np.array([1] + list(telluric_features) + [1])
    llskyo2 = np.array([1000] + list(llskyo2) + [15000])
    telluric_features_cut = np.interp(llsci, llskyo2, telluric_features)

    if _interactive is True:
        answ = 'y'
        while answ!='n':
            plt.clf()
            plt.plot(llsci,ffsci,'-r')
            plt.plot(llsci,ffsci/telluric_features_cut,'-b')
            plt.xlim(5000,10000)
            answ = kast.kastutil.ask('do you want to scale the telluric spectra manually  [y/n] [n]? ')
            if not answ: answ = 'n'
            if answ in ['yes','Y','y']:
                sh2o = kast.kastutil.ask(' h2o ? ' + str(sh2o) + ' ')
                sh2o = float(sh2o)
                so2 = kast.kastutil.ask(' o2 ? ' + str(so2) + ' ')
                so2 = float(so2)
                telluric_features = (ffskyh2o * sh2o) + 1 - sh2o + (ffskyo2 * so2) - so2
                telluric_features = np.array([1] + list(telluric_features) + [1])
                print(len(llskyo2),len(llsci),len(telluric_features))
                telluric_features_cut = np.interp(llsci, llskyo2, telluric_features)
    else:
        print('scale h2o file of ',sh2o)
        print('scale o2 file of ',so2)
                

    imge = re.sub('.fits','_e.fits',simg)
    shutil.copyfile(simg,imge)
    with fits.open(imge, mode='update') as hdus:
            hdus[0].data[0][0] = hdus[0].data[0][0]/telluric_features_cut
            hdus[0].data[1][0] = hdus[0].data[1][0]/telluric_features_cut
            hdus[0].data[2][0] = hdus[0].data[2][0]/telluric_features_cut
            hdus[0].data[3][0] = hdus[0].data[3][0]/telluric_features_cut
            hdus.close()  
            hdus.flush()  
                
    return llskyo2, telluric_features_cut

#################################################################
