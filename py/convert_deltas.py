#!/usr/bin/env python
import glob
import numpy as np 
import sys 
import pyfits
import fitsio 
import scipy 
import math 
import pylab
import argparse 

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-d','--indir', type = str, default = None, required=True,
                        help = 'deltas input directory')
parser.add_argument('-o', '--out',type = str, default = None, required=True,
                        help = 'output fits file ')
parser.add_argument('--subtract', action="store_true", required=False,
                    help = 'subtract first the mean delta as a function of wavelength, then subtract a and b*dlwave terme per quasar')
args = parser.parse_args()

in_dir=args.indir
out_file=args.out
print 'Opening '+in_dir+' ...'

nspec=0 
RA=[]
DEC=[]
Z=[]
PLATE=[]
MJD=[]
FIBER=[]
THING_ID=[]

nfile=0
plate_name=[]
fi = glob.glob(in_dir+"/*.fits.gz")
for i,f in enumerate(fi):
    nfile+=1
    plate_name.append(f)
plate_name=sorted(plate_name)
print "number of file = ",nfile

wave_list=[]
data_list=[]
ivar_list=[]
id_list=[]


for i in range(nfile):
    print "%2.2f %%, ndata  = %i"%(float(i)/float(nfile)*100,nspec)
    sys.stdout.flush()
    hdu=fitsio.FITS(plate_name[i])
    for d in hdu[1:]:
        nspec+=1
        wave_spec=10**d["LOGLAM"][:]
        ivar_spec=d["WEIGHT"][:]
        data_spec=d["DELTA"][:]
        data_spec*=(ivar_spec>0) #Avoid interpolation problem 
        id_spec=d.read_header()['THING_ID']

        
        wave_list.append(wave_spec)
        data_list.append(data_spec)
        ivar_list.append(ivar_spec)
        id_list.append(id_spec)

        RA.append(d.read_header()['RA'])
        DEC.append(d.read_header()['DEC'])
        Z.append(d.read_header()['Z'])
        PLATE.append(d.read_header()['PLATE'])
        MJD.append(d.read_header()['MJD'])
        FIBER.append(d.read_header()['FIBERID'])
        THING_ID.append(d.read_header()['THING_ID'])
    
    hdu.close()

nq   = len(data_list)
print 'nq = ',nq

lwmin=None
lwmax=None
lwstep=None
for q in range(nq) :
    if lwmin is None :
        lwmin = np.log10(wave_list[q][0])
        lwmax = np.log10(wave_list[q][-1])
        lwstep = np.log10(wave_list[q][1])-np.log10(wave_list[q][0])
    else :
        lwmin = min(lwmin,np.log10(wave_list[q][0]))
        lwmax = max(lwmax,np.log10(wave_list[q][-1]))
        # check
        tmp_lwstep = np.log10(wave_list[q][1])-np.log10(wave_list[q][0])
        if np.abs(lwstep-tmp_lwstep)>1.e-12 :
            print "error with lwstep"

print "meas. lwstep",lwstep
lwstep=0.0003
print "true lwstep",lwstep
nw=int((lwmax-lwmin)/lwstep)+1
print "delta ?= 0 = ",nw-(lwmax-lwmin)/lwstep-1
wave=10**(np.arange(nw)*lwstep+lwmin)
print np.log10(wave[-1])-lwmax

data = np.zeros((nq,wave.size))
ivar = np.zeros((nq,wave.size))

print 'data.shape = ',data.shape
for q in range(nq) :
    i0=int((np.log10(wave_list[q][0])-lwmin+0.001*lwstep)/lwstep)
    #print 'q = %i, i0:i0+data_list[q].size = %i:%i'%(q,i0,i0+data_list[q].size)
    data[q,i0:i0+data_list[q].size] = data_list[q]
    ivar[q,i0:i0+data_list[q].size] = ivar_list[q]

RA=np.array(RA)*180./math.pi 
DEC=np.array(DEC)*180./math.pi 
Z=np.array(Z)
PLATE=np.array(PLATE)
MJD=np.array(MJD)
FIBER=np.array(FIBER)
THING_ID=np.array(THING_ID)

nw=data.shape[1]
nspec=data.shape[0]
print 'nw = ',nw, 
print 'nspec = ',nspec 

print 'Erase quasars with large chi2'
ndata=np.sum(ivar>0,axis=1)
bad=np.where((ndata<20)&(ndata>0))[0]
if bad.size > 0 :
    data[bad] *= 0
    ivar[bad] *= 0
    print "erase %d QSOs with less than 20 points"%bad.size


if not args.subtract: 
    c1=pyfits.Column(name='RA',format = 'D',array=RA)
    c2=pyfits.Column(name = 'DEC', format = 'D',array=DEC)
    c3=pyfits.Column(name = 'Z', format = 'D',array=Z)
    c4=pyfits.Column(name = 'PLATE', format = 'J',array=PLATE)
    c5=pyfits.Column(name = 'MJD', format = 'J',array=MJD)
    c6=pyfits.Column(name = 'FIBER', format = 'J',array=FIBER)
    c7=pyfits.Column(name = 'THING_ID', format = 'J',array=THING_ID)
    table_hdu = pyfits.new_table([c1,c2,c3,c4,c5,c6,c7])
    print 'writing '+out_file
    dflux = pyfits.HDUList([pyfits.PrimaryHDU(data)])
    dflux.append(pyfits.ImageHDU(ivar, name = "IVAR"))
    dflux.append(pyfits.ImageHDU(wave ,name="WAVELENGTH"))     
    dflux.append(table_hdu)
    dflux.writeto(out_file,clobber=True)
else: 
    print 'subtract first the mean delta as a function of wavelength, then subtract a and b*dlwave terme per quasar'
    delta_tilde=data.copy()
    lwave=np.log10(wave)

    print "subtracting mean vs wave"
    sw=np.sum(ivar,axis=0)
    swx=np.sum(ivar*delta_tilde,axis=0)
    mean=swx/(sw+(sw==0))
    delta_tilde -= mean*(ivar>0)
    for spec in range(nspec): 
        if (spec%10000==0): print 'spec = %i, spec/npec =%2.2f'%(spec,float(spec)/float(nspec)) 
        if np.sum(ivar[spec])==0 : 
            #print 'skip empty spec %d'%spec
            continue

        dlwave=lwave-np.sum(ivar[spec]*lwave)/np.sum(ivar[spec])
        a = np.sum(ivar[spec]*delta_tilde[spec])/np.sum(ivar[spec]) # weighted mean of delta_tilde
        b = np.sum(ivar[spec]*delta_tilde[spec]*dlwave)/np.sum(ivar[spec]*dlwave**2) # slope parameter
        delta_tilde[spec,:] -= ( a*np.ones(nw) + b*dlwave )
        delta_tilde[spec]*=(ivar[spec]>0)
    
    print 'Erase quasars with large chi2'   
    chi2=np.sum(ivar*delta_tilde**2,axis=1)
    ndata=np.sum(ivar>0,axis=1)
    bad=np.where(chi2>25*ndata)[0]
    if bad.size > 0 :
        delta_tilde[bad] *= 0
        ivar[bad] *= 0
        print "erase %d QSOs with chi2>25*ndata"%bad.size
        print "first chi2s"
        for i in bad[:10] :
            print "id=",i,"chi2/data=",chi2[i]/ndata[i]

    c1=pyfits.Column(name='RA',format = 'D',array=RA)
    c2=pyfits.Column(name = 'DEC', format = 'D',array=DEC)
    c3=pyfits.Column(name = 'Z', format = 'D',array=Z)
    c4=pyfits.Column(name = 'PLATE', format = 'J',array=PLATE)
    c5=pyfits.Column(name = 'MJD', format = 'J',array=MJD)
    c6=pyfits.Column(name = 'FIBER', format = 'J',array=FIBER)
    c7=pyfits.Column(name = 'THING_ID', format = 'J',array=THING_ID)
    table_hdu = pyfits.new_table([c1,c2,c3,c4,c5,c6,c7])

    print 'writting %s ...'%out_file
    h=pyfits.HDUList([pyfits.PrimaryHDU(delta_tilde)])
    h.append(pyfits.ImageHDU(ivar,name ="IVAR"))
    h.append(pyfits.ImageHDU(wave,name="WAVELENGTH"))     
    h.append(table_hdu)
    h.writeto(out_file,clobber=True)

print 'done' 
