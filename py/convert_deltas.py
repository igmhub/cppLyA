#!/usr/bin/env python
import glob
import numpy as np
import sys
import fitsio
import math
import pylab
import argparse
import scipy as sp
import astropy.io.fits as pyfits

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-d','--indir', type = str, default = None, required=True,
                        help = 'deltas input directory')
parser.add_argument('-o', '--out',type = str, default = None, required=True,
                        help = 'output fits file ')
parser.add_argument('--nspec', type=int, required=False,
                    help = 'number of spec')
parser.add_argument('--substract', action="store_true", required=False,
                    help = 'substract first the mean delta as a function of wavelength, then substract a and b*dlwave terme per quasar')

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

lwave_list=[]
wave_list=[]
data_list=[]
ivar_list=[]
id_list=[]


for i in range(nfile):
    if (args.nspec):
        if (nspec>args.nspec): break
    print "%2.2f %%, nspec  = %i"%(float(i)/float(nfile)*100,nspec)
    sys.stdout.flush()
    hdu=fitsio.FITS(plate_name[i])
    for d in hdu[1:]:
        nspec+=1
        lwave_spec=d["LOGLAM"][:]
        wave_spec=10**d["LOGLAM"][:]
        ivar_spec=d["WEIGHT"][:]
        data_spec=d["DELTA"][:]
        data_spec*=(ivar_spec>0) #Avoid interpolation problem 
        id_spec=d.read_header()['THING_ID']

        lwave_list.append(lwave_spec)
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

if args.substract:
    print 'substract first the mean delta as a function of wavelength, then substract a and b*dlwave terme per quasar'
    for spec in range(nq):     
	mde = sp.average(data_list[spec],weights=ivar_list[spec])
        mll = sp.average(lwave_list[spec],weights=ivar_list[spec])
        mld = sp.sum(ivar_list[spec]*data_list[spec]*(lwave_list[spec]-mll))/sp.sum(ivar_list[spec]*(lwave_list[spec]-mll)**2)
        data_list[spec] -= mde + mld * (lwave_list[spec]-mll)
    print 'done with the substraction'

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

print 'lwmin = %f , lwmax = %f'%(10**lwmin,10**lwmax)
print "meas. lwstep",lwstep
lwstep=0.0003
print "true lwstep",lwstep
nw=int((lwmax-lwmin)/lwstep)+2
print 'nw = ',nw
print
wave=10**(np.arange(nw)*lwstep+lwmin)
print 'np.log10(wave[-1])-lwmax = ',np.log10(wave[-1])-lwmax
print
data = np.zeros((nq,wave.size))
ivar = np.zeros((nq,wave.size))
print 'data.shape = ',data.shape

"""
for q in range(nq) :
    i0=int((np.log10(wave_list[q][0])-lwmin+0.001*lwstep)/lwstep)
    i1=i0+data_list[q].size
    data[q,i0:i1] = data_list[q]
    ivar[q,i0:i1] = ivar_list[q]
"""
for q in range(nq) :
    i=(( np.log10(wave_list[q])-lwmin+0.001*lwstep) /lwstep).astype(int)
    data[q,i] = data_list[q]
    ivar[q,i] = ivar_list[q]

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

c1=pyfits.Column(name='RA',format = 'D',array=RA)
c2=pyfits.Column(name = 'DEC', format = 'D',array=DEC)
c3=pyfits.Column(name = 'Z', format = 'D',array=Z)
c4=pyfits.Column(name = 'PLATE', format = 'J',array=PLATE)
c5=pyfits.Column(name = 'MJD', format = 'J',array=MJD)
c6=pyfits.Column(name = 'FIBER', format = 'J',array=FIBER)
c7=pyfits.Column(name = 'THING_ID', format = 'J',array=THING_ID)
cols=pyfits.ColDefs([c1,c2,c3,c4,c5,c6,c7])
table_hdu = pyfits.BinTableHDU.from_columns(cols)
print 'writting %s ...'%out_file
h=pyfits.HDUList([pyfits.PrimaryHDU(data)])
h.append(pyfits.ImageHDU(ivar,name ="IVAR"))
h.append(pyfits.ImageHDU(wave,name="WAVELENGTH"))
h.append(table_hdu)
h.writeto(out_file,clobber=True)

print 'done'

        
