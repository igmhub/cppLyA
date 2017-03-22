#!/usr/bin/env python 
import pyfits 
import numpy as np 
import sys 
import argparse
import multiprocessing
import pylab 

def read_infile(infile): 
    h=pyfits.open(infile)
    dflux=h[0].data 
    ivar=h[1].data 
    wave=h[2].data 
    plate=h[3].data.PLATE
    return dflux, ivar, wave, plate 

def fill_histograms(spec): 
    hww=np.zeros((nbins))
    hwwdd=np.zeros((nbins))
    wave_selection1=np.where(ivar1[spec]>0)[0]
    if (wave_selection1.size != 0):
        for w1 in wave_selection1:
            RA=wave[w1]/wave
            selec=np.where((RA>Rmin)*(RA<Rmax)*(ivar2[spec]>0))[0]
            ind=((RA[selec]-Rmin)/rstep).astype(int)
            hww[ind]+=ivar1[spec][w1]*ivar2[spec][selec]
            hwwdd[ind]+=ivardelta1[spec][w1]*ivardelta2[spec][selec]
    return [hww,hwwdd]

def save(): 
    print "writing %s"%out_file
    corrfunc = pyfits.HDUList([pyfits.PrimaryHDU(XI)])
    corrfunc.append(pyfits.ImageHDU(Rbin, name = "R"))
    corrfunc.writeto(out_file,clobber=True)

######################## Parser ##############################

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f1','--input1', type = str, default = None, required=True,
                        help = 'path of input fits file 1')
parser.add_argument('-f2','--input2', type = str, default = None, required=True,
                        help = 'path of input fits file 2')
parser.add_argument('-o','--output', type = str, default = None, required=True,
                        help = 'path of correlation fits file')
parser.add_argument('--w1', type = str, default = False, required=False,
                        help = 'type of absorption in the 1st forest, by default Lya absorption')
parser.add_argument('--w2', type = str, default = False, required=False,
                        help = 'type of absorption in the 2nd forest, by default Lya absorption')
parser.add_argument('--ncpu', type = int, default = 1, required=False,
                        help = 'number of CPU for parallelization')
parser.add_argument('--nspec', type = int, default = 0, required=False,
                        help = 'max number of QSO spectra')
parser.add_argument('--nbins', type = int, default = 0, required=False,
                        help = 'number of bins')
parser.add_argument('--Rmin', type = float, default = 0, required=False,
                        help = 'minimum of wave1/wave2')
parser.add_argument('--Rmax', type = float, default = 0, required=False,
                        help = 'maximum of wave1/wave2')
args = parser.parse_args()

########################

in_file1=args.input1
in_file2=args.input2
out_file=args.output 

if args.w1: 
    forest1=args.w1
else: 
    forest1='lya'

if args.w2: 
    forest2=args.w2
else: 
    forest2='lya'

if args.ncpu:
    ncpu=args.ncpu
else: 
    ncpu=12

print 'reading file ',in_file1
dflux1,ivar1,wave1,plate1=read_infile(in_file1)
if in_file1 == in_file2: 
    dflux2=dflux1
    ivar2=ivar1
    wave2=wave1
    plate2=plate1
else: 
    print 'reading file ',in_file2
    dflux2,ivar2,wave2,plate2=read_infile(in_file2)
    if (np.abs(wave1[0]-wave2[0])>0.001) : 
        print 'cannot handdle different grids'
        sys.exit(12)

wave=wave1
nw = dflux1.shape[1]
if args.nspec: 
    nspec=args.nspec
else: 
    nspec=dflux1.shape[0]
print 'nspec = %i , nw = %i'%(nspec,nw)


if nspec != dflux1.shape[0]: 
    dflux1=dflux1[:nspec]
    ivar1=ivar1[:nspec]
    dflux2=dflux2[:nspec]
    ivar2=ivar2[:nspec]

if args.nbins:
    nbins=args.nbins
else: 
    nbins = 100
print 'nbins = ',nbins 

if args.Rmin: 
    Rmin=args.Rmin 
else: 
    Rmin=1.

if args.Rmax: 
    Rmax=args.Rmax
else: 
    Rmax=1.07
print 'Rmin = %2.2f , Rmax = %2.2f'%(Rmin,Rmax)

Rbin=np.linspace(Rmin,Rmax,nbins)
rstep=Rbin[1]-Rbin[0]
print 'rstep = ',rstep 

print 'preparation of ivar*delta ...'
ivardelta1=ivar1*dflux1
ivardelta2=ivar2*dflux2



print 'filling histograms...'
pool = multiprocessing.Pool(ncpu)
results = pool.map(fill_histograms, np.arange(nspec))



#res(spec,i,bin)
h_sum_w=np.zeros((nbins))
h_sum_wdd=np.zeros((nbins))
for spec in range(nspec): 
    h_sum_w+=results[spec][0]
    h_sum_wdd+=results[spec][1]

XI  = h_sum_wdd*(h_sum_w != 0)/(h_sum_w + (h_sum_w==0) )
print "done with computation" 

save()
print 'done' 

pylab.plot(Rbin,XI)
pylab.show()


sys.exit(12)



