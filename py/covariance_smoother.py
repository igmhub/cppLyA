#!/usr/bin/env python


import sys
import numpy as np
import pyfits
from math import sqrt
from math import factorial
import scipy.linalg

if len(sys.argv)<2 :
    print sys.argv[0],"input-covmat.fits (output-covmat.fits else overwrite)"
    sys.exit(12)

input_covmat_filename = sys.argv[1]
if len(sys.argv)>2 :
    output_covmat_filename = sys.argv[2]
else :
    output_covmat_filename = input_covmat_filename


# read
covmat=pyfits.open(input_covmat_filename)[0].data

# check if has nan
if np.isnan(np.sum(covmat)) :
    print "matrix has nan !!!"
    sys.exit(12)



print "estimate of variance of empirical covariance ..."
# the variance of an empirical variance based on N obs is var(var) = 2/N*var**2
#   IF we neglect the terms due to the covariance
nobs=2396 # number of plates
var=np.diag(covmat).copy()
var2=np.outer(var,var).ravel()
var_of_covmat=2*var2/nobs

print "set dummy values to null variances (can happen on rp,rt edges not used anyway) ..."
mvar=np.mean(var[var>0])
var[var==0]=mvar

print "indexing ..."
n2d=covmat.shape[0]
n1d=int(sqrt(n2d))
print "n2d=",n2d
print "n1d=",n1d

indices=np.arange(n2d)
rt=(indices%n1d)
rp=(indices/n1d)
rti=np.tile(rt,(n2d,1))
rtj=rti.T
rpi=np.tile(rp,(n2d,1))
rpj=rpi.T
drt=np.abs(rti-rtj)
drp=np.abs(rpi-rpj)
mrt=(rti+rtj)/2.
mrp=(rpi+rpj)/2.

#print "ndrt=",np.unique(drt).size,"drtmin=",np.min(drt),"drtmax=",np.max(drt)
#rr=np.sqrt(rt**2+rp**2)

print "compute correlation matrix ..."
corr=np.zeros(covmat.shape)
for i in range(n2d) :
    corr[i,:]=covmat[i,:]/(sqrt(var[i])*np.sqrt(var[:]))

print "compute mean correlation as a function of drt and drp ..."
sum_corr=np.zeros((n2d))

indices=drt.astype(int)+n1d*drp.astype(int)
bins=np.arange(n2d+1)-0.5
sum_n,junk    = np.histogram(indices,bins=bins)
sum_corr,junk = np.histogram(indices,bins=bins,weights=corr)
sum_corr2,junk = np.histogram(indices,bins=bins,weights=corr**2)
mean_corr = sum_corr/(sum_n+(sum_n==0))

print "compute smooth cov ..."
scov=np.zeros(covmat.shape)
for i in range(n2d) :
    scov[i,:]=(sqrt(var[i])*np.sqrt(var[:]))*mean_corr[indices[i,:]]

print "use original cov for drt=0 and drp=1 and mrp=0.5"
ok=np.where((drt.ravel()==0)&(drp.ravel()==1)&(mrp.ravel()<1))[0]
print "ok.size=",ok.size
scov=scov.ravel()
scov[ok]=covmat.ravel()[ok]
scov=scov.reshape(covmat.shape)



print "rm cov for r<10Mpc and r>180Mpc"
ri=np.sqrt( ((rti+0.5)*4.)**2 + ((rpi+0.5)*4.)**2 )
rj=np.sqrt( ((rtj+0.5)*4.)**2 + ((rpj+0.5)*4.)**2 )
mask=( ((ri<10.)|(ri>180.)|(rj<10.)|(rj>180.)) & ((drt>0)|(drp>0)) )
tscov=scov.copy()
tscov[mask]=0.

print "check definite positive or reduce cov ..."
while True :
    try :
        scipy.linalg.cho_factor(tscov, lower=True, overwrite_a=False)
    except :
        print "not pos def"
        print sys.exc_info()
        print "lower cov ..."
        j=np.arange(n2d)
        for i in range(n2d) :
           scov[i,j!=i] *= 0.95
           tscov[i,j!=i] *= 0.95
           
    break

print "save ..."
pyfits.writeto(output_covmat_filename,scov,clobber=True)
print "wrote",output_covmat_filename
