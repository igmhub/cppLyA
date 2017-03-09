import numpy as np
import astropy.io.fits as fits
import string
import os.path
import sys
from math import sqrt

def read_baofit_data(filename) :
    if filename.find(".data")<0 :
        data_filename="%s.data"%filename
    else :
        data_filename=filename
    print "reading data in %s"%data_filename

    if not os.path.isfile(data_filename) :
        print "error %s doesn't exist"%data_filename
        sys.exit(12)

    vals=np.loadtxt(data_filename).T
    return vals[1]

def read_baofit_data_sw(filename) :
    if filename.find(".data")<0 :
        data_filename="%s.data"%filename
    else :
        data_filename=filename
    print "reading data in %s"%data_filename

    if not os.path.isfile(data_filename) :
        print "error %s doesn't exist"%data_filename
        sys.exit(12)

    vals=np.loadtxt(data_filename).T
    return vals[1],vals[2]

def read_baofit_model(res_filename,n2d=2500) :
    if not os.path.isfile(res_filename) :
        print "error %s doesn't exist"%res_filename
        sys.exit(12)

    vals=np.loadtxt(res_filename).T
    i=vals[0].astype(int)
    ni=np.max(i)+1
    if ni>n2d :
        n2d=ni
    m=vals[7]
    res=np.zeros((n2d))
    res[i]=m
    return res


def read_baofit_cov(filename,n2d,convert=True) :
    print "reading cov in %s"%filename
    
    if filename.find(".fits")>0 :
        cov = fits.open(filename)[0].data
        if cov.shape != (n2d,n2d) :
            print "error, incorrect matrix size"
            sys.exit(12)
    
    if filename.find(".cov")>0 :
        cov_filename=filename
        base_filename=string.replace(filename,".cov","")
        fits_filename="%s-cov.fits"%base_filename
    else :
        base_filename=filename
        cov_filename="%s.cov"%base_filename
        fits_filename="%s-cov.fits"%base_filename
    
    if os.path.isfile(fits_filename) :
        # check date
        date_of_cov_filename  = os.path.getmtime(cov_filename)
        date_of_fits_filename =  os.path.getmtime(fits_filename)
        if date_of_fits_filename > date_of_cov_filename :
            print "using %s"%fits_filename
            return fits.open(fits_filename)[0].data
        else :
            print "%s exists, but use %s which is more recent"%(fits_filename,cov_filename)
    
    if not os.path.isfile(cov_filename) :
        print "warning, %s doesn't exist, returns null matrix"%cov_filename
        return np.zeros((n2d,n2d))

    print "reading %s"%cov_filename
    vals=np.loadtxt(cov_filename).T
    
    ii=vals[0].astype(int)
    jj=vals[1].astype(int)
    cc=vals[2]
    
    n=max(np.max(ii),np.max(jj))+1
    cov=np.zeros((n,n))
    if n != n2d :
        print "error, incorrect matrix size"
        sys.exit(12)
    
    for i,j,c in zip(ii,jj,cc) :
        cov[i,j]=c
        cov[j,i]=c
    
    if convert :
        fits.writeto(fits_filename,cov,clobber=True)
        print "wrote %s"%fits_filename
    
    return cov

