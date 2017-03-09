#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import string
import pyfits 

from lya_io import read_baofit_data
from lya_io import read_baofit_cov
from lya_io import read_baofit_model
from wedge import compute_wedge
from wedge import compute_wedge_with_ivar

plt.rcParams["font.family"]="serif"
plt.rcParams["font.size"]=16.0


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-d','--data', type = str, default = None, required=True, nargs='*', action='append',
                        help = 'baofit data')
parser.add_argument('--rpower', type = int, default = 2, required=False,
                        help = 'rpower')
parser.add_argument('--mu', type = str, default = None, required=False,
                        help = 'mu range for wedge of the form "mumin:mumax')
parser.add_argument('--rrange', type = str, default = "10:180", required=False,
                        help = 'r range for wedge of the form "rmin:rmax')
parser.add_argument('--rbin', type = float, default = 4.0, required=False,
                        help = 'r bin size')
parser.add_argument('--res', type = str, default = None, required=False, nargs='*', action='append', 
                    help = 'baofit residuals file to plot model')
parser.add_argument('--out', type = str, default = None, required=False,
		                        help = 'output prefix')
parser.add_argument('--chi2', action="store_true",
		help = 'compute chi2 of wedges (if only one data set and model)')

parser.add_argument('--noshow', action="store_true",
		            help = 'prevent the figure window from displaying')

parser.add_argument('--ivar_weight', action="store_true",
		            help = 'use inverse variance to combine the bins')

parser.add_argument('--same_cov', type = str, default = None, required=False,
		            help = 'use this covariance matrix for all the data sets')
args = parser.parse_args()

data=[]
cov=[]
title=[]
for d2 in args.data :
    for d1 in d2 :
        title.append(d1)
        data.append(read_baofit_data(d1))
        #cov_filename=d1+".cov"
        #cov.append(read_baofit_cov(cov_filename,n2d=data[0].size,convert=True))
        if args.same_cov: 
            tmp_cov = pyfits.open(args.same_cov)      
        else: 
            cov_filename=string.replace(d1,".data","-cov.fits")
            tmp_cov = pyfits.open(cov_filename)
        cov.append(tmp_cov[0].data)

models=[]
if args.res is not None :
    for res2 in args.res :
        for res in res2 :
            mod = read_baofit_model(res,n2d=data[0].size)
            if data[0].size != mod.size :
                print "error data and model don't have same size"
                sys.exit(12)
                models.append(mod)

rpower=2
if args.rpower is not None: 
    rpower = args.rpower

if args.mu :
    try :
        vals=string.split(args.mu,":")
        if len(vals)!=2 :
            print "incorrect format for mu range '%s', expect mumin:mumax"%args.mu
            sys.exit(12)
        mumin=string.atof(vals[0])
        mumax=string.atof(vals[1])
        wedges=[[mumin,mumax]]
    except ValueError,e:
        print e
        print "incorrect format for mu range '%s', expect mumin:mumax"%args.mu
        sys.exit(12)
else :
    wedges= [[0.8,1.0],[0.5,0.8],[0.0,0.5]]

if args.rrange :
    try :
        vals=string.split(args.rrange,":")
        if len(vals)!=2 :
            print "incorrect format for r range '%s', expect rmin:rmax"%args.rrange
            sys.exit(12)
        rmin=string.atof(vals[0])
        rmax=string.atof(vals[1])
        rrange=[rmin,rmax]
    except ValueError,e:
        print e
        print "incorrect format for r range '%s', expect rmin:rmax"%args.rrange
        sys.exit(12)
else :
    rrange=[10,180]


nw=len(wedges)
if nw > 1 :
    f, ax = plt.subplots(nw, sharex=True, sharey=False)
else :
    f=plt.figure()
    ax=[]
    ax.append(plt.subplot(1,1,1))

ax[nw/2].set_ylabel(r"$r^{%i} \xi(r)\mathrm{[h^{-2}Mpc^2]}$"%rpower)

data_colors=["b","r","g","k"]
model_colors=["r","k","k","k"]

for w,wedge in zip(range(nw),wedges) :
    print "plotting mu",wedge
    for d,C,c,t in zip(data,cov,data_colors,title) :
	if not args.ivar_weight: r,xidata,xierr,wedge_cov=compute_wedge(d,C,murange=wedge,rrange=rrange,rbin=args.rbin)      
	else: r,xidata,xierr,wedge_cov=compute_wedge_with_ivar(d,C,murange=wedge,rrange=rrange,rbin=args.rbin) 
        ax[w].errorbar(r,r**rpower*xidata,r**rpower*xierr,fmt="o",color=c,label=t)
	ax[w].grid(b=True)
    
    for model,C,c in zip(models,cov,model_colors)  :
	if not args.ivar_weight: r,ximod,junk,junk=compute_wedge(model,C,murange=wedge,rrange=rrange,rbin=args.rbin)
	else: r,ximod,junk,junk=compute_wedge_with_ivar(model,C,murange=wedge,rrange=rrange,rbin=args.rbin)
        ax[w].plot(r,r**rpower*ximod,"-",color=c)
    
    if args.chi2 and len(data)==1 and len(models)==1 :
        weight=np.linalg.inv(wedge_cov)
        res=xidata-ximod
        chi2=np.inner(res,weight.dot(res))
        ndata=res.size
        print "chi2/ndata=%f/%d=%f"%(chi2,ndata,chi2/ndata)
    
    ax[w].set_title(r"$%2.2f < \mu < %2.2f$"%(wedges[w][0],wedges[w][1]))
plt.legend(loc=3)
plt.xlabel(r"$r\mathrm{[h^{-1}Mpc]}$")
if not args.noshow: 
    plt.show()

if args.out != None:
	f.savefig(args.out+".png",bbox_inches="tight")

