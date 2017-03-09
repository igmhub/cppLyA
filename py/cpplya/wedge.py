import numpy as np
from math import sqrt
import sys

def block(covmat,indices) :
    res=np.zeros((indices.size,indices.size))
    for i in range(indices.size) :
        res[i,:]=covmat[indices[i],indices]
    return res

def compute_wedge(input_xi2d,input_cov,murange=[0.8,1.0],rrange=[10,180],rbin=4) : 
    
    # indexing
    n2d=input_xi2d.size
    n1d=np.sqrt(n2d).astype(int)
    rstep=4.
    rt=((np.arange(n2d)%n1d+0.5)*rstep).astype(float)
    rp=((np.arange(n2d)/n1d+0.5)*rstep).astype(float)
    rr=np.sqrt(rt**2+rp**2)
        
    rt_edges=np.zeros((rt.size,2,2))
    rp_edges=np.zeros((rp.size,2,2))
    for i in range(2) :
        for j in range(2) :
            rt_edges[:,i,j]=rt[:]-rstep/2.+i*rstep
            rp_edges[:,i,j]=rp[:]-rstep/2.+j*rstep
    rr_edges=np.sqrt(rt_edges**2+rp_edges**2)
    mu_edges=rp_edges/(rr_edges+(rr_edges==0))
    
    rr_min=np.min(np.min(rr_edges,axis=-1),axis=-1)
    rr_max=np.max(np.max(rr_edges,axis=-1),axis=-1)
    mu_min=np.min(np.min(mu_edges,axis=-1),axis=-1)
    mu_max=np.max(np.max(mu_edges,axis=-1),axis=-1)
    
    nr=int((rrange[1]-rrange[0])/rbin)
    r=rrange[0]+rbin/2+np.arange(nr)*rbin
    
    wedge_indices=np.where((mu_max>=murange[0])&(mu_min<=murange[1])&(rr_max>=rrange[0])&(rr_min<=rrange[1]))[0]
    wedge_data=input_xi2d[wedge_indices]
    rr=rr[wedge_indices]
    rt=rt[wedge_indices]
    rp=rp[wedge_indices]
    rr_edges=rr_edges[wedge_indices]
    mu_edges=mu_edges[wedge_indices]
    rr_min=rr_min[wedge_indices]
    rr_max=rr_max[wedge_indices]
    mu_min=mu_min[wedge_indices]
    mu_max=mu_max[wedge_indices]
        
    ndata=wedge_data.size
    wedge_cov=block(input_cov,wedge_indices)
    
    H=np.zeros((nr,ndata))
    for i in range(nr) :
        rmin=r[i]-rbin/2.
        rmax=r[i]+rbin/2.
        jj=np.where((mu_max>=murange[0])&(mu_min<=murange[1])&(rr_max>=rmin)&(rr_min<=rmax))[0]

        for j,jindex in zip(jj,np.arange(jj.size)) :
            # find fraction of each pixel in slice rmin,rmax,mu_min,mu_max with subsampling pixel
            n=7
            rtb=np.tile(np.linspace(rt[j]-rstep/2.+rstep/n/2,rt[j]+rstep/2.-rstep/n/2.,n),(n,1)).ravel()
            rpb=np.tile(np.linspace(rp[j]-rstep/2.+rstep/n/2,rp[j]+rstep/2.-rstep/n/2.,n),(n,1)).T.ravel()
            rrb=np.sqrt(rtb**2+rpb**2)
            mub=rpb/rrb
            frac=np.sum((mub>=murange[0])*(mub<=murange[1])*(rrb>=rmin)*(rrb<rmax))/float(n**2)
            H[i,j]=frac
        s=np.sum(H[i])
        if s>0 :
            H[i] /= s
    
    res=H.dot(wedge_data)
    cov=H.dot(wedge_cov.dot(H.T))
    return r,res,np.sqrt(np.diag(cov).copy()),cov



def compute_wedge_with_ivar(input_xi2d,input_cov,murange=[0.8,1.0],rrange=[10,180],rbin=4) : 
    
    # indexing
    n2d=input_xi2d.size
    n1d=np.sqrt(n2d).astype(int)
    rstep=4.
    rt=((np.arange(n2d)%n1d+0.5)*rstep).astype(float)
    rp=((np.arange(n2d)/n1d+0.5)*rstep).astype(float)
    rr=np.sqrt(rt**2+rp**2)
        
    rt_edges=np.zeros((rt.size,2,2))
    rp_edges=np.zeros((rp.size,2,2))
    for i in range(2) :
        for j in range(2) :
            rt_edges[:,i,j]=rt[:]-rstep/2.+i*rstep
            rp_edges[:,i,j]=rp[:]-rstep/2.+j*rstep
    rr_edges=np.sqrt(rt_edges**2+rp_edges**2)
    mu_edges=rp_edges/(rr_edges+(rr_edges==0))
    
    rr_min=np.min(np.min(rr_edges,axis=-1),axis=-1)
    rr_max=np.max(np.max(rr_edges,axis=-1),axis=-1)
    mu_min=np.min(np.min(mu_edges,axis=-1),axis=-1)
    mu_max=np.max(np.max(mu_edges,axis=-1),axis=-1)
    
    nr=int((rrange[1]-rrange[0])/rbin)
    r=rrange[0]+rbin/2+np.arange(nr)*rbin
    
    wedge_indices=np.where((mu_max>=murange[0])&(mu_min<=murange[1])&(rr_max>=rrange[0])&(rr_min<=rrange[1]))[0]
    wedge_data=input_xi2d[wedge_indices]
    rr=rr[wedge_indices]
    rt=rt[wedge_indices]
    rp=rp[wedge_indices]
    rr_edges=rr_edges[wedge_indices]
    mu_edges=mu_edges[wedge_indices]
    rr_min=rr_min[wedge_indices]
    rr_max=rr_max[wedge_indices]
    mu_min=mu_min[wedge_indices]
    mu_max=mu_max[wedge_indices]
        
    ndata=wedge_data.size
    wedge_cov=block(input_cov,wedge_indices)
    
    var=np.diag(wedge_cov)

    H=np.zeros((nr,ndata))
    for i in range(nr) :
        rmin=r[i]-rbin/2.
        rmax=r[i]+rbin/2.
        jj=np.where((mu_max>=murange[0])&(mu_min<=murange[1])&(rr_max>=rmin)&(rr_min<=rmax))[0]

        for j,jindex in zip(jj,np.arange(jj.size)) :
            # find fraction of each pixel in slice rmin,rmax,mu_min,mu_max with subsampling pixel
            n=7
            rtb=np.tile(np.linspace(rt[j]-rstep/2.+rstep/n/2,rt[j]+rstep/2.-rstep/n/2.,n),(n,1)).ravel()
            rpb=np.tile(np.linspace(rp[j]-rstep/2.+rstep/n/2,rp[j]+rstep/2.-rstep/n/2.,n),(n,1)).T.ravel()
            rrb=np.sqrt(rtb**2+rpb**2)
            mub=rpb/rrb
	    w = (mub>=murange[0]) * (mub<=murange[1]) * (rrb>=rmin) & (rrb<rmax)
            frac=np.sum(w)/var[j]
            H[i,j]=frac
        s=np.sum(H[i])
        if s>0 :
            H[i] /= s
    
    res=H.dot(wedge_data)
    cov=H.dot(wedge_cov.dot(H.T))
    return r,res,np.sqrt(np.diag(cov).copy()),cov
