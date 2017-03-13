#include <cmath>
#include <map>
#include <iostream>
#include <qso_spectrum.h>
#include <fiducial_cosmo.h>
#include <cosmology.h>
#include <matvect.h>
#include <correlation_tools.h>
#include <histo2d.h>
#include <fileutils.h>

using namespace std;

Vect compute_lya_distances(const Vect& wave,const string& line) {
  cout << "computing distances assuming it is " << line << endl;
  cout << " cosmo" << endl;

  GeneralCosmo cosmo(FIDUCIAL_OM,1.-FIDUCIAL_OM,-1.,0.);

  double rest_lambda=0;
  if(line=="lya")
    rest_lambda=lya_lambda;
  else if(line=="lyb")
    rest_lambda=lyb_lambda;
  else {
    cout << "don't know " << line << endl;
    exit(12);
  }
  double minz=wave(0)/rest_lambda-1;
  double maxz=wave(wave.size()-1)/rest_lambda-1;
  EqualStepBSpline1D sd;
  {
    int nz=100;
    double zbin=(maxz-minz)/(nz-1);
    Vect z(nz);
    Vect d(nz);
    for(int i=0;i<nz;i++) {
      z(i)=minz+zbin*i;
      d(i)=CSPEED/H0*(cosmo.Dm(z(i)));
    }
    sd.Fit(nz,z.Data(),d.Data());
  }
  Vect dist(wave.size());
  for(size_t i=0;i<wave.size();i++) {
    dist(i)=sd.Value(wave(i)/rest_lambda-1);
  }
  return dist;
}


void compute_qso_distances(vector<QSO>& qsos) {
  cout << "computing distances" << endl;
  cout << " cosmo" << endl;

  GeneralCosmo cosmo(FIDUCIAL_OM,1.-FIDUCIAL_OM,-1.,0.);

  double minz=1000.;
  double maxz=0;
  for(size_t k=0;k<qsos.size();k++) {
    minz=min(minz,qsos[k].z);
    maxz=max(maxz,qsos[k].z);
  }
  minz=max(minz,0.);
  
  EqualStepBSpline1D sd;
  {
    int nz=100;
    double zbin=(maxz-minz)/(nz-1);
    Vect z(nz);
    Vect d(nz);
    for(int i=0;i<nz;i++) {
      z(i)=minz+zbin*i;
      d(i)=CSPEED/H0*(cosmo.Dm(z(i)));
    }
    sd.Fit(nz,z.Data(),d.Data());
  }
  for(size_t k=0;k<qsos.size();k++) {
    if(qsos[k].z>0) 
      qsos[k].dist=sd.Value(qsos[k].z);
    else 
      qsos[k].dist=0.;
  }
}

void compute_qso_directions(vector<Spectrum>& spectra) {
  cout << "computing direction vectors ..." << endl;
  for(size_t k=0;k<spectra.size();k++) {
    if(!spectra[k].valid) continue;
    spectra[k].set_dir();
    if(k%10000==0) cout << "  done " << k << "/" << spectra.size() << endl;
  }  
}

void compute_qso_directions(vector<QSO>& qsos) {
  cout << "computing direction vectors ..." << endl;
  for(size_t k=0;k<qsos.size();k++) {    
    qsos[k].set_dir();
    if(k%10000==0) cout << "  done " << k << "/" << qsos.size() << endl;
  }  
}

void check_nan(const vector<Spectrum>& spectra) {
  cout << "look for nan ..." << endl;
  for(size_t k=0;k<spectra.size();k++) {
    const Spectrum& spec_k = spectra[k];
    const double *flux=spec_k.flux.Data();
    const double *weight=spec_k.weight.Data();
    for(size_t i=0;i<spec_k.flux.size();i++) {
      if(isnan(flux[i])) {
	cout << "NaN in spec " << k << " flux " << i << endl;
	exit(12);
      }
      if(isnan(weight[i])) {
	cout << "NaN in spec " << k << " flux " << i << endl;
	exit(12);
	}
    }
    if(k%10000==0) cout << "  done " << k << "/" << spectra.size() << endl;
  }
}
void compute_valid_range(const std::vector<Spectrum>& spectra,int *& wbegin,int *& wend) {
  cout << "computing valid ranges ..." << endl;
  wbegin=new int[spectra.size()];
  wend=new int[spectra.size()];
  // find valid ranges (to go faster)
  for(size_t k=0;k<spectra.size();k++) {
    
    const Spectrum& spec_k = spectra[k];
    if(!spec_k.valid) continue;

    int nw=spec_k.flux.size();
    {
      int& b=wbegin[k];
      wbegin[k]=0;
      while(spec_k.weight(b)==0) {
	b++;
	if(b==nw) break;
      };
    }
    {
      int& b=wend[k];
      wend[k]=nw;
      while(spec_k.weight(b-1)==0) {
	b--;
	if(b==0) break;
      };
    }
  }
}

void multiply_flux_by_weight(std::vector<Spectrum>& spectra) {
  cout << "multiply flux by weight ..." << endl;
  for(size_t k=0;k<spectra.size();k++) {
    Spectrum& spec_k = spectra[k];
    if(!spec_k.valid) continue;
    double *flux=spec_k.flux.NonConstData();
    const double *weight=spec_k.weight.Data();
    for(size_t i=0;i<spec_k.flux.size();i++)
      flux[i] *= weight[i];
    if(k%10000==0) cout << "  done " << k << "/" << spectra.size() << endl;
  }
}

void save_lya_correlation_results_fast(const string& filename, const map<int,Histo2d*>& h_sum_wdd, const map<int,Histo2d*>& h_sum_w,bool with_covmat) {
  cout << "computing results" << endl;

  map<int,Histo2d*>::const_iterator h_sum_wdd_it = h_sum_wdd.begin();
  map<int,Histo2d*>::const_iterator h_sum_w_it = h_sum_w.begin();

  int nplates = h_sum_wdd.size();
  int n2d = h_sum_wdd_it->second->nx * h_sum_wdd_it->second->ny;
  int npair = 0; 

  double swx[n2d];
  double sw[n2d]; 
  for(int i=0;i<n2d;i++) {
    swx[i]=0;
    sw[i]=0;
  }
  
  for (;h_sum_wdd_it != h_sum_wdd.end(); h_sum_wdd_it ++, h_sum_w_it ++) {
    const double* w  = h_sum_w_it->second->data;
    const double* wx = h_sum_wdd_it->second->data;
    for(int i=0;i<n2d;i++) {
      swx[i] += wx[i];
      sw[i]  += w[i];
      npair+=1; 
    }
  }
  
  double mx[n2d];
  for(int i=0;i<n2d;i++) {
    if(sw[i]>0)
      mx[i] = swx[i]/sw[i];
    else 
      mx[i] = 0;
  }
  cout << "writing " << filename << endl;
  FILE *file = fopen(filename.c_str(),"w");
  
  for(int i=0;i<n2d;i++) {

    fprintf(file,"%d %g %g \n",i ,mx[i], sw[i]);
  }
  fclose(file);
  
  if (with_covmat) {
    cout << "computing covmat ..." << endl;
  
    Mat cov(n2d,n2d);
    double wrx[n2d];
    double* cov_data = cov.NonConstData();
  
    h_sum_wdd_it = h_sum_wdd.begin();
    h_sum_w_it = h_sum_w.begin();
    int p=0;
    for (;h_sum_wdd_it != h_sum_wdd.end(); h_sum_wdd_it ++, h_sum_w_it ++, p++) {
      const double* w  = h_sum_w_it->second->data;
      const double* wx = h_sum_wdd_it->second->data;
      for(int i=0;i<n2d;i++) {
	wrx[i] = wx[i] - w[i]*mx[i];
	if(wrx[i]==0) continue;
	for(int j=0;j<=i;j++) {
	  cov_data[i+n2d*j] += (wrx[i]*wrx[j]);
	}
      }
      if(p%100==0)
	cout << "plate " << p << "/" << nplates << endl;
    }
    // normalize
    for(int i=0;i<n2d;i++) {
      for(int j=0;j<=i;j++) {
	if(sw[i]*sw[j]>0)
	  cov_data[i+n2d*j] /= (sw[i]*sw[j]);
      }
    }
    // symmetrize
    for(int i=0;i<n2d;i++) {
      for(int j=0;j<=i;j++) {
	cov_data[j+n2d*i] = cov_data[i+n2d*j];
      }
    }
      
    string covmat_baofit_ascii = CutExtension(filename)+".cov";
    string covmat_fits         = CutExtension(filename)+"-cov.fits";
    
    /*
    cout << "writing " << covmat_baofit_ascii << endl;
    file = fopen(covmat_baofit_ascii.c_str(),"w");
    for(int i=0;i<n2d;i++) {
      for(int j=0;j<n2d;j++) {
      fprintf(file,"%d %d %g\n",i,j,cov_data[i+n2d*j]);
      }
    }
    fclose(file);
    */
    
    cout << "writing " << covmat_fits << endl;
    cov.writeFits(covmat_fits);
    
  }
  cout << "npair = " << npair << endl; 
  cout << "done" << endl;
}

void save_lya_correlation_results(const string& filename, const map<int,Histo2d*>& h_sum_wdd, const map<int,Histo2d*>& h_sum_w,const map<int,Histo2d*>& h_sum_z, const map<int,Histo2d*>& h_sum_rp, const map<int,Histo2d*>& h_sum_rt, bool with_covmat) {
  cout << "computing results" << endl;

  map<int,Histo2d*>::const_iterator h_sum_wdd_it = h_sum_wdd.begin();
  map<int,Histo2d*>::const_iterator h_sum_w_it = h_sum_w.begin();
  map<int,Histo2d*>::const_iterator h_sum_z_it = h_sum_z.begin();
  map<int,Histo2d*>::const_iterator h_sum_rp_it = h_sum_rp.begin();
  map<int,Histo2d*>::const_iterator h_sum_rt_it = h_sum_rt.begin();

  int nplates = h_sum_wdd.size();
  int n2d = h_sum_wdd_it->second->nx * h_sum_wdd_it->second->ny;

  double swx[n2d];
  double sw[n2d]; 
  double sz[n2d];
  double srp[n2d]; 
  double srt[n2d];


  for(int i=0;i<n2d;i++) {
    swx[i]=0;
    sw[i]=0;
    sz[i]=0;
    srp[i]=0;
    srt[i]=0;
  }
  
  for (;h_sum_wdd_it != h_sum_wdd.end(); h_sum_wdd_it ++, h_sum_w_it ++,h_sum_z_it ++, h_sum_rp_it ++,h_sum_rt_it ++) {
    const double* w  = h_sum_w_it->second->data;
    const double* wx = h_sum_wdd_it->second->data;
    const double* z  = h_sum_z_it->second->data;
    const double* rp = h_sum_rp_it->second->data;
    const double* rt = h_sum_rt_it->second->data;

    for(int i=0;i<n2d;i++) {
      swx[i] += wx[i];
      sw[i]  += w[i];
      sz[i]  += z[i];
      srp[i] += rp[i];
      srt[i] += rt[i];
    }
  }
  
  double mx[n2d];
  double mz[n2d];
  double mrp[n2d];
  double mrt[n2d];
  for(int i=0;i<n2d;i++) {
    if(sw[i]>0){
      mx[i]  = swx[i]/sw[i];
      mz[i]  = sz[i]/sw[i];
      mrp[i] = srp[i]/sw[i];
      mrt[i] = srt[i]/sw[i];}
    else {
      mx[i]  = 0;
      mz[i]  = 0;
      mrp[i] = 0;
      mrt[i] = 0;}
  }
  cout << "writing " << filename << endl;
  FILE *file = fopen(filename.c_str(),"w");
  
  for(int i=0;i<n2d;i++) {

    fprintf(file,"%d %g %g %g %g %g \n",i ,mx[i], sw[i], mz[i], mrp[i], mrt[i]);
  }
  fclose(file);
  
  if (with_covmat) {
    cout << "computing covmat ..." << endl;
  
    Mat cov(n2d,n2d);
    double wrx[n2d];
    double* cov_data = cov.NonConstData();
  
    h_sum_wdd_it = h_sum_wdd.begin();
    h_sum_w_it = h_sum_w.begin();
    int p=0;
    for (;h_sum_wdd_it != h_sum_wdd.end(); h_sum_wdd_it ++, h_sum_w_it ++, p++) {
      const double* w  = h_sum_w_it->second->data;
      const double* wx = h_sum_wdd_it->second->data;
      for(int i=0;i<n2d;i++) {
	wrx[i] = wx[i] - w[i]*mx[i];
	if(wrx[i]==0) continue;
	for(int j=0;j<=i;j++) {
	  cov_data[i+n2d*j] += (wrx[i]*wrx[j]);
	}
      }
      if(p%100==0)
	cout << "plate " << p << "/" << nplates << endl;
    }
    // normalize
    for(int i=0;i<n2d;i++) {
      for(int j=0;j<=i;j++) {
	if(sw[i]*sw[j]>0)
	  cov_data[i+n2d*j] /= (sw[i]*sw[j]);
      }
    }
    // symmetrize
    for(int i=0;i<n2d;i++) {
      for(int j=0;j<=i;j++) {
	cov_data[j+n2d*i] = cov_data[i+n2d*j];
      }
    }
      
    string covmat_baofit_ascii = CutExtension(filename)+".cov";
    string covmat_fits         = CutExtension(filename)+"-cov.fits";
    
    /*
    cout << "writing " << covmat_baofit_ascii << endl;
    file = fopen(covmat_baofit_ascii.c_str(),"w");
    for(int i=0;i<n2d;i++) {
      for(int j=0;j<n2d;j++) {
      fprintf(file,"%d %d %g\n",i,j,cov_data[i+n2d*j]);
      }
    }
    fclose(file);
    */
    
    cout << "writing " << covmat_fits << endl;
    cov.writeFits(covmat_fits);
    
  }
  cout << "done" << endl;
}

void save_lya_correlation_results_with_npair(const string& filename, const map<int,Histo2d*>& h_sum_wdd, const map<int,Histo2d*>& h_sum_w,const map<int,Histo2d*>& h_npair,bool with_covmat) {
  cout << "computing results" << endl;

  map<int,Histo2d*>::const_iterator h_sum_wdd_it = h_sum_wdd.begin();
  map<int,Histo2d*>::const_iterator h_sum_w_it = h_sum_w.begin();
  map<int,Histo2d*>::const_iterator h_npair_it = h_npair.begin();

  int nplates = h_sum_wdd.size();
  int n2d = h_sum_wdd_it->second->nx * h_sum_wdd_it->second->ny;

  double swx[n2d];
  double sw[n2d]; 
  double sn[n2d];
  for(int i=0;i<n2d;i++) {
    swx[i]=0;
    sw[i]=0;
    sn[i]=0;
  }
  
  for (;h_sum_wdd_it != h_sum_wdd.end(); h_sum_wdd_it ++, h_sum_w_it ++, h_npair_it ++) {
    const double* w  = h_sum_w_it->second->data;
    const double* wx = h_sum_wdd_it->second->data;
    const double* n_p = h_npair_it->second->data;
    for(int i=0;i<n2d;i++) {
      swx[i] += wx[i];
      sw[i]  += w[i];
      sn[i]  += n_p[i];
    }
  }
  
  double mx[n2d];
  for(int i=0;i<n2d;i++) {
    if(sw[i]>0){
      mx[i] = swx[i]/sw[i];}
    else {
      mx[i] = 0;}
  
  }
  cout << "writing " << filename << endl;
  FILE *file = fopen(filename.c_str(),"w");
  
  for(int i=0;i<n2d;i++) {
    fprintf(file,"%d %g %g %g \n",i ,mx[i], sw[i], sn[i]);
  }
  fclose(file);
  if (with_covmat) {
    cout << "computing covmat ..." << endl;
  
    Mat cov(n2d,n2d);
    double wrx[n2d];
    double* cov_data = cov.NonConstData();
  
    h_sum_wdd_it = h_sum_wdd.begin();
    h_sum_w_it = h_sum_w.begin();
    int p=0;
    for (;h_sum_wdd_it != h_sum_wdd.end(); h_sum_wdd_it ++, h_sum_w_it ++, p++) {
      const double* w  = h_sum_w_it->second->data;
      const double* wx = h_sum_wdd_it->second->data;
      for(int i=0;i<n2d;i++) {
	wrx[i] = wx[i] - w[i]*mx[i];
	if(wrx[i]==0) continue;
	for(int j=0;j<=i;j++) {
	  cov_data[i+n2d*j] += (wrx[i]*wrx[j]);
	}
      }
      if(p%100==0)
	cout << "plate " << p << "/" << nplates << endl;
    }
    // normalize
    for(int i=0;i<n2d;i++) {
      for(int j=0;j<=i;j++) {
	if(sw[i]*sw[j]>0)
	  cov_data[i+n2d*j] /= (sw[i]*sw[j]);
      }
    }
    // symmetrize
    for(int i=0;i<n2d;i++) {
      for(int j=0;j<=i;j++) {
	cov_data[j+n2d*i] = cov_data[i+n2d*j];
      }
    }
      
    string covmat_baofit_ascii = CutExtension(filename)+".cov";
    string covmat_fits         = CutExtension(filename)+"-cov.fits";
    
    /*
    cout << "writing " << covmat_baofit_ascii << endl;
    file = fopen(covmat_baofit_ascii.c_str(),"w");
    for(int i=0;i<n2d;i++) {
      for(int j=0;j<n2d;j++) {
      fprintf(file,"%d %d %g\n",i,j,cov_data[i+n2d*j]);
      }
    }
    fclose(file);
    */
    
    cout << "writing " << covmat_fits << endl;
    cov.writeFits(covmat_fits);
    
  }
  cout << "done" << endl;
}  

void compute_indices_for_rmax(const Vect& dist, const double& rmax ,int *& begin_for_index,int *& end_for_index) {
  cout << "compute indices for rmax ..." << endl;
  int nw=dist.size();
  begin_for_index = new int[nw];
  begin_for_index[0]=0;
  for(int i=1;i<nw;i++) {
    begin_for_index[i]=0;
    for(int j=begin_for_index[i-1];j<=i;j++) {
      if(dist(i)-dist(j)<rmax) { 
	begin_for_index[i]=j;
	break;
      }
    }
  }
  end_for_index = new int[nw];
  end_for_index[nw-1]=nw;
  for(int i=nw-2;i>=0;i--) {
    end_for_index[i]=nw;
    for(int j=end_for_index[i+1]-1;j>i;j--) {
      if(dist(j)-dist(i)<rmax) { 
	end_for_index[i]=min(j+1,nw);
	break;
      }
    }
  }
  
  for(int i=0;i<nw;i++) {
    cout << i << " i_range=" << begin_for_index[i] << ":" << end_for_index[i] << " r_range=" << dist(begin_for_index[i])-dist(i) << ":" << dist(min(end_for_index[i],nw-1))-dist(i) << endl;
  }
}
