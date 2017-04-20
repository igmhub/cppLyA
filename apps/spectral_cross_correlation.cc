#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <time.h>

#include <matvect.h>
#include <cosmology.h>
#include <bspline1d.h>
#include <dictfile.h>
#include <histo1d.h>
#include <histo2d.h>
#include <fileutils.h>
#include <correlation_tools.h>

using namespace std;

#include <qso_spectrum.h>
#include <fiducial_cosmo.h>
#include <sampledfunction.h>

#ifdef USEOMP
#include <omp.h>
#endif

void usage(const string &pg) {
  cout << pg << " -f1 dflux1.fits -f2 dflux2.fits " << endl;
  cout << "options:" << endl;
  cout << " -F1 # (format of dflux1.fits, default 0)" << endl;
  cout << " -F2 # (format of dflux2.fits, default 0)" << endl;
  cout << " -o res.dat (output filename)" << endl;
  cout << " -w1 # (type of absorption 1, default lya)" << endl;
  cout << " -w2 # (type of absorption 2, default lya)" << endl;  
  cout << " -n (exclude pairs on same plate)" << endl;
  cout << " -B begin (start entry)" << endl;
  cout << " -E end (end entry, excluded)" << endl;
  cout << " -P (same half-plates)" << endl;
  cout << " -m # (mu min)" << endl;
  cout << " -W (wedges)" << endl;
  cout << " -z zref (default 2.25)" << endl;
  cout << " -a alpha (default 2.9)" << endl;  
  cout<<"spectral_cross_correlation computes the 2d auto or cross correlation function"<<endl; 
  exit(12);
}

#define _GNU_SOURCE 1
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <fenv.h>

int main(int argc, char** argv) {
  
  // to crash when NaN
  //feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);

  string fits_filename1="";
  string fits_filename2="";
  string line1="lya";
  string line2="lya";
  string output_filename="";
  bool not_same_plate = false;
  int format1 = 0;
  int format2 = 0;
  bool wedges=false;
  int begin_qso=0;
  int end_qso=-1;
  
  bool not_same_bundle=false;
  bool same_half_plates=false;

  double mu_min=0.;
  double zref=2.25; 
  double alpha=2.9;
  
  if(argc<2) {
    usage(argv[0]);
  }
  
  for (int i=1; i<argc; ++i) {
    char *arg = argv[i];
    if (arg[0] != '-') {
      usage(argv[0]);
    }
    switch (arg[1])
      {
      case 'h' : usage(argv[0]); break;
      case 'f' : {if(arg[2]=='1') fits_filename1=argv[++i]; else if(arg[2]=='2') fits_filename2=argv[++i]; else usage(argv[0]); } break;
      case 'F' : {if(arg[2]=='1') format1=atoi(argv[++i]); else if(arg[2]=='2') format2=atoi(argv[++i]); else usage(argv[0]); } break;
      case 'o' : output_filename=argv[++i]; break;
      case 'w' : {if(arg[2]=='1') line1=argv[++i]; else if(arg[2]=='2') line2=argv[++i]; else usage(argv[0]); } break;	
      case 'n' : not_same_plate=true; break;
      case 'b' : not_same_bundle=true; break;
      case 'N' : not_same_bundle=true; break;
      case 'B' : begin_qso=atoi(argv[++i]); break;
      case 'E' : end_qso=atoi(argv[++i]); break;
      case 'P' : same_half_plates = true; break;
      case 'm' : mu_min = atof(argv[++i]); break;
      case 'W' : wedges = true; break;
      case 'z' : zref = atof(argv[++i]); break;
      case 'a' : alpha = atof(argv[++i]); break;
      default : usage(argv[0]); break;
      }
  }

  if(not_same_bundle) {
    cout << "will ignore Lya pairs of QSO in same bundle" << endl;
  }
  if(output_filename=="") {
    usage(argv[0]);
  } 
  
  int ncpu = 1;
  char* OMP_NUM_THREADS = getenv("OMP_NUM_THREADS");
  if(OMP_NUM_THREADS) {
    ncpu = atoi(OMP_NUM_THREADS);
    cout << "will use " << ncpu << " CPU" << endl;
  }
  
  double ly_lambda1=lya_lambda; 
  double ly_lambda2=lya_lambda; 
  if (line1=="lyb") ly_lambda1=lyb_lambda;
  if (line2=="lyb") ly_lambda2=lyb_lambda;
  cout << line1<< " = "<< ly_lambda1<<endl; 
  cout << line2<< " = "<< ly_lambda2<<endl;

  cout << "reading data ..." << endl;
  Vect wave1;
  vector<Spectrum> spectra1;
  read_data(fits_filename1,wave1,spectra1,format1);
  Vect wave2;
  vector<Spectrum> spectra2;
  read_data(fits_filename2,wave2,spectra2,format2);  
  if(fabs(wave1(0)-wave2(0))>0.001) {
    cout << "cannot handle different wave grids" << endl;
    exit(12);
  }
  Vect& wave = wave1;
  
  GeneralCosmo cosmo(FIDUCIAL_OM,1.-FIDUCIAL_OM,-1.,0.);
  
  double rmin=0.;
  double rmax=200;
  double rstep=4;
  
  double deg2rad=M_PI/180.;
  double maxdist=rmax; // Mpc/h
  double minz=wave(0)/max(ly_lambda1,ly_lambda2)-1;
  double min_da=CSPEED/H0*(cosmo.Da(minz)*(1+minz)); // Mpc/h
  double max_angle=asin(maxdist/min_da);
  double cos_max_angle=cos(max_angle);
  int nbins1d=((rmax-rmin)/rstep);
  
  cout << "minimal redshift = " << minz << endl;
  cout << "minimal QSO dist = " << min_da << " Mpc/h" << endl;
  cout << "maximal angle    = " << max_angle/deg2rad << " deg" << endl;
  
  // create as many 2d correlation as there are plates
  // make plate indices 
  int twod=true;
  map<int,Histo2d*> h_sum_wdd;
  map<int,Histo2d*> h_sum_w;
  map<int,Histo2d*> h_sum_z;
  for(size_t s=0;s<spectra1.size();s++) {
    int plate=spectra1[s].plate;
    if(h_sum_wdd.find(plate)==h_sum_wdd.end()) {
      if(twod) {
	h_sum_wdd[plate] = new Histo2d(nbins1d,rmin,rmax,nbins1d,rmin,rmax);
	h_sum_w[plate]   = new Histo2d(nbins1d,rmin,rmax,nbins1d,rmin,rmax);
	h_sum_z[plate]   = new Histo2d(nbins1d,rmin,rmax,nbins1d,rmin,rmax);
      }
    }
  }
  
  if(false) check_nan(spectra1);
  if(false) check_nan(spectra2);
  Vect dist1=compute_lya_distances(wave,line1);
  Vect dist2=compute_lya_distances(wave,line2);

  Vect Z1=compute_redshift(wave,line1);
  multiply_weight_by_redshift_evolution(spectra1,zref,alpha,Z1);
  Vect Z2=compute_redshift(wave,line2);
  multiply_weight_by_redshift_evolution(spectra2,zref,alpha,Z2);
 
  compute_qso_directions(spectra1);
  compute_qso_directions(spectra2);

  multiply_flux_by_weight(spectra1);
  multiply_flux_by_weight(spectra2);

  int *wbegin1,*wend1; compute_valid_range(spectra1,wbegin1,wend1);
  int *wbegin2,*wend2; compute_valid_range(spectra2,wbegin2,wend2);  
  
  cout << "filling histograms ..." << endl;
  
  if(end_qso<0 || end_qso>int(spectra1.size()))
    { end_qso=spectra1.size(); }

  
  // LOOP ON QSO K //
  int cpu=0;
#ifdef USEOMP
#pragma omp parallel for 
#endif
  for(cpu=0;cpu<ncpu;cpu++) {
    int cpu_begin_qso = begin_qso + (cpu*(end_qso-begin_qso))/ncpu;
    int cpu_end_qso   = begin_qso + ((cpu+1)*(end_qso-begin_qso))/ncpu;
    if(cpu==(ncpu-1))
      cpu_end_qso = end_qso;
    
  for(int k=cpu_begin_qso;k<cpu_end_qso;k++) {
        
    Spectrum& spec_k = spectra1[k];
    if(!spec_k.valid) continue;
    if(k%1000==0) cout << "CPU #" << cpu << " done " << (100.*(k-cpu_begin_qso))/(cpu_end_qso-cpu_begin_qso) << "%" << endl;
    int plate_k=spec_k.plate;
    Histo2d* h_sum_wdd_plate = h_sum_wdd[plate_k];
    Histo2d* h_sum_w_plate = h_sum_w[plate_k];
    Histo2d* h_sum_z_plate = h_sum_z[plate_k];
      
    int spectro_k=0; if(spec_k.fiber>500) spectro_k=1;

    const double* wk=spec_k.weight.Data();
    const double* wfk=spec_k.flux.Data(); // flux have been multiplied by weights
    const double* dk=dist1.Data();
    
    int bundle_k=(spec_k.fiber-1)/20;
    
    for(size_t p=0;p<spectra2.size();p++) {
      
      Spectrum& spec_p = spectra2[p];
      if(!spec_p.valid) continue;
				   
      // selection 
      int plate_p=spec_p.plate;   
      int spectro_p=0; if(spec_p.fiber>500) spectro_p=1;
      if(not_same_bundle && (plate_k==plate_p) && (((spec_p.fiber-1)/20)==bundle_k)) continue;      
      if(not_same_plate && (plate_k==plate_p) && (spectro_k==spectro_p)) continue;
      if(same_half_plates && (! ( (plate_k==plate_p) && (spectro_k==spectro_p) ))) continue;   

      // compute angle between the two quasars
      double cos_angle=scalar_product(spec_k.dir,spec_p.dir);
      if(fabs(cos_angle-1.)<1.1e-11) continue; // no cross correlation (this can be a problem)
      if (cos_angle<=cos_max_angle) continue; 

      double angle=0;
      angle=acos(cos_angle); 
      bool in_same_half_plate = ( (plate_k==plate_p) && (spectro_k==spectro_p));

      const double* wp=spec_p.weight.Data();      
      const double* wfp=spec_p.flux.Data(); // flux have been multiplied by weights
      const double* dp=dist2.Data();
      
      for(int ik=wbegin1[k];ik<wend1[k];ik++) { // loop on waves
	
	const double& wkk=wk[ik];
	if(wkk==0) continue;
	const double& wfkk=wfk[ik];		
	const double& dkk=dk[ik];
	const double& zkk=(wave1(ik)/ly_lambda1 -1.)/2.;
	
	for(int ip=wbegin2[p];ip<wend2[p];ip++) { // loop on waves
	  
	  const double& wpp=wp[ip];
	  if(wpp==0) continue;
	  const double& wfpp=wfp[ip];	  
	  const double& dpp=dp[ip];
	  const double& zpp=(wave1(ip)/ly_lambda2 -1.)/2.;
	  
          double rpar=fabs(dkk-dpp)*cos(angle/2.); // Busca's definition 
	  if ((in_same_half_plate)&&(rpar<rstep)) continue; 
	  double rperp=(dkk+dpp)*sin(angle/2); // Busca's definition

	  if (rpar>=rmax) continue; 
	  if (rperp>=rmax) continue; 

	  h_sum_wdd_plate->Fill(rpar,rperp,wfkk*wfpp);
	  h_sum_w_plate->Fill(rpar,rperp,wkk*wpp);
	  h_sum_z_plate->Fill(rpar,rperp,wkk*wpp*(zkk+zpp));	    
	    
	} // end of loop on ip	  
      } // end of loop on ik 
    } // end of loop on p
    
  } // end of loop on k
  } // end of loop on cpu

  //save_lya_correlation_results_with_z(output_filename,h_sum_wdd,h_sum_w,h_sum_z);
  save_lya_correlation_results_fast(output_filename,h_sum_wdd,h_sum_w);
  
  cout << "spectral_cross_correlation exit successfully" << endl;
  return 0;
}
