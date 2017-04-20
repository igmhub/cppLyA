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
  cout << pg << " -f dflux.fits" << endl;
  cout << "options:" << endl;
  cout << " -F # (format of dflux.fits, default 0)" << endl;
  cout << " -o res.dat (output filename)" << endl;
  cout << " -n (exclude pairs on same plate)" << endl;
  cout << " -B begin (start entry)" << endl;
  cout << " -E end (end entry, excluded)" << endl;
  cout << " -P (same half-plates)" << endl;
  cout << " -m # (mu min)" << endl;
  cout << " -w (wedges)" << endl;
  cout << " -z zref (default 2.25)" << endl;
  cout << " -a alpha (default 2.9)" << endl;  
  exit(12);
}
#define _GNU_SOURCE 1
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <fenv.h>

int main(int argc, char** argv) {

  float computation_time = 0;
  clock_t k_start,k_stop;
  k_start=clock();

  // to crash when NaN
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);

  string fits_filename="";
  string output_filename="";
  bool not_same_plate = false;
  int format = 0;
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
      case 'f' : fits_filename=argv[++i]; break;
      case 'F' : format=atoi(argv[++i]); break;
      case 'o' : output_filename=argv[++i]; break;
      case 'n' : not_same_plate=true; break;
      case 'b' : not_same_bundle=true; break;
      case 'N' : not_same_bundle=true; break;
      case 'B' : begin_qso=atoi(argv[++i]); break;
      case 'E' : end_qso=atoi(argv[++i]); break;
      case 'P' : same_half_plates = true; break;
      case 'm' : mu_min = atof(argv[++i]); break;
      case 'w' : wedges = true; break;
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
  
  cout << "reading data ..." << endl;
  Vect wave;
  vector<Spectrum> spectra;
  read_data(fits_filename,wave,spectra,format);
  
  GeneralCosmo cosmo(FIDUCIAL_OM,1.-FIDUCIAL_OM,-1.,0.);
  
  double rmin=0.;
  double rmax=200;
  double rstep=4;
  
  double deg2rad=M_PI/180.;
  double maxdist=rmax; // Mpc/h
  double minz=wave(0)/lya_lambda-1;
  double min_da=CSPEED/H0*(cosmo.Da(minz)*(1+minz)); // Mpc/h
  double max_angle=asin(maxdist/min_da);
  double cos_max_angle=cos(maxdist/min_da);
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
  map<int,Histo2d*> h_sum_rp;
  map<int,Histo2d*> h_sum_rt;
  //map<int,Histo2d*> h_npair;
  for(size_t s=0;s<spectra.size();s++) {
    int plate=spectra[s].plate;
    if(h_sum_wdd.find(plate)==h_sum_wdd.end()) {
      if(twod) {
	h_sum_wdd[plate] = new Histo2d(nbins1d,rmin,rmax,nbins1d,rmin,rmax);
	h_sum_w[plate]   = new Histo2d(nbins1d,rmin,rmax,nbins1d,rmin,rmax);
	h_sum_z[plate]   = new Histo2d(nbins1d,rmin,rmax,nbins1d,rmin,rmax);
	h_sum_rp[plate]   = new Histo2d(nbins1d,rmin,rmax,nbins1d,rmin,rmax);
	h_sum_rt[plate]   = new Histo2d(nbins1d,rmin,rmax,nbins1d,rmin,rmax);
	//h_npair[plate]   = new Histo2d(nbins1d,rmin,rmax,nbins1d,rmin,rmax);
      }
    }
  }
  
  if(false) check_nan(spectra);
  Vect dist=compute_lya_distances(wave);
  int *begin_for_index,*end_for_index;
  compute_indices_for_rmax(dist,rmax,begin_for_index,end_for_index);

  Vect Z=compute_redshift(wave,"lya");
  multiply_weight_by_redshift_evolution(spectra,zref,alpha,Z);
  
  compute_qso_directions(spectra);
  multiply_flux_by_weight(spectra);

  int *wbegin,*wend;
  compute_valid_range(spectra,wbegin,wend);
  
  cout << "filling histograms ..." << endl;
  
  if(end_qso<0 || end_qso>int(spectra.size()))
    end_qso=spectra.size();
 
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
    Spectrum& spec_k = spectra[k];
    if(!spec_k.valid) continue;
    if(k%10==0) {
      cout << "CPU #" << cpu << " done " << (100.*(k-cpu_begin_qso))/(cpu_end_qso-cpu_begin_qso) << "%" << endl;
    }
    int plate_k=spec_k.plate;
    Histo2d* h_sum_wdd_plate = h_sum_wdd[plate_k];
    Histo2d* h_sum_w_plate = h_sum_w[plate_k];
    Histo2d* h_sum_z_plate = h_sum_z[plate_k];   
    Histo2d* h_sum_rp_plate = h_sum_rp[plate_k]; 
    Histo2d* h_sum_rt_plate = h_sum_rt[plate_k]; 
    //Histo2d* h_npair_plate = h_npair[plate_k];

    int spectro_k=0; if(spec_k.fiber>500) spectro_k=1;

    const double* wk=spec_k.weight.Data();
    const double* wfk=spec_k.flux.Data(); // flux have been multiplied by weights
    const double* dk=dist.Data();
    
    int bundle_k=(spec_k.fiber-1)/20;    
        
    for(size_t p=k+1;p<spectra.size();p++) {
      
      Spectrum& spec_p = spectra[p];
      if(!spec_p.valid) continue;
				   
      // selection 
      int plate_p=spec_p.plate;   
      int spectro_p=0; if(spec_p.fiber>500) spectro_p=1;
      if(not_same_bundle && (plate_k==plate_p) && (((spec_p.fiber-1)/20)==bundle_k)) continue;      
      if(not_same_plate && (plate_k==plate_p) && (spectro_k==spectro_p)) continue;
      if(same_half_plates && (! ( (plate_k==plate_p) && (spectro_k==spectro_p) ))) continue;
      
      bool in_same_half_plate = ( (plate_k==plate_p) && (spectro_k==spectro_p)); 
      // compute angle between the two quasars
      double cos_angle=scalar_product(spec_k.dir,spec_p.dir);
      if (cos_angle>1.) cos_angle=1.; 
      if (cos_angle<=cos_max_angle) continue; 

      double angle=0;
      angle=acos(cos_angle);


      const double* wp=spec_p.weight.Data();      
      const double* wfp=spec_p.flux.Data(); // flux have been multiplied by weights
      const double* dp=dist.Data();
            
      for(int ik=wbegin[k];ik<wend[k];ik++) { // loop on waves
	
	const double& wkk=wk[ik];
	if(wkk==0) continue;
	const double& wfkk=wfk[ik];		
	const double& dkk=dk[ik];
	const double& zkk=(wave(ik)/lya_lambda -1.)/2.;	
	
	int ipb=max(wbegin[p],begin_for_index[ik]);
	int ipe=min(wend[p],end_for_index[ik]);
	for(int ip=ipb;ip<ipe;ip++) { // loop on waves

	  const double& wpp=wp[ip];
	  if(wpp==0) continue;

	  const double& dpp=dp[ip];

	  double rpar=fabs(dkk-dpp)*cos(angle/2.); // Busca's definition 
	  if ((in_same_half_plate)&&(rpar<rstep)) continue;
	  double rperp=(dkk+dpp)*sin(angle/2); // Busca's definition

	  if (rpar>=rmax) continue; 
	  if (rperp>=rmax) continue; 

	  const double& wfpp=wfp[ip];	  
	  const double& zpp=(wave(ip)/lya_lambda -1.)/2.;

	  h_sum_wdd_plate->Fill(rpar,rperp,wfkk*wfpp);
	  h_sum_w_plate->Fill(rpar,rperp,wkk*wpp);
	  h_sum_z_plate->Fill(rpar,rperp,wkk*wpp*(zkk+zpp));	
	  h_sum_rp_plate->Fill(rpar,rperp,wkk*wpp*rpar);
	  h_sum_rt_plate->Fill(rpar,rperp,wkk*wpp*rperp);
	  //h_npair_plate->Fill(rpar,rperp,1);	    

	} // end of loop on ip	  
      } // end of loop on ik 
    } // end of loop on p
  } // end of loop on k
  } // end of loop on cpu

  
  //save_lya_correlation_results_fast(output_filename,h_sum_wdd,h_sum_w);
  save_lya_correlation_results(output_filename,h_sum_wdd,h_sum_w,h_sum_z,h_sum_rp,h_sum_rt);
  //save_lya_correlation_results_with_npair(output_filename,h_sum_wdd,h_sum_w,h_npair);   // We want to compute the number of pair per bin

  k_stop=clock();
  computation_time = float(k_stop-k_start) / CLOCKS_PER_SEC;
  cout << "time (h) =" << computation_time/3600. << " /" << ncpu << " = " << computation_time/3600./ncpu ; 
  cout << ' '; 
  cout << "lya_auto_correlation exit successfully" << endl;
  return 0;
}
