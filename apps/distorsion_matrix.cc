#include <iostream>
#include <string>
#include <cmath>
#include <map>

#include <matvect.h>
#include <cosmology.h>
#include <bspline1d.h>
#include <dictfile.h>
#include <correlation_tools.h>
#include <fileutils.h>

using namespace std;

#include <qso_spectrum.h>
#include <fiducial_cosmo.h>


#ifdef USEOMP
#include <omp.h>
#endif


void usage(const string &pg) {
  cout << pg << " -f1 dflux1.fits (-f2 dflux2.fits)" << endl;
  cout << "options:" << endl;
  cout << " -F1 # (format of dflux1.fits)" << endl;
  cout << " -F2 # (format of dflux2.fits)" << endl;
  cout << " -w1 (lya or lyb) (wavelength to determine redshift for 1st data set, default is lya)" << endl;
  cout << " -w2 (lya or lyb) (wavelength to determine redshift for 2nd data set, default is lya)" << endl;  
  cout << " -o res.dat (output filename)" << endl;
  cout << " -n (exclude pairs on same plate)" << endl;
  cout << " -r r_par_min (default 0)" << endl;
  cout << " -R r_par_max (default 200)" << endl;
  cout << " -P (same half-plates)" << endl;
  cout << " -O (same plate but other spectro)" << endl;
  cout << " -B1 begin1 (start entry reading dflux1.fits)" << endl;
  cout << " -E1 end1 (end entry reading dflux1.fits, excluded)" << endl;
  cout << " -s # (only data from this spectro)" << endl;
  cout << " -g # (change gamma value, default = 3.8)" << endl;
  cout << " -G apply gamma factor term so dmat applies to xi(z=zref)" << endl;
  cout << " -a (include same wavelength same plate pairs)" << endl;
  cout << "examples:" << endl;
  cout << "distorsion_matrix -f1 these/delta/delta_alpha.fits -f2 these/delta/delta_beta.fits -w1 lya -w2 lya -o dmat-B-aa.fits >& distorsion_matrix__B_aa.log &" << endl;
  cout << "for test :" << endl;
  cout << "distorsion_matrix -f1 these/delta/delta_beta.fits -w1 lya -o dmat-C-aa.fits -E1 1" << endl;
  cout << "choose number of CPUs with 'export OMP_NUM_THREADS=ncpu'" << endl;
  exit(12);
}



int main(int argc, char** argv) {
  
  string fits_filename1="";
  string fits_filename2="";
  string line1="lya";
  string line2="lya";
  string output_filename="";
  int format1 = 0;
  int format2 = 0;
  bool not_same_plate = false;
  bool other_spectro = false;
  double r_par_min=0.;
  double r_par_max=200;
  bool include_same_wavelength_same_plate_pairs = false;
  double gamma=3.8;
  double zref=2.25;
  bool apply_gamma_factor = false;
  int begin1=0;
  int end1=-1;

  bool same_half_plates=false;
  int only_spectro=0;
  
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
      case 'w' : {if(arg[2]=='1') line1=argv[++i]; else if(arg[2]=='2') line2=argv[++i]; else usage(argv[0]); } break;
      case 'o' :  output_filename=argv[++i]; break;
      case 'n' :  not_same_plate=true; break;
      case 'r' :  r_par_min=atof(argv[++i]); break;
      case 'R' :  r_par_max=atof(argv[++i]); break;
      case 'B' : {if(arg[2]=='1') begin1=atoi(argv[++i]); else usage(argv[0]);} break;
      case 'E' : {if(arg[2]=='1') end1=atoi(argv[++i]); else usage(argv[0]);} break;
      case 'P' : same_half_plates = true; break;
      case 'O' : other_spectro = true; break;
      case 's' : only_spectro = atoi(argv[++i]); break;
      case 'a' : include_same_wavelength_same_plate_pairs = true; break;
      case 'g' : gamma=atof(argv[++i]); break;
      case 'G' : apply_gamma_factor=true; break;
      default : usage(argv[0]); break; 
      }
  }
  
  bool r_par_abs=(r_par_min>=0);
  
  if(r_par_abs)
    cout << "will do |rp|" << endl;
  cout << "gamma=" << gamma << endl;
  if(apply_gamma_factor) {
    cout << "WARNING : will apply a gamma factor for the eta " << endl;
  }

  cout << "will read " << fits_filename1 << " " << format1 << endl;
  if (fits_filename2 != "") {
    cout << "will read " << fits_filename2 << " " << format2 << endl;
  }else{
    line2=line1;
  }
  cout << "wave1=" << line1 << " wave2=" << line2 << endl;
  
  if(output_filename=="") {
    usage(argv[0]);
  }
  
  int ncpu = 1;
  char* OMP_NUM_THREADS = getenv("OMP_NUM_THREADS");
  if(OMP_NUM_THREADS) {
    ncpu = atoi(OMP_NUM_THREADS);
    cout << "will use " << ncpu << " CPU" << endl;
  }
  
  
  
  Vect wave1; vector<Spectrum> spectra1;
  read_data(fits_filename1,wave1,spectra1,format1);  
  Vect wave2; vector<Spectrum> spectra2;
  
  bool autocorr = (fits_filename2 == "");
  if (! autocorr ){
    read_data(fits_filename2,wave2,spectra2,format2);
    if(fabs(wave1(0)-wave2(0))>0.001) {
      cout << "cannot handle different wave grids" << endl;
      exit(12);
    }
  }
    
  Vect& wave = wave1;
  double logwave[wave.size()];
  for(size_t i=0;i<wave.size();++i)
    logwave[i]=log10(wave(i));
  
  
  double half_gamma_factor1[wave.size()];
  double inv_half_gamma_factor1[wave.size()];
  double half_gamma_factor2[wave.size()];
  double inv_half_gamma_factor2[wave.size()];
  double waveref1,waveref2;
  if(line1=="lya")
    waveref1 = (1+zref)*lya_lambda;
  else 
    waveref1 = (1+zref)*lyb_lambda;
  if(line2=="lya")
    waveref2 = (1+zref)*lya_lambda;
  else 
    waveref2 = (1+zref)*lyb_lambda;
  
  for(size_t i=0;i<wave.size();++i) {
    
    half_gamma_factor1[i] = pow(wave(i)/waveref1,gamma/2.);
    half_gamma_factor2[i] = pow(wave(i)/waveref2,gamma/2.);
    inv_half_gamma_factor1[i] = 1./half_gamma_factor1[i];
    inv_half_gamma_factor2[i] = 1./half_gamma_factor2[i];
  }
  

  GeneralCosmo cosmo(FIDUCIAL_OM,1.-FIDUCIAL_OM,-1.,0.);

  double r_perp_min=0.;
  double r_perp_max=200.; 
  double rstep=4.;

  double deg2rad=M_PI/180.;
  double maxdist=180.; // Mpc/h 
  double rest_lambda=0;
  if((line1=="lya") && (line2=="lya"))
    rest_lambda=lya_lambda;
  else
    rest_lambda=lyb_lambda;
  double minz=wave(0)/rest_lambda-1;
  double min_da=CSPEED/H0*(cosmo.Da(minz)*(1+minz)); // Mpc/h
  double max_angle=asin(maxdist/min_da);
  //double min_cos=cos(maxdist/min_da);
  int nbins_par=((r_par_max-r_par_min)/rstep);
  int nbins_perp=((r_perp_max-r_perp_min)/rstep);

  cout << "rp nbins= " << nbins_par << " min= " << r_par_min << " max = " << r_par_max  << endl;
  cout << "rt nbins= " << nbins_perp << " min= " << r_perp_min << " max = " << r_perp_max  << endl;
  cout << "minimal redshift = " << minz << endl;
  cout << "maximal angle    = " << max_angle/deg2rad << " deg" << endl;
  
  Histo2d histo2d(nbins_par,r_par_min,r_par_max,nbins_perp,r_perp_min,r_perp_max);
  
  vector<Mat> WD_vector;
  vector<Vect> W_vector;
  for(int cpu=0;cpu<ncpu;cpu++) {
    WD_vector.push_back(Mat(nbins_par*nbins_perp,nbins_par*nbins_perp));
    W_vector.push_back(Vect(nbins_par*nbins_perp));
  }
  
  Vect dist1 = compute_lya_distances(wave,line1); 
  compute_qso_directions(spectra1);
  int *wbegin1,*wend1; compute_valid_range(spectra1,wbegin1,wend1);
  if(false) check_nan(spectra1);  
  
  Vect dist2;
  int *wbegin2,*wend2;
  if(autocorr) {
    dist2 = dist1;
    wbegin2=wbegin1;
    wend2=wend1;
  }else{
    dist2 = compute_lya_distances(wave,line2); 
    if(false) check_nan(spectra2);  
    compute_qso_directions(spectra2);
    compute_valid_range(spectra2,wbegin2,wend2);  
  }
  
  if(end1<0) end1=spectra1.size();

  
  

  size_t begin_k_given_i[wave.size()];
  size_t end_k_given_i[wave.size()];
  size_t begin_i_given_k[wave.size()];
  size_t end_i_given_k[wave.size()];
  // default
  for(size_t u=0; u<wave.size();u++) {
    begin_k_given_i[u]=wave.size();
    end_k_given_i[u]=0;
    begin_i_given_k[u]=wave.size();
    end_i_given_k[u]=0;
  }
  // not necessarily symmetric
  for(size_t i=0; i<wave.size();i++) {
    for(size_t k=0; k<wave.size();k++) {
      double d=dist1(i)-dist2(k);
      if(r_par_abs) d=fabs(d);
      if((d >= histo2d.minx) && (d <= (histo2d.minx + histo2d.nx/histo2d.scalex))) {
	begin_k_given_i[i]=min(k,begin_k_given_i[i]);
	end_k_given_i[i]=max(k,end_k_given_i[i]);
	begin_i_given_k[k]=min(i,begin_i_given_k[k]);
	end_i_given_k[k]=max(i,end_i_given_k[k]);	
      }
    }
  }
  for(size_t u=0; u<wave.size();u++) {
    end_k_given_i[u] +=1;
    end_i_given_k[u] +=1;
  }

  
  
  cout << "begin with loop ..." << endl; 

  /*
    delta_q' = dq' = dq - <dq> - l * <l dq> / <l**2>
    
    dq  == delta_q
    l  == log(lambda)-<log(lambda)>_q
    dqi == delta_{q,i}
    sqj  == sum_j (for the QSO q)
    
    dqi' = dqi - (sqj wqj*dqj)/(sqj wqj) - li*(sqj wqj*dqj*lj)/(sqj wqj*lj**2)
    
    another quasar :
    dpk'  = dpk - (spm wpm*dpm)/(spm wpm) - lk*(spm wpm*dpm*lm)/(spm wpm*lm**2) 
    
    
    dqi'*dpk' = dqi*dpk - (sqj wqj*(dqj*dpk))/(sqj wqj) - li*(sqj wqj*(dqj*dpk)*lj)/(sqj wqj*lj**2)
                                        - (spm wpm*(dpm*dqi))/(spm wpm) - lk*(spm wpm*((dpm*dqi)*lm)/(spm wpm*lm**2) 
					+ (sqj wqj*dqj)/(sqj wqj)*(spm wpm*dpm)/(spm wpm) 
					+ li*(sqj wqj*dqj*lj)/(sqj wqj*lj**2)*(spm wpm*dpm)/(spm wpm)
					+ lk*(spm wpm*dpm*lm)/(spm wpm*lm**2)*(sqj wqj*dqj)/(sqj wqj)
					+ li*(sqj wqj*dqj*lj)/(sqj wqj*lj**2)*lk*(spm wpm*dpm*lm)/(spm wpm*lm**2)
			
    swq=(sqj wqj)
    swl2q=(sqj wqj*lj**2)
    xabcd= dab*dcd
    
    
    dqi'*dpk' = dqi*dpk - (sqj wqj*(dqj*dpk))/swq - li*(sqj wqj*(dqj*dpk)*lj)/swl2q
                                        - (spm wpm*(dpm*dqi))/swp - lk*(spm wpm*((dpm*dqi)*lm)/swl2p
					+ (sqj wqj*dqj)/swq*(spm wpm*dpm)/swp
					+ li*(sqj wqj*dqj*lj)/swl2q*(spm wpm*dpm)/swp
					+ lk*(spm wpm*dpm*lm)/swl2p*(sqj wqj*dqj)/swq
					+ li*(sqj wqj*dqj*lj)/swl2q*lk*(spm wpm*dpm*lm)/swl2p
    dqi'*dpk' = dqi*dpk 
                                        - sqj (wqj/swq)*xqjpk
                                        - sqj (li*lj*wqj/swl2q)*xqjpk
                                        - spm (wpm/swp)*xpmqi
					- spm (lk*lm*wpm/swl2p)*xpmqi
					+ sqj spm (wqj*wpm/swp/swq)*xqjpm
					+ sqj spm (wpm*wqj*li*lj/swl2q/swp) xqjpm
					+ sqj spm (wpm*wqj*lk*lm/swl2p/swq) xqjpm
					+ sqj spm (wpm*wqj*li*lj*lk*lm/swl2q/swl2p) xqjpm
    dqi'*dpk' = dqi*dpk 
                                        - sqj [ (wqj/swq)+(li*lj*wqj/swl2q) ]*xqjpk
                                        - spm [ (wpm/swp)+(lk*lm*wpm/swl2p) ]*xpmqi				    
					+ sqj spm [ (wqj*wpm/swp/swq)+(wpm*wqj*li*lj/swl2q/swp)+(wpm*wqj*lk*lm/swl2p/swq)+(wpm*wqj*li*lj*lk*lm/swl2q/swl2p) ]*xqjpm
					
    in code, need to fill histograms of wqi*wpk*dqi'*dpk' and wqi*wpk
    correlation function wo distortion
    
    1) find rp,rt bin A of qi,pk
    2) fill hww[A] += wqi*wpk hwwxi[A] += wqi*wpk*dqi*dpk 
    
    
    correlation function with distortion
    1) find rp,rt bin A of qi,pk and rp,rt bin B of qj,pm , fill distort_hww[A,B] += 
    2) fill distort_hww[A,B] += ....
   */
  
  const double* dq = dist1.Data(); 
  const double* dp = dist2.Data(); 
  
  /* dry run */
  for(int cpu=0;cpu<ncpu;cpu++) {
    int cpu_begin1 = begin1 + (cpu*(end1-begin1))/ncpu;
    int cpu_end1   = begin1 + ((cpu+1)*(end1-begin1))/ncpu;
    if(cpu==(ncpu-1))
      cpu_end1 = end1;
    cout << "will do with CPU #" << cpu << " [" << cpu_begin1 << "," << cpu_end1 << "]" << endl;
  }
  
  
  int cpu=0;
#ifdef USEOMP
#pragma omp parallel for 
#endif
  for(cpu=0;cpu<ncpu;cpu++) {
    int cpu_begin1 = begin1 + (cpu*(end1-begin1))/ncpu;
    int cpu_end1   = begin1 + ((cpu+1)*(end1-begin1))/ncpu;
    if(cpu==(ncpu-1))
      cpu_end1 = end1;

    Mat& WD = WD_vector[cpu];
    Vect & W = W_vector[cpu];
    double *WDa = WD.NonConstData();
    int WDnx=WD.SizeX();
    double *Wa = W.NonConstData();
    
    // needed , and MUST instanciate per CPU
    double lq[wave.size()];
    double wqn[wave.size()];
    double wlqn[wave.size()];
    double lp[wave.size()];
    double wpn[wave.size()];
    double wlpn[wave.size()];
    double w2qipk_eqi[wave.size()];
    double epk[wave.size()];
    int *histo2d_indices_for_ik[wave.size()];
    for(size_t i=0; i<wave.size();i++) {
      histo2d_indices_for_ik[i]=new int[wave.size()];    
    }
    
  
    for(int q=cpu_begin1; q<cpu_end1;q++){ //  loop on spectra q of 1st sample
      cout << "CPU #" << cpu << " done " << (100.*(q-cpu_begin1))/(cpu_end1-cpu_begin1) << "% (doing q=" << q << ")" << endl;
      
      const Spectrum& spec_q = spectra1[q]; 
      if(!spec_q.valid) continue; 
      
      const double* wq = spec_q.weight.Data();
      size_t beginq=wbegin1[q];
      size_t endq=wend1[q];
      
      //cout << "q=" << q << " dm(end) - dm(begin) =" << dq[endq] << "-" << dq[beginq] << " = " << dq[endq]-dq[beginq] << endl; if(q>100) exit(12); continue;
      
      // compute sums for QSO q and lq
      double swq=0;
      double swl2q=0;
      
      {
	double mean_logwave=0;      
	for(size_t u=beginq;u<endq;u++) {
	  swq += wq[u]*inv_half_gamma_factor1[u];
	  mean_logwave += wq[u]*inv_half_gamma_factor1[u]*logwave[u];
	}
	if(swq==0) 
	  continue;
	mean_logwave /= swq;
	for(size_t u=beginq;u<endq;u++) {
	  lq[u] = logwave[u]-mean_logwave;
	  swl2q += wq[u]*inv_half_gamma_factor1[u]*lq[u]*lq[u];
	}
      }
      if(swl2q==0) continue; // can happen if only one pix
      
      // compute the normed quantities
      double swq_inv = 1./swq;
      double swl2q_inv = 1./swl2q;
      for(size_t u=beginq;u<endq;u++) {
	wqn[u]=wq[u]*inv_half_gamma_factor1[u]*swq_inv;
	wlqn[u]=wq[u]*inv_half_gamma_factor1[u]*lq[u]*swl2q_inv;
      }
      
      int begin_p,end_p;
      if(autocorr) {
	begin_p = q+1;
	end_p   = spectra1.size();
      }else{
	begin_p = 0;
	end_p   = spectra2.size();
      }
      
      for(int p=begin_p; p<end_p; ++p) { // spectra p 
	
	if(autocorr && p==q) continue;
	
	const Spectrum* spec_pointer=0;
	if(autocorr)
	  spec_pointer = & spectra1[p];
	else
	  spec_pointer = & spectra2[p];
	
	const Spectrum& spec_p = *spec_pointer;
	
	if(!spec_p.valid) continue;
	
	double cos_angle=scalar_product(spec_q.dir, spec_p.dir);
	if(fabs(cos_angle-1.)<1.1e-11) continue; // no cross correlation (this can be a problem)
	//if(abs(cos_angle)<min_cos) continue;
	if (abs(cos_angle)>max_angle) continue; 
	double angle=0;
	//if(abs(cos_angle)<0.999999) angle=acos(cos_angle);
	angle=acos(cos_angle); 
	double halfangle=0.5*angle;
	
	const double* wp = spec_p.weight.Data(); 
	//const double* fp = spec_p.flux.Data();
	size_t beginp=wbegin2[p];
	size_t endp=wend2[p];
	
	size_t beginqp=max(beginq,begin_i_given_k[beginp]);
	size_t endqp=min(endq,end_i_given_k[endp-1]);
	
	size_t beginpq=max(beginp,begin_k_given_i[beginq]);
	size_t endpq=min(endp,end_k_given_i[endq-1]);
	
	
	//#define CHECK_INDICES	

#ifdef CHECK_INDICES
	// by default, all = -2 to detect where we test
	for(int i = 0; i<wave.size(); ++i) {
	    for(int k = 0; k<wave.size(); ++k) {
	      histo2d_indices_for_ik[i][k]=-2;
	    }
	}
#endif
	// precompute all possible indices 
	if(r_par_abs) {
	  for(int i = beginqp; i<endqp; ++i) {
	    for(int k = beginpq; k<endpq; ++k) {
	      double rp=fabs(dq[i]-dp[k])*cos(halfangle);  
	      double rt=(dq[i]+dp[k])*sin(halfangle);
	      histo2d_indices_for_ik[i][k]=histo2d.index_1d(rp,rt); // rp , rt 
	    }
	  }
	}else{
	  for(int i = beginqp; i<endqp; ++i) {
	    for(int k = beginpq; k<endpq; ++k) {
	      double rp=(dq[i]-dp[k])*cos(halfangle);
              double rt=(dq[i]+dp[k])*sin(halfangle);
	      histo2d_indices_for_ik[i][k]=histo2d.index_1d(rp,rt); // rp , rt 
	    }
	  }
	}
	
	

	// compute sums for QSO p and lp
	double swp=0;
	double swl2p=0;
	
	{
	  double mean_logwave=0;      
	  for(size_t u=beginp;u<endp;u++) {
	    swp += wp[u]*inv_half_gamma_factor2[u];
	    mean_logwave += wp[u]*inv_half_gamma_factor2[u]*logwave[u];
	  }
	  if(swp==0) 
	    continue;
	  mean_logwave /= swp;	  
	  for(size_t u=beginp;u<endp;u++) {
	    lp[u] = logwave[u]-mean_logwave;
	    swl2p += wp[u]*inv_half_gamma_factor2[u]*lp[u]*lp[u];
	  }
	}
	if(swl2p==0) continue; // can happen if only one pix
	
	// compute the normed quantities
	double swp_inv = 1./swp;
	double swl2p_inv = 1./swl2p;
	for(size_t u=beginp;u<endp;u++) {
	  wpn[u]=wp[u]*inv_half_gamma_factor2[u]*swp_inv;
	  wlpn[u]=wp[u]*inv_half_gamma_factor2[u]*lp[u]*swl2p_inv;	  	 
	}
	
	for(int i = beginqp; i<endqp; ++i) {
	  const double& wqi = wq[i];
	  if (wqi == 0) continue; 
	  const double& li  = lq[i];
	  
	  for(int k = max(beginpq,begin_k_given_i[i]); k<min(endpq,end_k_given_i[i]); ++k) {

	    const double& wpk = wp[k];
	    if (wpk == 0) continue; 
	    const double& lk  = lp[k];
	    
	    int index_qipk = histo2d_indices_for_ik[i][k];
	    if(index_qipk<0) continue;
	    
	    double w2qipk = wqi*wpk;
	    // in matvect.h operator () (const size_t i, const size_t j) data[i+j*nx];

	    int index_qipk_n = index_qipk*WDnx;

	    // we transposed the matrix so that it is consistent with other implementation
	    // with this, if D is written to fits in this C++ code with D.writeFits("d.fits")
	    // and then read in python with D=pyfits.open("d.fits")[0].data
	    // the application of the distorsion matrix to a model xi (in the form of a vector)
	    // in python is : xi_distorted = D.dot(xi)
	    
	    WDa[index_qipk+index_qipk_n] += w2qipk;
	    Wa[index_qipk] += w2qipk;
	    	    
	    for(size_t j=beginqp;j<endqp;++j)
	      w2qipk_eqi[j]=w2qipk*(wqn[j]+li*wlqn[j]);
	    for(size_t m=beginpq; m<endpq; ++m)
	      epk[m]=wpn[m]+lk*wlpn[m];
	    
	    size_t beginj=max(beginqp,begin_i_given_k[k]);
	    size_t endj=min(endqp,end_i_given_k[k]);
	    size_t beginm=max(beginpq,begin_k_given_i[i]);
	    size_t endm=min(endpq,end_k_given_i[i]);
	    
	    // all of the computation time is below
	    
	    if(apply_gamma_factor) { // this is a special case for tests	      
	      for(size_t j=beginj;j<endj;++j) { // second loop on QSO q 	      
		int index_qjpk = histo2d_indices_for_ik[j][k];
		if(index_qjpk<0) continue;
		WDa[index_qjpk+index_qipk_n] -= inv_half_gamma_factor1[i]*half_gamma_factor1[j]*w2qipk_eqi[j];
	      }
	      for(size_t m=beginm; m<endm; ++m) { // second loop on QSO p	      
		int index_qipm  = histo2d_indices_for_ik[i][m];
		if(index_qipm<0) continue;
		WDa[index_qipm+index_qipk_n] -= inv_half_gamma_factor2[k]*half_gamma_factor2[m]*w2qipk*epk[m];
	      }
	      for(size_t j=beginqp; j<endqp; ++j) { // second loop on QSO q 
		if(wq[j]==0) continue;
		size_t beginm=max(beginpq,begin_k_given_i[j]);
		size_t endm=min(endpq,end_k_given_i[j]);
		for(size_t m=beginm; m<endm; ++m) { // second loop on QSO p
		  int index_qjpm = histo2d_indices_for_ik[j][m];
		  if(index_qjpm<0) continue;		  	     
		  WDa[index_qjpm+index_qipk_n] += inv_half_gamma_factor1[i]*half_gamma_factor1[j]*inv_half_gamma_factor2[k]*half_gamma_factor2[m]*w2qipk_eqi[j]*epk[m]; 
		}
	      }
	      
	    }else{ // standard calculation, with more comments and tests
	    
	      // term : - sqj [ (wqj/swq)+(li*lj*wqj/swl2q) ]*xqjpk
	      for(size_t j=beginj;j<endj;++j) { // second loop on QSO q 
		
		int index_qjpk = histo2d_indices_for_ik[j][k];
		if(index_qjpk<0) {
#ifdef CHECK_INDICES
		  if (index_qjpk==-2) {
		    cout << "unexpected index=-2 for j=" << j << "k=" << k << endl;
		    cout << "beginqp:endqp= " << beginqp << ":" << endqp << endl;
		    cout << "beginpq:endpq= " << beginpq << ":" << endpq << endl;
		    exit(12);
		  }
#endif
		  continue;
		}
		//WDa[index_qipk+index_qjpk*WDnx] += -w2qipk*(wq[j]/swq+li*lq[j]*wq[j]/swl2q);
		//WDa[index_qipk+index_qjpk*WDnx] -= (w2qipk_div_swq+w2qipk_li_div_swl2q*lq[j])*wq[j];
		//WDa[index_qipk+index_qjpk*WDnx] -= w2qipk*(wqn[j]+li*wlqn[j]);
		WDa[index_qjpk+index_qipk_n] -= w2qipk_eqi[j];
	      } // j
	      
	      // term : - spm [ (wpm/swp)+(lk*lm*wpm/swl2p) ]*xpmqi	
	      for(size_t m=beginm; m<endm; ++m) { // second loop on QSO p	      
		int index_qipm  = histo2d_indices_for_ik[i][m];
		if(index_qipm<0) {
#ifdef CHECK_INDICES
		  if(index_qipm==-2) {
		    cout << "unexpected index=-2 for i=" << i << "m=" << m << endl;
		    cout << "beginqp:endqp= " << beginqp << ":" << endqp << endl;
		    cout << "beginpq:endpq= " << beginpq << ":" << endpq << endl;		  
		    exit(12);
		  }
#endif
		  continue;
		}
		
		//WDa[index_qipk+index_qipm*WDnx] += -w2qipk*(wp[m]/swp+lk*lp[m]*wp[m]/swl2p);
		//WDa[index_qipk+index_qipm*WDnx] -= (w2qipk_div_swp+w2qipk_lk_div_swl2p*lp[m])*wp[m];
		//WDa[index_qipk+index_qipm*WDnx] -= w2qipk*(wpn[m]+lk*wlpn[m]);
		WDa[index_qipm+index_qipk_n] -= w2qipk*epk[m];
		
	      } // m
	      
	      // term : + sqj spm [ (wqj*wpm)* [ (1/swp/swq)+(li*lj/swl2q/swp)+(lk*lm/swl2p/swq)+(li*lj*lk*lm/swl2q/swl2p) ] ]*xqjpm
	      for(size_t j=beginqp; j<endqp; ++j) { // second loop on QSO q 
		if(wq[j]==0) continue;
		size_t beginm=max(beginpq,begin_k_given_i[j]);
		size_t endm=min(endpq,end_k_given_i[j]);
		for(size_t m=beginm; m<endm; ++m) { // second loop on QSO p
		  int index_qjpm = histo2d_indices_for_ik[j][m];
		  if(index_qjpm<0) {
#ifdef CHECK_INDICES
		    if(index_qjpm==-2) {
		      cout << "unexpected index=-2 for j=" << j << "m=" << m << endl;
		      cout << "beginqp:endqp= " << beginqp << ":" << endqp << endl;
		      cout << "beginpq:endpq= " << beginpq << ":" << endpq << endl;		  
		      exit(12);
		    }
#endif
		    continue;
		    
		  }
		  //WDa[index_qipk+index_qjpm*WDnx] += w2qipk*wq[j]*wp[m]*(1./swp/swq + li*lq[j]/swl2q/swp + lk*lp[m]/swl2p/swq + li*lq[j]*lk*lp[m]/swl2q/swl2p);
		  //WDa[index_qipk+index_qjpm*WDnx] += (w2qipk_div_swp_swq + w2qipk_li_div_swp_swl2q*lq[j] + w2qipk_lk_div_swl2p_swq*lp[m] + w2qipk_li_lk_div_swl2p_swl2q*lq[j]*lp[m])*wq[j]*wp[m];
		  //WDa[index_qipk+index_qjpm*WDnx] += w2qipk*(wqn[j]+li*wlqn[j])*(wpn[m]+lk*wlpn[m]);
		  WDa[index_qjpm+index_qipk_n] += w2qipk_eqi[j]*epk[m]; 
		} // m
	      } // j
	      
	    } // end of std case without gamma factor
	    
	  } // k
	} // i
      } // p
    } // q
  } // cpu
  
  cout << "stack WD and W of cpus ... " << endl;
  Mat WD(nbins_par*nbins_perp,nbins_par*nbins_perp);
  Vect W(nbins_par*nbins_perp);
  for(int cpu=0;cpu<ncpu;cpu++) {
    WD += WD_vector[cpu];
    W += W_vector[cpu];
  }
  cout << "D=WD/W  ... " << endl;
  for(size_t u=0;u<W.size();++u) {
    for(size_t v=0;v<W.size();++v) {
      if (W(v)>0) 
	WD(u,v) /= W(v); // it is index v for which we filled W
    }
  }
  cout << "Writing D in " << output_filename << " ... " << endl;
  WD.writeFits(output_filename);
  cout << "done with distorsion_matrix. " << endl; 
}
