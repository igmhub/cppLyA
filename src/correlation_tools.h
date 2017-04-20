#ifndef CORRELATION_TOOLS_H
#define CORRELATION_TOOLS_H
#include <map>
#include <qso_spectrum.h>
#include <histo2d.h>
#include <histo1d.h>
class Vect;
Vect compute_lya_distances(const Vect& wave,const string& line="lya");
Vect compute_redshift(const Vect& wave,const string& line="lya");
void compute_qso_distances(vector<QSO>& qsos);
void compute_qso_directions(vector<Spectrum>& spectra);
void compute_qso_directions(vector<QSO>& qsos);
void check_nan(const std::vector<Spectrum>& spectra);
void compute_valid_range(const std::vector<Spectrum>& spectra,int *& wbegin,int *& wend);
void compute_indices_for_rmax(const Vect& dist, const double& rmax ,int *& begin_for_index,int *& end_for_index);
void multiply_flux_by_weight(std::vector<Spectrum>& spectra);
void multiply_weight_by_redshift_evolution(std::vector<Spectrum>& spectra,const double& z_ref,const double& z_evol,const Vect& Z);
void save_lya_correlation_results_fast(const string& filename, const map<int,Histo2d*>& h_sum_wdd, const map<int,Histo2d*>& h_sum_w,bool with_covmat=true);
void save_lya_correlation_results(const string& filename, const map<int,Histo2d*>& h_sum_wdd, const map<int,Histo2d*>& h_sum_w, const map<int,Histo2d*>& h_sum_z, const map<int,Histo2d*>& h_sum_rp, const map<int,Histo2d*>& h_sum_rt, bool with_covmat=true);
void save_lya_correlation_results_with_npair(const string& filename, const map<int,Histo2d*>& h_sum_wdd, const map<int,Histo2d*>& h_sum_w, const map<int,Histo2d*>& h_npair,bool with_covmat=true);
void save_cross_correlation_results(const string& filename, double r_perp_min,double r_par_min, double r_perp_max, double rstep,const map<int,Histo2d*>& h_sum_wdd, const map<int,Histo2d*>& h_sum_w,const map<int,Histo2d*>& h_sum_wz, bool with_covmat=true);
#endif
