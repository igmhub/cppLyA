#ifndef QSO_SPECTRUM_H
#define QSO_SPECTRUM_H
#include <matvect.h>
#include <bspline1d.h>

#include <fitsio.h>
#define CHECKERROR if(status) {fits_report_error(stdout, status); cerr << "fits error" << endl; exit(12);}  

#define lyb_lambda 1025.728   
#define lya_lambda 1215.668

class PlateMJD {
 public :
  int plate;
  int mjd;
  PlateMJD() : plate(0), mjd(0) {}
  PlateMJD(int p, int m) : plate(p), mjd(m) {}
};

Vect compute_dir(const double& ra_deg, const double& dec_deg);

class QSO {
  public :
  
  int plate;
  int mjd;
  int fiber;
  double z;
  double ra; // deg
  double dec; // deg
  double xfocal;
  double yfocal;
  double dist;
  double weight;
  Vect dir;
  QSO() {
    weight=1;
    z=-1;
  };
  void set_dir();
};

class Spectrum {
public :
  int plate;
  int mjd;
  int fiber;
  double z;
  double ra; // deg
  double dec; // deg
  double xfocal;
  double yfocal;
  double xfocal_weight;
  Vect flux;
  Vect weight;
  Vect ivar;
  //Vect wave;
  Vect dir;
  //Vect zlya;
  //Vect dist;
  bool weights_have_been_modified;
  bool valid;
  Spectrum();
  void set_dir();
  //void compute_dist(const Vect& wave, const double& wave_lya, const EqualStepBSpline1D& spline);
};

using namespace std;

int get_column_number(const string& key, const map<string,int>& columns);
void read_data(const string& filename,Vect& wave,vector<Spectrum>& spectra,int format=1,int begin=0, int end=-1);
void read_dflux_lrg(const string& filename,Vect& wave,vector<Spectrum>& spectra, int format=0);
void read_dflux_lrg_quadrants(const string& filename,Vect& wave,vector<Spectrum>* spectra);
void rm_badplates(vector<Spectrum>& spectra,const vector<PlateMJD>& badplates);
void mask_CaII_HK_lines(const Vect& wave,vector<Spectrum>& spectra);
void select_plates(vector<Spectrum>& spectra, int plate_min, int plate_max);
void read_spall_or_drq(const string& filename,vector<QSO>& qsos, int table_hdu=2);
#endif
