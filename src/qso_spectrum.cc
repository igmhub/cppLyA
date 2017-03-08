#include <qso_spectrum.h>
#include <cmath>

//#include <TMath.h>

#include <fitsio.h>
#define CHECKERROR if(status) {fits_report_error(stdout, status); cerr << "fits error" << endl; exit(12);}  

Spectrum::Spectrum() {
  weights_have_been_modified=false;
}

Vect compute_dir(const double& ra_deg, const double& dec_deg) {
  double deg2rad=M_PI/180.;
  Vect dir(3); 
  dir(0)=cos(dec_deg*deg2rad)*cos(ra_deg*deg2rad);
  dir(1)=cos(dec_deg*deg2rad)*sin(ra_deg*deg2rad);
  dir(2)=sin(dec_deg*deg2rad);
  return dir;
}

void Spectrum::set_dir() {
  dir=compute_dir(ra,dec);
}

void QSO::set_dir() {
  dir=compute_dir(ra,dec);
}

int get_column_number(const string& key, const map<string,int>& columns) {
  if(columns.find(key)==columns.end()) {
    //cout << "ERROR no key '" << key << "' in table" << endl;
    //exit(12);
    cout << "WARNING no key '" << key << "' in table" << endl;
    return -1;
  }
  return columns.find(key)->second;
}
/*
FORMAT 0
Filename: dflux-all-qso.fits
No.    Name         Type      Cards   Dimensions   Format
0    PRIMARY     PrimaryHDU       6   (683, 96877)   float64   
1    IVAR        ImageHDU         8   (683, 96877)   float64   
2    WAVELENGTH  ImageHDU         7   (683,)       float64   
3                BinTableHDU     22   96877R x 7C   [D, D, D, J, J, J, J]   
4    COEFS       ImageHDU         8   (2, 96877)   float64   

FORMAT 1
Filename: selection4-fake-001.fits
No.    Name         Type      Cards   Dimensions   Format
0    PRIMARY     PrimaryHDU       7   (683, 11282)   float64   
1    FAKEIVAR    ImageHDU         8   (683, 11282)   float64   
2    QSOIVAR     ImageHDU         8   (683, 11282)   float64   
3    SKYIVAR     ImageHDU         8   (683, 11282)   float64   
4    WAVELENGTH  ImageHDU         7   (683,)       float64   
5                BinTableHDU     26   11282R x 9C   [D, D, D, D, D, J, J, J, J]

FORMAT 2
Filename: dflux-fake.fits
No.    Name         Type      Cards   Dimensions   Format
0    PRIMARY     PrimaryHDU       7   (683, 9350)   float64   
1    QSOIVAR     ImageHDU         8   (683, 9350)   float64   
2    SKYIVAR     ImageHDU         8   (683, 9350)   float64   
3    WAVELENGTH  ImageHDU         7   (683,)       float64   
4                BinTableHDU     26   9350R x 9C   [D, D, D, D, D, J, J, J, J]   
5    COEFS       ImageHDU         8   (2, 9350)    float64  

FORMAT 4
Filename: coadd.fits
No.    Name         Type      Cards   Dimensions   Format
0    PRIMARY     PrimaryHDU       7   (4128, 500)   float32   
1    IVAR        ImageHDU         9   (4128, 500)   float32   
2    WAVELENGTH  ImageHDU         9   (4128,)      float32   

*/   
void read_data(const string& filename,Vect& wave,vector<Spectrum>& spectra,int format,int begin, int end) {
  Mat mflux,mqsoivar,mskyivar,mwave;
  
  int flux_hdu,qsoivar_hdu,skyivar_hdu,wave_hdu,table_hdu;

  //bool angles_in_rad = false;

  if(format==0) {

    flux_hdu=1;
    qsoivar_hdu=2;
    skyivar_hdu=-1;
    wave_hdu=3;
    table_hdu=4;

  }else if(format==1) {
    
    flux_hdu=1;
    qsoivar_hdu=3;
    skyivar_hdu=4;
    wave_hdu=5;
    table_hdu=6;
    
  }else if(format==2) {
    
    flux_hdu=1;
    qsoivar_hdu=2;
    skyivar_hdu=3;
    wave_hdu=4;
    table_hdu=5;
    
  }else if(format==4) {
    
    flux_hdu=1;
    qsoivar_hdu=2;
    skyivar_hdu=-1;
    wave_hdu=3;
    table_hdu=-1;
    
  }else{
    cout << "unknown format" << endl;
    exit(12);
  }
  
  bool has_skyivar=(skyivar_hdu>0);

  cout << "mflux.readFits" << endl;
  mflux.readFits(filename,flux_hdu);
  cout << "mqsoivar.readFits" << endl;
  mqsoivar.readFits(filename,qsoivar_hdu);
  if(has_skyivar) {
    cout << "mskyivar.readFits" << endl;
    mskyivar.readFits(filename,skyivar_hdu);
  }
  cout << "mwave.readFits" << endl;
  mwave.readFits(filename,wave_hdu);
  cout << "ok, done reading images" << endl;
  
  cout << mflux.SizeX() << " " << mflux.SizeY() << endl;
  
  int spec_begin = begin;
  int spec_end   = end;
  
  int nspec_in_file=mflux.SizeY();
  if((spec_end<0) || (spec_end>nspec_in_file)) spec_end=nspec_in_file;
  
  int nwave=mflux.SizeX();
  wave.allocate(nwave);
  for(int i=0;i<nwave;i++) { 
    wave(i)=mwave(i,0);
  }

  


  int nempty=0;
  
  for(int s=spec_begin;s<spec_end;s++) {
    
    double sumw=0;
    for(int i=0;i<nwave;i++)  {
      sumw+=mqsoivar(i,s);
    }
    

    Spectrum spec;
    

    if(sumw==0) {
    
      nempty++;
      spec.valid=false;
      
    }else{
      spec.valid=true;
    
      spec.flux.allocate(nwave);
      spec.weight.allocate(nwave);
      if(has_skyivar)
	spec.ivar.allocate(nwave);
      
      for(int i=0;i<nwave;i++)  {
	spec.flux(i)=mflux(i,s);
	spec.weight(i)=mqsoivar(i,s);
	if(has_skyivar)
	  spec.ivar(i)=mskyivar(i,s);
      }
    }
    
    spectra.push_back(spec);
      
  }
  cout << "ntot=" << spectra.size() << " nempty=" << nempty << endl;
  
  
  // reading the table (harder)
  if(table_hdu>0) {
  fitsfile *fptr = 0;
  int status = 0;
  fits_open_file(&fptr, filename.c_str() , READONLY, &status);
  CHECKERROR;
  fits_movabs_hdu(fptr, table_hdu, NULL, &status);
  CHECKERROR;
  long nrows;
  int ncols;
  fits_get_num_rows(fptr, &nrows, &status); CHECKERROR; 
  fits_get_num_cols(fptr, &ncols, &status); CHECKERROR;

  

  // description
  vector<string> keys;
  map<string,int> columns;
  for(int c=0;c<ncols;c++) {
    char key[12];
    sprintf(key,"TTYPE%d", c+1);
    char a_C_string[80];
    fits_read_key(fptr, TSTRING, key, a_C_string, NULL, &status); CHECKERROR;
    keys.push_back(string(a_C_string));
    columns[string(a_C_string)]=c;
    
  }
  cout << "keys:";
  for(size_t k=0;k<keys.size();k++) {
    cout << " " << keys[k];
  }
  cout << endl;

  int ra_col=get_column_number("RA",columns);
  int dec_col=get_column_number("DEC",columns);
  //int xfocal_col=get_column_number("XFOCAL",columns);
  //int yfocal_col=get_column_number("YFOCAL",columns);
  int z_col=get_column_number("Z",columns);
  int plate_col=get_column_number("PLATE",columns);
  int fiber_col=get_column_number("FIBER",columns);
  int mjd_col=get_column_number("MJD",columns);
  
  // loop on rows
  double nullval = 0.;
  int anynul;
  int nvals=1;
  double value;
  //double ra2deg=180./3.14159265358979312;
  
  for(int r=spec_begin; r<min(nrows,(long)(spec_end)); r++) {
    
    Spectrum& spec = spectra[r-spec_begin];
    
    for(int c=0; c<ncols; c++) {
      fits_read_col(fptr, TDOUBLE, c+1, long(r+1), 1,nvals, &nullval, &value, &anynul, &status); CHECKERROR; 
      if(c==ra_col) {
	spec.ra=value;
	continue;
      }
      if(c==dec_col) {
	spec.dec=value;
	continue;
	}/*
      if(c==xfocal_col) {
	spec.xfocal=value;
	continue;
      }

      if(c==yfocal_col) {
    	spec.yfocal=value;
    	continue;
	}*/
      if(c==z_col) {
	spec.z=value;
	continue;
      }
      if(c==plate_col) {
	spec.plate=int(value);
	continue;
      }
      if(c==mjd_col) {
	spec.mjd=int(value);
	continue;
      }
      if(c==fiber_col) {
	spec.fiber=int(value);
	continue;
      }
    }
    
  }
  }else{
    // add fake z and plate
    for(size_t s=0;s<spectra.size();s++) {
      spectra[s].z=3.5;
      spectra[s].plate=s;
    }
  }
}

void read_dflux_lrg(const string& filename,Vect& wave,vector<Spectrum>& spectra, int format) {
  
  Mat m_plate,m_mjd,m_wave,m_dflux1,m_ivar1,m_dflux2,m_ivar2;
  if(format==0) {
    m_plate.readFits(filename,1);
    m_mjd.readFits(filename,2);
    m_wave.readFits(filename,3);
    m_dflux1.readFits(filename,4);
    m_ivar1.readFits(filename,5);
    // skip 6
    m_dflux2.readFits(filename,7);
    m_ivar2.readFits(filename,8);
  } else if(format==1) {
    m_plate.readFits(filename,1);
    m_mjd.readFits(filename,2);
    m_wave.readFits(filename,3);
    m_dflux1.readFits(filename,4);
    m_ivar1.readFits(filename,5);
    m_dflux2.readFits(filename,6);
    m_ivar2.readFits(filename,7);
  }else{
    cout << "cannot deal with format " << format << endl;
    exit(12);
  }
  
  int nspec=m_dflux1.SizeY();
  //int nwave_in_file=m_dflux1.SizeX();
  
  // rebin data to match a first wave at 3600
  // original bins are 0.0001 in log lambda
  // we want to go up to 1200*(1+3.44)=5328A
  // nbins=log10(5328/3600.)/0.0001=1702
  

  double wave0=3600;

  double log10_input_bin=0.0001;
  int rebin=3;
  double log10_output_bin=log10_input_bin*rebin;
  
  int nwave_rebinned=int(log10(5328/3600.)/log10_output_bin);
  wave.allocate(nwave_rebinned);
  for(int i=0;i<nwave_rebinned;i++) { 
    wave(i)=wave0*pow(10.,log10_output_bin*i);
  }
  
  // 
  for(int s=0;s<nspec;s++) {
    
    Spectrum spec;
    spec.valid=true;
    spec.flux.allocate(nwave_rebinned);
    spec.weight.allocate(nwave_rebinned);
    spec.ivar.allocate(nwave_rebinned);
    
    int index_in_rebinned=int(floor(log10(m_wave(0,s)/wave0)/log10_output_bin));
    int begin_in_rebinned=max(0,index_in_rebinned);
    int end_in_rebinned=nwave_rebinned;
    
    //cout << "DEBUG " << wave0 << " " << m_wave(0,s) << " "  << index_in_rebinned << " " << begin_in_rebinned << " " << begin_in_rebinned-index_in_rebinned << endl;

    for(int i=begin_in_rebinned;i<end_in_rebinned;i++)  {
      double w=0;
      double wf=0;
      for(int k=0;k<rebin;k++) {
	w +=  m_ivar1(i*rebin+k-index_in_rebinned*rebin,s);
	wf += m_dflux1(i*rebin+k-index_in_rebinned*rebin,s)*m_ivar1(i*rebin+k-index_in_rebinned*rebin,s);
      }
      if(w>0){
	spec.flux(i)=wf/w;
	spec.weight(i)=w;
	spec.ivar(i)=w;
      }
    }
    //spec.mjd=m_mjd(0,s);
    spec.plate=int(m_plate(s,0));
    spec.fiber=0; // first spectro
    spectra.push_back(spec);

    spec.flux*=0;
    spec.weight*=0;
    spec.ivar*=0;
    for(int i=begin_in_rebinned;i<end_in_rebinned;i++)  {
      double w=0;
      double wf=0;
      for(int k=0;k<rebin;k++) {
	w +=  m_ivar2(i*rebin+k-index_in_rebinned*rebin,s);
	wf += m_dflux2(i*rebin+k-index_in_rebinned*rebin,s)*m_ivar2(i*rebin+k-index_in_rebinned*rebin,s);
      }
      if(w>0){
	spec.flux(i)=wf/w;
	spec.weight(i)=w;
	spec.ivar(i)=w;
      }
    }
    
    spec.fiber=500; // sec. spectro
    //cout << "IGNORE DFLUX2" << endl;
    spectra.push_back(spec);
  }
  cout << "nok=" << spectra.size() << endl;
}

void read_dflux_lrg_quadrants(const string& filename,Vect& wave,vector<Spectrum>* spectra) {
  
  Mat m_plate,m_mjd,m_wave,m_dflux[4],m_ivar[4];
  
  m_plate.readFits(filename,1);
  m_mjd.readFits(filename,2);
  m_wave.readFits(filename,3);
  for(int q=0;q<4;q++) {
    m_dflux[q].readFits(filename,4+q*2);
    m_ivar[q].readFits(filename,4+q*2+1);
  }
  
  int nspec=m_dflux[0].SizeY();
  //int nwave_in_file=m_dflux1.SizeX();
  
  // rebin data to match a first wave at 3600
  // original bins are 0.0001 in log lambda
  // we want to go up to 1200*(1+3.44)=5328A
  // nbins=log10(5328/3600.)/0.0001=1702
  
  double wave0=3600;
  double log10_input_bin=0.0001;
  int rebin=1;
  double log10_output_bin=log10_input_bin*rebin;
  int nwave_rebinned=int(log10(5328/3600.)/log10_output_bin);
  wave.allocate(nwave_rebinned);
  for(int i=0;i<nwave_rebinned;i++) { 
    wave(i)=wave0*pow(10.,log10_output_bin*i);
  }
  
  // 
  for(int s=0;s<nspec;s++) {
    
    Spectrum spec;
    spec.valid=true;
    spec.flux.allocate(nwave_rebinned);
    spec.weight.allocate(nwave_rebinned);
    spec.ivar.allocate(nwave_rebinned);
    
    int index_in_rebinned=int(floor(log10(m_wave(0,s)/wave0)/log10_output_bin));
    int begin_in_rebinned=max(0,index_in_rebinned);
    int end_in_rebinned=nwave_rebinned;
    
    //cout << "DEBUG " << wave0 << " " << m_wave(0,s) << " "  << index_in_rebinned << " " << begin_in_rebinned << " " << begin_in_rebinned-index_in_rebinned << endl;
    spec.plate=int(m_plate(s,0));
    for(int q=0;q<4;q++) {
      spec.flux*=0;
      spec.weight*=0;
      spec.ivar*=0;
      for(int i=begin_in_rebinned;i<end_in_rebinned;i++)  {
	double w=0;
	double wf=0;
	for(int k=0;k<rebin;k++) {
	  w+=m_ivar[q](i*rebin+k-index_in_rebinned*rebin,s);
	  wf+=m_dflux[q](i*rebin+k-index_in_rebinned*rebin,s)*m_ivar[q](i*rebin+k-index_in_rebinned*rebin,s);
	}
	if(w>0) {
	  spec.flux(i)=(wf/w);
	  spec.weight(i)=w;
	  spec.ivar(i)=w;
	}
      }
      spec.fiber=250*q;
      spectra[q].push_back(spec);
    } 
  } // end of loop on plates
}

void rm_badplates(vector<Spectrum>& spectra,const vector<PlateMJD>& badplates) {
  
  int n_bad=0;
  cout << "start rm_badplates ..." << endl;
  for(size_t s=0;s<spectra.size();s++) {
    int plate=spectra[s].plate;
    for(size_t p=0;p<badplates.size();p++) {
      if(badplates[p].plate == plate) {
	spectra[s].valid=false;
	n_bad++;
	break;
      }
    }
  }
  cout << "removed " << n_bad << " quasars in bad plates" << endl;
 
}
void mask_CaII_HK_lines(const Vect& wave,vector<Spectrum>& spectra) {
  
  
  // CaII HK in air = 3933.66 3968.47 (NIST)
  // CaII HK in vacuum = 3934.77373115  3969.59281349 ( using specex_air_to_vacuum )
  
  int begin[2];
  int end[2];
  double line_wave[2];
  line_wave[0]=3934.77373115;
  line_wave[1]=3969.59281349;
  double hw=5; // mask 10A in total
  int nw=int(wave.size());
  for(int line=0;line<2;line++) {
    begin[line]=nw-1;
    end[line]=0;
    for(int i=0;i<nw;i++) {
      if(wave(i)>line_wave[line]-hw) {begin[line]=i; break;}
    }
    for(int i=nw-1;i>=0;i--) {
      if(wave(i)<line_wave[line]+hw) {end[line]=i+1; break;}
    }
    cout << line_wave[line] << " masking= " << begin[line] << " " << end[line] << endl;
  }
  
  for(size_t s=0;s<spectra.size();s++) {
    Spectrum& spec=spectra[s];
    if(spec.valid)
      for(int line=0;line<2;line++) {
	for(int i=begin[line];i<end[line];i++) {
	  spec.weight(i)=0;
	  spec.ivar(i)=0;
	}
      }
  }
  cout << "masked CaII H&K lines " << endl;
 
}

void select_plates(vector<Spectrum>& spectra, int plate_min, int plate_max) {
   for(size_t s=0;s<spectra.size();s++) {
     Spectrum& spec=spectra[s];
     if(spec.plate<plate_min || spec.plate>plate_max) {
       cout << "discarding QSO because plate = " << spec.plate << " not in range " << plate_min << " " << plate_max << endl;
       spec.valid=false;

     }
   }
}


void read_spall_or_drq(const string& filename,vector<QSO>& qsos, int table_hdu) {
  cout << "read_spall_or_drq in " << filename << endl;
  fitsfile *fptr = 0;
  int status = 0;
  fits_open_file(&fptr, filename.c_str() , READONLY, &status);
  CHECKERROR;
  fits_movabs_hdu(fptr, table_hdu, NULL, &status);
  CHECKERROR;
  long nrows;
  int ncols;
  fits_get_num_rows(fptr, &nrows, &status); CHECKERROR; 
  fits_get_num_cols(fptr, &ncols, &status); CHECKERROR;
  // description
  vector<string> keys;
  map<string,int> columns;
  for(int c=0;c<ncols;c++) {
    char key[12];
    sprintf(key,"TTYPE%d", c+1);
    char a_C_string[80];
    fits_read_key(fptr, TSTRING, key, a_C_string, NULL, &status); CHECKERROR;
    keys.push_back(string(a_C_string));
    columns[string(a_C_string)]=c;
    
  }
  cout << "keys:";
  for(size_t k=0;k<keys.size();k++) {
    cout << " " << keys[k];
  }
  cout << endl;

  int ra_col=get_column_number("RA",columns);
  int dec_col=get_column_number("DEC",columns);
  //int xfocal_col=get_column_number("XFOCAL",columns);
  //int yfocal_col=get_column_number("YFOCAL",columns);
  int z_col=get_column_number("Z_VI",columns);
  int plate_col=get_column_number("PLATE",columns);
  int fiber_col=get_column_number("FIBERID",columns);
  int mjd_col=get_column_number("MJD",columns);
  
  

  // loop on rows
  double nullval = 0.;
  int anynul;
  int nvals=1;
  double value;
  
  for(int r=0; r<nrows; r++) {
    
    QSO qso;
    
    for(int c=0; c<ncols; c++) {

      if( ! ( c==ra_col || c==dec_col 
	      || c==z_col 
	      || c==plate_col 
	      || c==mjd_col 
	      || c==fiber_col
	      //|| c==xfocal_col 
	      //|| c==yfocal_col
	      )) continue;
      
      fits_read_col(fptr, TDOUBLE, c+1, long(r+1), 1,nvals, &nullval, &value, &anynul, &status); CHECKERROR; 
      if(c==ra_col) {
	qso.ra=value;
	continue;
      }
      if(c==dec_col) {
	qso.dec=value;
	continue;
	}/*
      if(c==xfocal_col) {
	qso.xfocal=value;
	continue;
      }
      if(c==yfocal_col) {
    	qso.yfocal=value;
    	continue;
	}*/
      if(c==z_col) {
	qso.z=value;
	continue;
      }
      if(c==plate_col) {
	qso.plate=int(value);
	continue;
      }
      if(c==mjd_col) {
	qso.mjd=int(value);
	continue;
      }
      if(c==fiber_col) {
	qso.fiber=int(value);
	continue;
      }
    }
    qsos.push_back(qso);
  } // end of loop on rowss
  cout << "read " << qsos.size() << " qso coordinates in " << filename << endl;
}
