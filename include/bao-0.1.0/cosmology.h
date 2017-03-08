#ifndef COSMOLOGY__H
#define COSMOLOGY__H

#include <iostream>


//! interface class
class Cosmology {
 protected :
    double omegak;
    double sqrt_abs_omegak;
    double (*Sin)(double x); // sin  or sinh or identity (k = -1,0,1)
    double (*Cos)(double x); // cos  or cosh or 1.       (k = -1,0,1)
    /* Radiation density, assumed to be zero by default. It matters
       when computing "distance to CMB". It is neglected when computing 
       curvature. */
    double omegar; 

 public:

    // constructor
    Cosmology(const double &Omegak);

    //!
    double Omegak() const { return omegak;}

    //!-
    virtual double  Omegam() const = 0;
    
    
  //! the Luminosity distance (in units of c/H0)
  virtual double Dl(const double z, int Nsteps = 100) const;
  
  //! the Angular distance (in units of c/H0)
  double Da(const double z, int Nsteps = 100) const;

  //! Proper motion distance
  double Dm(const double z,int Nsteps = 100) const;
  
  //! H(z)/H(z=0)
  virtual double E(const double z) const = 0;


  //! EoS of Dark Energy
  virtual double W(const double &Z) const = 0;


    //! Differential volume
  double dVdz(const double z) const;

    //! volume from 0 to z
  double Volume(double z) const;

  //! Number of parameters (useful because there are derived classes)
  virtual int NPar() const = 0;

  //! The gradient of the (natural Log of) luminosity distance (w.r.t cosmological parameters)
  void GradLogDl(const double z, double *grad, int Nsteps = 100) const;
  

  virtual double integ_hinv(const double &ZMin, const double &ZMax, const int NStep=100) const = 0;

  virtual void dump(std::ostream &stream) const = 0;

  virtual Cosmology* Clone() const = 0; 
  

  //! set parameter number i
  virtual void SetParam(const int IPar, const double Val) = 0;

  virtual double GetParam(const int IPar) const = 0;


  void SetOmegar(const double &Val);

  double Omegar() const { return omegar;}

  //! enables cout << MyCosmology
  friend std::ostream & operator << ( std::ostream &stream, const Cosmology &This);  
  virtual ~Cosmology() {};
};


//! 4 params cosmology : omegam, omegax , w_0 and w_1
class GeneralCosmo : public Cosmology {
  private :
    double omegam;
    double omegax;
    double w0;
    double w1;
    friend class hinv; // the luminosity distance involves the integral of 1/H(z)

  public :
    //! The constructor :
    GeneralCosmo(const double Omegam, const double OmegaX, const double W0, const double W1);
  //!-
    double  Omegax() const { return omegax;}

    //!
    double  Omegam() const { return omegam;}
  //!-
    double  W0() const { return w0;}
  //!-
    double  W1() const { return w1;}

    double W(const double &Z) const { return w0+w1*Z;}

    void dump(std::ostream &stream) const;

    void SetParam(const int i, const double Val);
    double GetParam(const int IPar) const;

    double E(const double z) const;

    virtual int NPar() const {return 4;}

    double integ_hinv(const double &ZMin, const double &ZMax, const int NStep=100) const;

    virtual Cosmology* Clone() const {return new GeneralCosmo(*this);};

};


// ! classic one : omegam and omegal
class OmegamOmegalCosmo : public GeneralCosmo {
 public:
  //! Constructor
  OmegamOmegalCosmo();
  OmegamOmegalCosmo(const double omegam, const double omegal);
  
  //! -
  void SetParam(const int i, const double Val);
  double GetParam(const int IPar) const;
  int NPar() const {return 2;}
  void dump(std::ostream &stream) const;
  virtual Cosmology* Clone() const {return new OmegamOmegalCosmo(*this);};

};


//! flat LambdaCDM cosmo
class FlatLCDMCosmo : public GeneralCosmo {
 public:
  //! Constructor
  FlatLCDMCosmo();
  FlatLCDMCosmo(const double omegam);
  
  //! -
  void SetParam(const int i, const double Val);
  double GetParam(const int IPar) const;
  int NPar() const {return 1;}
  void dump(std::ostream &stream) const;
  virtual Cosmology* Clone() const {return new FlatLCDMCosmo(*this);};

};



//! less classic : flat with omegam and w0
class OmegamW0Cosmo : public GeneralCosmo {
  public : 
    //! Constructor.
    OmegamW0Cosmo(const double Omegam, const double W0);
    //! -
    void SetParam(const int i, const double Val);
    double GetParam(const int IPar) const;

    int NPar() const {return 2;}
    void dump(std::ostream  &stream) const;
    virtual Cosmology* Clone() const {return new OmegamW0Cosmo(*this);};
};



//! less classic : omegam, omegax and w0
class OmegamOmegaXW0Cosmo : public GeneralCosmo {
  public : 
    //! Constructor.
    OmegamOmegaXW0Cosmo(const double Omegam, const double OmegaX, const double W0);
    //! -
    void SetParam(const int i, const double Val);
    double GetParam(const int IPar) const;

    int NPar() const {return 3;}
    void dump(std::ostream  &stream) const;
    virtual Cosmology* Clone() const {return new OmegamOmegaXW0Cosmo(*this);};
};


//! flat cosmology with w(z) = w0+w1*z
class OmegamW0W1Cosmo : public GeneralCosmo {
  public : 
    //! Constructor.
    OmegamW0W1Cosmo(const double Omegam, const double OmegaX, const double W0);
    //! -
    void SetParam(const int i, const double Val);
    double GetParam(const int IPar) const;

    int NPar() const {return 3;}
    void dump(std::ostream  &stream) const;
    virtual Cosmology* Clone() const {return new OmegamW0W1Cosmo(*this);};
};



//! Om and Ox, with eos = w0+wa*z/(1+z) 
class OmegamOmegaxW0WaCosmo : public Cosmology {
  protected :

    friend class W0Wahinv ;
  double omegam;
  double omegax;
  double w0 ;
  double wa;

  public :
    
    OmegamOmegaxW0WaCosmo(const double OmegaM, const double OmegaX, 
	      const double W0, const double Wa) ;
  double Omegam() const { return omegam;}
  double Omegax() const { return omegax;}
  double W0() const { return w0;}
  double Wa() const { return wa;}
  double W(const double &Z) const { return w0+wa*Z/(1+Z);}

  void SetParam(const int i, const double Val);
  double GetParam(const int IPar) const;

  void dump(std::ostream &stream) const ;
  virtual double E(const double z) const;

  virtual int NPar() const {return 4;} // omegar is considered to be known.
  virtual double integ_hinv(const double &ZMin, const double &ZMax, 
			    const int NStep=100) const;
  Cosmology *Clone() const { return new OmegamOmegaxW0WaCosmo(*this);}



};


//! flat Om,OX ( 1-Om) with eos  = w0+wa*z/(1+z)
class OmegamW0WaCosmo : public OmegamOmegaxW0WaCosmo {
public :
  //! Constructor.
  OmegamW0WaCosmo(const double Omegam, const double W0, const double Wa) ;
  //! -
  void SetParam(const int i, const double Val);
  double GetParam(const int IPar) const;

  int NPar() const {return 3;}
  void dump(std::ostream  &stream) const;
  virtual Cosmology* Clone() const {return new OmegamW0WaCosmo(*this);}
};


/////////////////////////////////////////////////////////////////////////////////////////////
#ifdef STORAGE

class HuCosmo : public Cosmology {
  protected :
    friend class Huhinv ;
  double omegam ;
  double omegax ;
  double wn ;
  double wa ;
  double zn ;
  public :
    HuCosmo(const double OmegaM, const double OmegaX, const double Wn, const double Wa, const double Zn) ;
  double Omegam() const { return omegam;}
  double Omegax() const { return omegax;}
  double Wn() const { return wn;}
  double Wa() const { return wa;}
  double Zn() const { return zn;}
  void SetParam(const int i, const double Val);
  double GetParam(const int IPar) const;

  void dump(std::ostream &stream) const ;
  virtual double E(const double z) const;
  virtual int NPar() const {return 4;}
  virtual double integ_hinv(const double &ZMin, const double &ZMax, 
			    const int NStep=100) const;
  Cosmology *Clone() const { return new HuCosmo(*this);}
};

class OmegamWnWaCosmo : public HuCosmo {
public :
  //! Constructor.
  OmegamWnWaCosmo(const double Omegam, const double Wn, const double Wa, const double Zn) ;
  //! -
  void SetParam(const int i, const double Val);
  double GetParam(const int IPar) const;
  int NPar() const {return 3;}
  void dump(std::ostream  &stream) const;
  virtual Cosmology* Clone() const {return new OmegamWnWaCosmo(*this);}
};

#endif 


/////////////////////////////////////////////////////////////////////////////////////////////



#endif /* COSMOLOGY__H */
