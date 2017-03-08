#include <cmath>
#include "cosmology.h"
#include "fatalerror.h"
#include "snfitexception.h"
#include <stdlib.h>

using namespace std;

// Reynald appelle ca l'integrale de Riemann:
template <class Function>  double integ(const Function &F, const double Xmin, const double Xmax, const int Nsteps = 100)
{
  double sum = 0;
  double xstep = (Xmax-Xmin)/double(Nsteps);
  double x = Xmin + 0.5*xstep;
  for (int i=0; i< Nsteps; ++i)
    { 
      sum += F(x);
      x+= xstep;
    }
  return sum*xstep;
}



ostream& operator << ( ostream &stream, const Cosmology &This) 
{ This.dump(stream); return stream;}

static double ident(const double x) { return x;}
static double one(const double x) { return 1;}
static double cub(double x) { return x*x*x;}
static double sq(double x) { return x*x;}

Cosmology::Cosmology(const double &Omegak) : omegak(Omegak)
{
  sqrt_abs_omegak = sqrt(fabs(omegak));
  if (omegak > 0.0001)
    {
      Sin = &sinh;
      Cos = &cosh;
    }
  else if (omegak < -0.0001)
    {
      Sin = &sin;
      Cos = &cos;
    }
  else
    {
      Sin = &ident;
      Cos = &one;
      sqrt_abs_omegak = 1.;
    }
  omegar = 0;
}


double Cosmology::Dm(const double z,int Nsteps) const
{
  double tmp=sqrt_abs_omegak*integ_hinv(0,z,Nsteps);
  if(omegak > 0.0001 && tmp>100) tmp=100;  // it's crazy anyway
  return Sin(tmp)/sqrt_abs_omegak; 
}

double Cosmology::Dl(const double z,int Nsteps) const
{
  return (1+z)*Dm(z,Nsteps);
}

double Cosmology::Da(const double z,int Nsteps) const
{
  return Dm(z,Nsteps)/(1+z);
}

double Cosmology::dVdz(const double z) const
{
  /* still to be checked against CPT 92 curves */
  /* checked "experimentally" however that this routine is the derivative
     of Volume(z) */
  double integral = integ_hinv(0,z);
  double dm = Sin(sqrt_abs_omegak*integral)/sqrt_abs_omegak; 
  double ddmdz = Cos(sqrt_abs_omegak*integral) / E(z);
  return dm*dm*ddmdz/sqrt(1+omegak*sq(dm));
}
  
//! from 0 to z. it is in fact H_0^3 V
double Cosmology::Volume(double z) const
{
  double integral = integ_hinv(0,z);
  double dm = Sin(sqrt_abs_omegak*integral)/sqrt_abs_omegak; 
  double x = omegak*sq(dm);
  /* I checked experimentaly the continuity over omegak = 0 */
  if (fabs(x) > 1e-3)
    {
      // use exact (however unstable) expression (CPT eq 27.)
      return 0.5*(dm*sqrt(1+x) - integral)/omegak;
    }
  else
    {
      // use taylor expansion
      return cub(dm)*(1 - 0.3*x + 9.*sq(x)/56.- 5.*cub(x)/48.)/3.;
    }
}
  

/*
> a:=int(x^2/sqrt(1+ok*x^2),x=0..z);
                           2 1/2   3/2        1/2              2 1/2
                z (1 + ok z )    ok    - ln(ok    z + (1 + ok z )   ) ok
       a := 1/2 --------------------------------------------------------
                                           5/2
                                         ok
 
(Equivalent to Carrol Press Turner eq 27)

> series(a,ok=0);
                3         5            7   2          9   3       4
           1/3 z  - 1/10 z  ok + 3/56 z  ok  - 5/144 z  ok  + O(ok )
 
>                                                                               
*/



void Cosmology::GradLogDl(const double z, double *grad, 
			  int Nsteps) const
{
  Cosmology *copy = Clone();
  double dl0 = Dl(z, Nsteps);
  double eps = 0.001;
  for (int i=0; i< NPar(); ++i)
    {
      if (dl0 == 0) {grad[i] = 0; continue;}
      double pval = GetParam(i);
      copy->SetParam(i, pval+eps);
      double dl1 = copy->Dl(z,Nsteps);
      copy->SetParam(i, pval);
      grad[i] = (dl1-dl0)/(eps*dl0);
    }
  delete copy;
}

void Cosmology::SetOmegar(const double &Val)
{
  if (Val>0.01)
    throw(SnfitException("Cosmology::SetOmegaR : this class does not handle properly large radiation densities"));
  // should recompute omegak (i.e. go through the Cosmology constructor)
  omegar = Val;
}




GeneralCosmo::GeneralCosmo(const double OmegaM, const double OmegaX, 
			   const double W0, const double W1)
  : Cosmology(1-OmegaM- OmegaX), omegam(OmegaM), omegax (OmegaX), 
    w0(W0), w1(W1) 
{
}


class hinv {
public : 
  const GeneralCosmo &par;
  hinv(const GeneralCosmo &C) : par(C) {};
  double operator()(const double z) const
  {
    
    double h2 = par.omegam*cub(1.+z) + par.omegax*
    exp(3.*(par.w1*z + (1.+ par.w0-par.w1)*log(1.+z)))
      + sq(1+z)*par.omegak;
    return 1./sqrt(h2);
  }
};


double GeneralCosmo::integ_hinv(const double &ZMin, const double &ZMax, const int NStep) const
{
  // change z to u=1/sqrt(1+z)
  double umin = 1./sqrt(1+ZMax);
  double umax = 1./sqrt(1+ZMin);
  double step = (umax-umin)/NStep;
  double u = umin+0.5*step;
  double sum =0;
  for (int i=0; i<NStep; ++i, u+=step)
    {
      double u2 = sq(u);
      double fact_ox = exp(3*(-2*log(u)*(w0-w1)+w1*(1./u2-1)));
      double u3hinv = sqrt(omegam + omegax*fact_ox + omegak*u2 + omegar/u2);
      sum += 1./u3hinv;
    }
  return 2.*sum*step;
}
  

double GeneralCosmo::E(const double z) const
{
  // cannot use easily the routine above, because of the change of variable
  return sqrt(omegam*cub(1.+z) + omegax*
	      exp(3.*(w1*z + (1.+ w0-w1)*log(1.+z)))
	      + sq(1+z)*(omegak + sq(1+z)*omegar));
}

void GeneralCosmo::SetParam(const int i, const double Val)
{
  double toto = Omegar();
  if (i==0) *this = GeneralCosmo(Val, Omegax(), W0(),  W1()      );
  else if (i==1) *this = GeneralCosmo(Omegam(), Val,  W0(), W1()  );
  else if (i==2) this->w0 = Val;
  else if (i==3) this->w1 = Val;
  else FatalError(" wrong param number in SetParam");
  SetOmegar(toto);
}

double GeneralCosmo::GetParam(const int i) const
{
  if (i==0) return Omegam();
  else if (i==1) Omegax();
  else if (i==2) return this->w0;
  else if (i==3) return this->w1;
  else FatalError(" wrong param number in GetParam");
  return 42; // rever reached
}



void GeneralCosmo::dump(ostream &stream) const
{
  stream << " om-ox-w0-w1 cosmo : (" << omegam << ',' << omegax << ',' << w0 << ',' << w1 << ')' << endl;
}


//=============================== OmegamOmegalCosmo =========================
  
// Here, choose the cosmology 
OmegamOmegalCosmo::OmegamOmegalCosmo()
//  : GeneralCosmo(0.28,0.72, -1, 0) WMAP 
  : GeneralCosmo(0.3147,1.-0.3147, -1, 0) // Planck  
{}



OmegamOmegalCosmo::OmegamOmegalCosmo(const double omegam, const double omegal )
  : GeneralCosmo(omegam, omegal, -1, 0)
{}

void OmegamOmegalCosmo::dump(ostream &stream) const
{
  stream << " Om-Ol cosmo : (" << Omegam() << ',' << Omegax() << ')' << endl;
}


void OmegamOmegalCosmo::SetParam(const int i, const double Val)
{
  double toto = Omegar();
  if (i==0) *this = OmegamOmegalCosmo(Val, Omegax()      );
  else if (i==1) *this = OmegamOmegalCosmo(Omegam(), Val);  
  else FatalError(" wrong param number in SetParam");
  SetOmegar(toto);
}

double OmegamOmegalCosmo::GetParam(const int i) const
{
  if (i==0) return Omegam();
  else if (i==1) return Omegax();
  else FatalError(" wrong param number in GetParam");
  return 42; // rever reached
}

//=============================== FlatLCDMCosmo =========================
  
// Here, choose the cosmo 
FlatLCDMCosmo::FlatLCDMCosmo()
//  : GeneralCosmo(0.28,0.72, -1, 0) WMAP 
  : GeneralCosmo(0.3147,1.-0.3147, -1, 0) // Planck  
{}



FlatLCDMCosmo::FlatLCDMCosmo(const double omegam)
  : GeneralCosmo(omegam, 1-omegam, -1, 0)
{}

void FlatLCDMCosmo::dump(ostream &stream) const
{
  stream << " FLCDM cosmo : (" << Omegam() <<  ')' << endl;
}


void FlatLCDMCosmo::SetParam(const int i, const double Val)
{
  double toto = Omegar();
  if (i==0) *this = FlatLCDMCosmo(Val);
  else FatalError(" wrong param number in SetParam");
  SetOmegar(toto);
}

double FlatLCDMCosmo::GetParam(const int i) const
{
  if (i==0) return Omegam();
  else FatalError(" wrong param number in GetParam");
  return 42; // rever reached
}



//================================ OmegamW0 ===========================

OmegamW0Cosmo::OmegamW0Cosmo(const double Omegam, const double W0) 
  :GeneralCosmo(Omegam, 1.-Omegam, W0, 0.0)
{}


void OmegamW0Cosmo::dump(ostream &stream) const
{
  stream << " Om-w0 cosmo : (" << Omegam() << ',' << W0() << ')' << endl;
}

void OmegamW0Cosmo::SetParam(const int i, const double Val)
{
  double toto = Omegar();
  if (i==0) *this = OmegamW0Cosmo(Val,       W0());
  else if (i==1) *this = OmegamW0Cosmo(Omegam() , Val );
  else FatalError(" wrong param number in SetParam");
  SetOmegar(toto);
}

double OmegamW0Cosmo::GetParam(const int i) const
{
  if (i==0) return Omegam();
  else if (i==1) return W0();
  else FatalError(" wrong param number in GetParam");
  return 42; // rever reached
}



//================================ OmegamOmegaXW0 ===========================

OmegamOmegaXW0Cosmo::OmegamOmegaXW0Cosmo(const double Omegam, const double OmegaX, const double W0) 
  :GeneralCosmo(Omegam, OmegaX, W0, 0.0)
{}


void OmegamOmegaXW0Cosmo::dump(ostream &stream) const
{
  stream << " Om-Ox-w0 cosmo : (" << Omegam() << ',' << Omegax() << ',' << W0() << ')' << endl;
}

void OmegamOmegaXW0Cosmo::SetParam(const int i, const double Val)
{
  double toto = Omegar();
  if (i==0) *this = OmegamOmegaXW0Cosmo(Val, Omegax(),       W0()      );
  else if (i==1) *this = OmegamOmegaXW0Cosmo(Omegam(), Val, W0()  );
  else if (i==2) *this = OmegamOmegaXW0Cosmo(Omegam(),       Omegax(),  Val);
  else FatalError(" wrong param number in SetParam");
  SetOmegar(toto);
}

double OmegamOmegaXW0Cosmo::GetParam(const int i) const
{
  if (i==0) return Omegam();
  else if (i==1) return Omegax();
  else if (i==2) return W0();
  else FatalError(" wrong param number in GetParam");
  return 42; // rever reached
}


//================================ OmegamW0W1Cosmo ===========================
OmegamW0W1Cosmo::OmegamW0W1Cosmo(const double Omegam, const double W0, const double W1)
  : GeneralCosmo(Omegam, 1-Omegam, W0, W1)
{
}

void OmegamW0W1Cosmo::dump(ostream &stream) const
{
  stream << " Om-w0-w1 cosmo : (" << Omegam() << ',' << ',' << W0() << ',' << W1() << ')' << endl;
}

void OmegamW0W1Cosmo::SetParam(const int i, const double Val)
{
  double toto = Omegar();
  if (i==0) *this = OmegamW0W1Cosmo(Val, W0(), W1()      );
  else if (i==1) *this = OmegamW0W1Cosmo(Omegam(),  Val, W1()  );
  else if (i==2) *this = OmegamW0W1Cosmo(Omegam(),  W0(), Val);
  else FatalError("  wrong param number in SetParam");
  SetOmegar(toto);
}

double OmegamW0W1Cosmo::GetParam(const int i) const
{
  if (i==0) return Omegam();
  else if (i==1) return W0();
  else if (i==2) return W1();
  else FatalError(" wrong param number in GetParam");
  return 42; // rever reached
}


//================================ OmOlW0WaCosmo ===========================




OmegamOmegaxW0WaCosmo::OmegamOmegaxW0WaCosmo(const double OmegaM, const double OmegaX, 
			     const double W0, const double Wa)
  : Cosmology(1-OmegaM-OmegaX), omegam(OmegaM), omegax (OmegaX), w0(W0), wa(Wa)
{
}

#ifdef STORAGE
class W0Wahinv {
public :
  const OmegamOmegaxW0WaCosmo &par;
  W0Wahinv(const OmegamOmegaxW0WaCosmo &C) : par(C) { } ;
  double operator()(const double z) const 
  {
    /* the following formula assumes that w(z) = w0+wa*z/(1+z) */
    double h2 = par.omegam*cub(1.+z) 
      + par.omegax * exp(3.*(log(1+z)*(1+par.w0+par.wa)-par.wa*z/(1+z)))
      + sq(1+z)*par.omegak;
    return 1./sqrt(h2);
  }
};
#endif



double OmegamOmegaxW0WaCosmo::integ_hinv(const double &ZMin, const double &ZMax, 
			     const int NStep) const
{
  // change of variable u = 1./sqrt(1+z), because it gives more accurate results on all Z ranges.
  // Z = Z_CMB even works pretty well (6e-6, with 100 steps)
  double umax = 1./sqrt(1+ZMin);
  double umin = 1./sqrt(1+ZMax);
  double step = (umax-umin)/NStep;
  double u = umin+0.5*step;
  double sum =0;
  for (int i=0; i<NStep; ++i, u+=step)
    {
      double u2 = sq(u);
      double fact_ox = exp(3*(-2*log(u)*(w0+wa)-wa+wa*u2));
      double u3hinv = sqrt(omegam + omegax*fact_ox + omegak*u2 + omegar/u2);
      sum += 1./u3hinv;
    }
  return 2.*sum*step;
}

  
double OmegamOmegaxW0WaCosmo::E(const double z) const
{
  return sqrt(omegam*cub(1.+z) 
	      + omegax * exp(3.*(log(1+z)*(1+w0+wa)-wa*z/(1+z)))
	      + sq(1+z)*(omegak+sq(1+z)*omegar));
}

void OmegamOmegaxW0WaCosmo::SetParam(const int i, const double Val)
{
  double toto = Omegar();
  if (i==0) *this = OmegamOmegaxW0WaCosmo(Val, Omegax(), W0(),       wa      );
  else if (i==1) *this = OmegamOmegaxW0WaCosmo(Omegam(), Val,       W0(), wa  );
  else if (i==2) *this = OmegamOmegaxW0WaCosmo(Omegam(), Omegax(),Val, wa  );
  else if (i==3) *this = OmegamOmegaxW0WaCosmo(Omegam(), Omegax(),   W0(), Val);
  else FatalError(" ben merdealors");
  SetOmegar(toto);
}

double OmegamOmegaxW0WaCosmo::GetParam(const int i) const
{
  if (i==0) return Omegam();
  else if (i==1) return Omegax();
  else if (i==2) return W0();
  else if (i==3) return wa;
  else FatalError(" wrong param number in GetParam");
  return 42; // rever reached
}

void OmegamOmegaxW0WaCosmo::dump(ostream &stream) const
{
  stream << " om-ox-w0-wa cosmo : (" << omegam << ',' << omegax << ',' << w0 << ',' << wa << ')' << endl;
}


// ========= OmegamW0WaCosmo ====================

OmegamW0WaCosmo::OmegamW0WaCosmo(const double Omegam, const double W0, const double Wa)
  : OmegamOmegaxW0WaCosmo(Omegam, 1-Omegam, W0, Wa)
{
}

void OmegamW0WaCosmo::dump(ostream &stream) const
{
  stream << " Om-w0-wa cosmo : (" << Omegam() << ',' << ',' << W0() << ',' << wa << ')' << endl;
}

void OmegamW0WaCosmo::SetParam(const int i, const double Val)
{
  double toto = Omegar();
  if (i==0) *this = OmegamW0WaCosmo(Val, W0(),       wa      );
  else if (i==1) *this = OmegamW0WaCosmo(Omegam(),  Val, wa  );
  else if (i==2) *this = OmegamW0WaCosmo(Omegam(),       W0(), Val);
  else FatalError(" wrong param number in SetParam");
  SetOmegar(toto);
}

double OmegamW0WaCosmo::GetParam(const int i) const
{
  if (i==0) return Omegam();
  else if (i==1) return W0();
  else if (i==2) return wa;
  else FatalError(" wrong param number in GetParam");
  return 42; // rever reached
}


