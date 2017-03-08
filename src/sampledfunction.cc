#include <iostream>
#include <stdlib.h>

#include "sampledfunction.h"
#include "datafile.h"
#include <assert.h>
#include "snfitexception.h"

#include <vector>

#ifdef USE_CANVAS
// To allocate Tgraph from general1Dfunction
#include <TGraph.h>
#endif

#define INTERPOLATION_DEGREE 1

// on contourne ce qu'a fait Julien en particulier l'implementation de la bicubique --- les photoz n'ont pas été testés avec
//#define SAMPLED_OLD

double ZeroIfNegative(double x) 
{ if (x<0) return 0; else return x;}


ostream& operator << ( ostream &stream, const Sampled1DFunction &This) 
{ This.dump(stream); return stream;}



double Sampled1DFunction::SolveForX(const double &FVal, const double XLeft, const double XRight, double Prec) const
{
  double xleft = XLeft;
  double xright = XRight;
  double fleft = Value(xleft) - FVal;
  double fright = Value(xright) - FVal;
  if (fleft*fright > 0)
    {
      cerr << "Sampled1DFunction::SolveForX : bad guesses" << endl;
      return 0;
    }
  double xmid;
  while (xright - xleft > Prec)
    {
      xmid = (xright+xleft)*0.5;
      double fmid = Value(xmid) - FVal;
      if ( fmid * fleft < 0)
	{
	  xright = xmid;
	  fright = fmid;
	}
      else
	{
	  xleft = xmid;
	  fleft = fmid;
	}
    }
  return xmid;
}
  

#include <math.h> /* for floor() */


// implementation of Histo1D

Histo1D:: Histo1D(const string &FileName)
{
  DataFile file(FileName);
  xmin = file.Globals().getDoubleValue("XMIN");
  xmax = file.Globals().getDoubleValue("XMAX");
  nbin = int(file.Globals().getDoubleValue("NBIN"));

  std::cout << " Histo1D::Histo1D : untested constructor... remove this warning if test is succeesful... " << endl;
  init(nbin,xmin,xmax);
  double *xv = new double [nbin];
  double *hv = new double [nbin];
  file.ReadCols(xv,hv,nbin);
  for (int i=0; i< nbin; ++i)
    HFill(xv[i], hv[i]);
  delete [] xv;
  delete [] hv;
}

void  Histo1D::init(const int NBin, const double XMin, const double XMax)
{
  nbin = NBin;
  xmin = XMin;
  xmax = XMax;
  xstep = (xmax - xmin)/double(nbin);
  values = new double[nbin];
  memset(values, 0, nbin * sizeof(double)); 
  cumval = NULL;
}

Histo1D::Histo1D(const Histo1D &O) : Sampled1DFunction()
{
  init(O.nbin, O.xmin, O.xmax);
  memcpy(values,O.values,nbin*sizeof(double));
}


  
Histo1D::Histo1D(const int NBin, const double XMin, const double XMax)
{
  init(NBin, XMin, XMax);
}

void Histo1D::Cumulate()
{
  cumval = new double [nbin+1];
  cumval[0] = 0;
  for (int i = 0; i< nbin; ++i)
    {
      cumval[i+1] = cumval[i] + values[i];
    }
  double norm = 1./cumval[nbin];
  for (int i = 0; i<=nbin; ++i) cumval[i] *= norm;
}

int locate_val(const double *array, const int &length, const double &val)
{// assumes array is increasing
  int low = 0;
  int high = length-1;
  while (high-low > 1)
    {
      int mid = (low+high)/2;
      if (array[mid] > val)
	{
	  high = mid;
	}
      else 
	{
	  low = mid;
	}
    }
  return low;
}

double Histo1D::GetBinContent(const int I) const
{
  if (I<0 || I>=nbin) return -1;
  return values[I];
}

void Histo1D::SetBinContent(const int I, const double &Val)
{
  if (I<0 || I>=nbin) throw(SnfitException(" bad bin index in SetBinContent"));
  values[I] = Val;
}
  
void Histo1D::HFill(const double &X, const double &Weight=1.)
{
  int bin = int(floor((X - xmin)/xstep));
  if (bin<0 || bin>=nbin) return;
  values[bin] += Weight;
} 
 
double Histo1D::XVal(const double CumulFrac)
{
  if (!cumval) Cumulate();
  int ichannel = locate_val(cumval, nbin+1, CumulFrac);
  double dy = cumval[ichannel+1] - cumval[ichannel];
  if (dy > 0)
    {
      return xmin + (double(ichannel)+(CumulFrac-cumval[ichannel])/dy)*xstep;
    }
  else return xmin+xstep*ichannel;
}

double Histo1D::Value(const double &XWhere) const
{
  int bin = int(floor((XWhere - xmin)/xstep));
  if (bin<0 || bin>=nbin) {
    cerr << "  Histo1D:: XWhere out of range" << endl;
    return 0;
  }
  return values[bin];
}

void Histo1D::dump(ostream & stream) const
{
  stream  << " Histo1D::dump not implemented yet ";
}




Histo1D::~Histo1D() 
{ 
  if (values) delete [] values;
  if (cumval) delete [] cumval;
}
  

// implementation of EqualSteps1DFunction
#include <string> /* for memcpy */

EqualSteps1DFunction::EqualSteps1DFunction(const string &FileName)
{
  fval = 0;
  General1DFunction f(FileName);
    nstep = f.Nval();
    if (nstep == 0)
    {
      init(0,-1,0,NULL);
      return;
    }
  if (nstep == 1)
    {
      init(f.Xmin(), f.Xmax(), nstep, f.values());
      return;
    }
    // check that it is an "EqualStep"
  xstep = f.Xval(1)-f.Xval(0);
  xstepinv = 1/xstep;
  for (int i = 1; i<nstep; ++i)
    {
      if (fabs((f.Xval(i)-f.Xval(i-1))*xstepinv -1)>1e-6)
	{
	  cout << " File : " << FileName << " is NOT an EqualStep1DFunction " << endl;
	  init (0,-1,0,0);
	  return;
	}
    }
  init(f.Xmin(), f.Xmax(), nstep, f.values());
}



EqualSteps1DFunction::EqualSteps1DFunction(const EqualSteps1DFunction &Model) :
  Sampled1DFunction()
{
  fval = 0;
  init(Model.xmin, Model.xmax, Model.nstep, Model.fval);
}


EqualSteps1DFunction::EqualSteps1DFunction(const double Begin, const double End, const int Nval, const double *Values) :
  Sampled1DFunction()
{
  fval = 0;
  init(Begin, End, Nval, Values);
}

EqualSteps1DFunction::EqualSteps1DFunction(const double Begin, const double End, const double Step)
{
  int nstep = int (ceil((End-Begin)/Step)) + 1;
  fval = 0;
  init(Begin, End, nstep, 0);
}


EqualSteps1DFunction::EqualSteps1DFunction(const Sampled1DFunction &Input, int NStep)
{
  xmin = Input.Xmin();
  xmax = Input.Xmax();
  nstep = NStep;
  xstepinv = (nstep-1)/(xmax-xmin);
  fval = new double[nstep];
  xstep = 1./xstepinv;
  // fill
  for (int i=0; i<nstep; i++)
    {
      fval[i] = Input.Value(Xval(i));
    }
}


void EqualSteps1DFunction::Cumulate()
{
  //cumulate
  for (int i=2; i<nstep; ++i)
    {
      fval[i] += fval[i-1];
    }
  // normalize
  double norm = 1./fval[nstep-1];
  for (int i=0; i< nstep; ++i) fval[i] *= norm;
}  



EqualSteps1DFunction& EqualSteps1DFunction::operator = (const EqualSteps1DFunction &Model)
{
  init(Model.xmin, Model.xmax, Model.nstep, Model.fval);
  return *this;
}
  

void EqualSteps1DFunction::init(const double Begin, const double End, const int Nval, const double *Values)
{
  xmin = Begin;
  xmax = End;
  nstep = Nval;
  if (fval) delete [] fval;
  fval = (nstep)? new double [nstep] : NULL;
  if (Values)
    memcpy(fval, Values, sizeof(double)*nstep);
  else
    memset(fval, 0, sizeof(double)*nstep);
  xstep = (xmax-xmin)/(nstep-1);
  xstepinv = 1/xstep;
}


double EqualSteps1DFunction::Value(const double &XWhere) const
{
  double di = (XWhere-xmin)*xstepinv;
  int i = int(floor(di));
  if (i<-1) return 0;
  if (i>=nstep) return 0;
  if (i==-1) i=0;
  else if (i==nstep-1) i = nstep-2;
  double x = di-i;
  return (fval[i]*(1.-x) + fval[i+1]*x);
}

void EqualSteps1DFunction::operator *= (const double Factor)
{
  for (int i=0; i< nstep; ++i)   fval[i] *= Factor;
}

#include <fstream>
void EqualSteps1DFunction::Write(const string &FileName) const
{
  ofstream f(FileName.c_str());
  f << "# x : " << endl << "# val : " << endl << " #end " << endl;
  for (int k=0; k<nstep; ++k) f << xmin+k*xstep << ' ' << fval[k] << endl;
}


void EqualSteps1DFunction::dump(ostream & stream) const
{
  stream << "[" << xmin << ',' << xmax << ']' << " in " << nstep << " steps "  << endl; 
}

#ifdef USE_CANVAS
// who is going to release the memory allocated by new? root!
TGraph EqualSteps1DFunction::AllocateTGraph(const double & xxmin , 
					    const double & xxmax)
{
  double *xval = new double[nstep];
  double *val = new double[nstep];
  int j = 0;
  
  for (int i = 0; i < nstep; ++i)
    {
      if (fval[i] >-50 )

	{
	  double x = xmin + i * xstep;
	  if (x> xxmax || x < xxmin) continue; 
	  xval[j] = x;   
	  val[j] = fval[i];
	  //cout << j << " " <<  xval[j] << " " <<  val[j] << endl;
	  ++j;
	}
    }
  return(TGraph(j,xval,val));

}
#endif


// implementation of General1DFunction



General1DFunction::General1DFunction(const int Nval, const double *Values, const double *XValues)
{
  fval = 0; xval = 0;
  init(Nval, Values, XValues);
}

General1DFunction::General1DFunction(const General1DFunction &Model) :
  Sampled1DFunction()
{
  fval = 0; xval = 0;
  init(Model.nstep, Model.fval, Model.xval);
}


General1DFunction::General1DFunction()
{
  fval = 0; xval = 0;
  init (0,0,0);
}

/* 
   ranker stuff :
   found on :   http://sites.google.com/site/jivsoft/Home/
   author : Ken Wilder
   License : BSD
   doc : http://sites.google.com/site/jivsoft/Home/compute-ranks-of-elements-in-a-c---array-or-vector
   heavily simplified
*/

template <class T> class ranker
{
 private:
  const T* p;
  size_t sz;

 public:
  //  ranker(const vector<T>& v) : p(&v[0]), sz(v.size()) { }
  ranker(const T* tp, size_t s) : p(tp), sz(s) { }

  int operator()(size_t i1, size_t i2) const { return(p[i1] < p[i2]); }  
template <class S>
  void get_orders(vector<S>& w) const {
    w.resize(sz);
    w.front() = 0;
    for (typename vector<S>::iterator i = w.begin(); i != w.end() - 1; ++i)
      *(i + 1) = *i + 1;
    sort(w.begin(), w.end(), *this);
  }
};


static bool is_strictly_ordered(const double *xv, const unsigned ndat)
{
  if ( ndat < 2 )
    return true ;
  const double *pend = xv+ndat;
  const double *p = xv+1; 
  for ( ; p<pend; ++p)
    {
      if (*(p-1) > *p) break;
      if (*(p-1) == *p) {
	cerr << "ERROR in is_strictly_ordered, equal abcissa in sampled function at *p = " << *p << endl;
	
	throw(string(" equal abcissa in sampled function"));
      }
    }
  return (p == pend);
}



static bool reordered(double *xv, double *fv, const unsigned ndat)
{
  if (is_strictly_ordered(xv,ndat)) return false;
  ranker<double> rk(xv,ndat);
  vector<unsigned int> ranks;
  rk.get_orders(ranks);
  double *tmp = new double[ndat];
  memcpy(tmp,xv,sizeof(double)*ndat);
  for (unsigned k=0; k<ndat; ++k)    xv[k] = tmp[ranks[k]];
  memcpy(tmp,fv,sizeof(double)*ndat);
  for (unsigned k=0; k<ndat; ++k)    fv[k] = tmp[ranks[k]];
  delete [] tmp;
  // check if there are no equal abcissa
  is_strictly_ordered(xv, ndat);
  return true;
}

General1DFunction::General1DFunction(const string &Filename, const int Icol)
{
  Read(Filename, Icol);
}

void General1DFunction::Read(const string &Filename, const int Icol)
{
  DataFile file(Filename);
  int ndat = file.NDat();
  if (!ndat)
    cerr << " WARNING : no data in " << Filename << " ... " << endl;
  double *xv = new double[ndat];
  double *fv = new double[ndat];
  file.ReadCols(xv, fv, ndat, Icol);
  if (file.Globals().HasKey("WAVELENGTH_IN_NM"))
    for (int k=0; k< ndat; ++k) xv[k]*=10; // use angstroms internally
  // Check if data is in increasing x order
  try 
    {
      if (ndat && reordered(xv,fv,ndat))
	cout << "WARNING : data in "+Filename+" had to be sorted in ascending order" << endl;
    }
  catch (const string &message)
    {
      throw(SnfitException(message+" in file "+Filename));
    }	    
  fval =0; xval = 0;
  init(ndat,fv, xv);
  delete [] xv;
  delete [] fv;
}


void General1DFunction::init(const int Nval, 
			     const double *Values, const double *XValues)
{

  lasti = 0;
  nstep = Nval;
  if (Values)
    {
      if (fval) delete [] fval;
      fval = new double [nstep];
      memcpy(fval, Values, sizeof(double)*nstep);
    }
  if (XValues)
    {
      if (xval) delete [] xval;
      xval = new double [nstep];
      memcpy(xval, XValues, sizeof(double)*nstep);
      xmin = xval[0];
      xmax = xval[nstep-1];
    }
  else
    {
      xmin = xmax = 0;
    }
  
}

General1DFunction&   General1DFunction::operator = (const General1DFunction &Right)
{
  init(Right.nstep, Right.fval, Right.xval);
  return *this;
}


General1DFunction::~General1DFunction()
{ 
  if (fval) delete [] fval;
  if (xval) delete [] xval;
}



extern "C" {
  double divdifd_(const double *, const double *, const int *, const double *, const int*);
}

/*double General1DFunction::Xmax(const double Factor, const int Nval) const
  {
  double xmaxf = 0;
  for(int i = 0; i < Nval; ++i)
  if(Fval(i) == Factor)
  xmaxf = Xval(i);
  return xmaxf;
  }
*/

/* this routine assumes a sorted array */
static int LocateBelow(const double *array, int NValues, const double &Where)
{
  if (NValues == 0) return 0;
  int low =0;
  int high = NValues-1;
  if (Where > array[high]) return high;
  if (Where < array[low]) return -1;
  while (high - low > 1)
    {
      int mid = (low+high)/2;
      if (array[mid] > Where) high = mid;
      else low = mid;
    }
  return low;
}

#ifdef STORAGE /* actually there was one in the same file ... */
/* I did not find a routine that brackets a value in a sequence */
int index_lower_than(const double Val, const double *B, const double *E)
{
  const double *b = B;
  const double *e = E;
  if (Val >= *e ) return e-B;
  if (Val < *b) return -1;
  int dist = e-b;
  do
    {
      const double *mid = b+dist/2;
      if (Val > *mid) b = mid;
      else e=mid;
      dist = e-b;
    }
  while (dist >1);
  return b-B;
}
#endif

 int FastLocateInArray(const int lasti, const double *xval, 
			     const int nval, const double &XWhere)
{
  /* assumes that  XWhere is really between xval[0] and xval[nval-1]
     and the array is sorted. */
  int i = lasti;
  if (XWhere < xval[i])
    i = LocateBelow(xval, i+1, XWhere);
  else if (XWhere >= xval[i+1])
    {
      if (i <= nval-3 && XWhere < xval[i+2] ) i++;
      else 
	i = i+1+LocateBelow(xval+i+1, nval-i-1,XWhere);
    }
  // note that if ( xval[i] <= XWhere < xval[i+2]), there is no binary search.
  if (i==nval-1 && XWhere == xval[nval-1]) i--;
  return i;
}




double General1DFunction::Value(const double &XWhere) const
{
  #ifdef SAMPLED_OLD
    {
      if (XWhere < xmin || XWhere > xmax) return 0;
      int interp_degree = INTERPOLATION_DEGREE;
      return divdifd_(fval, xval, &nstep, &XWhere, &interp_degree);
    }
 #endif

  // ------------ ORIG

 if (XWhere < xmin || XWhere > xmax) return 0;
  /* This routine tries to speedup the computation of integrals
     by using the last index (lasti) of the interpolated value,
     assuming that the current request is likely satisfied by 
     interpolating within the same interval or the one just after.
     The binary search is only launched if both assumptions fail
  */
#if ((INTERPOLATION_DEGREE == 1))
  // lasti should always remain <= nstep-2
  assert (lasti<=nstep-2);
  int i = FastLocateInArray(lasti, xval, nstep, XWhere);
  if (i>nstep-2) 
    { 
      if (XWhere == xmax) return fval[nstep-1];
      cout << " General1DFunction::Value interpolation bug " << endl; 
      abort();
    }
  double x = /* return */ (XWhere-xval[i])/(xval[i+1]-xval[i]);
  assert(x>=0 && x<=1);
  (unsigned &) lasti = i; // this cast explicitly violates const'ness 
  return (fval[i]*(1.-x) + fval[i+1]*x);
#else
  int interp_degree = INTERPOLATION_DEGREE;
  return divdifd_(fval, xval, &nstep, &XWhere, &interp_degree);
#endif
#if 0
  assert((val1 ==0 && val2 == 0) || fabs((val1-val2)/val1)< 1e-6);
  return val1;
  if (val1 != val2)
    {
      cout << " val1 , val2 " << val1 << ' ' << val2 << endl;
      abort();
    }
#endif
}

void General1DFunction::RescaleX(const double Factor)
{
  for (int i=0; i< nstep; ++i)   xval[i] *= Factor;
  xmin *= Factor;
  xmax *= Factor;
}

void General1DFunction::ShiftX(const double & Shift)
{
  for (int i=0; i< nstep; ++i)   xval[i] += Shift;
  xmin += Shift;
  xmax += Shift;
}


void General1DFunction::operator *= (const double Factor)
{
  for (int i=0; i< nstep; ++i)   fval[i] *= Factor;
}



#ifdef USE_CANVAS
// these routine should NOT change the General1DFunction: redshifting the graph
// and redshifting the function are totally different things
TGraph General1DFunction::AllocateRedshifftedTGraph(const double & z)
{
  double factor = (1+z);
  
  for (int i=0; i< nstep; ++i)  {
    xval[i] *= factor;
    
  }
  return(TGraph(nstep,xval,fval));

}

TGraph General1DFunction::AllocateBlueshifftedTGraph(const double & z)
{
  // see comment above
  double factor = 1/(1+z);
  for (int i=0; i< nstep; ++i)   xval[i] *= factor;

  return(TGraph(nstep,xval,fval));

}

// this routine should be in the client code!
TGraph General1DFunction::AllocateTGraph() const
{
  
  return(TGraph(nstep,xval,fval));

}
#endif


void General1DFunction::dump(ostream & stream) const
{
  stream << "[" << xmin << ',' << xmax << ']' << " in " << nstep << " steps "  << endl; 
}


#include <fstream>
bool General1DFunction::Write(const string &Filename) const
{
  ofstream s(Filename.c_str());
  if (!s)
    {
      cerr << "General1DFunction::Write : could not open " << Filename << endl;
      return false;
    }
  return Write(s);
}

bool General1DFunction::Write(ostream &s) const
{
  for (int ix=0; ix<nstep; ix++)
    s << xval[ix] << ' ' << fval[ix] << endl;
  return true;
}


double General1DFunction::Integral() const
{
  double integral = 0;
  for (int i=1; i< nstep; ++i) 
    integral += (fval[i]+fval[i-1])*(xval[i]-xval[i-1]);
  return integral*0.5;
}

double General1DFunction::Mean() const
{
  double integral = 0;
  double acc =0;
  for (int i=1; i< nstep; ++i) 
    {  
      double toto = (fval[i]+fval[i-1])*(xval[i]-xval[i-1]);
      integral += toto;
      acc += (xval[i] + xval[i-1]) * toto;
    }
  return acc/integral/2;
}


double General1DFunction::BoundedIntegral(const double & MinBound, const double & MaxBound) const
{
  double integral = 0;
  for (int i=1; i< nstep; ++i) 
    if (xval[i] > MinBound && xval[i]<MaxBound )
     {
       if (xval[i-1] < MinBound)
	 integral += (fval[i]+fval[i-1])*(xval[i]-MinBound);
       else
	 integral += (fval[i]+fval[i-1])*(xval[i]-xval[i-1]);
       if (xval[i+1]> MaxBound) 
	 integral += (fval[i+1]+fval[i])*(MaxBound-xval[i]);
     }
  return integral*0.5;
}






#include <algorithm> /* for min and max */

// why do I return a General1DFunction rather than an EqualStep1DFunction ?? P.A CHECK 
// answer : because EqualSteps1DFunction misses write(FileName) !
General1DFunction Product(const double Step, const Sampled1DFunction &Left, const Sampled1DFunction &Right)
{

  double xmin = max(Left.Xmin(), Right.Xmin());
  double xmax = min(Left.Xmax(), Right.Xmax());
  int nstep = int (ceil((xmax-xmin)/Step)) + 1;
  double step = (xmax-xmin)/(nstep-1);
  double *values = new double [nstep];
  double *xval = new double[nstep];
  for (int i=0; i<nstep; ++i)
    {
      double x = xmin+i*step;
      xval[i] = x;
      values[i] = Left.Value(x)*Right.Value(x);
      //cout << x << ' ' << Left.Value(x) << ' ' << Right.Value(x) << endl;
    }
  General1DFunction res(nstep, values, xval);
  delete [] values;
  delete [] xval;
  return res;
}


#include <stdarg.h> // for va_list, va_start, va_end;

General1DFunction Product(const double Step,const int NTerms, ... )
{
  const Sampled1DFunction** terms = new const Sampled1DFunction*[NTerms];
  va_list arglist;
  va_start(arglist,NTerms);
  for (int i=0; i<NTerms; ++i) terms[i] = va_arg(arglist, Sampled1DFunction*);
  va_end(arglist);

  double xmin = terms[0]->Xmin();
  double xmax = terms[0]->Xmax();
  for (int i=1; i<NTerms; ++i)
    {
      xmin = max(xmin, terms[i]->Xmin());
      xmax = min(xmax, terms[i]->Xmax());
    }
  int nstep = int (ceil((xmax-xmin)/Step))+1;
  double step = (xmax-xmin)/(nstep-1);
  double *values = new double[nstep];
  double *xval= new double[nstep];
  for (int i=0; i<nstep; ++i)
    {
      double x = xmin+i*step;
      xval[i] = x;
      double value = terms[0]->Value(x);
      for (int j=1; j<NTerms; ++j) value *= terms[j]->Value(x);
      values[i] = value;
    }
  
  General1DFunction res(nstep, values, xval);
  delete [] values;
  delete [] xval;
  delete [] terms;
  return res;
}


General1DFunction Product(const double Step, vector<Sampled1DFunction*> terms)
{
  
  double xmin = terms[0]->Xmin();
  double xmax = terms[0]->Xmax();
  int NTerms = terms.size();
  for (int i=1; i<NTerms; ++i)
    {
      xmin = max(xmin, terms[i]->Xmin());
      xmax = min(xmax, terms[i]->Xmax());
    }
  int nstep = int (ceil((xmax-xmin)/Step))+1;
  double step = (xmax-xmin)/(nstep-1);
  double *values = new double[nstep];
  double *xval = new double[nstep];
  for (int i=0; i<nstep; ++i)
    {
      double x = xmin+i*step;
      xval[i] = x;
      double value = terms[0]->Value(x);
      for (int j=1; j<NTerms; ++j) value *= terms[j]->Value(x);
      values[i] = value;
    }
  General1DFunction res(nstep, values, xval);
  delete [] values;
  delete [] xval;
  return res;
}


//// ==================    implementation of General2DFunction

static double * new_and_copy( const double *Right, const int N)
{
  double *new_array = new double [N];
  memcpy(new_array, Right, N*sizeof(double));
  return new_array;
}

static bool not_identical(const double *array1, const double *array2, int count)
{
  for (int i=0; i< count; ++i) if (array1[i] != array2[i]) return true;
  return false;
}

/* assume a sorted array */
static int LocateNearest(const double *array, int NValues, const double &Where)
{
  if (NValues == 0) return 0;
  int low = LocateBelow(array, NValues, Where);
  if (low==NValues-1) return low;
  if (fabs(array[low]-Where) > fabs(array[low+1] - Where)) return low+1; 
  return low;
}




bool Grid2DFunction::add_slice(const double &X, const double *Fval, const double *Yval, const int Nval)
{
  if (xval == NULL)
    {
      assert(is_strictly_ordered(Yval, Nval));
      xval = new double [1];
      xval[0] = X;
      fval = new double* [1];
      fval[0] = new_and_copy(Fval,Nval);
      yval = new_and_copy(Yval,Nval);
      nstepx = 1;
      nstepy = Nval;
      ymin = yval[0];
      ymax = yval[nstepy-1];

    }
  else
    {
      // check that the number of values is the same as the first row
      if (Nval != nstepy || not_identical(yval, Yval, Nval))
	{
	  cerr << " Grid2DFunction :: cannot handle data with different y sizes or values for different x " << endl;
	  cerr << " this happened for x = " << X << endl;
	  return false;
	}
      int low = LocateBelow(xval,nstepx, X);
      double **newfval = new double* [nstepx+1];
      double *newxval = new double[nstepx+1];
      for (int i=0; i<=low; ++i) 
	{
	  newfval[i] = fval[i];
	  newxval[i] = xval[i];
	}
      newfval[low+1] = new_and_copy(Fval,Nval);
      newxval[low+1] = X;
      for (int i=low+1; i<nstepx; ++i) 
	{
	  newfval[i+1] = fval[i];
	  newxval[i+1] = xval[i];
	}
      delete [] fval; fval = newfval;
      delete [] xval; xval = newxval;
      nstepx++;
    }
  xmin = xval[0];
  xmax = xval[nstepx-1];
  return true;
}



Grid2DFunction::Grid2DFunction(const string &Filename) {
  Read(Filename);
}


bool Grid2DFunction::Read(const string &Filename){

  lastix=lastiy=0;

  DataFile file(Filename);
  if (!file.IsValid()) return false;
  int ndat = file.NDat();
  int ncol = file.NCol();
  file.Rewind();
  char line[1024];
  int count = 0;
  fval = 0; xval = 0; yval = 0;
  if (ncol != 3)
    {
      cerr << " cannot read " << Filename << " with Grid2DFunction::Grid2DFunction(const string &Filename)" << endl;
      return false;
    }
  double *xv = new double[ndat];
  double *yv = new double[ndat];
  double *fv = new double[ndat];
  double old_x;
  int county=0;
  while ((file.next_numerical_data(line,1024)))
    {
      if (count == ndat)
	{
	  count ++; // to trigger the error message
          break;
	}
      double x,y,f;
      if (sscanf(line, "%lf %lf %lf",&x,&y,&f) != 3)
	{
	  cerr << " Grid2DFunction : cannot read 3 double in line: " << endl
	       << line << endl
	       << " from file " << Filename << " (Stop reading there )" << endl;
	  break;
	}
      if (count != 0 && (x != old_x))
	{
	  add_slice(old_x, fv, yv, county);
	  county = 0;
	}
      yv[county] = y;
      fv[county] = f;
      county++;
      count++;
      old_x = x;
    }
  // fill the last slice
  add_slice(old_x, fv, yv, county);
  if (count != ndat)
    {
      cerr << " General1DFunction::General1DFunction : Inconsistency between DataFile::NDat() and " << endl 
	   << " ...... just counting the numerical data line in " << Filename << " end of file not read ... " << endl;
    }
  delete [] fv;
  delete [] xv;
  delete [] yv;
  return true;
}





//Grid2DFunction::Grid2DFunction(const int NXval, const int NYval, const double **Values, const double *XValues, const double *YValues)
//{
//  fval = 0; xval = 0; yval = 0;
//  init(NXval, NYval, Values, XValues, YValues);
// }

Grid2DFunction::Grid2DFunction(const Grid2DFunction &Right)
{
  fval = 0; xval = 0; yval = 0;
  init(Right.nstepx, Right.nstepy, Right.fval, Right.xval, Right.yval);
}

Grid2DFunction::Grid2DFunction(const int NXval, const int NYval, const double ** FValues, const double *XValues, const double *YValues) {
  fval = 0; xval = 0; yval = 0;
  init(NXval,NYval,FValues,XValues,YValues);
}

void Grid2DFunction::init(const int NXval, const int NYval, const double * const *FValues, const double *XValues, const double *YValues)
{
  assert(is_strictly_ordered(XValues,NXval));
  assert(is_strictly_ordered(YValues,NYval));
  lastix = lastiy = 0;
  nstepx = NXval;
  nstepy = NYval;
  double **newfval = new double* [nstepx];
  for (int i=0; i<nstepx; ++i)
      newfval[i] = new_and_copy(FValues[i], nstepy);
  // in case newfval and Fvalues are the same 
  if (fval) delete [] fval;
  fval = newfval;

  double *newxval = new_and_copy(XValues,nstepx);
  if (xval) delete [] xval;
  xval = newxval;
  xmin = xval[0];
  xmax = xval[nstepx-1];

  double *newyval = new_and_copy(YValues, nstepy);
  if (yval) delete [] yval;
  yval = newyval;
  ymin = yval[0];
  ymax = yval[nstepy-1];
}

Grid2DFunction::~Grid2DFunction()
{
  if (fval)
    { for (int i=0; i< nstepx; ++i) delete [] fval[i]; delete [] fval;}
  if (xval) delete [] xval;
  if (yval) delete [] yval;
} 

double Grid2DFunction::Value(const double &XWhere, const double &YWhere) const {
  //if(INTERPOLATION_DEGREE==1)
#ifdef SAMPLED_OLD
    return BilinearValue(XWhere,YWhere);
#endif

  return BiCubicConvolution(XWhere,YWhere); // ORIG
}

static double kernval(const double& xval) {
  double x = fabs(xval);
  if(x>2) return 0;
  double a = -0.5;
  if(x<1)
    return x*x*( (a+2)*x-(a+3) )+1;
  return a*( -4+x*(8+x*(-5+x)));
}

/*!
    bicubic convolution using kernel  
    
    x = fabs((distance of point to node)/(distance between nodes)) 

    W(x) = (a+2)*x**3-(a+3)*x**2+1 for x<=1
    W(x) = a( x**3-5*x**2+8*x-4) for 1<x<2
    W(x) = 0 for x>2
*/
double Grid2DFunction::BiCubicConvolution(const double &XWhere, const double &YWhere) const
{
  
  if (XWhere < xmin || XWhere > xmax || YWhere < ymin || YWhere > ymax) return 0.0;
  
  // find nearest node
  /* slow 
  int ix = LocateNearest(xval,nstepx, XWhere);
  int iy = LocateNearest(yval,nstepy, YWhere);
  if(xval[ix]>XWhere) { ix--;}
  if(yval[iy]>YWhere) { iy--;}
  */
  int ix = FastLocateInArray(lastix, xval, nstepx, XWhere);
  int iy = FastLocateInArray(lastiy, yval, nstepy, YWhere);
  (unsigned &) lastix = ix; // this cast explicitly violates const'ness
  (unsigned &) lastiy = iy; // this cast explicitly violates const'ness
  
  if(
     nstepx<3 || nstepy<3
     || ix==0 
     || iy==0 
     || ix>(nstepx-3) 
     || iy>(nstepy-3)
     ) 
    return BilinearValue(XWhere,YWhere);
  
  
  /*
  double invdx = 0;
  if(ix>0) invdx=1./(xval[ix]-xval[ix-1]); else invdx=1./(xval[ix+1]-xval[ix]);
  double invdy = 0;
  if(iy>0) invdy=1./(yval[iy]-yval[iy-1]); else invdy=1./(yval[iy+1]-yval[iy]);
  */
  
  double invdx=1./(xval[ix+1]-xval[ix]);
  double invdy=1./(yval[iy+1]-yval[iy]);
  
  double wy[4];
  double dy = (yval[iy]-YWhere)*invdy;
  for(int j=iy-1;j<iy+3;j++) {
    wy[j-(iy-1)]=kernval(dy+j-iy);
  }

  double v=0;
  double dx = (xval[ix]-XWhere)*invdx;
  for(int i=ix-1;i<ix+3;i++) {
    double wx = kernval(dx + i-ix);
    const double* wyv = wy;
    for(int j=iy-1;j<iy+3;j++,wyv++) {
      v += wx*(*wyv)*val(i,j);
    }
  }
  return v;
}
  

#ifdef STORAGE
double Grid2DFunction::BiCubicInterpolation(const double &XWhere, const double &YWhere) const
{
  if (XWhere < xmin || XWhere > xmax || YWhere < ymin || YWhere > ymax) return 0.0;
  int i_xl = LocateBelow(xval,nstepx, XWhere);
  int i_yl = LocateBelow(yval,nstepy, YWhere);
  int i_xu = i_xl+1;
  int i_yu = i_yl+1;
  double xl = xval[i_xl];
  double xu = xval[i_xl+1];
  double yl = yval[i_yl];
  double yu = yval[i_yl+1];
  
  double f[4];
  double dfdx[4];
  double dfdy[4];
  double d2fdxdy[4];
  double dx = xu-xl;
  double dy = yu-yl;

  
  /*
    ca va etre tres tres lent tout ca

    3--2
    |  |
    0--1
    
  
  
  f[0]=val(i_xl,i_yl); f[1]=val(i_xu,i_yl); f[2]=val(i_xu,i_yu); f[3]=val(i_xl,i_yu);  
  
  dfdx[0]=dfdx[1]=(f[1]-f[0])/dx;
  dfdx[2]=dfdx[3]=(f[2]-f[3])/dx;
  dfdy[0]=dfdy[3]=(f[3]-f[0])/dy;
  dfdy[1]=dfdy[2]=(f[2]-f[1])/dy;
  d2fdxdy[0]=d2fdxdy[3]=(dfdx[3]-dfdx[0])/dy;
  d2fdxdy[1]=d2fdxdy[2]=(dfdx[2]-dfdx[1])/dy;
  */

  /*
    1--2
    |  |
    0--3
    
  
  
  f[0]=val(i_xl,i_yl); f[1]=val(i_xl,i_yu); f[2]=val(i_xu,i_yu); f[3]=val(i_xu,i_yl);  
  
  dfdx[0]=dfdx[3]=(f[3]-f[0])/dx;
  dfdx[1]=dfdx[2]=(f[2]-f[1])/dx;
  dfdy[0]=dfdy[1]=(f[1]-f[0])/dy;
  dfdy[2]=dfdy[3]=(f[2]-f[3])/dy;
  d2fdxdy[0]=d2fdxdy[1]=(dfdx[1]-dfdx[0])/dy;
  d2fdxdy[2]=d2fdxdy[3]=(dfdx[2]-dfdx[3])/dy;
  */

  /*

  avec ca les valeurs aux noeuds sont ok :

 y^    0--3
  |    |  |
   ->  1--2
    x

  */
  
  f[0]=val(i_xl,i_yu); f[1]=val(i_xl,i_yl); f[2]=val(i_xu,i_yl); f[3]=val(i_xu,i_yu);  
  
  dfdx[0]=dfdx[3]=(f[3]-f[0])/dx;
  dfdx[1]=dfdx[2]=(f[2]-f[1])/dx;
  dfdy[0]=dfdy[1]=(f[0]-f[1])/dy;
  dfdy[2]=dfdy[3]=(f[3]-f[2])/dy;
  d2fdxdy[0]=d2fdxdy[1]=-(dfdx[0]-dfdx[1])/dy;
  d2fdxdy[2]=d2fdxdy[3]=-(dfdx[3]-dfdx[2])/dy;
  

  double res_f,res_dfdx,res_dfdy;
  bcuint(f,dfdx,dfdy,d2fdxdy,xl,xu,yl,yu,XWhere,YWhere,&res_f,&res_dfdx,&res_dfdy);
  return res_f; // ouf (un peu complique tout ca ...)
}
#endif

// bi-linear interpolation is probably a problem for Minuit 
double Grid2DFunction::BilinearValue(const double &XWhere, const double &YWhere) const
{
  if (XWhere < xmin || XWhere > xmax || YWhere < ymin || YWhere > ymax) return 0.0;

  int ix = FastLocateInArray(lastix, xval, nstepx, XWhere);
  int iy = FastLocateInArray(lastiy, yval, nstepy, YWhere);
  (unsigned &) lastix = ix; // this cast explicitly violates const'ness
  (unsigned &) lastiy = iy; // this cast explicitly violates const'ness
  /*
  int ix = LocateBelow(xval, nstepx, XWhere);
  int iy = LocateBelow(yval, nstepy, YWhere);
  */

  //tentative checks
  //if the grid has only 1 point 
  if ( nstepx==1 && nstepy==1 )
    return(val(ix,iy));


  if ( ix == nstepx-1 || iy == nstepy-1 )
    {
      cerr << "Grid2DFunction::Value() : something went wrong in  LocateBelow()" 
	   << endl;
      abort();
    }
  //  int ixnext = (nstepx) ? 0 : ix+1;
  //  int iynext = (nstepy) ? 0 : iy+1;
  double ax = (XWhere - xval[ix])/(xval[ix+1] - xval[ix]);
  double ay = (YWhere - yval[iy])/(yval[iy+1] - yval[iy]);
  //tentative checks
  if (ax<0 || ax >1 || ay<0 || ay>1)
    {
      cerr << "Grid2DFunction::Value() : something went wrong in  LocateBelow()" << endl;
      abort();
    }
 
  return (  (val(ix,iy)  *(1-ax) + val(ix+1,iy)    *ax)*(1-ay)
	   +(val(ix,iy+1)*(1-ax) + val(ix+1,iy+1)*ax)* ay );
}



void Grid2DFunction::dump(ostream & stream ) const
{
  stream << " nstepx " << nstepx << " nstepy " << nstepy << endl;
}

ostream & operator << ( ostream &stream, const Grid2DFunction &This)
{
  This.dump(stream); return stream;
}

Grid2DFunction&  Grid2DFunction::operator = (const Grid2DFunction &Right)
{
  init(Right.nstepx, Right.nstepy, Right.fval, Right.xval, Right.yval);
  return *this;
}

Grid2DXSlice Grid2DFunction::ConstXSlice(const double XVal) const
{ 
  return Grid2DXSlice(*this, XVal);
}  

void Grid2DFunction::RescaleXSlice(const double Xval, const double Factor, const bool Verbose)
{
    int ix = LocateNearest(xval,nstepx, Xval);
    if (Verbose)
      {
	cout << " Rescaling slice at X = " << Xval << " (i="<< ix << ") by " << Factor << endl;
      }
    for (int iy=0; iy< nstepy; iy++) val(ix,iy)*= Factor;
}

bool Grid2DFunction::AddSlice(const double &XWhere, const Sampled1DFunction &NewSlice)
{
  double *fv = new double[nstepy];
  for (int i=0; i< nstepy; ++i)
    {
      fv[i] = NewSlice.Value(yval[i]);
    }
  return (add_slice(XWhere, fv, yval, nstepy));
}

#include <stdio.h>

bool Grid2DFunction::Write(const string &Filename) const
{
  FILE *f = fopen(Filename.c_str(),"w");
  if (!f)
    {
      cerr << "Grid2DFunction::Write : could not open " << Filename << endl;
      return false;
    }
  for (int ix=0; ix<nstepx; ix++)
    {
      for (int iy = 0; iy < nstepy; iy++)
	{
	  fprintf(f, " %f  %f  %12.10g\n", xval[ix], yval[iy], val(ix,iy));
	}
    }
  fclose(f);
  return true;
}



#ifdef STORAGE
General1DFunction Grid2DFunction::ConstXSlice(const double XVal) const
{
  double *values = new double[nstepy];
  for (int iy =0; iy< nstepy; ++iy)
    {
      values[iy] = Value(XVal, yval[iy]);
    }
  return General1DFunction(nstepy, values, yval);
}
#endif


double Constant1DFunction::Value(const double &XWhere) const { return value;}
