#ifndef SAMPLEDFUNCTION__H
#define SAMPLEDFUNCTION__H

//! to handle arrays of values from which we get interpolations.

#include <iostream>
#include <string>
#include <vector>
// To allocate Tgraph from general1Dfunction
// this is a non-sense to do that in sampledfunction.cc


#include "globalval.h"


using namespace std;

double ZeroIfNegative(double x);

class EqualSteps1DFunction;

//! interface class
class Sampled1DFunction {

  public:

  virtual double Value(const double &XWhere) const = 0;
  virtual void dump(ostream & stream = cout) const = 0;
  friend ostream & operator << (ostream &stream, const Sampled1DFunction &This);
		
  virtual double Xmin() const =0;
  virtual double Xmax() const =0;
  //! rough solver for F(x) = FVal
  virtual double SolveForX(const double &FVal, const double XLeft, 
			   const double XRight, double Prec) const;
  virtual ~Sampled1DFunction() {};
};


class Histo1D : public Sampled1DFunction {
 private:
  double *values;
  int nbin;
  double xmin,xmax,xstep;
  void init(const int NBin, const double XMin, const double XMax);
  double *cumval;
  void Cumulate();

  public :
    Histo1D(const string & FileName);
    Histo1D(const int NBin, const double XMin, const double XMax);
    Histo1D(const Histo1D &O);

    void HFill(const double &X, const double &Weight);
    double Xmin() const {return xmin;};
    double Xmax() const {return xmax;};
    int  NBin() const { return nbin;}
    double GetBinContent(const int I) const;
    void SetBinContent(const int I, const double &Val);

    double Value(const double &XWhere) const;
   //! returns the abcissa value up to which the histogram cumulates do the argument.
   double XVal(const double CumulFrac);
   void dump(ostream & = cout) const;
   ~Histo1D();
};


class EqualSteps1DFunction : public Sampled1DFunction{

private:
  double xmin;
  double xmax;
  int nstep;
  double xstepinv;
  double xstep;
  double *fval;
  void init(const double Begin, const double End, const int Nval, const double *Values);

public:
  double Xmin() const { return xmin;}
  double Xmax() const { return xmax;}
  int Nval() const { return nstep;};

  //! reads a file eligible for a General1DFunction and checks that the steps are equal
  EqualSteps1DFunction(const string &FileName);

  EqualSteps1DFunction() {
    nstep=0;
    fval=0;
  }

  EqualSteps1DFunction(const double Begin, const double End, const int Nval, 
		       const double *Values);

  EqualSteps1DFunction(const double Begin, const double End, 
		       const double Step);

  EqualSteps1DFunction(const Sampled1DFunction &Input, int nstep);
  //! transforms into the repartition function assuming it is a p.d.f.
  void Cumulate();
  //! abcissa value of bin number i (first bin is 0)
  double Xval(const int i) const {return xmin+i*xstep;};
  //! function value at bin number i (first bin is 0)
  double Fval(const int i) const {return fval[i];};
  double &Fval(const int i) { return fval[i];};

  double Step() const {return xstep;}

  const double *values() const { return fval;}
  
  //! Rescales the function values
  void operator *= (const double Factor);
  
  EqualSteps1DFunction(const EqualSteps1DFunction&);
  EqualSteps1DFunction& operator = (const EqualSteps1DFunction &);
  double Value(const double &XWhere) const;
  void dump(ostream & stream = cout) const;
  void Write(const string &FileName) const;
  virtual ~EqualSteps1DFunction() { if (fval) delete [] fval;};

};


class General1DFunction : public Sampled1DFunction {

private:
  double xmin;
  double xmax;
  int nstep;
  double *xval;
  double *fval;
  // index in xval&fval of last X requested  to the interpolator.
  int lasti; 
  void init(const int Nval, const double *Values, const double *XValues);

public:
  double Xmin() const { return xmin;}
  double Xmax() const { return xmax;}
  //double Xmax(const double Factor, const int Nval) const;
  double Value(const double &XWhere) const;
  int Nval() const { return nstep;}
  double Xval(const int i) const {return xval[i];};
  double Fval(const int i) const {return fval[i];};
  double& Fval(const int i) { return fval[i];};

  void SetXval(const int i, const double& x) {xval[i]=x;};

  const double *values() const { return fval;}
  const double *xvals() const { return xval;}


  //! the second parameter defines the column in which the function value is to be read.
  General1DFunction(const string &Filename, const int ICol = 2);

  General1DFunction(const int Nval, const double *Values, const double *XValues);
  General1DFunction(const General1DFunction&);
  General1DFunction();
  
  //! rescales the x axis
  void RescaleX(const double Factor);
  
  //! Shift the function 
  void ShiftX(const double & Shift);
  
  //! Rescales the function values
  void operator *= (const double Factor);

  //! Multiply the function by something (that has an operator ())
  template<class Fun> void operator *= (Fun &F)
      {  
	  for (int i=0; i< nstep; ++i) fval[i] *=  F(xval[i]);
      }
	  
  //! replace all function values by their image through fun(val)
  template <class Fun> void ApplyFun( Fun &F)
      {
	  for (int i=0; i<nstep; ++i) fval[i] = F(fval[i]);
      }


  //! the integral of the function
  double Integral() const;

  //! the mean of the function
  double Mean() const;

  //! Integral of the function between MinBound and MaxBound
  double BoundedIntegral(const double & MinBound, const double & MaxBound) const;

  void dump(ostream & stream = cout) const;
  General1DFunction&  operator = (const General1DFunction &Right);

  void Read(const string &Filename, const int ICol = 2);
  bool Write(std::ostream &s) const;
  bool Write(const std::string &Filename) const;

  virtual ~General1DFunction();
  /* do not write operator = (const SampledFuction &Right) because 
     SampledFunction(const &SampledFunction) uses the fact that there is none  */
};

//! computes a product of 2 (1-dim)  terms
General1DFunction Product(const double Step, 
			  const Sampled1DFunction &Fun1, 
			  const Sampled1DFunction &Fun2);

//! computes a product of N (1-dim)  terms. pointers to Sampled1DFunction's expected! 
General1DFunction Product(const double Step,const int NTerms, ... );

//! Compute a product of N (1-dim) terms
General1DFunction Product(const double Step, vector<Sampled1DFunction*> terms);


class Constant1DFunction : public Sampled1DFunction 
{
 private:
  double value;
 public:
  Constant1DFunction(const double Value) : value(Value) {};
  double Value(const double &XWhere) const; //  { return value;} in .cc
  double Xmin() const { return -1e30;}
  double Xmax() const { return  1e30;}
  void dump(ostream & stream = cout) const { stream << " Constant1DFunction : value = " 
						    << value << ' ';}
};

class TopHat1DFunction : public Sampled1DFunction 
{
 private:
  double xmin;
  double xmax;
  double value;
 public:
  TopHat1DFunction(const double input_xmin, const double input_xmax, const double input_value) : 
    xmin(input_xmin),
    xmax(input_xmax),
    value(input_value) {};
  double Value(const double &XWhere) const { if(XWhere>=xmin && XWhere<xmax) return value; else return 0;}
  double Xmin() const { return xmin;}
  double Xmax() const { return xmax;}
  void dump(ostream & stream = cout) const { stream << " TopHat1DFunction : value = " 
						    << value << ' ';}
};



typedef double (AnalyticFunction)(const double &, const void*);

class Computed1DFunction : public Sampled1DFunction 
{
 private:
  AnalyticFunction *f;
  const void *someData;
  double x_min;
  double x_max;
 public:
  Computed1DFunction(AnalyticFunction *AFunction, const void* p=NULL) : 
    f(AFunction), someData(p),x_min(-1e30), x_max(1.e30)  {}
  double Value(const double &XWhere) const { return f(XWhere, someData);}
  double Xmin() const { return x_min;}
  double Xmax() const { return x_max;}
  void SetXmin(const double &x) { x_min=x;}
  void SetXmax(const double &x) { x_max=x;}
  
  void dump(ostream & stream = cout) const 
    { 
      stream << " Computed1DFunction, can only provide the address of the code :" << f 
	     << endl;
    }
};

int FastLocateInArray(const int lasti, const double *xval, 
		      const int nval, const double &XWhere);

class Grid2DXSlice;

class Grid2DFunction {

private:
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  int nstepx;
  int nstepy;
  int lastix,lastiy;
  double *xval;
  double *yval;
  double **fval;
  void init(const int NXval, const int NYval, const double * const *FValues, 
	    const double *XValues, const double *YValues);

  

public:

    double &val(const int ix, const int iy) {return fval[ix][iy];};
  double val(const int ix, const int iy) const {return fval[ix][iy];};
  
  bool add_slice(const double &X, const double *Fval, const double *Yval, const int NVal);
  
  double Xmin() const { return xmin;}
  double Xmax() const { return xmax;}
  double Ymin() const { return ymin;}
  double Ymax() const { return ymax;}
  //  double Xstep() const { return (Xmax()-Xmin())/(nstepx-1);}
  // double Ystep() const { return (Ymax()-Ymin())/(nstepy-1);}
  int NstepX() const             { return nstepx;}
  double Xval(const int i) const {return xval[i];}
  int NstepY() const             { return nstepy;}
  double Yval(const int i) const {return yval[i];}
  
  const double *xvals() const {return xval;}
  const double *yvals() const {return yval;}
  double **fvals() {return fval;} // cannot pass this as const, don't know why

  double Value(const double &XWhere, const double &YWhere) const;
  double BilinearValue(const double &XWhere, const double &YWhere) const;
  //double BiCubicInterpolation(const double &XWhere, const double &YWhere) const;
  double BiCubicConvolution(const double &XWhere, const double &YWhere) const;
  bool AddSlice(const double &XWhere, const Sampled1DFunction &NewSlice);

  Grid2DFunction() {fval = 0; xval = 0; yval = 0; nstepx=0;nstepy=0;lastix=0; lastiy=0;};
  Grid2DFunction(const string &Filename);
  Grid2DFunction(const int NXval, const int NYval, const double **Values, 
		 const double *XValues, const double *YValues);
  // General1DFunction ConstXSlice(const double XVal) const;
  void dump(ostream & stream = cout) const;
  
  Grid2DFunction(const Grid2DFunction&);
  bool Read(const string &Filename);
  Grid2DFunction&  operator = (const Grid2DFunction &Right);
  friend ostream & operator << ( ostream &stream, const Grid2DFunction &This);
  Grid2DXSlice ConstXSlice(const double XVal) const;
  void RescaleXSlice(const double Xval, const double Factor, const bool Verbose =false);
  bool Write(const string &Filename) const;
  virtual ~Grid2DFunction();
  double **GetNonConstValues(){return fval;};
};


class Grid2DXSlice : public Sampled1DFunction 
{
 private:
  const Grid2DFunction &Parent2DFunction;
  double xVal;
  
 public:
  Grid2DXSlice(const Grid2DFunction &F, const double X) : Parent2DFunction(F), xVal(X) {};
  double Value(const double &XWhere) const { return Parent2DFunction.Value(xVal, XWhere);}
  void dump (ostream & stream = cout) const {stream << "dump not implemented Grid2DXSlice";};
  double Xmin() const { return Parent2DFunction.Ymin();};
  double Xmax() const { return Parent2DFunction.Ymax();};
};

#endif /* SAMPLEDFUNCTION__H */
