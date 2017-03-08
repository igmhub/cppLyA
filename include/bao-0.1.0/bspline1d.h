
#ifndef BSPLINE1D_H 
#define BSPLINE1D_H

#include <sampledfunction.h>
#include <matvect.h>

class EqualStepBSpline1D {
private :
  
  double xmin;
  double xmax;
  double scale;
  double offset;
  
public :
  Vect params;
  EqualStepBSpline1D(){   
    xmin=0;
    xmax=1;
    scale=1;
    offset=0;
  }
  void SetNPar(int nparams) {
    if(xmax>xmin)
      scale = (nparams-2)/(xmax-xmin);
    offset = -xmin*scale+2;
    if(nparams != int(params.size()))
      params.allocate(nparams);
  }
  void SetRange(const double& i_xmin, const double& i_xmax) {
    xmin = i_xmin;
    xmax = i_xmax;
    SetNPar(params.size());
  }
  
  EqualStepBSpline1D(const double& i_xmin, const double& i_xmax, int nparams)
  {
    SetRange(i_xmin,i_xmax);
    SetNPar(nparams);
  } 
  virtual double Value(const double& x, bool zero_outside_range = false) const;
  virtual Vect ParamDerivatives(const double& x) const;
  
  double Xmin() const{return xmin;}
  double Xmax() const{return xmax;}
  virtual bool Fit(size_t n, const double*x, const double* y , bool adapt = true);
  virtual bool FitWithWeights(size_t n, const double*x, const double* y, const double* w, bool adapt = true);
  virtual bool Fit(const General1DFunction &func, bool adapt = true);
  EqualStepBSpline1D(const General1DFunction &func){ 
    if( ! Fit(func)) abort();
  }  
  
};

class BSpline1D : public EqualStepBSpline1D {

 private : 
  
  Vect original_x_nodes;
  Vect remapped_x_nodes;
  

public :
  
 BSpline1D() : EqualStepBSpline1D() {};
  double Xmin() const{return original_x_nodes(0);}
  double Xmax() const{return original_x_nodes(original_x_nodes.size()-1);}
  virtual double Value(const double& x, bool zero_outside_range = false) const;
  virtual bool Fit(size_t n, const double*x, const double* y, bool adapt = true);
  double RemappedX(const double& x) const;
  void Init(size_t n, const double* x);
  void Init(const Vect& X) {Init(X.size(),X.Data());}
  virtual bool Fit(const General1DFunction &func, bool adapt = true);
  BSpline1D(const General1DFunction &func){ 
    if( ! Fit(func)) abort();
  }  
};
#endif
