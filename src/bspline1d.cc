
#include <bspline1d.h>
#include <bspline.c>


double EqualStepBSpline1D::Value(const double& x, bool zero_outside_range) const {
  if(zero_outside_range) {
    if(x<xmin) return 0;
    if(x>xmax) return 0;  
  }
  double rx = scale*x+offset;
  size_t n = params.size();
  size_t imin = (rx>2) ? rx-2 : 0;
  size_t imax = (imin+3<n) ? imin+3 : n;
  double v=0;
  const double *p = params.Data()+imin;
  for(size_t i=imin;i<imax;i++,p++)
    v += Bspline3(rx,i)*(*p);
  return v;
}

Vect EqualStepBSpline1D::ParamDerivatives(const double &x) const {
  size_t n = params.size();
  Vect H(n);
  if(x<xmin) return H;
  if(x>xmax) return H;
  double rx = scale*x+offset;
  size_t imin = (rx>2) ? rx-2 : 0;
  size_t imax = (imin+3<n) ? imin+3 : n;
  for(size_t i=imin;i<imax;i++)
    H(i) = Bspline3(rx,i);
  return H;
}

bool EqualStepBSpline1D::Fit(size_t n, const double*x, const double* y , bool adapt) { // assumes it is ordered with equal x steps
  double i_xmin = x[0];
  double i_ymin = x[n-1];
  if(adapt) {
    SetRange(i_xmin,i_ymin);
    SetNPar(n); // pas sur
  }
  size_t npar = params.size();
  Mat A(npar,npar);
  Vect B(npar);
  Vect H(npar);
  
  double rx;
  
  
  for(size_t i=0;i<n;i++,x++,y++) {

    if(*x<xmin) continue;
    if(*x>xmax) continue;
    

    rx = scale*(*x)+offset;
    
    size_t jmin = (rx>2) ? rx-2 : 0;
    size_t jmax = (jmin+3<npar) ? jmin+3 : npar;
    H.Zero();
    for(size_t j=jmin;j<jmax;j++) { // loop on params
      H(j) = Bspline(3,rx,j);
    }
    
    const double *hj = H.Data();
    double *a = A.NonConstData();
    double *b = B.NonConstData();
    for(size_t j=0;j<npar;j++,hj++,b++) {
      if( ! *hj) {a+=npar; continue;}
      *b += (*y)*(*hj);
	const double *hi = hj;
	a+=j;
	for(size_t i=j;i<npar;i++,hi++,a++)
	  *a += (*hi)*(*hj);//  Mat(H)*H.transposed();
    }
  }
  
  //A.writeFits("A.fits");
  int status = cholesky_solve(A,B,"L");
  if(status != 0) {
    cout << "ERROR EqualStepBSpline1D::Fit cannot solve system" << endl;
    return false;
  } 
  params=B;
  
  return true;
}
bool EqualStepBSpline1D::FitWithWeights(size_t n, const double*x, const double* y , const double* w , bool adapt) { // assumes it is ordered with equal x steps
  double i_xmin = x[0];
  double i_ymin = x[n-1];
  if(adapt) {
    SetRange(i_xmin,i_ymin);
    SetNPar(n); // pas sur
  }
  size_t npar = params.size();
  Mat A(npar,npar);
  Vect B(npar);
  Vect H(npar);
  
  double rx;
  
  
  for(size_t i=0;i<n;i++,x++,y++,w++) {

    if(*x<xmin) continue;
    if(*x>xmax) continue;
    

    rx = scale*(*x)+offset;
    
    size_t jmin = (rx>2) ? rx-2 : 0;
    size_t jmax = (jmin+3<npar) ? jmin+3 : npar;
    H.Zero();
    for(size_t j=jmin;j<jmax;j++) { // loop on params
      H(j) = Bspline(3,rx,j);
    }
    
    const double *hj = H.Data();
    double *a = A.NonConstData();
    double *b = B.NonConstData();
    for(size_t j=0;j<npar;j++,hj++,b++) {
      if( ! *hj) {a+=npar; continue;}
      *b += (*w)*(*y)*(*hj);
	const double *hi = hj;
	a+=j;
	for(size_t i=j;i<npar;i++,hi++,a++)
	  *a += (*w)*(*hi)*(*hj);//  Mat(H)*H.transposed();
    }
  }
  
  //A.writeFits("A.fits");
  int status = cholesky_solve(A,B,"L");
  if(status != 0) {
    cout << "ERROR EqualStepBSpline1D::Fit cannot solve system" << endl;
    return false;
  } 
  params=B;
  
  return true;
}

bool EqualStepBSpline1D::Fit(const General1DFunction &func, bool adapt) {
  return Fit(func.Nval(),func.xvals(),func.values(),adapt);
}

////////

extern "C" {
  double divdifd_(const double *, const double *, const int *, const double *, const int*);
}

double BSpline1D::RemappedX(const double& x) const {
  int n = original_x_nodes.size();
  if(n<2) {
    cout <<" ERROR in BSpline1D::RemappedX n=" << n << endl;
    abort();
  }
  int deg = 2;
  return divdifd_(remapped_x_nodes.Data(), original_x_nodes.Data(), &n, &x, &deg);
}

void BSpline1D::Init(size_t n, const double* x) {
  original_x_nodes.allocate(n);
  remapped_x_nodes.allocate(n);
  double *o = original_x_nodes.NonConstData();
  double *r = remapped_x_nodes.NonConstData();
  for(size_t i=0;i<n;i++,x++,o++,r++) {
    *o = *x;
    *r = i;
  }
  SetRange(remapped_x_nodes(0),remapped_x_nodes(n-1));
  SetNPar(n);
  //cout << "BSpline1D::Init with n = " << n << " xmin xmax = " << remapped_x_nodes(0) << " " << remapped_x_nodes(n-1) << endl;
}

bool BSpline1D::Fit(size_t n, const double* x, const double* y, bool adapt) { // assumes it is ordered with equal x steps
  if(! adapt) {
    cout << "WARNING  BSpline1D::Fit adapt anyway" << endl;
  }  
 
  Init(n,x);
  
  bool status = EqualStepBSpline1D::Fit(n,remapped_x_nodes.Data(),y);
  if( ! status) cout << "ERROR BSpline1D::Fit, EqualStepBSpline1D::Fit failed" << endl;
  return status;
}

double BSpline1D::Value(const double& x, bool zero_outside_range) const {
  return EqualStepBSpline1D::Value(RemappedX(x),zero_outside_range);
}

bool BSpline1D::Fit(const General1DFunction &func, bool adapt) {
  return Fit(func.Nval(),func.xvals(),func.values(),adapt);
}
