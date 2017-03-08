#include <iostream>
#include "histo2d.h"
#include <math.h> /* for floor */
#include <string.h> /* for memset*/

using namespace std;

Histo2d::Histo2d(int nnx, const double& mminx, const double& mmaxx, int nny, const double& mminy, const double& mmaxy) {
  nx = nnx;
  ny = nny;
  minx = mminx;
  miny = mminy;
  if (mmaxx!= mminx) 
    scalex = nx/(mmaxx-mminx); 
  else 
    {
      cerr << " Histo2d: minx = maxx requested" << endl;
      scalex = 1.0;
    }
  if (mmaxy != mminy)
    scaley = ny/(mmaxy-mminy);
  else
    {
      cerr << " Histo2d : maxy = miny requested" << endl;
      scaley = 1.0;
    }
  data = new double[nx*ny];
  memset(data, 0, nx*ny*sizeof(double));
}

Histo2d::Histo2d(const Histo2d &Other) {
  memcpy(this, &Other, sizeof(Histo2d));
  data = new double[nx*ny];
  memcpy(this->data, Other.data, nx*ny*sizeof(double));
}

bool Histo2d::indices(const double &X, const double &Y, int &ix, int &iy) const {
  ix = (int) floor(( X - minx)*scalex);
  if (ix <0 || ix >= nx) return false;
  iy = (int) floor((Y - miny)*scaley);
  return (iy >=0 && iy < ny);
}

double Histo2d::BinContent(const double &X, const double &Y) const {
  int ix, iy;
  if (indices(X,Y,ix,iy)) return data[iy + ny*ix];
  std::cout << " Histo2D::BinContent outside limits requested " << std::endl;
  return -1e30;
}
