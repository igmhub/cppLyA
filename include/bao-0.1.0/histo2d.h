#ifndef HISTO2D__H
#define HISTO2D__H

#include <cmath>

class Histo2d {
  
  public:
  
  double *data;
  int nx,ny;
  double minx,miny;
  double scalex, scaley;

  bool indices(const double &X, const double &Y, int &ix, int &iy) const;
  Histo2d() {data=NULL;}
  Histo2d(int nx, const double& minx, const double& maxx, int ny, const double& miny, const double& maxy);
  Histo2d(const Histo2d &Other);
  inline void Fill(const double& x, const double& y, const double& weight=1.) {
    int ix = (int) floor(( x - minx)*scalex);
    if (ix <0 || ix >= nx) return;
    int iy = (int) floor((y - miny)*scaley);
    if (iy <0 || iy >= ny) return;
    data[iy + ny*ix] += weight;
  }
  inline int index_1d(const double &x, const double &y) const {
    int ix = (int) floor(( x - minx)*scalex);
    if (ix <0 || ix >= nx) return -1;
    int iy = (int) floor((y - miny)*scaley);
    if (iy <0 || iy >= ny) return -1;
    return iy + ny*ix;
  }
  void BinWidth(double &Hdx, double &Hdy) const { Hdx = 1./scalex; Hdy = 1./scaley;}
  double BinContent(const double &X, const double &Y) const;
  ~Histo2d() { if (data) delete [] data;}

 private:
  void operator = (const Histo2d &Right);
  
};

#endif /* HISTO2D__H */
