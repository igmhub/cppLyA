#ifndef HISTO1D__H
#define HISTO1D__H

#include <iostream>
#include <math.h>


class Histo1d {

 public:
  double *data;
  int nx;
  double minx;
  double scalex;

 public:
  Histo1d() { data=NULL;}
  Histo1d(int Nx, double Minx, double Maxx);
  void Fill(const double& X, const double& weight=1.)
      {int bin = int(floor((X - minx)*scalex)); 
       if (bin<0 || bin>=nx) return; 
       data[bin] += weight;
      }
  int BinAt(const double& X) const;
  double BinCenter(int bin) const;
  double BinWidth() const { return (1./scalex);}
 void dump(std::ostream &stream) const;
  friend std::ostream& operator << (std::ostream &stream, const Histo1d &h)
    { h.dump(stream); return stream;}
  
  ~Histo1d() { if (data) delete [] data;}
};

#endif /* HISTO1D__H */
