#include <iostream>
#include "histo1d.h"
#include <string.h> /* for memset */

using namespace std;

Histo1d::Histo1d(int Nx, double Minx, double Maxx) : nx(Nx), minx(Minx)
{
if (Maxx != minx)
  scalex = Nx/(Maxx - minx);
else 
  {
  cerr << " Histo1d :: minx = maxx requested" << endl;
  scalex = 1;
  }
data = new double[nx];
memset(data,0,nx*sizeof(double));
}

int Histo1d::BinAt(const double& X) const { 
  int bin = int(floor((X - minx)*scalex));
  if (bin<0 || bin>=nx) return -1;
  return bin;
} 

double Histo1d::BinCenter(int bin) const{
  if (bin<0 || bin>=nx) {
    cerr << "Histo1d::BinCenter ERROR out of range" << endl;
    return 0;
  }
  return (minx + (bin+0.5)/scalex);
}

void Histo1d::dump(ostream &stream) const {
  stream << nx << " " << minx << " " << minx+nx/scalex << endl;
  for(int bin=0;bin<nx;++bin)
    stream << bin << " " << BinCenter(bin) << " " << data[bin] << endl;
}
