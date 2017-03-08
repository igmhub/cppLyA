#ifndef DIST_MAT__H
#define DIST_MAT__H

#include <iostream>

class Vect; 
void save_comp_mat_results(const string& filename, double D) 
double alpha(int spec1, int spec2, int w1, int w2);
double beta(int spec1, int spec2, int w1, int w2);
double alpha_beta(int spec1, int spec2, int w1, int w2);
double nu(int spec1, int spec2, int w1, int w2);
#endif
