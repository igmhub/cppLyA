#ifndef VUTILS__H
#define VUTILS__H


//! return the median of the array. The array is mixed up. 
/*!  Use next function if you want to keep your array untouched  */
double DArrayMedian(double *array, const int size);

//! same as above but does not mix up the array. (calls previous one on a copy of its input ... ). 
double DConstArrayMedian(const double *array, const int size);

//!
float FArrayMedian(float *array, const int size);

//! returns mean, median and rms of an array.
void Dmean_median_sigma(double *values, const int nval, double &mean,  double &median, double &sigma);

//! returns mean, median and rms of an array
void Fmean_median_sigma(float *values, const int nval, float &mean, float &median, float &sigma);
//! returns median and rms of an array
float Fmedian_sigma(float *values, const int nval, float &sigma);

//! computes the clipped-mean and sigma with cutting at k-sigma
double clipmean(double *values, int &nval, double &sigma, const double &k=3.5, const int niter=4);

#endif
