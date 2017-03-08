#ifndef MATVECT__H
#define MATVECT__H

#include <iostream>
#include <string>

#include <cstdlib> // for abort
#include <cstring> // for memset

class Vect;
class Mat;

// if you link with lapack uncomment this
//#define USELAPACK
// if you link with cfitsio uncomment this
//#define USECFITSIO
// if you link with cernlib uncomment this
#define USECERNLIB




#define MATVECT_CHECK_BOUNDS

// Routines for solving linear systems + inversion 
//==================================================================

#ifdef USELAPACK

// solving linear system A.X = B 
// Uses lapack dposv_
// Matrix A is assumed to be symmetric (you'll get a core dump if it is not
// (actually this is not compulsory with dposv (see 'man dposv') but we do it
// for consistency).
// You just need to fill half (n*(n+1)/2 parameters) of the matrix
// if you have filled matrix M parameters for which x>=y (with M(x,y)=...), use UorL="L"
// else use UorL="U"
// Matrix A is modified in this routine (also B which is at the end the solution X)
int cholesky_solve(Mat &A, Vect &B, const char* UorL); // "L" <=> A(i,j)!=0 for i>=j 

// for symmetric matrices
int general_solve(Mat &A, Vect &B, bool invert_A, const char* UorL);

// general (not symmetric)
int truly_general_solve(Mat &A, Vect &B);

double determinant_of_symmetric_matrix(const Mat &A); // assume A "U" triangle is filled
double determinant_using_cholesky_decomposition(const Mat &A, const char* UorL = "L"); 
double log_of_determinant_using_cholesky_decomposition(const Mat &A, const char* UorL = "L"); 

// Inverts matrix A using the factorization done in cholesky_solve
// Uses lapack dpotri_
// Matrix A is assumed to be symmetric (you'll get a core dump if it is not
// (actually this is not compulsory (see 'man dptri') but we do it
// for consistency).
// This routine must be called after cholesky_solve, make sure the value of UorL
// is the same as that used with dposv
int cholesky_invert(Mat &A, const char* UorL = "L"); // when cholesky_solve is called first, "L" <=> A(i,j)!=0 for i>=j 

//! calls cholesky_solve with a dummy RHS and then cholesky_invert.
int posdef_invert(Mat &A, const char* UorL);

/* Routine to "solve" a **posdef** system with constraints.
   The found solution minimizes:
   X**t A X  -  B**t X   under the constraints : C**t X = Y

   - dimensions : A(n,n), B(n), C(n,m), Y(m)
   - return value : 0 == OK
   - On exit the solution is in B and the Lagrange multiplier values in Y

   The technique used assumes that the unconstrained problem is posdef
   (because Cholesky factorization is used).  In particular, if the
   constraints are mandatory for the problem to have a solution at
   all, this is the wrong routine.  You then have e.g. to build your
   constraints into the matrix and call GeneralSolve.
*/
int cholesky_solve_with_constraints(Mat &A, Vect &B, 
				    Mat &Constraints, Vect& Y,
				    const char* UorL);

//! Diagonalize (symetric) matrix A using the dsysevd_ lapack routine. Assumes by default a "U" filled matrix.
int diagonalize(Mat const& A, Vect& eigenvals, Mat& eigenvects, const char* UorL="U");

/* computes  C = alpha*B*A + beta*C  */
void DGEMM(Mat &C, const Mat& A, const Mat& B, const double& alpha, const double& beta);
void DGEMM(Vect &C, const Vect& A, const Mat& B, const double& alpha, const double& beta);

#endif


/*
==================================================================
 returns A^t B C, usefull when there are lots of zeros in A
 - assumes B is symmetric
 - if B.SizeX()=1, assumes this is a diagonal matrix
==================================================================
*/
double scalar_product(const Vect& A, const Mat& B, const Vect& C);
double scalar_product(const Vect& V1,const Vect& V2);
double norm2(const Vect& V);
Vect vect_product(const Vect& V1,const Vect& V2);
Vect vect_product(const Vect& V1,const Vect& V2,const Vect& V3);

/*
==================================================================
Fast routine to compute largest eignevalue of a matrix
==================================================================
*/
double largest_eigenvalue(const Mat& M, Vect& eigenvector);


// Mat and Vect classes
//====================================================================


class Mat {
  friend class Vect;
 protected:
  double *data;
  size_t nx,ny;
  
  
 public:
  
  Mat() : data(NULL), nx(0), ny(0) {};  
  Mat(const size_t NX, const size_t NY);
  Mat(const Mat& other);
  ~Mat() { if(data) delete [] data;}
  
  void allocate(const size_t NX, const size_t NY);
  
  inline double& operator () (const size_t i, const size_t j) const {
#ifdef MATVECT_CHECK_BOUNDS
    if (i>=nx || j>=ny) { 
      std::cout << "Mat::operator () overflow i,j nx,ny " 
	 << i << "," << j << " "
	   << nx << "," << ny << " "
	   << std::endl;
      abort();
    }
#endif
    return data[i+j*nx];
  }
  
  
  size_t SizeX() const { return nx;}
  size_t SizeY() const { return ny;}
  
  void Zero() {memset(data, 0, nx*ny*sizeof(double));};
  void Identity();

  const double* Data() const {return data;};
  double* NonConstData() {return data;};

  double sum() const;

  
  friend std::ostream& operator << (std::ostream &stream, const Mat &m)
    { m.writeASCII(stream); return stream;}
  
  // get a block of this matrix as a new matrix
  // size (x_max-x_min+1)*(y_max-y_min+1) (both min and max included)
  // remember  0 <= x < nx ,  0 <= y < ny
  Mat SubBlock
    (size_t x_min,size_t x_max,size_t y_min,size_t y_max) const;
  
  Mat WithoutRows(size_t y_min,size_t y_max) const;
  Mat WithoutColumns(size_t x_min,size_t x_max) const;

#ifdef USECFITSIO
  // i/o in fits for matrices
  int readFits(const std::string &FitsName,int hdu=1);
  int writeFits(const std::string &FitsName) const;
#endif
  
  int readASCII(std::istream& Stream);
  int readASCII(const std::string &FileName);
  int writeASCII(std::ostream& Stream) const;
  int writeASCII(const std::string &FileName) const;
  
  void Symmetrize(const char* UorL); // "L" <=> (*this)(i,j)!=0 for i>=j 
  
#ifdef USECERNLIB
  // inverts a symetric posdef matrix using DSINV  CERNLIB's routine
  int SymMatInvert();
#endif
  // inverts a posdef matrix using Cholesky factorization from lapack.
  int CholeskyInvert(const char *UorL, bool abort_if_failed=true);
  
  // operators
  Mat operator +(const Mat& Right) const;
  Mat operator -(const Mat& Right) const;
  Mat operator *(const Mat& Right) const;
  Vect operator *(const Vect& Right) const;
  Mat & operator =(const Mat& Right);
  
  void operator +=(const Mat& Right);
  void operator -=(const Mat& Right);
  void operator *=(const Mat& Right);
 
  
  Mat operator *(const double Right) const;
  friend Mat operator *(const double Left, const Mat &Right);
  void operator *=(const double Right);

  operator double() const;
  operator Vect() const;
  Mat transposed() const;

};

#ifdef USECFITSIO
void AddFitsKey(const std::string &FileName, 
		const std::string &KeyName, 
		const double &KeyVal);
double  ReadFitsKey(const std::string &FileName, 
		    const std::string &KeyName);

#endif


class Vect {

  friend class Mat;

 protected:
  
  double *data;
  size_t n;
  
  
 public:
  
  Vect() : data(NULL), n(0) {};
  Vect(const size_t N);
  Vect(const Vect&);
  ~Vect() { if(data) delete [] data;}

  void allocate(const size_t N);
  
  inline double operator () (const size_t i) const {
#ifdef MATVECT_CHECK_BOUNDS
    if (i>=n ) {
      std::cout << "Vec::operator () overflow i,n " 
	   << i << "," << n << std::endl;
      abort();
    }
#endif
    return data[i];
  }
  
  inline double& operator () (const size_t i) {
#ifdef MATVECT_CHECK_BOUNDS
    if (i>=n) {
      std::cout << "Vec::operator () overflow i,n " 
	   << i << "," << n << std::endl;
      abort();
    }
#endif
    return data[i];
  }
  
  
  size_t Size() const { return n;}
  size_t size() const { return n;}
  bool empty() const {return n==0;}
  void Zero() {memset(data, 0, n*sizeof(double));};

  const double* Data() const {return data;};
  double* NonConstData() {return data;};

  int writeASCII(std::ostream& Stream) const;
  int readASCII(std::istream& Stream);
  
  friend std::ostream& operator << (std::ostream &stream, const Vect &v)
    { v.writeASCII(stream); return stream;}

  // operators
  Vect operator +(const Vect& Right) const;
  Vect operator -(const Vect& Right) const;
  Vect & operator =(const Vect& Right);
  void   operator =(const double &Val);
  double operator *(const Vect& Right) const; // scalar product
  
  void operator +=(const Vect& Right);
  void operator +=(double Right);
  void operator -=(const Vect& Right);
  
  Vect operator *(const double Right) const;
  friend Vect operator *(const double Left, const Vect &Right);
  void operator *=(const double Right);

  double sum() const;

  Mat transposed() const;
  operator Mat() const; // useful for diadic products
  operator double() const;
};

typedef float MFLOAT;

// to save memory
class SymmetricFloatMatContainer {

 protected:
  
  MFLOAT  *data;
  size_t n;
  
  
 public:
  SymmetricFloatMatContainer() : data(NULL), n(0) {};
  SymmetricFloatMatContainer(const Mat& other,const char UoL='L'); // "L" <=> (other(i,j)!=0) if i>=j 
  SymmetricFloatMatContainer(const SymmetricFloatMatContainer&);
  ~SymmetricFloatMatContainer();
  operator Mat() const;
  void Load(const Mat& other,const char UoL='L');
};
#endif /*MATVECT__H */
