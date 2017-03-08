#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <matvect.h>
#include <fileutils.h>

#ifdef USECFITSIO
#include <fitsio.h>
#endif


using namespace std; 

#ifdef USECERNLIB
#define dsinv dsinv_
#define dfact dfact_
#define dfinv dfinv_
#define dfeqn dfeqn_
#define eisrs1 eisrs1_


// using cernstuff (from cernlib)
extern "C" 
{
  void dsinv(int *N, double *A, int *IDIM, int *IFAIL);
  void dfact(int *N, double *A, int *idim, double *r, int *ifail, double *Det, int *jfail);
  void dfinv(int *n, double *A, int *idim, double *r);
  void dfeqn(int *n, double *a, int *idim, double *r, int *k, double *b);
  void eisrs1(int *NM,int *N,double *AR,double *WR,double *ZR,int *IERR,double *WORK);
}
#endif

#ifdef USELAPACK
// using lapack. 
// Signatures (with a fortran accent) can be obtained on many linux systems from e.g. "man dposv" 
extern "C" {
  void dposv_(const char *, int *, int *, double *, int *, double *, int *, int *);
  void dpotri_(const char *, int *, double *, int *, int *);
  void dsysv_(const char*, int* n, int* nrhs, 
	      double* a, int* lda, int* ipiv, 
	      double* b, int* ldb, 
	      double* work, int* lwork, int* info);
  void dsyev_(const char* jobz,   /* 'N': eigenval only. 'V': eigenval + eigenvec */
	      const char* uplow,  /* 'U': uppertriangle of A. 'L': lower triangle of A */
	      int* n,             /* order of A */
	      double* a,          /* the matrix */
	      int* lda,           /* leading dimension of A */
	      double* w,          /* eigenvalues: dim = n*/
	      double* workspace,  /* workspace. Allocated by the user. */
	      int* lwork,         /* dimension of the workspace */
	      int* info);         /* successful if zero. ==-i : the i-th element had an illegal val.
				      >0 failed to converge */
  void dsyevd_(const char* jobz,   /* 'N': eigenval only. 'V': eigenval + eigenvec */
	      const char* uplow,  /* 'U': uppertriangle of A. 'L': lower triangle of A */
	      int* n,             /* order of A */
	      double* a,          /* the matrix */
	      int* lda,           /* leading dimension of A */
	      double* w,          /* eigenvalues: dim = n*/
	      double* workspace,  /* workspace. Allocated by the user. */
	      int* lwork,         /* dimension of the workspace */
	       int *iworkspace,   /* another workspace user allocated */
	       int *liwork,       /* its size */
	      int* info);         /* successful if zero. ==-i : the i-th element had an illegal val.
				      >0 failed to converge */
  void dpotrs_(char *UPLO, int *N, int *NRHS, double *A, int *LDA, 
	       double *B, int *LDB, int *INFO);
  void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFO);
  void dsytri_(const char*, int* n, double* a , int* lda, int* ipiv, double* work, int* info);
  

  /*
    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

*  DGESV computes the solution to a real system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  The LU decomposition with partial pivoting and row interchanges is
*  used to factor A as
*     A = P * L * U,
*  where P is a permutation matrix, L is unit lower triangular, and U is
*  upper triangular.  The factored form of A is then used to solve the
*  system of equations A * X = B.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N coefficient matrix A.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          The pivot indices that define the permutation matrix P;
*          row i of the matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS matrix of right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
*                has been completed, but the factor U is exactly
*                singular, so the solution could not be computed.
*
*  =====================================================================
  */
  void dgesv_(int *n, int *nrhs, 
	 double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

  /*
    SUBROUTINE DGETRF( M,	N, A, LDA, IPIV, INFO )
    
    INTEGER	     INFO, LDA,	M, N    
    INTEGER	     IPIV( * )
    DOUBLE	     PRECISION A( LDA, * )
  */
  void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
  
  /*
    SUBROUTINE DGETRS( TRANS, N, NRHS, A,	LDA, IPIV, B, LDB, INFO	)

    CHARACTER	     TRANS
    INTEGER	     INFO, LDA,	LDB, N,	NRHS
    INTEGER	     IPIV( * )
    DOUBLE	     PRECISION A( LDA, * ), B( LDB, * )

    TRANS	  (input) CHARACTER*1
    Specifies the	form of	the system of equations:
    = 'N':  A * X	= B  (No transpose)
    = 'T':  A'* X	= B  (Transpose)
    = 'C':  A'* X	= B  (Conjugate	transpose = Transpose)
    
  */

  
  void dgetrs_(char *NTC, int* n, int *nrhs, double *a, int* lda, int* ipiv, double *b, int *ldb, int *info);


  /*
    SUBROUTINE DGERFS( TRANS, N, NRHS, A,	LDA, AF, LDAF, IPIV, B,	LDB, X,	LDX,
    FERR, BERR, WORK, IWORK, INFO )
    
    CHARACTER	     TRANS
    INTEGER	     INFO, LDA,	LDAF, LDB, LDX,	N, NRHS
    INTEGER	     IPIV( * ),	IWORK( * )
    DOUBLE	     PRECISION A( LDA, * ), AF(	LDAF, *	), B( LDB, * ),	BERR(
     * ), FERR(	* ), WORK( * ),	X( LDX,	* )

     PURPOSE
     DGERFS improves the computed solution	to a system of linear equations	and
     provides error bounds	and backward error estimates for the solution.
     
     ARGUMENTS
     
     TRANS	  (input) CHARACTER*1
     Specifies the	form of	the system of equations:
     = 'N':  A * X	= B	(No transpose)
     = 'T':  A**T * X = B	(Transpose)
     = 'C':  A**H * X = B	(Conjugate transpose = Transpose)
     
     N	  (input) INTEGER
     The order of the matrix A.  N	>= 0.
     
     NRHS	  (input) INTEGER
     The number of	right hand sides, i.e.,	the number of columns of the
     matrices B and X.  NRHS >= 0.
     
     A	  (input) DOUBLE PRECISION array, dimension (LDA,N)
     The original N-by-N matrix A.
     
     LDA	  (input) INTEGER
     The leading dimension	of the array A.	 LDA >=	max(1,N).
     
     AF	  (input) DOUBLE PRECISION array, dimension (LDAF,N)
     The factors L	and U from the factorization A = P*L*U as computed by
     DGETRF.
     
     LDAF	  (input) INTEGER
     The leading dimension	of the array AF.  LDAF >= max(1,N).
     
     IPIV	  (input) INTEGER array, dimension (N)
     The pivot indices from DGETRF; for 1<=i<=N, row i of the matrix was
     interchanged with row	IPIV(i).
     
     B	  (input) DOUBLE PRECISION array, dimension (LDB,NRHS)
     The right hand side matrix B.
     
     LDB	  (input) INTEGER
     The leading dimension	of the array B.	 LDB >=	max(1,N).
     
     X	  (input/output) DOUBLE	PRECISION array, dimension (LDX,NRHS)
     On entry, the	solution matrix	X, as computed by DGETRS.  On exit,
     the improved solution	matrix X.
     
     LDX	  (input) INTEGER
     The leading dimension	of the array X.	 LDX >=	max(1,N).
     
     FERR	  (output) DOUBLE PRECISION array, dimension (NRHS)
     The estimated	forward	error bounds for each solution vector X(j)
     (the j-th column of the solution matrix X).  If XTRUE	is the true
     solution, FERR(j) bounds the magnitude of the	largest	entry in
     (X(j)	- XTRUE) divided by the	magnitude of the largest entry in
     X(j).	 The quality of	the error bound	depends	on the quality of the
     estimate of norm(inv(A)) computed in the code; if the	estimate of
     norm(inv(A)) is accurate, the	error bound is guaranteed.
     
     BERR	  (output) DOUBLE PRECISION array, dimension (NRHS)
     The componentwise relative backward error of each solution vector
     X(j) (i.e., the smallest relative change in any entry	of A or	B
     that makes X(j) an exact solution).

     WORK	  (workspace) DOUBLE PRECISION array, dimension	(3*N)
     
     IWORK	  (workspace) INTEGER array, dimension (N)
     
     INFO	  (output) INTEGER
     = 0:	successful exit
     < 0:	if INFO	= -i, the i-th argument	had an illegal value

     PARAMETERS
     
     ITMAX	is the maximum number of steps of iterative refinement.
  */
  
  
  void dgerfs_(char *trans, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, double *b,	int *ldb, double *x, int *ldx , double *ferr, double *berr, double *work, int *iwork, int *info);

  /*
   *  Purpose
   *  =======
   *
   *  DSPSV computes the solution to a real system of linear equations
   *     A * X = B,
   *  where A is an N-by-N symmetric matrix stored in packed format and X
   *  and B are N-by-NRHS matrices.
   *
   *  The diagonal pivoting method is used to factor A as
   *     A = U * D * U**T,  if UPLO = 'U', or
   *     A = L * D * L**T,  if UPLO = 'L',
   *  where U (or L) is a product of permutation and unit upper (lower)
   *  triangular matrices, D is symmetric and block diagonal with 1-by-1
   *  and 2-by-2 diagonal blocks.  The factored form of A is then used to
   *  solve the system of equations A * X = B.
   *
   *  Arguments
   *  =========
   *
   *  UPLO    (input) CHARACTER*1
   *          = 'U':  Upper triangle of A is stored;
   *          = 'L':  Lower triangle of A is stored.
   *
   *  N       (input) INTEGER
   *          The number of linear equations, i.e., the order of the
   *          matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
   *          On entry, the upper or lower triangle of the symmetric matrix
   *          A, packed columnwise in a linear array.  The j-th column of A
   *          is stored in the array AP as follows:
   *          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
   *          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
   *          See below for further details.
   *
   *          On exit, the block diagonal matrix D and the multipliers used
   *          to obtain the factor U or L from the factorization
   *          A = U*D*U**T or A = L*D*L**T as computed by DSPTRF, stored as
   *          a packed triangular matrix in the same storage format as A.
   *
   *  IPIV    (output) INTEGER array, dimension (N)
   *          Details of the interchanges and the block structure of D, as
   *          determined by DSPTRF.  If IPIV(k) > 0, then rows and columns
   *          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1
   *          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,
   *          then rows and columns k-1 and -IPIV(k) were interchanged and
   *          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and
   *          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and
   *          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2
   *          diagonal block.
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the N-by-NRHS right hand side matrix B.
   *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
   *                has been completed, but the block diagonal matrix D is
   *                exactly singular, so the solution could not be
   *                computed.
   */
  void dspsv_(const char *, int *, int *, double*, int *, double *, int *, int *);
}

/*
The Lapack routine dspsv solves the linear system of equations Ax=b, where A is a symmetric matrix in packed storage format. However, there appear to be no Lapack functions that compute the determinant of such a matrix. We need to compute the determinant, for instance, in order to compute the multivariate normal density function. The dspsv function performs the factorization A=UDU’, where U is a unitriangular matrix and D is a block diagonal matrix where the blocks are of dimension 1×1 or 2×2. In addition to the solution for x, the dspsv function also returns the matrices U and D. The matrix D may then be used to compute the determinant of A. Recall from linear algebra that det(A) = det(UDU’) = det(U)det(D)det(U’). Since U is unitriangular, det(U) = 1. The determinant of D is the product of the determinants of the diagonal blocks. If a diagonal block is of dimension 1×1, then the determinant of the block is simply the value of the single element in the block. If the diagonal block is of dimension 2×2 then the determinant of the block may be computed according to the well-known formula b11*b22-b12*b21, where bij is the value in the i’th row and j’th column of the block. The following C code snip demonstrates the procedure.
*/


/*
** A and D are upper triangular matrices in packed storage
** A[] = a00, a01, a11, a02, a12, a22, a03, a13, a23, a33, ...
** use the following macro to address the element in the
** i'th row and j'th column for i <= j
*/
#define UMAT(i, j) (i + j * ( j + 1 ) / 2)


double determinant_using_cholesky_decomposition(const Mat &A, const char* UorL) {
  Mat M=A;
  Vect V(M.SizeX());
  if (cholesky_solve(M,V,UorL) != 0) {
    cout << "error in determinant_using_cholesky_decomposition" << endl;
    abort();
  }
  // now M is such that MM^* = A
  double det=1;
  for(size_t k =0; k<M.SizeX(); ++k) det *= M(k,k);
  det *= det; // square for determinant of A
  return det;
}
double log_of_determinant_using_cholesky_decomposition(const Mat &A, const char* UorL) {
  Mat M=A;
  Vect V(M.SizeX());
  if (cholesky_solve(M,V,UorL) != 0) {
    cout << "error in determinant_using_cholesky_decomposition" << endl;
    abort();
  }
  // now M is such that MM^* = A
  double log_det=0;
  for(size_t k =0; k<M.SizeX(); ++k) log_det += log(M(k,k));
  log_det *= 2; // det(A)=det(M)^2
  return log_det;
}

double determinant_of_symmetric_matrix(const Mat &A) {
  
  int n = A.SizeX();
  double *a = new double[(n*(n+1))/2];
  double *b = new double[n];
  int *ipiv = new int[n];

  const char* UorL = "U";
  
  // need to pack the data
  int k=0;
  for(int i=0;i<n;i++) 
    for(int j=i;j<n;j++,k++)
      a[ UMAT(i, j) ] = A(i,j);
  
  int info;
  int one = 1;
  
  
  /*
  ** solve Ax=B using A=UDU' factorization, D is placed in A
  ** where A represents a qxq matrix, b a 1xq vector
  */
  dspsv_(UorL, &n, &one, a, ipiv, b, &n, &info);
  if( info > 0 ) { /*issue warning, system is singular*/ }
  if( info < 0 ) { /*issue error, invalid argument*/ }
  
  /*
  ** compute the determinant det = det(A)
  ** if ipiv[i] > 0, then D(i,i) is a 1x1 block diagonal
  ** if ipiv[i] = ipiv[i-1] < 0, then D(i-1,i-1),
  ** D(i-1,i), and D(i,i) form a 2x2 block diagonal
  */
  
  double *D = a;
  double det = 1.0;
  for(int  i = 0; i < n; i++ ) {
    if( ipiv[ i ] > 0 ) {
      det *= D[ UMAT(i,i) ];
    } else if( i > 0 && ipiv[ i ] < 0 && ipiv[ i-1 ] == ipiv[ i ] ) {
      det *= D[ UMAT(i,i) ] * D[ UMAT(i-1,i-1) ] - D[ UMAT(i-1,i) ] * D[ UMAT(i-1,i) ];
    }
  }
  
  delete [] ipiv;
  delete [] a;
  delete [] b;
  return det;
}

int cholesky_solve(Mat &A, Vect &B, const char* UorL)
{
#ifdef FNAME
  cout << " >  cholesky_solve" << endl;
#endif  

  if(A.SizeX() != A.SizeY() || A.SizeX()<=0) {
    cout << "error in matvect.cc, cholesky_solve Matrix A must be symmetric and you have nx,ny = "
	 << A.SizeX() << "," << A.SizeY() << endl;
    abort();
  }
  if(B.Size() != A.SizeX()) {
    cout << "error in matvect.cc, cholesky_solve Vector B must have a dimention B.Size()=A.SizeY() and you have B.Size(),A.SizeY() = "
	 << B.Size() << "," << A.SizeY() << endl;
    abort();
  }

  if (A.SizeX() == 0)
    {
      cout << " error in matvect.cc, cholesky_solve : Matrix dimension is 0 "<< endl;
      abort();
    }

  double *a = A.NonConstData();
  double *b = B.NonConstData();
  
  int n = B.Size();
  int nhrs = 1, info = 0;

  dposv_(UorL, &n, &nhrs, a, &n, b, &n, &info);

  if (info != 0) 
    cerr << " cholesky_solve(" << a << "," << b << "," << n
	 << ") : Warning: Cholesky factorization failure . info =" 
	 << info <<  " (>0 is not pos.def)" << endl;

  return info;
}

int cholesky_invert(Mat &A, const char* UorL)
{  
  if(A.SizeX() != A.SizeY() || A.SizeX()<=0) {
    cout << "error in matvect.cc, cholesky_invert Matrix A must be symmetric and you have nx,ny = "
	 << A.SizeX() << "," << A.SizeY() << endl;
    abort();
  }
  
  int info = 0;
  double *a = A.NonConstData();
  int n = A.SizeX();

  //  Now invert using the factorization done in dposv_

  dpotri_(UorL, &n, a, &n, &info);

  if (info != 0) 
    cerr << " cholesky_invert(" << a << "," << n
	 << ") : Warning: inversion failure . info =" 
	 << info <<  " (>0 is not pos.def)" << endl;

  return info;
}

int posdef_invert(Mat &A, const char* UorL)
{
  Vect B(A.SizeX());
  int info = cholesky_solve(A, B, UorL);
  if (info == 0)
    info = cholesky_invert(A, UorL);
  return info;
}

#define DO10(A) A;A;A;A;A;A;A;A;A;A;

static double fast_scal_prod(double *x, double *y, const unsigned size)
{
  int nblock = int(size)/10;
  int remainder = size - nblock*10;
  double sum = 0;
  for (int i=0; i<nblock; ++i) 
    {
      DO10( sum+= (*x)*(*y); ++x; ++y;)
    }
  for (int i=0; i<remainder; ++i) {sum+= (*x)*(*y); ++x; ++y;}
  return sum;
}

int cholesky_solve_with_constraints(Mat &A, Vect &B, 
				    Mat &Constraints, Vect& ConstraintsValues,
				    const char* UorL)
{
#ifdef FNAME
  cout << " >  cholesky_solve_with_constraints" << endl;
#endif  

  if(A.SizeX() != A.SizeY() || A.SizeX()<=0) {
    cout << "error in matvect.cc, cholesky_solve_with_constraints Matrix A must be symmetric and you have nx,ny = "
	 << A.SizeX() << "," << A.SizeY() << endl;
    abort();
  }
  if(B.Size() != A.SizeX()) {
    cout << "error in matvect.cc, cholesky_solve_with_constraints " << endl 
	 << " Vector B must have a dimension B.Size()=A.SizeY() and you have B.Size(),A.SizeY() = "
	 << B.Size() << "," << A.SizeY() << endl;
    abort();
  }

  if (Constraints.SizeX() != A.SizeX())
    {
      cout << "error in matvect.cc, cholesky_solve_with_constraints " << endl
	 << " Constraints.SizeY() != A.Size  : "
	 << Constraints.SizeY() << ' ' << A.SizeX() << endl;
    abort();
    }

  if (Constraints.SizeY() != ConstraintsValues.Size())
    {
      cout << "error in matvect.cc, cholesky_solve_with_constraints " << endl
	   << " Constraints.SizeY() != ConstraintsValues.Size() : "
	   << Constraints.SizeY() << ' ' << ConstraintsValues.Size() << endl;
      abort();
    }

  double *a = A.NonConstData();
  double *b = B.NonConstData();
  
  int n = B.Size();
  int nrhs = 1, info = 0;
  char uorl[6]; strncpy(uorl,UorL,6);


  dposv_(uorl, &n, &nrhs, a, &n, b, &n, &info);

  if (info != 0) 
    {
      cerr << " cholesky_solve_with_constraints(" << a << "," << b << "," << n
	   << ") : Warning: Cholesky factorization failure . info =" 
	   << info <<  " (>0 is not pos.def)" << endl;
      
      return info;
    }
  

  size_t nconst = Constraints.SizeY();
  //need to keep a copy of the constraints
  Mat am1c(Constraints);
  double *c = am1c.NonConstData();
  nrhs = nconst;
  dpotrs_(uorl, &n, &nrhs, a, &n, c, &n, &info);

  Mat smalla(nconst,nconst);
  Vect smallb(nconst);
  for (size_t j=0; j <nconst; ++j)
    {
      for (size_t i = 0; i<=j ; ++i)
	{
	  smalla(i,j) = fast_scal_prod(&Constraints(0,i), &am1c(0,j), n);
	}
      smallb(j) = fast_scal_prod(&Constraints(0,j),&B(0), n) - ConstraintsValues(j);
    }

  // solve the small system
  if (cholesky_solve(smalla, smallb,"U") != 0)
    {
      cout << " cholesky_solve_with_constraints small sub-system could not be solved " << endl;
      return info;
    }

  // small b now contains the values of Lagrange multiplier

  B -=  am1c.transposed()*smallb;

  cout << " cholesky solve with constraints : residuals " << Constraints*B << endl;
  return 0;

}


// calls dsytri after the use of dsysv
int general_solve(Mat& A, Vect& B, bool invert_A, const char* UorL)
{
  int nrhs=1;
  int n = B.Size();
  int lwork=n;
  int info;
    
  double* a = A.NonConstData();
  double* b = B.NonConstData();
  int*    ipiv = new int[n];
  double* work = new double[n];
  
  dsysv_(UorL, &n, &nrhs, 
	 a, &n, ipiv, b, &n, work, &lwork, &info);
  if(info!=0) {
    std::cout << " Pb solving sys: info="
	      << info << std::endl;
    delete[] ipiv;
    delete[] work;
    return info;
  }
  
  if(! invert_A) {
    delete[] ipiv;
    delete[] work;
    return info;
  }

  dsytri_(UorL, &n, a , &n, ipiv, work, &info);

  if(info!=0) {
    std::cout << " Pb inverting the matrix: info="
	      << info << std::endl;
  }
  
  delete[] ipiv;
  delete[] work;
  return info; 
}

// calls dgesv
int truly_general_solve(Mat& A, Vect& B)
{
  Mat  Af=A;
  Vect X=B;
  


  int m      = X.Size();
  int n      = m;
  double* af = Af.NonConstData();
  int  lda   = n;
  int* ipiv  = new int[n];
  int info;

  dgetrf_(&m,&n,af,&lda,ipiv,&info);
  if(info!=0) {
    std::cout << " Pb with dgetrf_ info=" << info << std::endl;
    delete [] ipiv;
    return info;
  }
  
  int nrhs  = 1;
  double* x = X.NonConstData();
  int ldb   = n;
  char NTC = 'T';
  
  dgetrs_(&NTC,&n,&nrhs,af,&lda,ipiv,x,&ldb,&info);
  
  if(info!=0) {
    std::cout << " Pb with dgetrs_: info="
	      << info << std::endl;
    delete [] ipiv;
    return info;
  }
  
  double* a = A.NonConstData();
  double* b = B.NonConstData();
  
  double ferr,berr;
  double *work = new double[3*n];
  int *iwork = new int[n];
  int ldx  = n;
  int ldaf = n;
  dgerfs_(&NTC,&n,&nrhs,a,&lda,af,&ldaf,ipiv,b,&ldb,x,&ldx,&ferr,&berr,work,iwork,&info);
  
   
  if(info!=0) {
    std::cout << " Pb with dgerfs_: info="
	      << info << std::endl;
  }
  
  B=X;
  delete[] work;
  delete[] iwork;
  delete[] ipiv;
  return info;
}

/* calls dgesv : does not work
   int truly_general_solve(Mat& A, Vect& B)
   {
   
   int n     = B.Size();
   int nrhs  = 1;
   double* a = A.NonConstData();
   int  lda  = n;
   int* ipiv = new int[n];
   double* b = B.NonConstData();
   int ldb   = n;
   int info;
   
   dgesv_(&n, &nrhs, 
   a, &lda, ipiv, b, &ldb, &info);
   
   if(info!=0)
   std::cout << " Pb solving sys: info="
   << info << std::endl;
   
   delete[] ipiv;
   return info;
   }
*/



/* From the next routine, the eigenvector with lowest eigenvalue is
   eigenvals(0,*) */

int diagonalize(Mat const& A, 
		Vect& eigenvals, Mat& eigenvects, 
		const char* UorL)
{
  //cout << "[diagonalize]" << endl;
    
    Mat A_TMP(A);
    double* a  = A_TMP.NonConstData();
    char jobz  = 'V';
    char uplow = 'U';
    if (strlen(UorL) >=1) uplow = UorL[0];
    else
      cerr << "diagonalize() : UorL argument should be at least 1 char, assuming \'U\'" << endl;
    int  n     = A.SizeX();
    int info = 0;
    
    eigenvals.allocate(n);
    eigenvects.allocate(n,n);
    /* 
       There are several routines in lapack to diagonalize real
       symetric matrices. dsyev is said on the web to be slower than
       dsyevd. This is perhaps true at the 20 to 30 % level.  For an
       unknown reason both are *way* slower than the cernlib
       eis... package (translated to double precision). Maybe
       precision issues?
    */
#if 0 /* use dsyev (lapack) */
    // the function needs some workspace, this first call is meant to ask it "how much?"
    int lwork = -1;
    double how_much;
    dsyev_(&jobz, &uplow, &n, a, &n, eigenvals.NonConstData(),
	   &how_much, &lwork, &info);
    
    lwork = int(how_much);

    // reserve the required workspace, and call the function
    Vect  work(lwork);
    double* workspace = work.NonConstData();
    dsyev_(&jobz, &uplow, &n, a, &n, eigenvals.NonConstData(),
	   workspace, &lwork, &info);
#else /* use dsyevd (lapack) */
    // the function needs some workspace, this first call is meant to ask it "how much?"
    int lwork = -1;
    double how_much_d;
    int liwork = -1;
    int how_much_i;
    dsyevd_(&jobz, &uplow, &n, a, &n, eigenvals.NonConstData(),
	   &how_much_d, &lwork, &how_much_i, &liwork, &info);
    
    lwork = int(how_much_d);
    liwork = how_much_i;

    // reserve the required workspaces, and call the function
    Vect  work(lwork);
    double* workspace = work.NonConstData();
    int * iwork = new int [liwork];
    dsyevd_(&jobz, &uplow, &n, a, &n, eigenvals.NonConstData(),
	   workspace, &lwork, iwork, &liwork, &info);
    delete [] iwork;
#endif


    if(info != 0) {
	cerr << "[diagonalize] ERROR, failed to diagonalize matrix, info="
	     << info << endl;
	return info;
    }
    else {
      //cout << "[diagonalize] OK, info=" << info << endl;
    }
    
    
    // now, we want to update the eigenstuff
    //    memcpy(eigenvects.NonConstData(), 
    //	   A_TMP.NonConstData(),
    //	   n*n*sizeof(double));
    
    int i,j;
    for(i=0;i<n;i++)
    	for(j=0;j<n;j++)
   	    eigenvects(i,j) = A_TMP(j,i);
    
    return info;
}




#endif




//==================================================================




Mat::Mat(const size_t NX, const size_t NY) { 
 data=NULL;
 nx=ny=0;
 allocate(NX,NY);
}

Mat::Mat(const Mat& other) {
  data=NULL;
  nx=ny=0;
  allocate(other.nx,other.ny);
  memcpy(data,other.data,sizeof(double)*nx*ny); 
}

void Mat::allocate(const size_t NX, const size_t NY) {
  if(nx!=NX || ny!=NY) {
    nx=NX; 
    ny=NY; 
    if (data) 
      delete [] data;
    data = new double[nx*ny];
  } 
  Zero();
}





int Mat::writeASCII(const std::string &FileName) const
{
  std::ofstream S(FileName.c_str());
  if (!S)
    {
      cout << " Mat::writeASCII() : cannot open " << FileName << endl;
      return 0;
    }
  int status = writeASCII(S);
  S.close();
  return status;
}

int Mat::writeASCII(ostream& Stream) const {
  Stream << nx << " " << ny << endl;
  int oldprec = Stream.precision();
  Stream << setprecision(10);
  for(size_t j=0;j<ny;j++) {
    //Stream << "0.." << nx-1 << "," << j << ": ";
    for(size_t i=0;i<nx;i++) {
      Stream << " " << (*this)(i,j);
    }
    Stream << endl;
  }
  Stream << setprecision(oldprec);
  return 0;
}


int Mat::readASCII(const std::string &FileName)
{
  std::ifstream S(FileName.c_str());
  if (!S)
    {
      cout << " Mat::readASCII() : cannot open " << FileName << endl;
      return -1;
    }
  int status = readASCII(S);
  S.close();
  return status;
}


int Mat::readASCII(istream& Stream) {
  size_t  fnx,fny;
  
  Stream >> fnx >> fny;
  allocate(fnx,fny);
  
  double val;
  for(size_t j=0;j<ny;j++) {
    for(size_t i=0;i<nx;i++) {
      if(Stream >> val)
	(*this)(i,j)=val;
      else{
	cerr << "error reading MATRIX" << endl;
	abort();
      }
    }
  }
  return 1;
}


void Mat::Identity() {
  if(nx!=ny) {
    cout << "Mat::Identity ERROR nx!=ny" <<endl;
    abort();
  }
  Zero();
  for(size_t i=0;i<nx;++i)
    (*this)(i,i)=1.;
}

double Mat::sum() const
{
  double res =0;
  const double *a = data;
  for (size_t i=nx*ny; i; i--,a++) res += *a;
  return res;
}



static bool same_size(const Mat& m1, const Mat& m2)
{
  if ((m1.SizeX() == m2.SizeX()) && (m1.SizeY() == m2.SizeY())) return true;
  cout << " matrices have different sizes" << endl;
  abort(); // we should in fact throw an exception.
  return false;
}

static bool same_size(const Vect& v1, const Vect& v2)
{
  if (v1.Size() == v2.Size()) return true;
  cout << " vectors have different sizes" << endl;
  abort(); // we should in fact throw an exception.
  return false;
}

Mat Mat::operator +(const Mat& Right) const
{
  same_size(*this,Right);
  Mat res = (*this);
  res += Right;
  return res;
}

Mat Mat::operator -(const Mat& Right) const
{
  same_size(*this,Right);
  Mat res = (*this);
  res -= Right;
  return res;
}

void Mat::operator +=(const Mat& Right)
{
  same_size(*this,Right);
  double *a = data;
  const double *b = Right.data;
  size_t size = nx*ny;
  for(size_t i=0;i<size;++i, ++a, ++b)
    *a += *b;
}

void Mat::operator -=(const Mat& Right)
{
  same_size(*this,Right);
  double *a = data;
  const double *b = Right.data;
  size_t size = nx*ny;
  for(size_t i=0;i<size;++i, ++a, ++b)
    *a -= *b;
}

Mat Mat::operator *(const double Right) const 
{
  Mat res = (*this);
  res *= Right;
  return res;
}

Mat operator *(const double Left, const Mat &Right)
{
  Mat res = Right;
  res *= Left;
  return res;
}
 
void Mat::operator *=(const double Right)
{
  size_t size = nx*ny;
  double *a = data;
  for(size_t i=0;i<size;++i, ++a)
    *a *= Right;
}
/* OLD VERY SLOW
Mat Mat::operator *(const Mat& Right) const
{
  if(nx != Right.SizeY()) {
    cout << "Mat::operator *= ERROR nx != Right.SizeY()" << endl;
    abort();
  }
  Mat res(Right.nx,ny);
  
  for(size_t j=0;j<res.ny;++j) {
    for(size_t i=0;i<res.nx;++i) {  
      for(size_t k=0;k<nx;++k) {
	res(i,j) += (*this)(k,j)*Right(i,k);
      }
    }
  }
  return res;
}
*/


extern "C" {
void dgemm_(char &transa,
	    char &transb,
	    int &m, // number of rows of op(A)
	    int &n, // number of cols of op(B)
	    int &k, // number of columns of op(a) and numb of rows of op(B)
	    double &alpha,
	    const double *A, // const invented from the doc...
	    int &lda, // leading dim of A
	    const double *B, // const invented form the doc...
	    int &ldb,
	    double &beta,
	    double *C,
	    int &ldc);
	    
  // returns C = alpha A*B + beta C
  
  // note : if A(x,y)  in blas parlance :
  // number of rows = nx, number of columns = ny
  // which means that the "printout" reads (!= ours)
  // for (int i=0; i<nx; ++i)
  //   {
  //   for (int j=0;j<ny;++j)
  //       cout << A(i,j) << ' ';
  //   cout << endl;
  //   }

}


Mat Mat::operator *(const Mat& Right) const
{
  if(nx != Right.ny) {
    cout << "Mat::operator *= ERROR nx != Right.ny" << endl;
    abort();
  }
  Mat ResultMat(Right.nx,ny);

#define USE_DGEMM
#ifdef USE_DGEMM
  /* The definition chosen originally conflicts with usual definitions
      for(unsigned int k=0;k<nx;++k) {
	res(i,j) += (*this)(k,j)*Right(i,k);

	we will get the expected behaviour from DGEMM by swapping A and B
  */

  char transa='N';
  char transb='N';
  int m = Right.nx; // 
  int n = ny; //
  int k = nx; // = Right.ny 
  double alpha = 1;
  int lda = Right.nx;
  int ldb = nx;
  double beta = 0;
  int ldc = ResultMat.nx;
  
  if (Right.nx == 0) // we should throw an exception
    {
      std::cerr << " in matrix multiplication : Right.nx == 0 cannot work " << std::endl;
      abort ();
    }
	

  dgemm_(transa,transb,m,n,k,
	 alpha, Right.data,lda, 
	 data, ldb, 
	 beta, ResultMat.data, ldc);

#else
  double *res = ResultMat.data;
  for(size_t j=0;j<ResultMat.ny;++j) {
    for(size_t i=0;i<ResultMat.nx;++i,++res) {  
      double *left  = & data[j*nx];
      double *right = & Right.data[i];
      for(size_t k=0;k<nx;++k,++left,right+=Right.nx) {
	*res += (*left)*(*right);
      }
    }
  }
#endif /*USE_DGEMM  */
  return ResultMat;
}




Vect Mat::operator *(const Vect& Right) const
{
  if (Right.size() != nx)
    {
      cout << "Mat*Vect operator ERROR Mat.nx != Vect.size()" << endl;
      abort();
    }
     
  Vect res(ny);
  for (unsigned k =0; k<ny; ++k)
    {
      res(k) = fast_scal_prod(&(*this)(0,k), Right.data, nx);
    }
  return res;
}

void Mat::operator *=(const Mat& Right)
{
  Mat res = (*this)*Right;
  (*this) = res;
}


Mat & Mat::operator =(const Mat& Right){
  allocate(Right.nx,Right.ny);
  memcpy(data,Right.data,nx*ny*sizeof(double));
  return (*this);
}

Mat::operator double() const
{
  if(nx!=1 || ny !=1) {
    cout << "Mat::operator double() error, nx=ny=1 needed, you have nx=" 
	 << nx <<" ny=" << ny << endl;
    abort();
  }
  return (*this)(0,0);
}

Mat::operator Vect() const
{
  if(nx!=1) {
    cout << "Mat::operator Vect() error, nx=1 needed, you have nx=" 
	 << nx << endl;
    abort();
  }
  Vect res(ny);
  memcpy(res.data,data,ny*sizeof(double));
  return res;
}

Mat Mat::transposed() const {
  Mat res(ny,nx);
  
  for(size_t i=0;i<nx;i++) {
    for(size_t j=0;j<ny;j++) {
      res(j,i)=(*this)(i,j);
    }
  }
  return res;
}

Mat Mat::SubBlock
(size_t x_min,size_t x_max,size_t y_min,size_t y_max) const {
  if( x_max >= nx || y_max>= ny ) {
    cout << "Mat::SubBlockFromIndexes ERROR, trying to get a sub-matrix with indices" << endl;
    cout << "x_min,x_max,y_min,y_max = "
	 << x_min << ","
	 << x_max << ","
    	 << y_min << ","
	 << y_max << endl;
    cout << "nx,ny = "<< nx << "," << ny << endl;
    abort();
  }
  size_t nx_new = (x_max-x_min+1);
  size_t ny_new = (y_max-y_min+1);
  Mat res(nx_new,ny_new);
  for(size_t j=0;j<ny_new;++j)
    for(size_t i=0;i<nx_new;++i)
      res(i,j) = (*this)(i+x_min,j+y_min);
  return res;
}

Mat Mat::WithoutRows(size_t y_min,size_t y_max) const {
  if( y_max < y_min || y_max >= ny ) {
    cout << "Mat::WithoutRows ERROR y_min,y_max,ny = "
	 << y_min << ","
	 << y_max << ","
	 << ny
	 << endl;
    abort();
  } 
  size_t nrows = y_max-y_min+1;
  Mat res(nx,ny-nrows);
  for(size_t j = 0 ;j<y_min ;j++)
    for(size_t i=0;i<nx;i++)
      res(i,j)=(*this)(i,j);
  for(size_t j = y_max+1 ;j<ny ;j++)
    for(size_t i=0;i<nx;i++)
      res(i,j-nrows)=(*this)(i,j);
  return res;
}

Mat Mat::WithoutColumns(size_t x_min,size_t x_max) const {
  if( x_max < x_min || x_max >= nx ) {
    cout << "Mat::WithoutColumns ERROR " << endl;
    abort();
  } 
  size_t ncols = x_max-x_min+1;
  Mat res(nx-ncols,ny);
  for(size_t i=0;i<x_min;i++)
    for(size_t j = 0 ;j<ny ;j++)
      res(i,j)=(*this)(i,j);
  for(size_t i=x_max+1;i<nx;i++)
    for(size_t j = 0 ;j<ny ;j++)
      res(i-ncols,j)=(*this)(i,j);
  return res;
}




#ifdef USECFITSIO
int Mat::readFits(const string &FitsName, int hdu) {

  
  if( ! FileExists(FitsName) ) {
    cerr << "Mat::readFits : cannot open " << FitsName << endl;
    exit(-12);
  }
  
  


  int status = 0;
  fitsfile *fptr = 0;

#pragma omp critical // apparently cfitsio does not like multi thread
  {
    fits_open_file(&fptr,FitsName.c_str(),0,&status);
  }
  
  if (status)
    {
      cerr << " when opening file : " << FitsName ;
      cerr << " we got these successive errors:" << endl;
      
      
#pragma omp critical
      {
	fits_report_error(stderr, status);
      }
      
      return status;
    }

  // move hdu
  
#pragma omp critical // apparently cfitsio does not like multi thread
  {
    int exttype;
    fits_movabs_hdu(fptr,hdu,&exttype,&status);
  }
  
  if (status)
    {
      cerr << " when opening file : " << FitsName ;
      cerr << " we got these successive errors:" << endl;
      
      
#pragma omp critical
      {
	fits_report_error(stderr, status);
      }
      
      return status;
    }
  



  // first get the size of the image/matrix
  status=0;
  char value[256];
  char keyname[16] = "NAXIS1";
  
#pragma omp critical // apparently cfitsio does not like multi thread
  {
    fits_read_key(fptr, TSTRING, keyname, &value, NULL, &status);
  }
  
  if (status)
    {
      cerr << " when reading content of : " << FitsName ;
      cerr << " we got these successive errors:" << endl;
      
#pragma omp critical // apparently cfitsio does not like multi thread
      {
	fits_report_error(stderr, status);
      }
      
      return status;
    }
  int n1 = atoi(value);
  strcpy(keyname,"NAXIS2");
  
#pragma omp critical // apparently cfitsio does not like multi thread
  {
    fits_read_key(fptr, TSTRING, keyname, &value, NULL, &status);
  }
  
  int n2=1;
  
  if (status==0)
    {
      n2 = atoi(value);
    }
  
  if(n1<=0 || n2<=0) {
    cout << "Mat::readFits error NAXIS1,NAXIS2 = " << n1 << "," << n2 << endl;
    return -1;
  }
  allocate(n1,n2);
  
  status = 0;
  float nullval = 0;
  int anynull;
  
#pragma omp critical // apparently cfitsio does not like multi thread
  {
    fits_read_img(fptr, TDOUBLE, 1, nx*ny, &nullval,  data, &anynull, &status);
  }
  
  if (status)
    {
      cerr << " when reading content of : " << FitsName ;
      cerr << " we got these successive errors:" << endl;
#pragma omp critical
      {
	fits_report_error(stderr, status);
      }
      return status;
    }
  
  status = 0;
#pragma omp critical // apparently cfitsio does not like multi thread
  {
    
    fits_close_file(fptr, &status);
  }

  if (status)
    {
     cerr << " when closing file : " << FitsName ;
     cerr << " we got these successive errors:" << endl;
#pragma omp critical // apparently cfitsio does not like multi thread
     {
       fits_report_error(stderr, status);
     }
    }
  return status;  
}


int Mat::writeFits(const string &FitsName) const {
 

// we cannot use fitsimage cause we want to write it in double precision
  
  int status = 0;
  fitsfile *fptr = 0;
  
  remove(FitsName.c_str());
  
  // enum FitsFileMode {RO = 0, RW = 1};
#pragma omp critical // apparently cfitsio does not like multi thread
  {
    fits_create_file(&fptr,FitsName.c_str(),  &status);
  }
  if (status)
    {
      cerr << " when creating file : " << FitsName ;
      cerr << " we got these successive errors:" << endl;
#pragma omp critical // apparently cfitsio does not like multi thread
      {
	fits_report_error(stderr, status);
      }
      return status;
    }
  // set a minimal header
  status = 0;
  long naxes[2];
  naxes[0]=nx;
  naxes[1]=ny;
#pragma omp critical // apparently cfitsio does not like multi thread
  {
    fits_write_imghdr(fptr,-64,2,naxes,&status);
  }
  if (status)
    {
      cerr << " when writing minial header  ";
      cerr << " we got these successive errors:" << endl;
#pragma omp critical // apparently cfitsio does not like multi thread
      {
	fits_report_error(stderr, status);
      }
      return status;
    }
  
  // say to cfitsio to take into account this new BITPIX
  status = 0;
#pragma omp critical // apparently cfitsio does not like multi thread
  {
    fits_flush_file(fptr,&status);
  }
  if (status)
    {
      cerr << " when flushing  ";
      cerr << " we got these successive errors:" << endl;
#pragma omp critical // apparently cfitsio does not like multi thread
      {
	fits_report_error(stderr, status);
      }
      return status;
    }
  status = 0;
#pragma omp critical // apparently cfitsio does not like multi thread
  {
    fits_write_img(fptr, TDOUBLE, 1, nx*ny, data, &status);
  }
  if (status)
    {
      cerr << " when writing data ";
      cerr << " we got these successive errors:" << endl;
#pragma omp critical // apparently cfitsio does not like multi thread
      {
	fits_report_error(stderr, status);
      }
      return status;
    }
  status = 0;
#pragma omp critical // apparently cfitsio does not like multi thread
  {
    fits_close_file(fptr, &status);
  }
  if (status)
    {
      cerr << " when closing file : " << FitsName ;
      cerr << " we got these successive errors:" << endl;
#pragma omp critical // apparently cfitsio does not like multi thread
      {
	fits_report_error(stderr, status);
      }
    }
  return status;
}

void AddFitsKey(const std::string &FileName, 
		const string &KeyName, 
		const double &KeyVal)
{
  fitsfile *fptr = 0;
  int status=0;
  fits_open_file(&fptr,FileName.c_str(), 1, &status);
  char key_name[80];
  sprintf(key_name,"%s",KeyName.c_str());
  double value = KeyVal;
  fits_write_key(fptr, TDOUBLE, key_name, &value, NULL, &status);
  fits_close_file(fptr, &status);
  if (status)
    fits_report_error(stderr, status);
}

double  ReadFitsKey(const std::string &FileName, 
		    const string &KeyName)
{
  fitsfile *fptr = 0;
  int status=0;
  fits_open_file(&fptr,FileName.c_str(), 0, &status);
  double value;
  char key_name[80];
  sprintf(key_name,"%s",KeyName.c_str());
  fits_read_key(fptr, TDOUBLE, key_name, &value, NULL, &status);
  fits_close_file(fptr, &status);
  if (status)
    fits_report_error(stderr, status);
  return value;
}


#endif


void Mat::Symmetrize(const char* UorL) {
  
  if(nx!=ny) {
    cout << "Mat::Symmetrize ERROR nx!=ny nx,ny = " << nx << "," << ny << endl;
    abort();
  }
  
  
  for(size_t j=0;j<ny;j++)
    for(size_t i=j+1;i<nx;i++)
      if(UorL[0]=='L') { // x >= y
	(*this)(j,i)=(*this)(i,j);
      }else{
	(*this)(i,j)=(*this)(j,i);
      }
} 

#ifdef USECERNLIB
int Mat::SymMatInvert() {
  if(nx!=ny) {
    cout << "Mat::Symmetrize ERROR nx!=ny nx,ny = " << nx << "," << ny << endl;
    abort();
  }
  double *A = data;
  int n = nx;
   int ierr;
  dsinv(&n,A,&n,&ierr);
  if(ierr!=0) {
    cerr << "fatal error in Mat::SymMatInvert() ierr=" << ierr << endl;
    exit(-12);
  }
  return (!ierr);
}
#endif

int Mat::CholeskyInvert(const char *UorLorS, bool abort_if_failed) {
  if(nx!=ny || nx == 0) {
    cout << "Mat::CholeskyInvert ERROR (bad dimensions) nx,ny = " << nx << "," << ny << endl;
    abort();
  }

  char uorl[6];
  bool is_symmetric = (UorLorS[0]=='S');
  if(is_symmetric)
    uorl[0]='L'; // ou U, on s'en fout
  else
    strncpy(uorl,UorLorS,6);
  
  double *a = NonConstData();
  int n = SizeX();
  int info = 0;

  dpotrf_(uorl, &n, a,  &n, &info);
  if (info != 0) 
    {
      
      cerr << "fatal error in Mat::CholeskyInvert : could not factorize, info = " << info << endl;
      if(abort_if_failed) 
	abort();
      else 
	return info;
      //exit(info);
    }
  
  dpotri_(uorl, &n, a, &n, &info);
  if (info != 0) 
    {
      cerr << "Mat::CholeskyInvert : could not invert, info = " << info << endl;
      if(abort_if_failed) 
	abort();
      else 
	return info;
    }
  
  if(is_symmetric) {
    this->Symmetrize(uorl);
  }
  
  return info;
}


//=================================================================

Vect::Vect(const size_t N) {
  data = 0; 
  n=0;
  if(N<=0) {
    //cout << "Vect::Vect ERROR N = " << N <<  endl;
  }else{
    allocate(N);
  }
}

Vect::Vect(const Vect& other) {
  data = 0; 
  n=0;
  allocate(other.n);
  memcpy(data,other.data,sizeof(double)*n); 
}

void Vect::allocate(const size_t N) {
  if(n!=N) {
    n=N;
    if (data) {
	delete [] data;
    }
    data = new double[n];
  }
  Zero();
}




int Vect::writeASCII(ostream& Stream) const {
  Stream << n << endl;
  for(size_t i=0;i<n;i++) {
    Stream << (*this)(i) << endl;
  }
  return 0;
}


int Vect::readASCII(istream& Stream)  {
  size_t fn;
  Stream >> fn;
  allocate(fn);
  double val;
  for(size_t i=0;i<n;i++) {
    Stream >> val;
    (*this)(i) = val;
  }
  return 0;
}

Vect Vect::operator *(const double Right) const 
{
  Vect res = (*this);
  res *= Right;
  return res;
}

Vect operator *(const double Left, const Vect &Right)
{
  Vect res = Right;
  res *= Left;
  return res;
}
 
void Vect::operator *=(const double Right)
{
  double *a = data;
  for(size_t i=0;i<n;++i, ++a)
    *a *= Right;
}


Vect Vect::operator +(const Vect& Right) const
{
  same_size(*this,Right);
  Vect res = (*this);
  res += Right;
  return res;
}

Vect Vect::operator -(const Vect& Right) const
{
  same_size(*this,Right);
  Vect res = (*this);
  res -= Right;
  return res;
}

void Vect::operator +=(const Vect& Right)
{
  same_size(*this,Right);
  double *a = data;
  const double *b = Right.data;
  for(size_t i=0;i<n;++i, ++a, ++b)
    *a += *b;
}

void Vect::operator +=(double Right)
{
  double *a = data;
  for(size_t i=0;i<n;++i)
    a[i] += Right;
}

void Vect::operator -=(const Vect& Right)
{
  same_size(*this,Right);
  double *a = data;
  const double *b = Right.data;
  for(size_t i=0;i<n;++i, ++a, ++b)
    *a -= *b;
}

double Vect::operator *(const Vect& Right) const
{
  same_size(*this,Right);
  double res = 0;
  const double *a = data;
  const double *b = Right.data;
  for(size_t i=0;i<n;++i, ++a, ++b)
    res += (*a)*(*b);
  return res;
}

Vect & Vect::operator =(const Vect& Right){
  allocate(Right.Size());
  memcpy(data,Right.data,n*sizeof(double));
  return (*this);
}

void Vect::operator =(const double &Val){
  double *pend = data+n;
  for (double *p=data; p<pend; ++p) *p = Val; 
}

double Vect::sum() const
{
  double res =0;
  const double *a = data;
  for (size_t i=n; i; i--,a++) res += *a;
  return res;
}


Mat Vect::transposed() const {
  Mat trans(n,1);
  memcpy(trans.data,data,n*sizeof(double));
  return trans;
}

Vect::operator double() const
{
  if(n!=1) {
    cout << "Vect::operator double() error, n=1 needed, you have n=" 
	 << n << endl;
    abort();
  }
  return (*this)(0);
}


Vect::operator Mat() const {
  Mat mat(1,n);
  memcpy(mat.data,data,n*sizeof(double));
  return mat;
}

//================================================

SymmetricFloatMatContainer::SymmetricFloatMatContainer(const Mat& other,const char UoL)
{
  n=0;
  data = NULL;
  Load(other,UoL);
}

void SymmetricFloatMatContainer::Load(const Mat& other,const char UoL){
  if(other.SizeX()!=other.SizeY()) {
    cerr << "SymmetricFloatMatContainer::SymmetricFloatMatContainer fatal error : mat is no symmetric" << endl;
    abort();
  }
  if(data)
    delete [] data;

  n = other.SizeX();
  data = new MFLOAT[n*(n+1)/2];
  
  switch(UoL) {
  case 'L':
    {
      MFLOAT *data_pointer = data;
      for(size_t i=0;i<n;i++)
	for(size_t j=0;j<=i;++j,++data_pointer)
	  *data_pointer = other(i,j);
    }
    break;
  case 'U':
    {
      MFLOAT *data_pointer = data;
      for(size_t i=0;i<n;i++)
	for(size_t j=0;j<=i;++j,++data_pointer)
	  *data_pointer = other(j,i);
    }
    break;
  default :
    cerr << "SymmetricFloatMatContainer::SymmetricFloatMatContainer fatal error : U or L and not " << UoL << endl;
    abort();
  }
  
}


SymmetricFloatMatContainer::operator Mat() const {
  
  Mat A(n,n);
  
  const MFLOAT *data_pointer = data;
  for(size_t i=0;i<n;i++)
    for(size_t j=0;j<=i;++j,++data_pointer) {
      A(i,j) = *data_pointer;
      A(j,i) = *data_pointer;
    }
  
  return A;
}

SymmetricFloatMatContainer::~SymmetricFloatMatContainer() {

  if(data) 
    delete [] data;
}
  
SymmetricFloatMatContainer::SymmetricFloatMatContainer(const SymmetricFloatMatContainer& other) {
  n = other.n;
  data = 0;
  if(n>0) {
    size_t size = n*(n+1)/2;
    data = new MFLOAT[size];
    memcpy(data,other.data,sizeof(MFLOAT)*size); 
  }
}



/*
==================================================================
 returns A^t B C, usefull when there are lots of zeros in A
 - assumes B is symmetric
 - if B.SizeX()=1, assumes this is a diagonal matrix
==================================================================
*/
double scalar_product(const Vect& A, const Mat& B, const Vect& C) {
  
  if(B.SizeX()==1) {
    
    const double* a = A.Data();
    const double* b = B.Data();
    const double* c = C.Data();
    size_t n = A.Size();
    size_t nblock = n/10;
    size_t remainder = n - nblock*10;
    double res=0;
    
    for (size_t i=0; i<nblock; ++i) {
      DO10( res += (*a)*(*b)*(*c); ++a; ++b; ++c;)
	}
    
    for(size_t i=0;i<remainder;++i,++a,++b,++c)
      res += (*a)*(*b)*(*c);
    
    return res;
  }
    
  double res=0;
  size_t size = A.Size(); // size is not very large here cause light curve points
  const double* a = A.Data();
  const double* b = B.Data();
  const double* c0 = C.Data();
  const double* c;
  double bc;
  size_t i;
  
  for(size_t j=0;j<size;j++,a++) {
    if(*a) {
      c = c0;
      bc = 0;
      for(i=0;i<size;++i,++b,++c)
	bc += (*b)*(*c);
      res += (*a)*bc; // r += A(j) sum_i B(i,j) C(i)
    }else{
      b += size;
    }
  }
  return res; // = sum_j A(j) sum_i B(i,j) C(i)
}


double scalar_product(const Vect& V1,const Vect& V2) {
  
  const double* v1 = V1.Data();
  const double* v2 = V2.Data();  
  size_t n = V1.Size();
  size_t nblock = n/10;
  size_t remainder = n - nblock*10;
  double res=0;

  for (size_t i=0; i<nblock; ++i) {
    DO10( res+= (*v1)*(*v2); ++v1; ++v2;)
  }

  for(size_t i=0;i<remainder;++i,++v1,++v2)
    res += (*v1)*(*v2);

  return res;
}


double norm2(const Vect& V) {
  
  const double* v = V.Data();
  size_t n = V.Size();
  size_t nblock = n/10;
  size_t remainder = n - nblock*10;
  double res=0;

  for (size_t i=0; i<nblock; ++i) {
    DO10( res+= (*v)*(*v); ++v;)
  }

  for(size_t i=0;i<remainder;++i,++v)
    res += (*v)*(*v);

  return res;
}

/*
==================================================================
Fast routine to save CPU
==================================================================
*/
Vect vect_product(const Vect& V1,const Vect& V2) {
  
  Vect Res = V1;
  double* res = Res.NonConstData();
  const double* v2 = V2.Data();  
  size_t n = Res.Size();
  size_t nblock = n/10;
  size_t remainder = n - nblock*10;

  for (size_t i=0; i<nblock; ++i) {
    DO10( *res *= (*v2); ++res; ++v2;)
  }

  for(size_t i=0;i<remainder;++i,++res,++v2)
    *res *= (*v2);

  return Res;
}
/*
==================================================================
Fast routine to save CPU
==================================================================
*/
Vect vect_product(const Vect& V1,const Vect& V2,const Vect& V3) {
  
  Vect Res = V1;
  double* res = Res.NonConstData();
  const double* v2 = V2.Data();
  const double* v3 = V3.Data();   
  size_t n = Res.Size();
  size_t nblock = n/10;
  size_t remainder = n - nblock*10;

  for (size_t i=0; i<nblock; ++i) {
    DO10( *res *= (*v2)*(*v3); ++res; ++v2; ++v3;)
  }

  for(size_t i=0;i<remainder;++i,++res,++v2,++v3)
    *res *= (*v2)*(*v3);

  return Res;
}

/*
==================================================================
Fast routine to compute largest eignevalue of a matrix
==================================================================
*/
double largest_eigenvalue(const Mat& M, Vect& eigenvector) {
  
  size_t s = M.SizeX();
  // init
  eigenvector.allocate(s);
  double* val = eigenvector.NonConstData();
  for(size_t i=0;i<s;i++,val++) 
    *val = 1;
  Vect eigenvector_tmp;
  
  double s1,s2;
  int count = 0;
  while(count <30) {
    eigenvector_tmp = M*eigenvector;
    s1 = scalar_product(eigenvector_tmp,eigenvector);
    s2 = norm2(eigenvector_tmp);
    // if the same, sqr(s1)=s2;
    if( fabs(s1*s1/s2-1.)<1.e-8 ) 
      break;
    eigenvector = (1./sqrt(s2))*eigenvector_tmp;
    count ++;
  }
  cout << "convergence after " << count << " iterations" << endl;
  eigenvector_tmp = M*eigenvector;
  return sqrt(norm2(eigenvector_tmp));
}



/*
  DGEMM  performs one of the matrix-matrix operations
  C = alpha*B*A+beta*C,
*/
void DGEMM(Mat &C, const Mat& A, const Mat& B, const double& i_alpha, const double& i_beta) {

  if(A.SizeX() != C.SizeX()) {cout << "ERROR A.SizeX() != C.SizeX()" << endl; abort();}
  if(B.SizeY() != C.SizeY()) {cout << "ERROR B.SizeY() != C.SizeY()" << endl; abort();}
  if(A.SizeY() != B.SizeX()) {cout << "ERROR A.SizeY() != B.SizeX()" << endl; abort();}
  
  char transa='N';
  char transb='N';
  int m = A.SizeX(); // M specifies  the number  of rows  of the  matrix op( A )  and of the  matrix  C
  int n = B.SizeY(); // N  specifies the number  of columns of the matrix op( B ) and the number of columns of the matrix C. N must be at least zero.
  int k = A.SizeY(); // K  specifies  the number of columns of the matrix op( A ) and the number of rows of the matrix op( B )
  int lda = m;
  int ldb = k;
  int ldc = m;
  double alpha = i_alpha;
  double beta = i_beta;
  
  /* SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */
  dgemm_(transa,transb,m,n,k,
	 alpha, A.Data(),lda, 
	 B.Data(), ldb, 
	 beta, C.NonConstData(), ldc);
}

/*
  DGEMM  performs one of the matrix-matrix operations
  C = alpha*B*A+beta*C,
*/
void DGEMM(Vect &C, const Vect& A, const Mat& B, const double& i_alpha, const double& i_beta) {

  if(B.SizeY() != C.Size()) {cout << "ERROR B.SizeY() != C.SizeY()" << endl; abort();}
  if(A.Size() != B.SizeX()) {cout << "ERROR A.SizeY() != B.SizeX()" << endl; abort();}
  
  char transa='N';
  char transb='N';
  int m = 1; // M specifies  the number  of rows  of the  matrix op( A )  and of the  matrix  C
  int n = B.SizeY(); // N  specifies the number  of columns of the matrix op( B ) and the number of columns of the matrix C. N must be at least zero.
  int k = A.Size(); // K  specifies  the number of columns of the matrix op( A ) and the number of rows of the matrix op( B )
  int lda = m;
  int ldb = k;
  int ldc = m;
  double alpha = i_alpha;
  double beta = i_beta;
  
  /* SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */
  dgemm_(transa,transb,m,n,k,
	 alpha, A.Data(),lda, 
	 B.Data(), ldb, 
	 beta, C.NonConstData(), ldc);
}

