/************************************************************/
/*                Utility functions                         */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/

/**
 * @file utilities.c
 * @brief Utility functions for simple arithmetic calculation,
 * matrix computation and numerical derivatives. 
 * @author Wayne Zhang                         
 */

#include "common.h"
#include "utilities.h"
/** constant used in numerical differentiation */
#define DIFF_EPS 0.001


/************************************************************/
/*               Simple arithmetic utility                  */
/************************************************************/

/**
 * Compute sample variance 
 *
 * @param x samples in long vector 
 * @param n number of samples
 *
 * @return sample variance
 */
double var(double *x, int n){
  double ans = 0.0, m = mean(x, n) ;
  for (int i = 0; i < n; i++)
    ans += (x[i] - m) * (x[i] - m);
  ans /= n - 1.0  ;
  return ans;
}

/**
 * Compute sample covariance matrix 
 *
 * @param n number of samples
 * @param p number of variables (columns), p >2
 * @param x samples in long vector 
 * @param ans vector to store computed covariance matrix
 *
 */
void cov(int n, int p, double *x, double *ans){
  double *one = Calloc(n * n, double),
    *x2 = Calloc(n * p, double),
    *x3 = Calloc(n * p, double);
  double alpha = -1.0 / n, beta = 1.0, beta2 = 0.0;

  // subtract mean    
  for (int i = 0; i < n * n; i++) one[i] = 1.0 ;
  Memcpy(x2, x, n * p) ;
  Memcpy(x3, x, n * p); 
  F77_CALL(dgemm)("N", "N", &n, &p, &n, &alpha, one,
		  &n, x2, &n, &beta, x3, &n);
  Memcpy(x2, x3, n * p) ;
  AZERO(ans, p * p) ;
    
  // compute covariance 
  F77_CALL(dgemm)("T", "N", &p, &p, &n, &beta, x2,
		  &n, x3, &n, &beta2, ans, &p);
  for (int i = 0; i < p * p; i++) ans[i] /= (n - 1) * 1.0 ;
  Free(one) ;
  Free(x2) ;
  Free(x3) ;
}

/************************************************************/
/*                 Matrix computations                      */
/************************************************************/


/**
 * Multiply a matrix and a vector 
 *
 * @param trans transpose of matrix?
 * @param m row count of matrix
 * @param n column count of matrix
 * @param A input matrix
 * @param x input vector
 * @param out output vector 
 *
 */
void mult_mv(char *trans, int m, int n, double *A,
             double *x, double *out){
  double one = 1.0, zero = 0.0 ;
  int incx = 1;
  F77_CALL(dgemv)(trans, &m, &n, &one, A, &m, x, &incx,
		  &zero, out, &incx) ;
}


/**
 * compute t(x) * x
 *
 * @param m row dimension of the matrix
 * @param n column dimension of the matrix
 * @param x the input matrix  
 * @param out output results
 *
 */

void mult_xtx(int m, int n, double *x, double *out){
  double alpha = 1.0, beta = 0.0, *x2 = Calloc(m * n, double);
  Memcpy(x2, x, m * n) ;
  F77_CALL(dgemm)("T", "N", &n, &n, &m, &alpha, x2, &m,
		  x, &m, &beta, out, &n) ;
  Free(x2) ;
}

/**
 * compute the lower cholesky factor
 *
 * @param d dimension of the matrix
 * @param v input matrix
 * @param iv output cholesky factor
 *
 */
void chol(int d, double *v, double *iv){    
  int info = 0;
  // cholesky factor of v
  Memcpy(iv, v, d * d) ;   
  F77_CALL(dpotrf)("L", &d, iv, &d, &info) ;
  if (info)  error(_("Error %d in Cholesky decomposition."), info) ;   
}

/**
 * invert a positive symmetric matrix 
 *
 * @param d dimension of the matrix
 * @param v input matrix
 * @param iv output inverse of the matrix
 *
 */
void solve_po(int d, double *v, double *iv){    
  int info = 0;
  // cholesky factor of v
  chol(d, v, iv) ;
  // compute inverse    
  F77_CALL(dpotri)("L", &d, iv, &d, &info) ;    
  if (info) error(_("Error %d in inverting matrix."), info) ;
  // fill upper triangle 
  for (int i = 0; i < d - 1; i++){
    for (int j = i + 1; j < d; j++)
      iv[j * d + i] = iv[i * d + j] ;
  }    
}


/************************************************************/
/*                Compute numerical derivatives             */
/************************************************************/


/**
 * Compute numerical gradient 
 *
 * @param n length of parmaters
 * @param x values at which to evaluate the gradient
 * @param myfunc user specified function 
 * @param data struct used in myfunc
 * @param ans vector to store the gradient 
 *
 */

void grad(int n, double *x,  double (*myfunc)(double *x, void *data), 
          void *data, double *ans){
  double y1, y2 ;
  for (int i = 0; i < n; i++){
    x[i] += DIFF_EPS ;
    y1 = myfunc(x, data) ;
    x[i] -= 2 * DIFF_EPS ;        
    y2 = myfunc(x, data) ;
    ans[i] = (y1 - y2) / DIFF_EPS * 0.5 ;
    x[i] += DIFF_EPS ;
  }
}


/**
 * Compute numerical hessian matrix  
 *
 * @param n length of parmaters
 * @param x values at which to evaluate the hessian
 * @param myfunc user specified function 
 * @param data struct used in myfunc
 * @param ans n*n vector to store the hessian matrix 
 *
 */
void hess(int n, double *x, double (*myfunc)(double *x, void *data), 
          void *data, double *ans){
  double *y1 = Calloc(n, double),
    *y2 = Calloc(n, double)  ;
  for (int i = 0; i < n; i++){
    x[i] += DIFF_EPS ;
    grad(n, x, myfunc, data, y1) ;
    x[i] -= 2 * DIFF_EPS ;
    grad(n, x, myfunc, data, y2) ;
    for (int j = 0; j < n; j++)
      ans[j + i * n] = (y1[j] - y2[j]) / DIFF_EPS * 0.5 ;
    x[i] += DIFF_EPS ;
  }
  Free(y1) ; Free(y2) ;
}


/************************************************************/
/*                Simulate a Wishart variable               */
/************************************************************/

/**
 * Simulate the Cholesky factor of a standardized Wishart variate with
 * dimension p and nu degrees of freedom.
 *
 * @param nu degrees of freedom
 * @param p dimension of the Wishart distribution
 * @param upper if 0 the result is lower triangular, otherwise upper
                triangular
 * @param ans array of size p * p to hold the result
 *
 * @return ans
 */
static double *std_rWishart_factor(double nu, int p, int upper, double ans[])
{
    int pp1 = p + 1;

    if (nu < (double) p || p <= 0)
      error(_("inconsistent degrees of freedom and dimension"));

    AZERO(ans, p * p);
    for (int j = 0; j < p; j++) {	/* jth column */
	ans[j * pp1] = sqrt(rchisq(nu - (double) j));
	for (int i = 0; i < j; i++) {
	    int uind = i + j * p, /* upper triangle index */
		lind = j + i * p; /* lower triangle index */
	    ans[(upper ? uind : lind)] = norm_rand();
	    ans[(upper ? lind : uind)] = 0;
	}
    }
    return ans;
}

/**
 * Simulate a sample of random matrix from a Wishart distribution
 *
 * @param d row (=column) dimension of the matrix 
 * @param nu Degrees of freedom
 * @param scal Positive-definite scale matrix
 * @param out simulated matrix (d*d)
 *
 */
void rwishart(int d, double nu, double *scal, double *out)
{
    int  info,  psqr;
    double *scCp, *tmp, one = 1, zero = 0;

    psqr = d*d;
    tmp = Calloc(psqr, double);
    scCp = Calloc(psqr, double);

    Memcpy(scCp, scal, psqr);
    AZERO(tmp, psqr);
    F77_CALL(dpotrf)("U", &d, scCp, &d, &info);
    if (info)
	error(_("scale matrix is not positive-definite"));
    GetRNGstate();    
    std_rWishart_factor(nu, d, 1, tmp);
    F77_CALL(dtrmm)("R", "U", "N", "N", &d, &d,
			&one, scCp, &d, tmp, &d);
    F77_CALL(dsyrk)("U", "T", &d, &d, &one, tmp, &d,
			&zero, out, &d);
    for (int i = 1; i < d; i++){
        for (int k = 0; k < i; k++)
            out[i + k * d] = out[k + i * d];
    }
    PutRNGstate();
    Free(tmp) ;
    Free(scCp) ;
}
