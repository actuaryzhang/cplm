/**
 * @file mh.c
 * @brief Functions for the random work Metropolis-Hastings algorithm
 * @author Wayne Zhang                        
 *
 */

#include <R.h>            /* for R related */
#include <Rmath.h>        /* for random number simulators */ 
#include <R_ext/BLAS.h>   /* for dtrmv etc */
#include "utilities.h"    /* for chol, mult_mv related */
#include <Rinternals.h>   /* for definition of R callables */

/********************************************************************/
/*       One dimensional truncated random walk M-H                  */
/********************************************************************/

/**
 * Simulate a truncated Normal random variable 
 *
 * @param m mean of the untruncated normal
 * @param sd standard deviation of the untruncated normal
 * @param lb left bound of the truncated normal
 * @param rb right bound of the truncated normal
 *
 * @return one simulated truncated normal 
 */
static R_INLINE double rtnorm(double m, double sd, double lb, double rb){
  double u = runif(R_FINITE(lb) ? pnorm(lb, m, sd, 1, 0) : 0.0,
		   R_FINITE(rb) ? pnorm(rb, m, sd, 1, 0) : 1.0);
  return qnorm(u, m, sd, 1, 0);
}


/**
 * compute the log density of a truncated normal
 *
 * @param x the point at which the log density is computed
 * @param m the mean of the untruncated normal
 * @param sd the standard deviation of the untruncated normal
 * @param lb the left bound of the truncated normal
 * @param rb the right bound of the truncated normal
 *
 * @return the log density at the point x
 */
static R_INLINE double dtnorm(double x, double m, double sd, double lb, double rb){
  double c =  (R_FINITE(rb) ? pnorm(rb, m, sd, 1, 0) : 1.0) -
    (R_FINITE(lb) ? pnorm(lb, m, sd, 1, 0) : 0.0) ; 
  return dnorm(x, m, sd, 1) - log(c) ;
}


/**
 * RW Metropolis update using the truncated Normal 
 *
 * @param m the mean of the untruncated normal
 * @param sd the standard deviation of the untruncated normal
 * @param lb the left bound of the truncated normal
 * @param rb the right bound of the truncated normal
 * @param sn pointer to store simulated value
 * @param myfunc user specified function to compute the log posterior 
 * @param data the struct used in myfunc
 *
 * @return  a 0-1 integer: 0 means not accepted and 1 accepted
 */
int metrop_tnorm_rw(double m, double sd, double lb, double rb, double *sn, 
		     double (*myfunc)(double x, void *data), void *data){
  *sn = rtnorm(m, sd, lb, rb);
  double C = (lb == R_NegInf && rb == R_PosInf) ? 0.0 :  
    (dtnorm(m, *sn, sd, lb, rb) - dtnorm(*sn, m, sd, lb, rb));
  /* determine whether to accept the sample */
  double A = exp(myfunc(*sn, data) - myfunc(m, data) + C) ;  
  if (A < 1 && runif(0, 1) >= A){ 
    *sn = m ;
    return 0 ;
  }
  else return 1 ;  
}


/********************************************************************/
/*       Multivariate (untruncated) random walk M-H                 */
/********************************************************************/

/**
 * simulation of a multivariate normal
 *
 * @param d the dimension
 * @param m the mean vector
 * @param v the positive-definite covariance matrix
 * @param s the vector to store the simulation
 *
 */
static void rmvnorm(int d, double *m, double *v, double *s){
  int incx = 1;
  double *lv = Calloc(d * d, double) ;
  /* simulate d univariate normal r.v.s */
  for (int i = 0; i < d; i++)  s[i] = rnorm(0, 1) ;    
  /* cholesky factor of v */
  chol(d, v, lv) ;
  /* scale and shift univariate normal r.v.s */
  F77_CALL(dtrmv)("L", "N", "N", &d, lv, &d, s, &incx) ;
  for (int i = 0; i < d; i++)  s[i] += m[i] ;    
  Free(lv) ;    
}

/**
 * return the exponent of a multivariate normal
 *
 * @param d dimension of the matrix
 * @param x sample vector 
 * @param m mean vector
 * @param iv inverse of the covariance matrix 
 *
 * @return exponent of MVN
 */
double dmvnorm(int d, double *x, double *m, double *iv){
  double ep = 0.0, *dx = Calloc(d, double), *tmp = Calloc(d, double) ;
  for (int i = 0; i < d; i++)
    dx[i] = m ? (x[i] - m[i]) : x[i];
  mult_mv("N", d, d, iv, dx, tmp) ;
  for (int i = 0; i < d; i++)
    ep += dx[i] * tmp[i] ;
  ep *= -0.5 ;
  Free(dx); Free(tmp);
  return ep ;
}

/**
 * random walk metropolis sampling for a vector of parameters 
 * of length d using multivariate normal proposal
 *
 * @param d the dimension of the parameter
 * @param m current values of the parameter (also mean vector
 *      in the multivariate Normal)
 * @param v covariance matrix in the proposal distribution
 * @param sn simulated new vector 
 * @param myfunc user specified function to compute the log posterior 
 * @param data the struct used in myfunc
 *
 * @return  a 0-1 integer: 0 means not accepted and 1 accepted
 *
 */
int metrop_mvnorm_rw(int d, double *m, double *v, double *sn, 
		     double (*myfunc)(double *x, void *data), void *data){
  rmvnorm(d, m, v, sn) ;
  /* determine whether to accept the sample */
  double A = exp(myfunc(sn, data) - myfunc(m, data) ) ;
  if (A < 1 && runif(0, 1) >= A){
    Memcpy(sn, m, d) ;
    return 0 ;
  }
  else return 1 ;
  
}

/********************************************************************/
/*                      R Callable functions                        */
/********************************************************************/


/**
 * struct used in the general M-H algorithm admitting R functions.
 * This struct is passed to the metrop_tnorm_rw function.
 *
 */
typedef struct MH_STRUCT{
  SEXP R_fcall;    /**<a user-defined R function */
  SEXP R_env;      /**<where to evaluate the function calls */
} mh_str;

/**
 * evaluate an arbitrary univaraite function stored in mh_str
 *
 * @param x the point at which to evaluate a function
 * @param data a void struct to be converted to mh_str
 *
 * @return function evaluation result
*/

static double R_fun(double x, void *data){
  mh_str *da = data ;
  SEXP R_x, s ;
  PROTECT_INDEX ipx;  
  PROTECT(R_x = allocVector(REALSXP, 1));
  REAL(R_x)[0] = x ;
  SETCADR(da->R_fcall, R_x);          /* assign the argument */
                                      /* evaluate function calls */
  PROTECT_WITH_INDEX(s = eval(da->R_fcall, da->R_env), &ipx);
  REPROTECT(s = coerceVector(s, REALSXP), ipx);
  if (LENGTH(s) != 1)
    error(("objective function evaluates to length %d not 1"), LENGTH(s));
  if (!R_FINITE(REAL(s)[0]) || R_IsNaN(REAL(s)[0]) || R_IsNA(REAL(s)[0])) 
    error("objective funtion evaluates to Inf, NaN or NA");
  UNPROTECT(2);
  return REAL(s)[0];
}

/**
 * generic random walk Metropolis algorithms for an arbitrary unvariate
 * R function
 *
 * @param n number of samples 
 * @param m the starting value 
 * @param sd the proposal standard deviation 
 * @param lb left bound
 * @param rb right bound 
 * @param fun the name of the R function
 * @param rho the environment to evaluate the function
 *
 * @return simulated results
*/

SEXP bcplm_metrop_rw(SEXP n, SEXP m, SEXP sd, SEXP lb, SEXP rb, 
		     SEXP fun, SEXP rho){
  mh_str *da;
  SEXP ans, acc;  
  double sm;                /* the mean used in each iteration */
  int ns = INTEGER(n)[0];
  if (!isFunction(fun)) error(("'fun' is not a function"));
  if (!isEnvironment(rho)) error(("'rho'is not an environment"));

  /* construct the mh_str object */
  da = (mh_str *) R_alloc(1, sizeof(mh_str));
  PROTECT(da->R_fcall = lang2(fun, R_NilValue));
  da->R_env = rho;

  /* run the random walk metropolis algorithm */
  PROTECT(ans = allocVector(REALSXP, ns));
  PROTECT(acc = allocVector(INTSXP, 1));
  INTEGER(acc)[0] = 0;
  GetRNGstate();
  for (int i = 0; i < ns; i++){
    sm = (i) ? REAL(ans)[i - 1] : REAL(m)[0];
    INTEGER(acc)[0] += metrop_tnorm_rw(sm, REAL(sd)[0], REAL(lb)[0], REAL(rb)[0], 
		    REAL(ans) + i, R_fun, da);
  }
  setAttrib(ans, install("accept"), acc);
  UNPROTECT(3);
  PutRNGstate();
  return ans;
}
