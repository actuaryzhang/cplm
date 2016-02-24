/**
 * @file utilities.h
 * @brief header files for utility functions 
 * @author Wayne Zhang                         
*/

#ifndef CPLM_UTILS_H
#define CPLM_UTILS_H

/* inline functions */

/**
 * cumulative sum of a vector of double 
 *
 * @param x double vector to be summed
 * @param n number of elements
 *
 * @return cumulative sum
 */
static R_INLINE double dcumsum(double *x, int n){
  double sm = 0.0 ;
  for (int i = 0; i < n; i++)  sm += x[i] ;
  return sm ;
}

/**
 * get the maximum value of a double vector
 *
 * @param x a double vector
 * @param n length of the vector
 *
 * @return the maximum value
 */
static R_INLINE double dmax (double *x, int n){
  double s = x[0] ;
  for (int i = 1; i < n; i++)
      if (x[i] > s) s = x[i] ; 
  return s ;
}

/**
 * get the minmimum value of a double vector
 *
 * @param x a double vector
 * @param n length of the vector
 *
 * @return the minimum value
 */
static R_INLINE double dmin(double *x, int n){
    double s = x[0];
    for (int i = 1; i < n; i++){
	if (x[i] < s) s = x[i];
    }
    return s;
}

/**
 * get the maximum value of an int vector
 *
 * @param x an int vector
 * @param n length of the vector
 *
 * @return the maximum value
 */
static R_INLINE int imax (int *x, int n){
  int s = x[0] ;
  for (int i = 1; i < n; i++)
      if (x[i] > s) s = x[i] ; 
  return s ;
}


/**
 * Compute sample mean
 *
 * @param n number of samples
 * @param x samples in long vector 
 *
 * @return mean 
 */
static R_INLINE double mean(double *x, int n){
  return dcumsum(x, n) / n ;
}

/**
 * inverse link function
 *
 * @param eta linear predictors 
 * @param lp  link power
 *
 * @return mean
 */ 
static R_INLINE double link_inv(double eta, double lp){
    return (lp == 0) ? exp(eta) : pow(eta, 1.0 / lp);
}

/**
 * derivative d(mu) / d(eta)
 *
 * @param eta linear predictors 
 * @param lp  link power
 *
 * @return d(mu)/d(eta)
 */ 
static R_INLINE double mu_eta(double eta, double lp){
    return (lp == 0) ? exp(eta) : pow(eta, 1.0 / lp - 1.0) / lp ;
}

/**
 * Return the sum of squares of the first n elements of x
 *
 * @param n
 * @param x
 *
 * @return sum of squares
 */
static R_INLINE double sqr_length(double *x, int n)
{
    double ans = 0;
    for (int i = 0; i < n; i++) ans += x[i] * x[i];
    return ans;
}

/**
 * Return the index of the term associated with parameter index ind
 *
 * @param ind an index in [0, Gp[nt] - 1]
 * @param nt total number of terms
 * @param Gp group pointers, a vector of length nt+1 with Gp[0] = 0
 *
 * @return sum of squares
 */
static R_INLINE int Gp_grp(int ind, int nt, const int *Gp)
{
    for (int i = 0; i < nt; i++) if (ind < Gp[i + 1]) return i;
    error(("invalid row index %d (max is %d)"), ind, Gp[nt]);
    return -1;                  /* -Wall */
}



/* prototypes defined in utilities.c */

double var(double *x, int n) ;
void cov(int n, int p, double *x, double *ans);

// matrix computation
void solve_po(int d, double *v, double *iv) ;
void mult_xtx(int m, int n, double *x, double *out) ;
void mult_mv(char *trans, int m, int n, double *A,
             double *x, double *out) ;
void chol(int d, double *v, double *iv);
void solve_po(int d, double *v, double *iv);

// numerical derivatives
void grad(int n, double *x, 
          double (*myfunc)(double *x, void *data), 
          void *data, double *ans) ;
void hess(int n, double *x, 
          double (*myfunc)(double *x, void *data), 
          void *data, double *ans) ;

// wishart simulation
void rwishart(int d, double nu, double *scal, double *out);


#endif
