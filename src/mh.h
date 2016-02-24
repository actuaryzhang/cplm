/**
 * @file mh.h
 * @brief header files for Metropolis-Hastings algorithms
 * @author Wayne Zhang                         
*/

#ifndef CPLM_MH_H
#define CPLM_MH_H

/* one-dimensional update with a truncated normal proposal */
int metrop_tnorm_rw(double m, double sd, double lb, double rb, double *sn, 
		     double (*myfunc)(double x, void *data), void *data);
/* R callable 1-d random walk metropolis */
SEXP bcplm_metrop_rw(SEXP n, SEXP m, SEXP sd, SEXP lb, SEXP rb, 
		     SEXP fun, SEXP rho);

/* multivariate update with a multivariate normal proposal */
double dmvnorm(int d, double *x, double *m, double *iv) ;
int metrop_mvnorm_rw(int d, double *m, double *v, double *sn, 
		     double (*myfunc)(double *x, void *data), void *data);


#endif 
