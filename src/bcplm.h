/**
 * @file bcplm.h
 * @brief header files for Bayesian models
 * @author Wayne Zhang                         
*/

#ifndef BCPLM_H
#define BCPLM_H

/* bcplm */
SEXP bcplm_mcmc(SEXP x);
SEXP bcplm_update_mu(SEXP da);
SEXP bcplm_post_p(SEXP x, SEXP da);
SEXP bcplm_post_phi(SEXP x, SEXP da);
SEXP bcplm_post_betak(SEXP x, SEXP da);
SEXP bcplm_post_uk(SEXP x, SEXP da);


#endif

