/**
 * @file cpglmm.h
 * @brief header files for mixed models
 * @author Wayne Zhang                         
*/

#ifndef CPGLMM_H
#define CPGLMM_H

/* cpglmm */
SEXP cpglmm_optimize(SEXP x) ;
SEXP cpglmm_update_mu(SEXP x) ;
SEXP cpglmm_update_u(SEXP x) ;
SEXP cpglmm_update_L(SEXP x) ;
SEXP cpglmm_update_dev(SEXP x, SEXP pm) ;
SEXP cpglmm_setPars(SEXP x, SEXP pm) ;
SEXP cpglmm_ST_getPars(SEXP x);
SEXP cpglmm_ST_chol(SEXP x);
SEXP cpglmm_update_ranef(SEXP x);
SEXP cpglmm_update_RX(SEXP x);
SEXP cpglmm_ghq(SEXP np);
SEXP mer_ST_initialize(SEXP ST, SEXP Gpp, SEXP Zt);
SEXP mer_create_L(SEXP CmP);
SEXP mer_postVar(SEXP x, SEXP which);
SEXP mer_update_projection(SEXP x);
#endif
