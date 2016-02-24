/**
 * @file init.c
 * @brief Register native routines and  initialize
 * cholmod_common struct used in mixed effect models     
 * @author Wayne Zhang                          
*/

#include "common.h"
#include "cpglmm.h"
#include "bcplm.h"
#include "mh.h"
#include <R_ext/Rdynload.h>

/** utitlity macro in registering native routines */
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

/** function that registers native routines  */
static R_CallMethodDef CallEntries[] = {

    CALLDEF(cpglmm_optimize, 1),
    CALLDEF(cpglmm_update_mu, 1),
    CALLDEF(cpglmm_update_u, 1),
    CALLDEF(cpglmm_update_L, 1),
    CALLDEF(cpglmm_update_dev, 2),
    CALLDEF(cpglmm_setPars, 2),    
    CALLDEF(bcplm_mcmc, 1),
    CALLDEF(bcplm_update_mu, 1),   
    CALLDEF(bcplm_post_p, 2),   
    CALLDEF(bcplm_post_phi, 2),
    CALLDEF(bcplm_post_betak, 2),
    CALLDEF(bcplm_post_uk, 2),
    CALLDEF(bcplm_metrop_rw, 7),
    CALLDEF(mer_ST_initialize, 3),
    CALLDEF(mer_create_L, 1),
    CALLDEF(mer_postVar, 2),
    CALLDEF(mer_update_projection, 1),
    CALLDEF(cpglmm_ghq, 1),
    CALLDEF(cpglmm_update_RX, 1),
    CALLDEF(cpglmm_update_ranef, 1),
    CALLDEF(cpglmm_ST_chol, 1),
    CALLDEF(cpglmm_ST_getPars, 1),
    {NULL, NULL, 0}
};

/** cholmod_common struct local to the cplm package */
cholmod_common c;

/** This is the CHOLMOD error handler from lme4*/
void attribute_hidden
cplm_R_cholmod_error(int status, const char *file, int line, const char *message)
{
    if(status < 0) {
#ifdef Failure_in_matching_Matrix
/* This fails unexpectedly with
 *  function 'cholmod_l_defaults' not provided by package 'Matrix'
 * from ../tests/lmer-1.R 's  (l.68)  lmer(y ~ habitat + (1|habitat*lagoon)
 */
	M_cholmod_defaults(&c);/* <--- restore defaults,
				* as we will not be able to .. */
	c.final_ll = 1;	    /* LL' form of simplicial factorization */
#endif

	error(_("Cholmod error '%s' at file:%s, line %d"), message, file, line);
    }
    else
	warning("Cholmod warning '%s' at file:%s, line %d",
		message, file, line);
}

/** Initializer for cplm, called upon loading the package.
 *
 *  Initialize CHOLMOD and require the LL' form of the factorization.
 *  Install the symbols to be used by functions in the package.
 */
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
void R_init_cplm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

    M_R_cholmod_start(&c);
    c.final_ll = 1;	    /* LL' form of simplicial factorization */

    /* need own error handler, that resets  final_ll (after *_defaults()) : */
    c.error_handler = cplm_R_cholmod_error;

}

/** Finalizer for cplm called upon unloading the package.
 *
 */
void R_unload_cplm(DllInfo *dll){
    M_cholmod_finish(&c);
}
