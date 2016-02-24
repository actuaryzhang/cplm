/************************************************************/
/*         Common header file for cplm functions            */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/

/**
 * @file common.h
 * @brief header files to be included in the cplm C functions
 * @author Wayne Zhang                         
*/

#ifndef CPLM_COMMON_H 
#define CPLM_COMMON_H 

/* common headers for R related definitions */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>        /* for BLAS and Lapack related */
#include "Matrix.h"		 /* for cholmod functions and S4 structures (GET_SLOT)*/

/** zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

/* When appropriate, alloca is cleaner than malloc/free.  The storage
 * is freed automatically on return from a function. When using gcc the
 * builtin version is much faster. */
#ifdef __GNUC__
#undef alloca
#define alloca(x) __builtin_alloca((x))
#else
/* this is necessary (and sufficient) for Solaris 10: */
#ifdef __sun
#include <alloca.h>
#endif
#endif

/** alloca n elements of type t */
#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

/** Allow for translation of error messages */
#ifdef ENABLE_NLS		
#include <libintl.h>
#define _(String) dgettext ("cplm", String)
#else
#define _(String) (String)
#endif


/**
 * Extract the slot named str from the object obj and return a null pointer
 * if the slot has length 0 or a pointer to the REAL contents.
 *
 * @param obj pointer to an S4 object
 * @param str pointer to a symbol naming the slot to extract
 *
 * @return pointer to the REAL contents, if nonzero length, otherwise
 * a NULL pointer
 *
 */
static R_INLINE double *SLOT_REAL_NULL(SEXP obj, char *str)
{
    SEXP pt = GET_SLOT(obj, install(str));
    return LENGTH(pt) ? REAL(pt) : (double*) NULL; 
}

/************************************************/
/*           Slots used in cpglmm               */  
/************************************************/

/** Return the double pointer to the X slot */
#define X_SLOT(x) SLOT_REAL_NULL(x, "X")

/** Return the double pointer to the y slot */
#define Y_SLOT(x) SLOT_REAL_NULL(x, "y")

/** Allocate (alloca) a cholmod_sparse struct, populate it with values
 * from the Zt slot and return the pointer. */
#define Zt_SLOT(x) AS_CHM_SP(GET_SLOT(x, install("Zt")))

/** Return the double pointer to the offset slot or (double*) NULL if
 * offset has length 0) */
#define OFFSET_SLOT(x) SLOT_REAL_NULL(x, "offset")

/** Return the double pointer to the pWt slot or (double*) NULL if
 * pWt has length 0) */
#define PWT_SLOT(x) SLOT_REAL_NULL(x, "pWt")

/** Return the integer pointer to the dims slot */
#define DIMS_SLOT(x) INTEGER(GET_SLOT(x, install("dims")))

/** Return the double pointer to the fixef slot */
#define FIXEF_SLOT(x) SLOT_REAL_NULL(x, "fixef")

/** Return the double pointer to the u slot */
#define U_SLOT(x) SLOT_REAL_NULL(x, "u")

/** Return the double pointer to the eta slot */
#define ETA_SLOT(x) SLOT_REAL_NULL(x, "eta")

/** Return the double pointer to the mu slot */
#define MU_SLOT(x) SLOT_REAL_NULL(x, "mu")

/** Return the double pointer to the muEta slot or (double*) NULL if
 * muEta has length 0) */
#define MUETA_SLOT(x) SLOT_REAL_NULL(x, "muEta")

/** Return the double pointer to the var slot or (double*) NULL if
 * var has length 0) */
#define VAR_SLOT(x) SLOT_REAL_NULL(x, "var")

/** Return the double pointer to the resid slot */
#define RESID_SLOT(x) SLOT_REAL_NULL(x, "resid")

/** Allocate (alloca) a cholmod_sparse struct, populate it with values
 * from the A slot and return the pointer. */
#define A_SLOT(x) AS_CHM_SP(GET_SLOT(x, install("A")))

/** Allocate (alloca) a cholmod_factor struct, populate it with values
 * from the L slot and return the pointer. */
#define L_SLOT(x) AS_CHM_FR(GET_SLOT(x, install("L")))

/** Return the integer pointer to the Gp slot */
#define Gp_SLOT(x) INTEGER(GET_SLOT(x, install("Gp")))

/** Return the double pointer to the Cx slot or (double*) NULL if
 * Cx has length 0) */
#define Cx_SLOT(x) SLOT_REAL_NULL(x, "Cx")

/** Return the double pointer to the deviance slot */
#define DEV_SLOT(x) SLOT_REAL_NULL(x, "deviance")

/** Return the double pointer to the sqrtrWt slot or (double*) NULL if
 *  sqrtrWt has length 0) */
#define SRWT_SLOT(x) SLOT_REAL_NULL(x, "sqrtrWt")

/** Return the double pointer to the sqrtXWt slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define SXWT_SLOT(x) SLOT_REAL_NULL(x, "sqrtXWt")

/** Return the double pointer to the p slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define P_SLOT(x) SLOT_REAL_NULL(x, "p")

/** Return the double pointer to the phi slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define PHI_SLOT(x) SLOT_REAL_NULL(x, "phi")

/** Return the double pointer to the link.power slot or (double*) NULL if
 *  sqrtXWt has length 0) */
#define LKP_SLOT(x) SLOT_REAL_NULL(x, "link.power")

/** Return the double pointer to the bound.p slot  */
#define BDP_SLOT(x) SLOT_REAL_NULL(x, "bound.p")

/** Return the integer pointer to the permutation vector in the L slot */
#define PERM_VEC(x) INTEGER(GET_SLOT(GET_SLOT(x, install("L")), install("perm")))

/** Return the double pointer to the ranef slot or (double*) NULL if
 *  ranef has length 0) */
#define RANEF_SLOT(x) SLOT_REAL_NULL(x, "ranef")

/** Return the double pointer to the ghw slot */
#define GHW_SLOT(x) SLOT_REAL_NULL(x, "ghw")

/** Return the double pointer to the ghx slot */
#define GHX_SLOT(x) SLOT_REAL_NULL(x, "ghx")

/** Return the double pointer to the RX slot */
#define RX_SLOT(x) SLOT_REAL_NULL(x, "RX")

/** Return the double pointer to the RZX slot */
#define RZX_SLOT(x) SLOT_REAL_NULL(x, "RZX")

/** Allocate (alloca) a cholmod_sparse struct, populate it with values
 * from the Cm slot and return the pointer. */
#define Cm_SLOT(x) AS_CHM_SP(GET_SLOT(x, install("Cm")))


/************************************************/
/*       Additional slots used in bcplm         */  
/************************************************/

/** Return the integer pointer to the ygt0 slot or (double*) NULL if
 * pWt has length 0) */
#define YPO_SLOT(x) INTEGER(GET_SLOT(x, install("ygt0")))

/** Return the double pointer to the bound.phi slot  */
#define BDPHI_SLOT(x) SLOT_REAL_NULL(x, "bound.phi")

/** Return the double pointer to the pbeta.mean slot  */
#define PBM_SLOT(x) SLOT_REAL_NULL(x, "pbeta.mean")

/** Return the double pointer to the pbeta.var slot  */
#define PBV_SLOT(x) SLOT_REAL_NULL(x, "pbeta.var")

/** Return the double pointer to the mh.sd slot  */
#define MHSD_SLOT(x) SLOT_REAL_NULL(x, "mh.sd")

/** Return the integer pointer to the k slot */
#define K_SLOT(x) INTEGER(GET_SLOT(x, install("k")))

/** Return the double pointer to the cllik slot  */
#define CLLIK_SLOT(x) SLOT_REAL_NULL(x, "cllik")

/** Return the double pointer to the Xb slot  */
#define XB_SLOT(x) SLOT_REAL_NULL(x, "Xb")

/** Return the double pointer to the Zu slot  */
#define ZU_SLOT(x) SLOT_REAL_NULL(x, "Zu")

/** Return the integer pointer to the ncol slot */
#define NCOL_SLOT(x) INTEGER(GET_SLOT(x, install("ncol")))

/** Return the integer pointer to the nlev slot */
#define NLEV_SLOT(x) INTEGER(GET_SLOT(x, install("nlev")))

/** Return the double pointer to the accept slot  */
#define ACC_SLOT(x) SLOT_REAL_NULL(x, "accept")

#endif


