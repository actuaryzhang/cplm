/************************************************************/
/*   Function for fitting the Compound Poisson Generalized  */
/*    Linear Mixed models using the Laplace approximation   */
/*    and the adaptive Gaussian-Hermite quadrature method.  */
/*    This is based on the code from the lme4 package.      */
/*             Orignal authors (lme4.0):                    */
/*         Bates D, Maechler M, Bolker B and Walker S       */
/*                                                          */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/


/**
 * @file cpglmm.c
 * @brief Functions that provide the Laplacian approximation  
 * and the adaptive Gaussian-Hermite quadrature methods for the         
 * Compound Poisson Generalized Linear Mixed Model.        
 * The program is based on the C code from lme4, where    
 * I removed the unnecessary pieces for linear mixed      
 * models and nonlinear models, and changed the code for    
 * updating the mean and the deviance and optimization. 
 * The compound Poisson density approximation method is used 
 * for computing the loglikelihood.      
 * @author Wayne Zhang                        
*/

#include "common.h"                /* for common include headers */ 
#include "tweedie.h"               /* for evaluating tweedie densities */
#include "utilities.h"             /* for utility functions */
#include "cpglmm.h"                  /* prototypes */
#include <R_ext/stats_package.h>   /* for S_nlminb_iterate declarations */

/** cholmod_common struct initialized  */
extern cholmod_common c;

/* Constants used in Penalized Least Squares*/
/** Maximum number of iterations in update_u */
#define CM_MAXITER  300
/** Tolerance level for convergence criterion in update_u */
#define CM_TOL      1e-10
/** Minimum step factor in update_u */
#define CM_SMIN     1e-5
/** precision and maximum number of iterations used in GHQ */
#define GHQ_EPS    1e-15
#define GHQ_MAXIT  40

/** positions in the deviance vector */
enum devP {
    ML_POS = 0,			/**<Maximum likelihood estimation criterion  */
    REML_POS,			/**<REML criterion */
    ldL2_POS,			/**<2*log-determinant of L */
    ldRX2_POS,			/**<2*log-determinant of RX */
    sigmaML_POS,		/**<current ML estimate of sigma */
    sigmaREML_POS,		/**<current REML estimate of sigma */
    pwrss_POS,			/**<penalized weighted residual sum of squares */
    disc_POS,			/**<discrepancy */
    usqr_POS,			/**<squared length of u */
    wrss_POS,			/**<weighted residual sum of squares  */
    dev_POS,			/**<deviance - defined for quasi families  */
    llik_POS,			/**<log-likelihood - undefined for quasi families  */
    NULLdev_POS			/**<null deviance */
};

/** positions in the dims vector */
enum dimP {
    nt_POS = 0,			/**<number of terms in random effects */
    n_POS,			/**<number of observations */
    p_POS,			/**<number of fixed-effects parameters */
    q_POS,			/**<number of random effects */
    s_POS,			/**<number of variables in h (1 unless nonlinear) */
    np_POS,			/**<total number of parameters for T and S */
    LMM_POS,			/**<is the model a linear mixed model? */
    isREML_POS,			/**<indicator of REML estimation */
    fTyp_POS,			/**<family type for generalized model */
    lTyp_POS,			/**<link type for generalized model */
    vTyp_POS,			/**<variance type for generalized model */
    nest_POS,			/**<indicator of nested grouping factors */
    useSc_POS,			/**<does the family use a separate scale parameter */
    nAGQ_POS,			/**<number of adaptive Gauss-Hermite quadrature pts */
    verb_POS,			/**<verbose output in mer_optimize? */
    mxit_POS,			/**<maximum # of iterations in mer_optimize */
    mxfn_POS,			/**<maximum # of function evaluations in mer_optimize */
    cvg_POS			/**<convergence indictor from port optimization  */
};

/**
 * Permute the vector src according to perm into dest
 *
 * @param dest destination
 * @param src source
 * @param perm NULL or 0-based permutation of length n
 * @param n length of src, dest and perm
 *
 * @return dest
 *
 * \note If perm is NULL the first n elements of src are copied to dest.
 */
static R_INLINE double*
apply_perm(double *dest, const double *src, const int *perm, int n)
{
    for (int i = 0; i < n; i++) dest[i] = src[perm ? perm[i] : i];
    return dest;
}

/**
 * Extract the parameters from ST list
 *
 * @param x a list object
 * @param pars vector of the appropriate length
 *
 * @return pointer to the parameter vector
 */
static double *ST_getPars(SEXP x, double *pars)
{
    SEXP ST = GET_SLOT(x, install("ST"));
    int nt = LENGTH(ST), pos = 0;
    for (int i = 0; i < nt; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int ncp1 = nci + 1;

	for (int j = 0; j < nci; j++)
	    pars[pos++] = st[j * ncp1];  /**<get parm along the diagonal (T) */
	for (int j = 0; j < (nci - 1); j++) //get S
	    for (int k = j + 1; k < nci; k++)
		pars[pos++] = st[k + j * nci];
    }
    return pars;
}

/**
 * Populate the st, nc and nlev arrays.  Return the maximum element of nc.
 *
 * @param ST pointer to a list (length nt) of matrices
 * @param Gp group pointers (length nt + 1)
 * @param st length nt array of (double*) pointers to be filled with
 * pointers to the contents of the matrices in ST.  Not used if NULL.
 * @param nc length nt array to be filled with the number of columns
 * @param nlev length nt array to be filled with the number of
 *        levels of the grouping factor for each term
 *
 * @return maximum element of nc
 */
/* populate the st, nc and nlev arrays */
int  ST_nc_nlev(const SEXP ST, const int *Gp, double **st, int *nc, int *nlev)
{
    int ans = 0, nt = LENGTH(ST);

    for (int i = 0; i < nt; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	int nci = *INTEGER(getAttrib(STi, R_DimSymbol));

	if (nci > ans) ans = nci;
	if (st) st[i] = REAL(STi);
	nc[i] = nci;
	nlev[i] = (Gp[i + 1] - Gp[i]) / nci;
    }
    return ans;
}

/**
 * Multiply A on the left by T'
 *
 * @param A sparse model matrix
 * @param Gp group pointers
 * @param nc number of columns per term
 * @param nlev number of levels per term
 * @param st ST arrays for each term
 * @param nt number of terms
 *
 */
static void Tt_Zt(CHM_SP A, int *Gp, int *nc, int *nlev, double **st, int nt)
{
    int *ai = (int*)(A->i), *ap = (int *)(A->p);
    double *ax = (double*)(A->x), one[] = {1,0};

    for (int j = 0; j < A->ncol; j++) /* multiply column j by T' */
	for (int p = ap[j]; p < ap[j + 1];) {
	    int i = Gp_grp(ai[p], nt, Gp);

	    if (nc[i] <= 1) p++;
	    else {
		int nr = p;	/* number of rows in `B' in dtrmm call */
		while ((ai[nr] - Gp[i]) < nlev[i]) nr++;
		nr -= p;	/* nr == 1 except in models with carry-over */
		F77_CALL(dtrmm)("R", "L", "N", "U", &nr, nc + i,
				one, st[i], nc + i, ax + p, &nr);
		p += (nr * nc[i]);
	    }
	}
}

/**
 * Evaluate the sparse model matrix A from the Zt and ST slots
 * A = S'T'Z', the transpose of the design matrix with unit random effects 
 *
 * @param x a list object
 */
static void update_A(SEXP x)
{
    CHM_SP A = A_SLOT(x), Zt = Zt_SLOT(x);
    int *Gp = Gp_SLOT(x), *ai = (int*)(A->i), *ap = (int*)(A->p),
	 *zp = (int*)(Zt->p), ncmax,
	nt = DIMS_SLOT(x)[nt_POS];
    int annz = ap[A->ncol], znnz = zp[Zt->ncol];
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    double **st = Alloca(nt, double*), *ax = (double*)(A->x),
	*zx = (double*)(Zt->x);
    R_CheckStack();

    ncmax = ST_nc_nlev(GET_SLOT(x, install("ST")), Gp, st, nc, nlev);

     /* Copy Z' to A unless A has new nonzeros */
    Memcpy(ax, zx, znnz);
    /* When T != I multiply A on the left by T' */
    if (ncmax > 1) Tt_Zt(A, Gp, nc, nlev, st, nt);
				/* Multiply A on the left by S */
    for (int p = 0; p < annz; p++) {
	int i = Gp_grp(ai[p], nt, Gp);
	ax[p] *= st[i][((ai[p] - Gp[i]) / nlev[i]) * (nc[i] + 1)];
    }
}


/**
 * Update the L, sqrtrWt and resid slots.  It is assumed that
 * update_mu has already been called at the current values of u and
 * the model parameters and that A has been updated.
 *
 * @param x pointer to an mer object
 *
 * @return penalized weighted residual sum of squares
 *
 */
static double cp_update_L(SEXP x)
{
    int *dims = DIMS_SLOT(x);
    int n = dims[n_POS] ;
    double  *cx = Cx_SLOT(x), *d = DEV_SLOT(x),
	*res = RESID_SLOT(x), *mu = MU_SLOT(x), *muEta = MUETA_SLOT(x),
	*pwt = PWT_SLOT(x), *sXwt = SXWT_SLOT(x), *srwt = SRWT_SLOT(x),
	*var =  VAR_SLOT(x), *y = Y_SLOT(x), one[] = {1, 0};
    CHM_SP A = A_SLOT(x);
    CHM_FR L = L_SLOT(x);
    R_CheckStack();

    // Update srwt and res. Reevaluate wrss. 
    d[wrss_POS] = 0;
    for (int j = 0; j < n; j++) {
        srwt[j] = sqrt((pwt ? pwt[j] : 1.0) / (var ? var[j] : 1.0)) ;
        res[j] = srwt[j] * (y[j] - mu[j]);
        d[wrss_POS] += res[j] * res[j];
    }
    int  *ap = (int*)A->p;
    double *ax = (double*)(A->x);

    // sXwt is the sqrt of the weight in the penalized regression
    for (int i = 0; i < n; i++) 
          sXwt[i] = (srwt ? srwt[i] : 1) * (muEta ? muEta[i] : 1) ;
    
    //  a scaled version of A: A * sqrt(W) 
    for (int j = 0; j < n; j++)
        for (int p = ap[j]; p < ap[j + 1]; p++)
            cx[p] = ax[p] * sXwt[j] ;

    A->x = (void*)cx;
    
    // compute the cholesky factor of AA'+I with permutation
    if (!M_cholmod_factorize_p(A, one, (int*)NULL, 0, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);

    
    d[ldL2_POS] = M_chm_factor_ldetL2(L);
    d[pwrss_POS] = d[usqr_POS] + d[wrss_POS];
    
    return d[pwrss_POS];
}


/**
 * Update the eta, v, mu, resid and var slots according to the current
 * values of the parameters and u.  Also evaluate d[wrss_POS] using
 * the current sqrtrWt slot.  The sqrtrWt slot is changed in update_L.
 * This function is different from that in lme4, and hence the prefix
 * "cp" is added. Specifically, the update of mu, muEta and variance
 * is from a Tweedie compound Poisson model.
 *
 * @param x pointer to a list object
 *
 * @return penalized, weighted residual sum of squares
 */
static double cp_update_mu(SEXP x)
{
    int *dims = DIMS_SLOT(x);
    int i1 = 1, n = dims[n_POS], p = dims[p_POS];
    double *d = DEV_SLOT(x), *eta = ETA_SLOT(x),
        *mu = MU_SLOT(x),
	*muEta = MUETA_SLOT(x), *offset = OFFSET_SLOT(x),
	*srwt = SRWT_SLOT(x), *res = RESID_SLOT(x),
        *lkp = LKP_SLOT(x), *vp = P_SLOT(x),
        *var = VAR_SLOT(x), 
        *y = Y_SLOT(x), one[] = {1, 0};
    CHM_FR L = L_SLOT(x);
    CHM_SP A = A_SLOT(x);
    CHM_DN Ptu, ceta, cu = AS_CHM_DN(GET_SLOT(x, install("u")));
    R_CheckStack();

    /* eta := offset (offset initialized to zero if NULL) */
    Memcpy(eta, offset, n);
    /* eta := eta + X \beta */
    F77_CALL(dgemv)("N", &n, &p, one, X_SLOT(x), &n,
		    FIXEF_SLOT(x), &i1, one, eta, &i1);
    /* eta := eta + C' P' u */
    Ptu = M_cholmod_solve(CHOLMOD_Pt, L, cu, &c);
    ceta = N_AS_CHM_DN(eta, n, 1);
    R_CheckStack();
    if (!M_cholmod_sdmult(A, 1 , one, one, Ptu, ceta, &c))
	error(_("cholmod_sdmult error returned"));
    M_cholmod_free_dense(&Ptu, &c);

    /* update mu, muEta and var */
    for (int i = 0; i < n; i++){
        mu[i] = link_inv(eta[i], *lkp);
        muEta[i] = mu_eta(eta[i], *lkp);
        var[i] = pow(mu[i], *vp); 
    }
    
    /* update resid slot and d[wrss_POS] */
    d[wrss_POS] = 0;		
    for (int i = 0; i < n; i++) {
	res[i] = (y[i] - mu[i]) * (srwt ? srwt[i] : 1);
	d[wrss_POS] += res[i] * res[i];
    }
    /* store u'u */
    d[usqr_POS] = sqr_length((double*)(cu->x), dims[q_POS]);
    d[pwrss_POS] = d[usqr_POS] + d[wrss_POS];
    return d[pwrss_POS];
}



/**
 * Iterate to determine the conditional modes of the random effects.
 * This is the same as that in lme4. 
 * @param x pointer to a list object
 *
 * @return number of iterations to convergence (0 for non-convergence)
 */

static int cp_update_u(SEXP x)
{
    int *dims = DIMS_SLOT(x);
    int i = 0, n = dims[n_POS], q = dims[q_POS], verb = dims[verb_POS];
    double *Cx = Cx_SLOT(x), 
	*res = RESID_SLOT(x), *u = U_SLOT(x),
	cfac = ((double)n) / ((double)q),
	crit, pwrss, pwrss_old, step;
    double *tmp = Calloc(q, double), *tmp1 = Calloc(q, double),
	*uold = Calloc(q, double), one[] = {1,0}, zero[] = {0,0};
    CHM_FR L = L_SLOT(x);
    CHM_DN cres = N_AS_CHM_DN(res, n, 1),
	ctmp = N_AS_CHM_DN(tmp, q, 1), sol;
    
    R_CheckStack();

    if (!(L->is_ll)) error(_("L must be LL', not LDL'"));
    CHM_SP C  = A_SLOT(x);
    R_CheckStack();
    C->x = (void*)Cx;
    
    // update mu related parms
    AZERO(u, q);
    cp_update_mu(x);
 
    for (i = 0; ; i++) {
	Memcpy(uold, u, q);
	pwrss_old = cp_update_L(x);
        // tmp := PC %*% wtdResid 
	M_cholmod_sdmult(C, 0 , one, zero, cres, ctmp, &c);
	Memcpy(tmp1, tmp, q);
	apply_perm(tmp, tmp1, (int*)L->Perm, q);
				// tmp := tmp - u 
	for (int j = 0; j < q; j++) tmp[j] -= u[j];
				// solve L %*% sol = tmp 
	if (!(sol = M_cholmod_solve(CHOLMOD_L, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_L) failed"));
	Memcpy(tmp, (double*)(sol->x), q);
	M_cholmod_free_dense(&sol, &c);
				// evaluate convergence criterion 
	crit = cfac * sqr_length(tmp, q) / pwrss_old;
	if (crit < CM_TOL) break; // don't do needless evaluations 
				// solve t(L) %*% sol = tmp 
	if (!(sol = M_cholmod_solve(CHOLMOD_Lt, L, ctmp, &c)))
	    error(_("cholmod_solve (CHOLMOD_Lt) failed"));
	Memcpy(tmp, (double*)(sol->x), q);
	M_cholmod_free_dense(&sol, &c);
	for (step = 1; step > CM_SMIN; step /= 2) { // step halving 
	    for (int j = 0; j < q; j++) u[j] = uold[j] + step * tmp[j];
	    pwrss = cp_update_mu(x);
	    if (verb < 0)
		Rprintf("%2d,%8.6f,%12.4g: %15.6g %15.6g %15.6g %15.6g\n",
			i, step, crit, pwrss, pwrss_old, u[1], u[2]);
	    if (pwrss < pwrss_old) {
		pwrss_old = pwrss;
		break;
	    }
	}
	if (step <= CM_SMIN || i > CM_MAXITER) return 0;
    }

    Free(tmp); Free(tmp1); Free(uold);
    return i;
}

/**
 * Update the ST and A slots of an mer object.
 *
 * @param x an mer object
 * @param pars double vector of the appropriate length
 * @return updated deviance
 *
 */
static void
ST_setPars(SEXP x, const double *pars)
{
    int *Gp = Gp_SLOT(x), nt = DIMS_SLOT(x)[nt_POS], pos = 0;
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    double **st = Alloca(nt, double*);
    R_CheckStack();

    ST_nc_nlev(GET_SLOT(x, install("ST")), Gp, st, nc, nlev);
				/* install the parameters in the ST slot */
    for (int i = 0; i < nt; i++) {
	int nci = nc[i], ncp1 = nc[i] + 1;
	double *sti = st[i];

	for (int j = 0; j < nci; j++)
	    sti[j * ncp1] = pars[pos++];
	for (int j = 0; j < (nci - 1); j++)
	    for (int k = j + 1; k < nci; k++)
		sti[k + j * nci] = pars[pos++];
    }
    update_A(x);
}


/**
 * Update the ST, A, fixef, phi and p slots 
 *
 * @param x an cpglmm object
 * @param pars double vector of the appropriate length: theta, beta, log(phi) and p in order
 * @return updated deviance
 *
 */
static void cp_setPars(SEXP x, double *pars){
    int *dims = DIMS_SLOT(x) ;
    int phi_POS = dims[np_POS] + dims[p_POS] ;    
    double  *phi = PHI_SLOT(x),
      *p = P_SLOT(x), *fixef = FIXEF_SLOT(x) ;
    // update ST from parm
    ST_setPars(x, pars);
    // update fixef, phi and p
    Memcpy(fixef, pars + dims[np_POS], dims[p_POS]);
    phi[0] = exp(pars[phi_POS])  ;
    p[0] = pars[phi_POS + 1] ;       
}

/**
 * Compute twice negative marginal logliklihood to be used in optimization. 
 * This is different from that in lme4, so I added the prefix "cp".
 *
 * @param x a cpglmm object
 * @param parm vector of parameters: theta, beta, log(phi) and p in order
 *            
 * @return twice negative marginal logliklihood
 *
 */
static double cp_update_dev(SEXP x, double *parm)
{
    int *dims = DIMS_SLOT(x) ;
    int n = dims[n_POS], q = dims[q_POS], nAGQ = dims[nAGQ_POS] ;
    SEXP flistP = GET_SLOT(x, install("flist"));
    double *d = DEV_SLOT(x), *y = Y_SLOT(x),
      *mu = MU_SLOT(x), *phi = PHI_SLOT(x),
      *p = P_SLOT(x), *pwt = PWT_SLOT(x),
      *u = U_SLOT(x);
    CHM_FR L = L_SLOT(x);
    
    if (parm) cp_setPars(x, parm) ;
    // find conditional mode
    cp_update_u(x);

    if (nAGQ < 1) error(_("nAGQ must be positive"));
    if ((nAGQ > 1) & (LENGTH(flistP) != 1))
      error(_("AGQ method requires a single grouping factor"));
    d[ML_POS] = d[ldL2_POS];
    // Laplace Approximation 
    if (nAGQ == 1) {
      d[disc_POS] = dl2tweedie(n, y, mu, phi[0], p[0], pwt);
      d[ML_POS] += d[disc_POS] + d[usqr_POS] / phi[0] ; 
      return d[ML_POS]; 
    } else {
      // Adaptive Gauss-Hermite quadrature 
      const int nl = nlevels(VECTOR_ELT(flistP, 0));
      const int nre = q / nl; 	             /* number of random effects per level */
      int pos = 0, nt = 1, *fl0 = INTEGER(VECTOR_ELT(flistP, 0)), 
	*pointer = Alloca(nre, int);
      for (int i = 0; i < nre; i++) nt *= nAGQ;
	
      double *ghw = GHW_SLOT(x), *ghx = GHX_SLOT(x),          
	*uold = Memcpy(Calloc(q, double), u, q), 
	*z = Calloc(q, double),              /* current abscissas */
	*ans_tmp = Calloc(n, double),        /* log data density */
	*ans = Calloc(nl, double);           /* loglikelihood by group */
      double mm, W = 0.0;                    /* values needed in AGQ evaluation */
      double sigma = sqrt(phi[0]);
      double **tmp = Calloc(nl, double*);    /* matrix of exponents as in log(exp(x1) + ...) */
      for (int i = 0; i < nl; i++)  tmp[i] = Calloc(nt, double);
      R_CheckStack();
      AZERO(pointer, nre);

      while(pointer[nre - 1] < nAGQ){                 /* the Cartesian product */
	/* update abscissas and weights */
	for(int i = 0; i < nre; ++i){
	  for(int j = 0; j < nl; ++j)
	    z[j + i * nl] = ghx[pointer[i]];          /* zeros corresp to each u */ 
	  W += log(ghw[pointer[i]]) + ghx[pointer[i]] * ghx[pointer[i]]; 
	  /* accumulate zeros and weights (const used in AGQ)*/
	}
	
	/* scale the zeros: u = u_old + sigma * L^{-1} * z */
	CHM_DN cz = N_AS_CHM_DN(z, q, 1), sol;
	if(!(sol = M_cholmod_solve(CHOLMOD_L, L, cz, &c)))
	  error(_("cholmod_solve(CHOLMOD_L) failed"));
	Memcpy(z, (double *)sol->x, q);
	M_cholmod_free_dense(&sol, &c);   
	for(int i = 0; i < q; ++i) 
	  u[i] = uold[i] + M_SQRT2 * sigma * z[i];   
	/* get data likelihood */
	cp_update_mu(x);      
        dtweedie(n, y, mu, phi[0], p[0], pwt, ans_tmp);
	AZERO(ans, nl);
	for (int i = 0; i < n; i++) 
	  ans[fl0[i] - 1] +=  ans_tmp[i] ;        

	/* add the contribution from u */
	for(int i = 0; i < nre; ++i)
	  for(int j = 0; j < nl; ++j)
	    ans[j] -= 0.5 * u[j + i * nl] * u[j + i * nl] / phi[0];
	
	/* the exponent in the AGQ objective */
	for(int i = 0; i < nl; ++i)
	  tmp[i][pos] = ans[i] + W;
	pos++ ;

	/* move pointer to next combination of weights and zeros */
	int count = 0;
	pointer[count]++;
	while(pointer[count] == nAGQ && count < nre - 1){
	  pointer[count] = 0;
	  pointer[++count]++;
	}
	W = 0;
      }

      /* compute log(exp(x1) + exp(x2) + ...) */
      AZERO(ans, nl);
      for (int i = 0; i < nl; i++){
	mm = dmax(tmp[i], nt);                    /* max of the exponents */
	for (int k = 0; k < nt; k++)
	  ans[i] += exp(tmp[i][k] - mm);          /* subtract max before exponentiation */
	d[ML_POS] -= 2 * (log(ans[i]) + mm);      /* -2 loglik */
      }

      /* restore values */
      Memcpy(u, uold, q);
      cp_update_mu(x);
      Free(uold);
      Free(z);
      Free(ans_tmp); 
      Free(ans);
      for (int i = 0; i < nl; i++) Free(tmp[i]);
      Free(tmp);
      return d[ML_POS] ;
    }
}


/**
 * Optimize the  Laplace approximation to the deviance
 * of a cpglmm object. The parameter set includes
 * theta (variance component), beta (fixed effects),
 * phi (scale parameter) and p (index parameter) 
 *
 * @param x pointer to an mer object
 *
 * @return R_NilValue
 */
SEXP cpglmm_optimize(SEXP x)
{
    SEXP ST = GET_SLOT(x, install("ST"));
    int *dims = DIMS_SLOT(x);
    int  nt = dims[nt_POS];
    int nv1 = dims[np_POS] +  dims[p_POS], verb = dims[verb_POS];
    int nv = nv1 +2 ;
    int liv = S_iv_length(OPT, nv), lv = S_v_length(OPT, nv);
    double *g = (double*)NULL, *h = (double*)NULL, fx = R_PosInf;
    double *fixef = FIXEF_SLOT(x);
    int *iv = Alloca(liv, int);
    double *b = Alloca(2 * nv, double), *d = Alloca(nv, double),
	*v = Alloca(lv, double), *xv = Alloca(nv, double);
    R_CheckStack();

    // retrieve parameter values from x 
    ST_getPars(x, xv); /* update xv for theta */
    Memcpy(xv + dims[np_POS], fixef, dims[p_POS]); /*update xv for beta*/
    xv[nv1] = log(PHI_SLOT(x)[0]) ;
    xv[nv1+1] = P_SLOT(x)[0] ;
    
    double eta = 1.e-5; /* estimated rel. error on computed lpdisc */
    v[31] = eta;		/* RFCTOL */
    v[36] = eta;		/* SCTOL */
    v[41] = eta;		/* ETA0 */
				/* initialize the state vectors v and iv */
    S_Rf_divset(OPT, iv, liv, lv, v);
    iv[OUTLEV] = (verb < 0) ? -verb : verb;
    iv[MXFCAL] = dims[mxfn_POS];
    iv[MXITER] = dims[mxit_POS];
				/* set the bounds to plus/minus Infty  */
    for (int i = 0; i < nv; i++) {
	b[2 * i] = R_NegInf; b[2 * i + 1] = R_PosInf; d[i] = 1;
    }
				/* reset lower bounds on elements of theta_S */
    for (int i = 0, pos = 0; i < nt; i++) {
	int nc = *INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol));
	for (int j = 0; j < nc; j++) b[pos + 2 * j] = 0;
	pos += nc * (nc + 1);
    }
    /* reset bound for p */
    b[(nv1 + 1) * 2] = BDP_SLOT(x)[0] ;
    b[nv * 2 - 1] = BDP_SLOT(x)[1] ;

    /* run optimization */
    do {
        fx = cp_update_dev(x, xv);        
	S_nlminb_iterate(b, d, fx, g, h, iv, liv, lv, nv, v, xv);
    } while (iv[0] == 1 || iv[0] == 2);
    fx = cp_update_dev(x, xv) ; //update slots using the values at the minima
    dims[cvg_POS] = iv[0];
    return R_NilValue;
}


/**
 * Create PAX in dest.
 *
 * @param dest values to be calculated
 * @param perm NULL or a 0-based permutation vector defining P
 * @param A sparse matrix
 * @param X dense matrix
 * @param nc number of columns in X
 *
 */
static void
P_sdmult(double *dest, const int *perm, const CHM_SP A,
	 const double *X, int nc)
{
    int *ai = (int*)(A->i), *ap = (int*)(A->p), m = A->nrow, n = A->ncol;
    double *ax = (double*)(A->x), *tmp = Calloc(m, double);
    R_CheckStack();

    for (int k = 0; k < nc; k++) {
	AZERO(tmp, m);
	for (int j = 0; j < n; j++) {
	    for (int p = ap[j]; p < ap[j + 1]; p++)
		tmp[ai[p]] += X[j + k * n] * ax[p];
	}
	apply_perm(dest + k * m, tmp, perm, m);
    }
    Free(tmp);
}


/**
 * Update the RZX and RX slots in an mer object. update_L should be
 * called before update_RX
 *
 * @param x pointer to an mer object
 *
 * @return profiled deviance or REML deviance
 */
static double update_RX(SEXP x)
{
    int *dims = DIMS_SLOT(x), info;
    int n = dims[n_POS], p = dims[p_POS], q = dims[q_POS], s = dims[s_POS];
    double *cx = Cx_SLOT(x), *d = DEV_SLOT(x),
	*RZX = RZX_SLOT(x), *RX = RX_SLOT(x), *sXwt = SXWT_SLOT(x),
	*WX = (double*) NULL, *X = X_SLOT(x),
	mone[] = {-1,0}, one[] = {1,0}, zero[] = {0,0};
    CHM_SP A = A_SLOT(x);
    CHM_FR L = L_SLOT(x);
    CHM_DN cRZX = N_AS_CHM_DN(RZX, q, p), ans;
    R_CheckStack();

    if (sXwt) {			/* Create W^{1/2}GHX in WX */
	WX = Calloc(n * p, double);

	AZERO(WX, n * p);
	for (int j = 0; j < p; j++)
	    for (int k = 0; k < s; k++)
		for (int i = 0; i < n; i++)
		    WX[i + j * n] +=
			sXwt[i + k * n] * X[i + n * (k + j * s)];
	X = WX;
	/* Replace A by C, either just the x component or the entire matrix */
	if (cx) A->x = (void*)cx;
	else {
	    A = Cm_SLOT(x);
	    R_CheckStack();
	}
    }
				/* solve L %*% RZX = PAW^{1/2}GHX */
    P_sdmult(RZX, (int*)L->Perm, A, X, p); /* right-hand side */
    ans = M_cholmod_solve(CHOLMOD_L, L, cRZX, &c); /* solution */
    Memcpy(RZX, (double*)(ans->x), q * p);
    M_cholmod_free_dense(&ans, &c);
    				/* downdate X'X and factor  */
    F77_CALL(dsyrk)("U", "T", &p, &n, one, X, &n, zero, RX, &p); /* X'X */
    F77_CALL(dsyrk)("U", "T", &p, &q, mone, RZX, &q, one, RX, &p);
    F77_CALL(dpotrf)("U", &p, RX, &p, &info);
    if (info)
	error(_("Downdated X'X is not positive definite, %d."), info);
				/* accumulate log(det(RX)^2)  */
    d[ldRX2_POS] = 0;
    for (int j = 0; j < p; j++) d[ldRX2_POS] += 2 * log(RX[j * (p + 1)]);

    if (WX) Free(WX);
    return d[ML_POS];
}



/**
 * Generate zeros and weights of Hermite polynomial of order N, for AGQ method
 *
 * changed from fortran in package 'glmmML'
 * @param N order of the Hermite polynomial
 * @param x zeros of the polynomial, abscissas for AGQ
 * @param w weights used in AGQ
 *
 */

static void internal_ghq(int N, double *x, double *w)
{
    int NR, IT, I, K, J;
    double Z = 0, HF = 0, HD = 0;
    double Z0, F0, F1, P, FD, Q, WP, GD, R, R1, R2;
    double HN = 1/(double)N;
    double *X = Calloc(N + 1, double), *W = Calloc(N + 1, double);

    for(NR = 1; NR <= N / 2; NR++){
	if(NR == 1)
	    Z = -1.1611 + 1.46 * sqrt((double)N);
	else
	    Z -= HN * (N/2 + 1 - NR);
	for (IT = 0; IT <= GHQ_MAXIT; IT++) {
	    Z0 = Z;
	    F0 = 1.0;
	    F1 = 2.0 * Z;
	    for(K = 2; K <= N; ++K){
		HF = 2.0 * Z * F1 - 2.0 * (double)(K - 1.0) * F0;
		HD = 2.0 * K * F1;
		F0 = F1;
		F1 = HF;
	    }
	    P = 1.0;
	    for(I = 1; I <= NR-1; ++I){
		P *= (Z - X[I]);
	    }
	    FD = HF / P;
	    Q = 0.0;
	    for(I = 1; I <= NR - 1; ++I){
		WP = 1.0;
		for(J = 1; J <= NR - 1; ++J){
		    if(J != I) WP *= ( Z - X[J] );
		}
		Q += WP;
	    }
	    GD = (HD-Q*FD)/P;
	    Z -= (FD/GD);
	    if (fabs((Z - Z0) / Z) < GHQ_EPS) break;
	}

	X[NR] = Z;
	X[N+1-NR] = -Z;
	R=1.0;
	for(K = 1; K <= N; ++K){
	    R *= (2.0 * (double)K );
	}
	W[N+1-NR] = W[NR] = 3.544907701811 * R / (HD*HD);
    }

    if( N % 2 ){
	R1=1.0;
	R2=1.0;
	for(J = 1; J <= N; ++J){
	    R1=2.0*R1*J;
	    if(J>=(N+1)/2) R2 *= J;
	}
	W[N/2+1]=0.88622692545276*R1/(R2*R2);
	X[N/2+1]=0.0;
    }

    Memcpy(x, X + 1, N);
    Memcpy(w, W + 1, N);

    if(X) Free(X);
    if(W) Free(W);
}


/**
 * Update the contents of the ranef slot in an mer object using the
 * current contents of the u and ST slots.
 *
 * b = T  %*% S %*% t(P) %*% u
 *
 * @param x an mer object
 */
static void update_ranef(SEXP x)
{
    int *Gp = Gp_SLOT(x), *dims = DIMS_SLOT(x), *perm = PERM_VEC(x);
    int nt = dims[nt_POS], q = dims[q_POS];
    double *b = RANEF_SLOT(x), *u = U_SLOT(x), one[] = {1,0};
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    double **st = Alloca(nt, double*);
    R_CheckStack();

    ST_nc_nlev(GET_SLOT(x, install("ST")), Gp, st, nc, nlev);
				/* inverse permutation */
    for (int i = 0; i < q; i++) b[perm[i]] = u[i];
    for (int i = 0; i < nt; i++) {
	for (int k = 0; k < nc[i]; k++) { /* multiply by \tilde{S}_i */
	    double dd = st[i][k * (nc[i] + 1)];
	    int base = Gp[i] + k * nlev[i];
	    for (int kk = 0; kk < nlev[i]; kk++) b[base + kk] *= dd;
	}
	if (nc[i] > 1) {	/* multiply by \tilde{T}_i */
	    F77_CALL(dtrmm)("R", "L", "T", "U", nlev + i, nc + i, one,
			    st[i], nc + i, b + Gp[i], nlev + i);
	}
    }
}


/*******************************************************
 *            functions callable from R                *
********************************************************/

/**
 * R callable function to update L
 *
 * @param x a R list object
 *
 * @return penalized, weighted residual sum of squares
*/
SEXP cpglmm_update_L(SEXP x){
    return ScalarReal(cp_update_L(x));
}

/**
 * R callable function to update mu
 *
 * @param x a R list object
 *
 * @return penalized, weighted residual sum of squares
*/
SEXP cpglmm_update_mu(SEXP x){
    return ScalarReal(cp_update_mu(x));
}


/**
 * R callable function to update u
 *
 * @param x a R list object
 *
 * @return number of iterations to convergence (0 for non-convergence)
*/
SEXP cpglmm_update_u(SEXP x){
    return ScalarInteger(cp_update_u(x));
}

/**
 * R callable function to update deviance
 *
 * @param x a cpglmm object
 * @param pm vector of parameters: theta, beta, log(phi) and p in order
 *
 * @return twice negative loglikelihood
*/

SEXP cpglmm_update_dev(SEXP x, SEXP pm){  
  return ScalarReal(cp_update_dev(x, (pm != R_NilValue) ? REAL(pm) : (double*) NULL));            
}


/**
 * R callable function to update parameters
 *
 * @param x a cpglmm object
 * @param pm vector of parameters: theta, beta, log(phi) and p in order
 *
*/
SEXP cpglmm_setPars(SEXP x, SEXP pm){  
    cp_setPars(x, REAL(pm)) ;
    return R_NilValue ;      
}

/**
 * Extract the parameters from the ST slot of an mer object
 *
 * @param x an mer object
 *
 * @return pointer to a REAL vector
 */
SEXP cpglmm_ST_getPars(SEXP x)
{
    SEXP ans = PROTECT(allocVector(REALSXP, DIMS_SLOT(x)[np_POS]));
    ST_getPars(x, REAL(ans));

    UNPROTECT(1);
    return ans;
}


/**
 * Return a list of (upper) Cholesky factors from the ST list
 *
 * @param x an mer object
 *
 * @return a list of upper Cholesky factors
 */
SEXP cpglmm_ST_chol(SEXP x)
{
  SEXP ans = PROTECT(duplicate(GET_SLOT(x, install("ST"))));
    int nt = DIMS_SLOT(x)[nt_POS];
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    double **st = Alloca(nt, double*);
    R_CheckStack();

    ST_nc_nlev(ans, Gp_SLOT(x), st, nc, nlev);
    for (int k = 0; k < nt; k++) {
	if (nc[k] > 1) {	/* nothing to do for nc[k] == 1 */
	    int nck = nc[k], nckp1 = nc[k] + 1;
	    double *stk = st[k];

	    for (int j = 0; j < nck; j++) {
		double dd = stk[j * nckp1]; /* diagonal el */
		for (int i = j + 1; i < nck; i++) {
		    stk[j + i * nck] = dd * stk[i + j * nck];
		    stk[i + j * nck] = 0;
		}
	    }
	}
    }

    UNPROTECT(1);
    return ans;
}



/**
 * Externally callable update_ranef.
 * Update the contents of the ranef slot in an mer object.  For a
 * linear mixed model the conditional estimates of the fixed effects
 * and the conditional mode of u are evaluated first.
 *
 * @param x an mer object
 *
 * @return R_NilValue
 */
SEXP cpglmm_update_ranef(SEXP x)
{
    update_ranef(x);
    return R_NilValue;
}


/**
 * Externally callable update_RX
 *
 * @param x pointer to an mer object
 *
 * @return profiled deviance or REML deviance
 */
SEXP cpglmm_update_RX(SEXP x)
{
    return ScalarReal(update_RX(x));
}

/**
 * Return zeros and weights of Hermite polynomial of order n as a list
 *
 * @param np pointer to a scalar integer SEXP
 * @return a list with two components, the abscissas and the weights.
 *
 */
SEXP cpglmm_ghq(SEXP np)
{
    int n = asInteger(np);
    SEXP ans = PROTECT(allocVector(VECSXP, 2));

    if (n < 1) n = 1;
    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, n ));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, n ));

    internal_ghq(n, REAL(VECTOR_ELT(ans, 0)), REAL(VECTOR_ELT(ans, 1)));
    UNPROTECT(1);
    return ans;
}



/************************************************************
*  The following functions are simply copied from
*  lme4.0 since it is replaced by a new implementation
*  that is no longer compatible with cplm
************************************************************/

/**
 * Evaluate starting estimates for the elements of ST
 *
 * @param ST pointers to the nt ST factorizations of the diagonal
 *     elements of Sigma
 * @param Gpp length nt+1 vector of group pointers for the rows of Zt
 * @param Zt transpose of Z matrix
 *
 */
SEXP mer_ST_initialize(SEXP ST, SEXP Gpp, SEXP Zt)
{
    int *Gp = INTEGER(Gpp),
	*Zdims = INTEGER(GET_SLOT(Zt, install("Dim"))),
	*zi = INTEGER(GET_SLOT(Zt, install("i"))), nt = LENGTH(ST);
    int *nc = Alloca(nt, int), *nlev = Alloca(nt, int),
	nnz = INTEGER(GET_SLOT(Zt, install("p")))[Zdims[1]];
    double *rowsqr = Calloc(Zdims[0], double),
	**st = Alloca(nt, double*),
	*zx = REAL(GET_SLOT(Zt, install("x")));
    R_CheckStack();

    ST_nc_nlev(ST, Gp, st, nc, nlev);
    AZERO(rowsqr, Zdims[0]);
    for (int i = 0; i < nnz; i++) rowsqr[zi[i]] += zx[i] * zx[i];
    for (int i = 0; i < nt; i++) {
	AZERO(st[i], nc[i] * nc[i]);
	for (int j = 0; j < nc[i]; j++) {
	    double *stij = st[i] + j * (nc[i] + 1);
	    for (int k = 0; k < nlev[i]; k++)
		*stij += rowsqr[Gp[i] + j * nlev[i] + k];
	    *stij = sqrt(nlev[i]/(0.375 * *stij));
	}
    }
    Free(rowsqr);
    return R_NilValue;
}




/**
 * Create and initialize L
 *
 * @param CmP pointer to the model matrix for the orthogonal random
 * effects (transposed)
 *
 * @return L
 */
SEXP mer_create_L(SEXP CmP)
{
    double one[] = {1, 0};
    CHM_SP Cm = AS_CHM_SP(CmP);
    CHM_FR L;
    R_CheckStack();

    L = M_cholmod_analyze(Cm, &c);
    if (!M_cholmod_factorize_p(Cm, one, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);

    return M_chm_factor_to_SEXP(L, 1);
}



/**
 * Extract the conditional variances of the random effects in an mer
 * object.  Some people called these posterior variances, hence the name.
 *
 * @param x pointer to an mer object
 * @param which pointer to a logical vector
 *
 * @return pointer to a list of arrays
 */
SEXP mer_postVar(SEXP x, SEXP which)
{
    int *Gp = Gp_SLOT(x), *dims = DIMS_SLOT(x), *ww;
    SEXP ans, flistP = GET_SLOT(x, install("flist"));
    const int nf = LENGTH(flistP), nt = dims[nt_POS], q = dims[q_POS];
    int nr = 0, pos = 0;
    /* int *asgn = INTEGER(getAttrib(flistP, install("assign"))); */
    double *vv, one[] = {1,0}, sc;
    CHM_SP sm1, sm2;
    CHM_DN dm1;
    CHM_FR L = L_SLOT(x);
    int *Perm = (int*)(L->Perm), *iperm = Calloc(q, int),
	*nc = Alloca(nt, int), *nlev = Alloca(nt, int);
    double **st = Alloca(nt, double*);
    R_CheckStack();

/* FIXME: Write the code for nt != nf */
    if (nt != nf) error(_("Code not written yet"));
				/* determine length of list to return */
    if (!isLogical(which) || LENGTH(which) != nf)
	error(_("which must be a logical vector of length %d"), nf);
    ww = LOGICAL(which);
    for (int i = 0; i < nt; i++) if (ww[i]) nr++;
    if (!nr) return(allocVector(VECSXP, 0));
    ans = PROTECT(allocVector(VECSXP, nr));

    ST_nc_nlev(GET_SLOT(x, install("ST")), Gp, st, nc, nlev);
    for (int j = 0; j < q; j++) iperm[Perm[j]] = j; /* inverse permutation */
    sc = dims[useSc_POS] ?
	(DEV_SLOT(x)[dims[isREML_POS] ? sigmaREML_POS : sigmaML_POS]) : 1;
    for (int i = 0; i < nt; i++) {
	if (ww[i]) {
	    const int nci = nc[i];
	    const int ncisqr = nci * nci;
	    CHM_SP rhs = M_cholmod_allocate_sparse(q, nci, nci,
						   1/*sorted*/, 1/*packed*/,
						   0/*stype*/, CHOLMOD_REAL, &c);

	    SET_VECTOR_ELT(ans, pos, alloc3DArray(REALSXP, nci, nci, nlev[i]));
	    vv = REAL(VECTOR_ELT(ans, pos));
	    pos++;
	    for (int j = 0; j <= nci; j++) ((int *)(rhs->p))[j] = j;
	    for (int j = 0; j < nci; j++)
		((double *)(rhs->x))[j] = st[i][j * (nci + 1)] * sc;
	    for (int k = 0; k < nlev[i]; k++) {
		double *vvk = vv + k * ncisqr;
		for (int j = 0; j < nci; j++)
		    ((int*)(rhs->i))[j] = iperm[Gp[i] + k + j * nlev[i]];
		sm1 = M_cholmod_spsolve(CHOLMOD_L, L, rhs, &c);
		sm2 = M_cholmod_transpose(sm1, 1 /*values*/, &c);
		M_cholmod_free_sparse(&sm1, &c);
		sm1 = M_cholmod_aat(sm2, (int*)NULL, (size_t)0, 1 /*mode*/, &c);
		dm1 = M_cholmod_sparse_to_dense(sm1, &c);
		M_cholmod_free_sparse(&sm1, &c); M_cholmod_free_sparse(&sm2, &c);
		Memcpy(vvk, (double*)(dm1->x), ncisqr);
		M_cholmod_free_dense(&dm1, &c);
		if (nci > 1) {
		    F77_CALL(dtrmm)("L", "L", "N", "U", nc + i, nc + i,
				    one, st[i], nc + i, vvk, nc + i);
		    F77_CALL(dtrmm)("R", "L", "T", "U", nc + i, nc + i,
				    one, st[i], nc + i, vvk, nc + i);
		}
	    }
	    M_cholmod_free_sparse(&rhs, &c);
	}
    }
    Free(iperm);
    UNPROTECT(1);
    return ans;
}



/**
 * Update the projections of the response vector onto the column
 * spaces of the random effects and the fixed effects.  This function
 * is needed separately for the one-argument form of the anova function.
 *
 * @param x an mer object
 * @param pb position to store the random-effects projection
 * @param pbeta position to store the fixed-effects projection
 */
static void lmm_update_projection(SEXP x, double *pb, double *pbeta)
{
    int *dims = DIMS_SLOT(x), i1 = 1;
    int n = dims[n_POS], p = dims[p_POS], q = dims[q_POS];
    double *WX = (double*) NULL, *X = X_SLOT(x),
	*cx = Cx_SLOT(x), *d = DEV_SLOT(x),
	*off = OFFSET_SLOT(x), *RZX = RZX_SLOT(x),
	*RX = RX_SLOT(x), *sXwt = SXWT_SLOT(x),
	*wy = (double*)NULL, *y = Y_SLOT(x),
	mone[] = {-1,0}, one[] = {1,0}, zero[] = {0,0};
    CHM_SP A = A_SLOT(x);
    CHM_FR L = L_SLOT(x);
    CHM_DN cpb = N_AS_CHM_DN(pb, q, 1), sol;
    R_CheckStack();

    wy = Calloc(n, double);
    for (int i = 0; i < n; i++) wy[i] = y[i] - (off ? off[i] : 0);
    if (sXwt) {		     /* Replace X by weighted X and weight wy */
	if (!cx) error(_("Cx slot has zero length when sXwt does not."));

	A->x = (void*)cx;
	WX = Calloc(n * p, double);

	for (int i = 0; i < n; i++) {
	    wy[i] *= sXwt[i];
	    for (int j = 0; j < p; j++)
		WX[i + j * n] = sXwt[i] * X[i + j * n];
	}
	X = WX;
    }
				/* solve L del1 = PAy */
    P_sdmult(pb, (int*)L->Perm, A, wy, 1);
    sol = M_cholmod_solve(CHOLMOD_L, L, cpb, &c);
    Memcpy(pb, (double*)sol->x, q);
    M_cholmod_free_dense(&sol, &c);
				/* solve RX' del2 = X'y - RZX'del1 */
    F77_CALL(dgemv)("T", &n, &p, one, X, &n,
		    wy, &i1, zero, pbeta, &i1);
    F77_CALL(dgemv)("T", &q, &p, mone, RZX, &q,
		    pb, &i1, one, pbeta, &i1);
    F77_CALL(dtrsv)("U", "T", "N", &p, RX, &p, pbeta, &i1);
    d[pwrss_POS] = sqr_length(wy, n)
	- (sqr_length(pbeta, p) + sqr_length(pb, q));
    if (d[pwrss_POS] < 0)
	error(_("Calculated PWRSS for a LMM is negative"));
    Free(wy);
    if (WX) Free(WX);
}



/**
 * Externally callable lmm_update_projection.
 * Create the projections onto the column spaces of the random effects
 * and the fixed effects.
 *
 * @param x an mer object
 *
 * @return a list with two elements, both REAL vectors
 */
SEXP mer_update_projection(SEXP x)
{
    SEXP ans = PROTECT(allocVector(VECSXP, 2));
    int *dims = DIMS_SLOT(x);

    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, dims[q_POS]));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, dims[p_POS]));
    lmm_update_projection(x,
			  /* pb = */
			  REAL(VECTOR_ELT(ans, 0)),
			  /* pbeta = */
			  REAL(VECTOR_ELT(ans, 1)));
    UNPROTECT(1);
    return ans;
}
