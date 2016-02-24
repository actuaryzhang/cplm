/************************************************************/
/*    Markov Chain Monte Carlo algorithms for the Bayesian  */
/*  Compound Poisson Linear Model via density evaluation    */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/

/**
 * @file bcplm.c
 * @brief Functions for implementing the MCMC algorithm for 
 * Compound Poisson Linear Models using direct tweedie 
 * density evaluation
 * @author Wayne Zhang                         
 */

/** common R headers and macros that extract slots */ 
#include "common.h"
/** tweedie density evaluation functions */
#include "tweedie.h"
/** univariate/multivariate random walk Metropolis-Hastings algorithms */
#include "mh.h"
/** utility and inline functions */
#include "utilities.h"
/** declarations for R callable functions */
#include "bcplm.h"

/** cholmod_common struct initialized  */
extern cholmod_common c;

/** target acceptance rate for tuning the univariate scale parameter */
#define MH_ACC_TAR 0.5
/** threshold in specifying the interval of allowed acceptance rates */
#define MH_ACC_TOL 0.1
/** the number of mark times required for a parameter to be deemed as fully tuned */
#define MH_ACC_MARK 3
/** the shape parameter of the inverse-Gamma prior */
#define IG_SHAPE 0.001
/** the scale parameter of the inverse-Gamma prior */
#define IG_SCALE 0.001

/** positions in the dims vector in Bayesian models */
enum dimB {
  nO_POS = 0,			/**<number of observations */
  nB_POS,			/**<number of fixed effects */
  nP_POS,			/**<number of positive observations */
  nT_POS,                       /**<number of terms (grouping factors) in bcpglmm */
  nU_POS,                       /**<number of random coefficients */
  nA_POS,                       /**<number of parameters in Gibbs */
  chn_POS,			/**<number of chains */
  itr_POS,			/**<number of iterations */
  bun_POS,			/**<number of burn-in */
  thn_POS,			/**<number of thinning */
  kp_POS,			/**<number of simulations to keep in each chain */
  sims_POS,			/**<total number of simulations */
  rpt_POS,			/**<report frequency */
  tnit_POS,			/**<number of iterations used for tuning */
  ntn_POS,			/**<number of tuning loops */
  nmh_POS                       /**<number of M-H proposal variances to be tuned **/
};


/************************************************/
/*   Level 0 utility functions for printing    */  
/************************************************/

/**
 * utility function to print out acceptance rates
 *
 * @param n number of iterations
 * @param p the length of acc
 * @param acc a vector that stores acceptance times or percentages
 * @param pct indicating whether acc is the acceptance percentage or the unscaled acceptance times
 *
 */
static R_INLINE void print_acc(int n, int p, double *acc, int pct){
  double C = (pct) ? 100 : (100.0/n);
  Rprintf(_("Acceptance rate: min(%4.2f%%), mean(%4.2f%%), max(%4.2f%%)\n"),
	  dmin(acc, p) * C, mean(acc, p) * C, dmax(acc, p) * C);     
}

/**
 * utility function to print out division lines 
 */
static R_INLINE void print_line(){
  Rprintf("-----------------------------------------\n");
}


/************************************************/
/*   Level 1 utility functions for computing    */  
/*        full conditionals in bcplm            */
/************************************************/

/**
 * update the eta and mu slots in cpglm according to x
 *
 * @param x pointer to the vector of values for beta
 * @param da pointer to an SEXP struct
 *
 */

static void cpglm_fitted(double *x, SEXP da){
  int *dm = DIMS_SLOT(da) ;
  int nO = dm[nO_POS], nB = dm[nB_POS];
  double *X = X_SLOT(da),
    *beta = FIXEF_SLOT(da), *eta = ETA_SLOT(da),
    *mu = MU_SLOT(da), *offset = OFFSET_SLOT(da), 
    lp = LKP_SLOT(da)[0] ;
  if (x)  beta = x ;   /* point beta to x if x is not NULL */
  /* eta = X %*% beta */
  mult_mv("N", nO, nB, X, beta, eta) ;
  for (int i = 0; i < nO; i++){
    eta[i] += offset[i] ;
    mu[i] = link_inv(eta[i], lp);
  }
}

/**
 * Update the Xb, Zu, eta and mu slots in cpglmm according to x
 *
 * @param x pointer to the vector of values for beta or u
 * @param is_beta indicates whether x contains the values for beta or u. 
 *   1: x contains values of beta, and Zu is not updated. If x is null, the fixef slot is used; 
 *   0: x contains values of u, and Xb is not updated. If x is null, the u slot is used;
 *  -1: x is ignored, and the fixef and u slots are used.
 * @param da an SEXP object
 *
 */
static void cpglmm_fitted(double *x, int is_beta, SEXP da){    
  int *dm = DIMS_SLOT(da) ;
  int nO = dm[nO_POS], nB = dm[nB_POS], nU = dm[nU_POS];
  double  *X = X_SLOT(da), *eta = ETA_SLOT(da), 
    *mu = MU_SLOT(da), *beta = FIXEF_SLOT(da),
    *u = U_SLOT(da), *offset= OFFSET_SLOT(da), 
    *Xb = XB_SLOT(da), *Zu = ZU_SLOT(da), 
    lp = LKP_SLOT(da)[0], one[] = {1, 0}, zero[] = {0, 0};

  if (is_beta == -1) x = NULL ;            /* update from the fixef and u slots */

  /* update Xb */
  if (is_beta == 1 || is_beta == -1){      /* beta is updated */
    if (x)  beta = x ;                     /* point beta to x if x is not NULL */
    mult_mv("N", nO, nB, X, beta, Xb) ;    /* Xb = x * beta */
  }
  /* update Zu */
  if (is_beta == 0 || is_beta == -1){      /* u is updated */ 
    SEXP tmp;                              /* create an SEXP object to be coerced to CHM_DN */
    PROTECT(tmp = allocVector(REALSXP, nU));
    if (x) Memcpy(REAL(tmp), x, nU);      
    else Memcpy(REAL(tmp), u, nU);
    CHM_DN ceta, us = AS_CHM_DN(tmp);
    CHM_SP Zt = Zt_SLOT(da);
    R_CheckStack();
    ceta = N_AS_CHM_DN(Zu, nO, 1);          /* update Zu */
    R_CheckStack();                         /* Y = alpha * A * X + beta * Y */
    if (!M_cholmod_sdmult(Zt, 1 , one, zero, us, ceta, &c))
      error(_("cholmod_sdmult error returned"));
    UNPROTECT(1) ;
  }
  for (int i = 0; i < nO; i++){             /* update mu */
    eta[i] = Xb[i] + Zu[i] + offset[i];     /* eta = Xb + Z * u  + offset*/
    mu[i] = link_inv(eta[i], lp);
  }
}

/**
 * Update the Xb, Zu, eta and mu slots in bcplm using FIXEF and U 
 *
 * @param da a list object
 *
 */
static void update_mu(SEXP da){
  if (DIMS_SLOT(da)[nU_POS]) 
    cpglmm_fitted((double *) NULL, -1, da);
  else 
    cpglm_fitted((double *) NULL,  da);
}


/**
 * the data loglikelihood related to mu, assuming all but 
 * the mean parameters are constants
 *
 * @param da an SEXP struct
 *
 * @return the loglikelihood related to the mean
 */

static double llik_mu(SEXP da){  
  int *dm = DIMS_SLOT(da), *ygt0 = YPO_SLOT(da), k = 0;
  double *Y = Y_SLOT(da), *mu = MU_SLOT(da), *pwt = PWT_SLOT(da), 
    p = P_SLOT(da)[0], phi = PHI_SLOT(da)[0] ;
  double ld = 0.0, p2 = 2.0 - p, p1 = p - 1.0 ;

  for (int i = 0; i < dm[nO_POS]; i++)
    ld += pow(mu[i], p2) * pwt[i];
  ld /= - phi * p2 ;
  for (int i = 0; i < dm[nP_POS]; i++){
    k = ygt0[i] ;
    ld += - Y[k] * pow(mu[k], -p1) * pwt[k] / (phi * p1);
  }
  return ld ;
}

/**
 * compute the log prior distribution for a specified group  
 *
 * @param gn the group number
 * @param da an SEXP struct
 *
 * @return the log prior contribution of the specified group 
 */

static double prior_u_Gp(int gn, SEXP da){
  SEXP V = GET_SLOT(da, install("Sigma")); 
  int *Gp = Gp_SLOT(da), *nc = NCOL_SLOT(da), *nlev = NLEV_SLOT(da) ;
  double *v = REAL(VECTOR_ELT(V, gn)), *u = U_SLOT(da), ans = 0.0; 
  if (nc[gn] == 1) {                 /* univariate normal */
    for (int j = 0; j < nlev[gn]; j++)
      ans += -0.5 * u[Gp[gn] + j] * u[Gp[gn] + j] / v[0];
    return ans ;
  }
  else {                             /* multivariate normal */
    double *xv = Alloca(nc[gn], double), 
      *iv = Alloca(nc[gn] * nc[gn], double);
    R_CheckStack() ;    
    solve_po(nc[gn], v, iv) ;
    for (int j = 0; j < nlev[gn]; j++){
      for (int i = 0; i < nc[gn]; i++)
	xv[i] =  u[Gp[gn] + i * nlev[gn] + j] ;  
      ans += dmvnorm(nc[gn], xv, (double*) NULL, iv) ;
    }
    return ans ;
  }
}

/**
 * compute the contribution of the prior specification for u_k
 *
 * @param x vector of values for u_k
 * @param da an SEXP struct
 *
 * @return the loglikelihood related to the mean
 */

static double prior_uk(double x, SEXP da){
  int *dm = DIMS_SLOT(da), *Gp = Gp_SLOT(da), k = K_SLOT(da)[0];
  int gn = Gp_grp(k, dm[nT_POS], Gp); /* group number of u_k */
  double *u = U_SLOT(da), tmp = U_SLOT(da)[k], ans; 
  u[k] = x ;
  ans = prior_u_Gp(gn, da);           /* compute log prior for group gn */              
  u[k] = tmp ;
  return ans;
}

/************************************************/
/*   Level 2 utility functions to compute full  */
/*      conditionals used in the simulation     */  
/************************************************/

/**
 * the log posterior density of the index parameter p
 * (this is the same in both bcpglm and bcpglmm)
 *
 * @param x the value of p at which the log density is to be calculated
 * @param data a void struct, cocerced to SEXP internally
 *
 * @return log posterior density for p
 */
static double post_p(double x, void *data){
  SEXP da = data;
  double *Y = Y_SLOT(da), *mu = MU_SLOT(da), phi = PHI_SLOT(da)[0],
    *pwt = PWT_SLOT(da);
  return -0.5 * dl2tweedie(DIMS_SLOT(da)[nO_POS], Y, mu, phi, x, pwt) ;
}


/**
 * the log posterior density of the dispersion parameter phi
 * (this is the same in both bcpglm and bcpglmm)
 *
 * @param x the value of phi at which the log density is to be calculated
 * @param data a void struct, cocerced to SEXP internally
 *
 * @return log posterior density for phi
 */
static double post_phi(double x, void *data){
  SEXP da = data ;
  double *Y = Y_SLOT(da), *mu = MU_SLOT(da), p = P_SLOT(da)[0],
    *pwt = PWT_SLOT(da) ;
  return -0.5 * dl2tweedie(DIMS_SLOT(da)[nO_POS], Y, mu, x, p, pwt)  ;
}

/**
 * the log posterior density of beta_k, assuming Zu has been updated in cpglmm
 *
 * @param x the value of beta_k
 * @param data a void struct that is coerced to SEXP
 *
 * @return log posterior density for beta_k
 *
 */
static double post_betak(double x,  void *data){
  SEXP da = data ;
  int k = K_SLOT(da)[0], nU = DIMS_SLOT(da)[nU_POS]; 
  double pm = PBM_SLOT(da)[k], pv = PBV_SLOT(da)[k], 
    *l = CLLIK_SLOT(da), *beta = FIXEF_SLOT(da);
  double tmp = beta[k];
  beta[k] = x ;       /* update beta to compute mu */
  if (nU) cpglmm_fitted(beta, 1, da);
  else cpglm_fitted(beta, da) ;
  beta[k] = tmp ;     /* restore old beta values */
  *l = llik_mu(da) ;  /* this is stored and reused to speed up the simulation */
  return  *l - 0.5 * (x - pm) * (x - pm) / pv ;
}

/**
 * the log posterior density of u_k, assuming Xb has been updated
 *
 * @param x vector of values for u_k
 * @param data a void struct that is coerced to SEXP
 *
 * @return the log posterior density for u_k
 * 
 */
static double post_uk(double x,  void *data){
  SEXP da = data ;
  int k = K_SLOT(da)[0];
  double *u = U_SLOT(da),  *l = CLLIK_SLOT(da), tmp = U_SLOT(da)[k];    
  u[k] = x ;          /* update u to compute mu */
  cpglmm_fitted(u, 0, da) ;
  u[k] = tmp ;        /* restore old u values */
  *l = llik_mu(da);  
  return *l + prior_uk(x, da);
}


/************************************************/
/*   Level 3 utility functions for simulating   */  
/*          the parameters in bcplm             */
/************************************************/

/**
 * Simulate beta using the naive Gibbs update
 *
 * @param da an SEXP struct
 *
 */
static void sim_beta(SEXP da){
  int *dm = DIMS_SLOT(da), *k = K_SLOT(da);
  int nB = dm[nB_POS];
  double *beta = FIXEF_SLOT(da), *mh_sd = MHSD_SLOT(da), *l = CLLIK_SLOT(da), 
    *pm = PBM_SLOT(da), *pv = PBV_SLOT(da), *acc = ACC_SLOT(da);
  double xo, xn, l1, l2, A;

  /* initialize llik_mu*/
  *l = llik_mu(da);
  for (int j = 0; j < nB; j++){
    *k = j;
    xo = beta[j];
    xn = rnorm(xo, mh_sd[j]);
    l1 = *l;
    l2 = post_betak(xn, da);
    A = exp(l2 - l1 + 0.5 * (xo - pm[j]) * (xo - pm[j]) / pv[j]);
    /* determine whether to accept the sample */
    if (A < 1 && runif(0, 1) >= A){ /* not accepted */
      *l = l1;       /* revert the likelihood (this is updated in post_betak) */
    }
    else {
      beta[j] = xn;
      acc[j]++;    
    }
  }                  /* update the mean using the new beta */                    
  if (dm[nU_POS]) cpglmm_fitted(beta, 1, da);
  else cpglm_fitted(beta, da);  
}

/**
 * Simulate phi and p using the univariate RW M-H update.
 *
 * @param da a bcplm_input object
 *
 */

static void sim_phi_p(SEXP da){
  int *dm = DIMS_SLOT(da);
  int nB = dm[nB_POS]; 
  double *p = P_SLOT(da), *phi = PHI_SLOT(da), 
    *mh_sd = MHSD_SLOT(da), *acc = ACC_SLOT(da), 
    xn = 0.0;
  /* update phi */
  acc[nB] += metrop_tnorm_rw(*phi, mh_sd[nB], 0, BDPHI_SLOT(da)[0], 
			&xn, post_phi, (void *) da);
  *phi = xn ;
  /* update p */
  acc[nB + 1] += metrop_tnorm_rw(*p, mh_sd[nB + 1], BDP_SLOT(da)[0], 
			BDP_SLOT(da)[1], &xn, post_p, (void *) da);	
  *p = xn ;
}


/**
 * Simulate u using the naive Gibbs method. Block M-H update
 * generally does not work well and hence is not implemented.
 *
 * @param da a bcplm_input object
 *
 */

static void sim_u(SEXP da){
  int *dm = DIMS_SLOT(da), *k = K_SLOT(da);
  int nB = dm[nB_POS], nU = dm[nU_POS];
  double *u = U_SLOT(da), *l = CLLIK_SLOT(da), 
    *mh_sd = MHSD_SLOT(da) + nB + 2, /* shift the proposal variance pointer */
    *acc = ACC_SLOT(da) + nB + 2;    /* shift the acc pointer */
  double xo, xn, l1, l2, A;

  /* initialize llik_mu*/
  *l = llik_mu(da);
  for (int j = 0; j < nU; j++){
    *k = j ;
    xo = u[j];
    xn = rnorm(xo, mh_sd[j]);
    l1 = *l;
    l2 = post_uk(xn, da);
    A = exp(l2 - (l1 + prior_uk(xo, da)));  
    /* determine whether to accept the sample */
    if (A < 1 && runif(0, 1) >= A){ 
      *l = l1;  /* revert llik_mu (this is updated in post_uk) */
    }
    else{
      u[j] = xn;
      acc[j]++;    
    }
  }
  cpglmm_fitted(u, 0, da) ;  /* update the mean using the new u */
}

/**
 * Simulate the variance components. If there is only one random effect 
 * per level, the posterior is the inverse Gamma. Otherwise, the posterior 
 * is the inverse Wishart.  
 *
 * @param da a bcplm_input object
 *
 */

static void sim_Sigma(SEXP da){
  SEXP V = GET_SLOT(da, install("Sigma")) ;
  int *dm = DIMS_SLOT(da), *Gp = Gp_SLOT(da),  
    *nc = NCOL_SLOT(da), *nlev = NLEV_SLOT(da); 
  int nT = dm[nT_POS], mc = imax(nc, nT);
  double *v, su, *u = U_SLOT(da), 
    *scl = Alloca(mc * mc, double);
  R_CheckStack();

  for (int i = 0; i < nT; i++){
    v = REAL(VECTOR_ELT(V, i));
    if (nc[i] == 1){         /* simulate from the inverse-Gamma */
      su = sqr_length(u + Gp[i], nlev[i]);                    
      v[0] = 1/rgamma(0.5 * nlev[i] + IG_SHAPE, 1.0/(su * 0.5 + IG_SCALE));      
    }
    else {                   /* simulate from the inverse-Wishart */
      mult_xtx(nlev[i], nc[i], u + Gp[i], scl);            /* t(x) * (x) */
      for (int j = 0; j < nc[i]; j++) scl[j * j] += 1.0;   /* add prior (identity) scale matrix  */
      solve_po(nc[i], scl, v);
      rwishart(nc[i], (double) (nlev[i] + nc[i]), v, scl);
      solve_po(nc[i], scl, v);                  
    }
  }
}

/************************************************/
/* Additional utility functions for implementing*/  
/*        the MCMC sampler in bcplm             */
/************************************************/

/**
 * Set parameters to the k_th initial values provided in the inits slot
 *
 * @param da a bcplm_input object
 * @param k indicates the k_th set of initial values
 *
 * \note the order of the parameters in inits is: beta, phi, p, u and Sigma
 */

static void get_init(SEXP da, int k){
  
  SEXP inits = GET_SLOT(da, install("inits"));  /* inits is a list */
  int *dm = DIMS_SLOT(da);
  int nB = dm[nB_POS], nU = dm[nU_POS], nT = dm[nT_POS];
  double *init = REAL(VECTOR_ELT(inits, k));

  /* set beta, phi, p*/
  Memcpy(FIXEF_SLOT(da), init, nB) ;
  *PHI_SLOT(da) = init[nB] ;
  *P_SLOT(da) = init[nB + 1] ;

  /* set U and Sigma */
  if (nU) {
    SEXP V = GET_SLOT(da, install("Sigma"));
    int pos = 0, st = nB + 2, *nc = NCOL_SLOT(da) ;
    double *v ;
    Memcpy(U_SLOT(da), init + st, nU);
    /* set Sigma */
    for (int i = 0; i < nT; i++){
      v = REAL(VECTOR_ELT(V, i)) ;
      Memcpy(v, init + st + nU + pos, nc[i] * nc[i]) ;
      pos += nc[i] * nc[i] ;
    }    
  }
}

/**
 * Save parameters to the ns_th row of the simulation results
 *
 * @param da a bcplm_input object 
 * @param ns indicates the ns_th row
 * @param nS number of rows in sims 
 * @param sims a long vector to store simulations results 
 *
 */
static void set_sims(SEXP da, int ns, int nS, double *sims){

  int *dm = DIMS_SLOT(da) ;
  int pos = 0, nB = dm[nB_POS], nU = dm[nU_POS],
    nT = dm[nT_POS];
  double *beta = FIXEF_SLOT(da), *u = U_SLOT(da);

  for (int j = 0; j < nB; j++)  
    sims[j * nS + ns] = beta[j] ;
  sims[nB * nS + ns] = *PHI_SLOT(da) ;
  sims[(nB + 1) * nS + ns] = *P_SLOT(da) ;

  /* set U and Sigma */
  if (nU) {
    SEXP V = GET_SLOT(da, install("Sigma"));
    int *nc = NCOL_SLOT(da), st = nB + 2;
    double *v; 
    for (int j = 0; j < nU; j++)
      sims[(j + st) * nS + ns] = u[j] ; 
    for (int i = 0; i < nT; i++){
      v = REAL(VECTOR_ELT(V, i));
      for (int j = 0; j < nc[i] * nc[i]; j++)
	sims[(st + nU + pos + j) * nS + ns] = v[j] ;
      pos += nc[i] * nc[i] ;
    }
  } 
}


/**
 * update the proposal standard deviations 
 *
 * @param p the number of parameters to be tuned
 * @param acc pointer to the acceptance rate in the M-H update
 * @param mh_sd pointer to the vector of proposal standard deviations
 * @param mark pointer to the vector of marks
 *
 */
static R_INLINE void tune_var(int p, double *acc, double *mh_sd, int *mark) {
  double acc_bd;
  for (int j = 0; j < p; j++){
    acc_bd = fmin2(fmax2(acc[j], 0.01), 0.99); /* bound the empirical acceptance */
    if (acc[j] < (MH_ACC_TAR - MH_ACC_TOL))
      mh_sd[j] /= 2 - acc_bd/MH_ACC_TAR;
    else if (acc[j] > (MH_ACC_TAR - MH_ACC_TOL))	    
      mh_sd[j] *= 2 - (1 - acc_bd)/(1 - MH_ACC_TAR);      
    else mark[j]++;      
  }  
}



/************************************************/
/*         MCMC simulations for cplm            */
/************************************************/

/**
 * main workhorse for the MCMC simulation
 *
 * @param da a bcplm_input object
 * @param nit number of iterations
 * @param nbn number of burn-in
 * @param nth thinning rate
 * @param nS the number of rows of sims
 * @param nR report interval
 * @param sims a 2d array to store simulation results
 *
 * \note It's assumed that eta and mu are already updated. Also,
 * nS is passed as an argument rather than determined from da. 
 * This is because the tuning process uses a different number of 
 * simulations than the n.keep element in the dims slot. 
 *
 */

static void do_mcmc(SEXP da, int nit, int nbn, int nth, int nS, 
		    int nR, double *sims){

  int *dm = DIMS_SLOT(da);
  int nU = dm[nU_POS], nmh = dm[nmh_POS],
    ns = 0, do_print = 0;
  /* initialize acc */
  double *acc = ACC_SLOT(da);
  AZERO(acc, nmh);

  /* run MCMC simulatons */
  GetRNGstate();
  for (int iter = 0; iter < nit; iter++){
    do_print = (nR > 0 && (iter + 1) % nR == 0);
    if (do_print) Rprintf(_("Iteration: %d \n "), iter + 1);

    /* update parameters */
    sim_beta(da);
    sim_phi_p(da);
    if (nU){
      sim_u(da);
      sim_Sigma(da);
    }
        
    /* store results  */
    if (iter >= nbn &&  (iter + 1 - nbn) % nth == 0 ){
      ns = (iter + 1 - nbn) / nth - 1;
      set_sims(da, ns, nS, sims);
    } 

    /* print out acceptance rate if necessary */
    if (do_print) print_acc(iter + 1, nmh, acc, 0);
    R_CheckUserInterrupt();
  }
  PutRNGstate();
  /* compute acceptance percentage */
  for (int i = 0; i < nmh; i++) acc[i] /= nit ;
}


/**
 * tune the proposal variances 
 *
 * @param da an input list object
 *
 */
static void tune_mcmc(SEXP da){
  int *dm = DIMS_SLOT(da);
  int nmh = dm[nmh_POS],
    etn = ceil(dm[tnit_POS] * 1.0 / dm[ntn_POS]) ;  // # iters per tuning loop;
  double *mh_sd = MHSD_SLOT(da), *acc = ACC_SLOT(da), 
    *sims = Calloc(etn * dm[nA_POS], double);
  int nmark = 0, *mark = Calloc(nmh, int);
  AZERO(mark, nmh);

  /* run MCMC and tune parameters */
  if (dm[rpt_POS]) Rprintf(_("Tuning phase...\n"));
  for (int i = 0; i < dm[ntn_POS]; i++) {
    do_mcmc(da, etn, 0, 1, etn, 0, sims);      /* run mcmc */    
    tune_var(nmh, acc, mh_sd, mark);           /* adjust proposal sd's */
    /* determine whether the parameters are fully tuned */
    nmark = 0;
    for (int j = 0; j < nmh; j++)
      if (mark[j] >= 3) nmark++;
    if (nmark == nmh) break;
  }
  if (dm[rpt_POS]){
    print_acc(1, nmh, acc, 1);
    print_line();
  }
  Free(sims);
  Free(mark);
}

/**
 * implement MCMC for compound Poisson linear models
 *
 * @param da an bcplm_input object
 *
 * @return the simulated values
 *
 */

SEXP bcplm_mcmc (SEXP da){
  /* get dimensions */
  int *dm = DIMS_SLOT(da) ;
  int nR = dm[rpt_POS];
  SEXP ans, ans_tmp;

  /* tune the scale parameter for M-H update */
  if (dm[tnit_POS]) {
    update_mu(da);                    /* update eta mu*/
    tune_mcmc(da);
  }
    
  /* run Markov chains */
  PROTECT(ans = allocVector(VECSXP, dm[chn_POS])) ;

  for (int k = 0; k < dm[chn_POS]; k++){
    if (nR) Rprintf(_("Start Markov chain %d\n"), k + 1);   
    get_init(da, k) ;                /* initialize the chain */
    update_mu(da) ;                  /* update eta and mu */ 
  
    /* run MCMC and store result */
    PROTECT(ans_tmp = allocMatrix(REALSXP, dm[kp_POS], dm[nA_POS]));
    do_mcmc(da, dm[itr_POS], dm[bun_POS], dm[thn_POS], 
	    dm[kp_POS], nR, REAL(ans_tmp));
    SET_VECTOR_ELT(ans, k, ans_tmp);
    UNPROTECT(1) ;
    if (nR) print_line();
  }
  UNPROTECT(1) ;
  if (nR)  Rprintf(_("Markov Chain Monte Carlo ends!\n"));
  return ans ;
    
}


/************************************************/
/*      Export R Callable functions             */
/************************************************/


/**
 * R callable function to update the mean in a bcplm object
 *
 * @param da  an bcplm_input object
 *
 */
SEXP bcplm_update_mu(SEXP da){
  update_mu(da);
  return R_NilValue;
}


/**
 * R callable function to compute the conditional posterior of p
 *
 * @param x value at which to evaluate the density
 * @param da  an bcplm_input object
 *
 * @return the log of the conditional posterior of p 
 */
SEXP bcplm_post_p(SEXP x, SEXP da){
  return ScalarReal(post_p(REAL(x)[0], (void *) da));
}


/**
 * R callable function to compute the conditional posterior of phi
 *
 * @param x value at which to evaluate the density
 * @param da  an bcplm_input object
 *
 * @return the log of the conditional posterior of phi
 */
SEXP bcplm_post_phi(SEXP x, SEXP da){
  return ScalarReal(post_phi(REAL(x)[0], (void *) da));
}


/**
 * R callable function to compute the conditional posterior of betak
 *
 * @param x value at which to evaluate the density
 * @param da  an bcplm_input object
 *
 * @return the log of the conditional posterior of betak
 */
SEXP bcplm_post_betak(SEXP x, SEXP da){
  return ScalarReal(post_betak(REAL(x)[0], (void *) da));
}


/**
 * R callable function to compute the conditional posterior of uk
 *
 * @param x value at which to evaluate the density
 * @param da  an bcplm_input object
 *
 * @return the log of the conditional posterior of uk
 */
SEXP bcplm_post_uk(SEXP x, SEXP da){
  return ScalarReal(post_uk(REAL(x)[0], (void *) da));
}
