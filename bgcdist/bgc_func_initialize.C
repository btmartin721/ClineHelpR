#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include <float.h>

#include <math.h>

#include "bgc.h"

using namespace std;

/* -------------------------------------------------------------- */
/*           parameter initialization functions                   */
/* -------------------------------------------------------------- */

/* main function for parameter initialization */
void initParams(paramcont * param, datacont * data, datadimcont * datadim,
		auxcont * auxvar, int ngs, int linkMod){
  int i, k, n;
  int nk;
  double a, b;

  // for ngs, p0cnt and p1cnt are the unobserved genotypes, which must be initialized
  if(ngs == 1){
    initPcnts(data,datadim,auxvar->vecAllelen1,auxvar->vecAllele_int);
  }

  // init pi, parental allele frequencies by sampling directly from posterior
  for(i=0; i<datadim->nloci; i++){
    nk = gsl_vector_int_get(datadim->nallele, i);
    // pi0
    for(k=0; k<nk; k++){
      gsl_vector_set(auxvar->vecAllelen1, k, gsl_matrix_int_get(data->p0cnt, i, k) + 1); // 1 added for dirichlet prior
    }
    gsl_ran_dirichlet(r, nk, auxvar->vecAllelen1->data, auxvar->vecAllelen2->data);
    for(k=0; k<nk; k++){
      gsl_matrix_set(param->pi0Cur, i, k, gsl_vector_get(auxvar->vecAllelen2, k));
    }
    // pi1
    for(k=0; k<nk; k++){
      gsl_vector_set(auxvar->vecAllelen1, k, gsl_matrix_int_get(data->p1cnt, i, k) + 1); // 1 added for dirichlet prior
    }
    gsl_ran_dirichlet(r, nk, auxvar->vecAllelen1->data, auxvar->vecAllelen2->data);
    for(k=0; k<nk; k++){
      gsl_matrix_set(param->pi1Cur, i, k, gsl_vector_get(auxvar->vecAllelen2, k));
    }
  }

  // init z, ancestry of admixed individuals at random
  for(i=0; i<datadim->nloci; i++){
    for(n=0; n<datadim->nind; n++){
      for(k=0; k<PLOIDY; k++){
	gsl_matrix_int_set(param->zCur[i].mat, n, k, gsl_ran_bernoulli(r, 0.5)); // 0.5 gives equal prob. of ancestry
      }
    }
  }

  // init h, hybrid index using uniform
  for(n=0; n<datadim->nind; n++){
    gsl_matrix_set(param->hi, n, 0, gsl_ran_flat(r,0,1));
  }

  // init fis, set to zero
  for(n=0; n<datadim->nind; n++){
    gsl_vector_set(param->fis, n, 0);
  }
  // init gamma, locus effect for cline center (alpha)
  initRandomEffect(param->gamma, datadim->nloci, 1, INITSD);
  
  // init eta, nested population effect for cline center (alpha)
  if (datadim->npop > 1){
    initRandomEffect(param->etaCur, datadim->nloci, datadim->npop, INITSDP);
  }
  
  // init zeta, locus effect for cline width (beta)
  initRandomEffect(param->zeta, datadim->nloci, 1, INITSD);

  // init kappa, nested population effect for cline width (beta)
  if (datadim->npop > 1){
    initRandomEffect(param->kappaCur, datadim->nloci, datadim->npop, INITSDP);
  }

  // init tau_alpha, precision for gamma random effects, uses expected precision
  a = datadim->nloci / 2;
  b = 0.5 * gsl_pow_2(INITSD);
  param->tauAlphaCur = gsl_ran_gamma(r,a,b);
  // init tau_beta, precision for zeta random effects, uses expected precision
  param->tauBetaCur = gsl_ran_gamma(r,a,b);

  // init nu and omega, precision for eta and kappa respectively, uses expected precision
  a = datadim->npop / 2;
  b = 0.5 * gsl_pow_2(INITSDP);
  for(i=0; i<datadim->nloci; i++){
    gsl_matrix_set(param->nu, i, 0, gsl_ran_gamma(r,a,b));
    gsl_matrix_set(param->omega, i, 0, gsl_ran_gamma(r,a,b));    
  }
  
  // init rho for linkage model
  if (linkMod == 1){
    param->rhoCur = gsl_ran_flat(r,0,1);
  }
  // init alphas, based on gamma and eta
  initClineParameter(param->alphaCur,param->gamma,param->etaCur,datadim->nloci,datadim->npop);

  // init betas, based on zeta and kappa
  initClineParameter(param->betaCur,param->zeta,param->kappaCur,datadim->nloci,datadim->npop);
}
/* alternative function for parameter initialization that uses information from data */
void initParams2(paramcont * param, datacont * data, datadimcont * datadim,
		auxcont * auxvar, int ngs, int linkMod){
  int i, k, n;
  int nk;
  double a, b;
  double hctr, tctr;

  // for ngs, p0cnt and p1cnt are the unobserved genotypes, which must be initialized
  if(ngs == 1){
    initPcnts(data,datadim,auxvar->vecAllelen1,auxvar->vecAllele_int);
  }

  // init pi, parental allele frequencies by sampling directly from posterior
  for(i=0; i<datadim->nloci; i++){
    nk = gsl_vector_int_get(datadim->nallele, i);
    // pi0
    for(k=0; k<nk; k++){
      gsl_vector_set(auxvar->vecAllelen1, k, gsl_matrix_int_get(data->p0cnt, i, k) + 1); // 1 added for dirichlet prior
    }
    gsl_ran_dirichlet(r, nk, auxvar->vecAllelen1->data, auxvar->vecAllelen2->data);
    for(k=0; k<nk; k++){
      gsl_matrix_set(param->pi0Cur, i, k, gsl_vector_get(auxvar->vecAllelen2, k));
    }
    // pi1
    for(k=0; k<nk; k++){
      gsl_vector_set(auxvar->vecAllelen1, k, gsl_matrix_int_get(data->p1cnt, i, k) + 1); // 1 added for dirichlet prior
    }
    gsl_ran_dirichlet(r, nk, auxvar->vecAllelen1->data, auxvar->vecAllelen2->data);
    for(k=0; k<nk; k++){
      gsl_matrix_set(param->pi1Cur, i, k, gsl_vector_get(auxvar->vecAllelen2, k));
    }
  }

  // init z, ancestry of admixed individuals based on data and parental allele frequencies 
  if (ngs == 1){
    initZngs(param->zCur, data->phcnt, param->pi0Cur, param->pi1Cur, datadim->nallele, 
	     datadim->nloci, datadim->nind);
  }
  else {
    cerr << "This initialization scheme only works with ngs data" << endl;
    exit(1);
  }

  // init h, hybrid index from beta based on estimated ancestry
  for (n=0; n<datadim->nind; n++){
    hctr = 1;
    tctr = 1;
    for (i=0; i<datadim->nloci; i++){
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	hctr++;
      }
      else {
	tctr++;
      }
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 1) == 1){
	hctr++;
      }
      else {
	tctr++;
      }
    }
    // cerr << hctr << " : " << tctr << endl;
    gsl_matrix_set(param->hi, n, 0, gsl_ran_beta(r, hctr, tctr));
  }
 
  // init fis, set to zero
  for(n=0; n<datadim->nind; n++){
    gsl_vector_set(param->fis, n, 0);
  }

  // init gamma, locus effect for cline center (alpha)
  initRandomEffect(param->gamma, datadim->nloci, 1, INITSD);
  
  // init eta, nested population effect for cline center (alpha)
  if (datadim->npop > 1){
    initRandomEffect(param->etaCur, datadim->nloci, datadim->npop, INITSDP);
  }
  
  // init zeta, locus effect for cline width (beta)
  initRandomEffect(param->zeta, datadim->nloci, 1, INITSD);

  // init kappa, nested population effect for cline width (beta)
  if (datadim->npop > 1){
    initRandomEffect(param->kappaCur, datadim->nloci, datadim->npop, INITSDP);
  }

  // init tau_alpha, precision for gamma random effects, uses expected precision
  a = datadim->nloci / 2;
  b = 0.5 * gsl_pow_2(INITSD);
  param->tauAlphaCur = gsl_ran_gamma(r,a,b);
  // init tau_beta, precision for zeta random effects, uses expected precision
  param->tauBetaCur = gsl_ran_gamma(r,a,b);

  // init nu and omega, precision for eta and kappa respectively, uses expected precision
  a = datadim->npop / 2;
  b = 0.5 * gsl_pow_2(INITSDP);
  for(i=0; i<datadim->nloci; i++){
    gsl_matrix_set(param->nu, i, 0, gsl_ran_gamma(r,a,b));
    gsl_matrix_set(param->omega, i, 0, gsl_ran_gamma(r,a,b));    
  }
  
  // init rho for linkage model
  if (linkMod == 1){
    param->rhoCur = gsl_ran_flat(r,0,1);
  }
  // init alphas, based on gamma and eta
  initClineParameter(param->alphaCur,param->gamma,param->etaCur,datadim->nloci,datadim->npop);

  // init betas, based on zeta and kappa
  initClineParameter(param->betaCur,param->zeta,param->kappaCur,datadim->nloci,datadim->npop);
}

/* general function for initializing random effect parameters with mean = 0 */
void initRandomEffect(gsl_matrix * reParam, int dim1, int dim2, double sd){
  int i, j;
  
  for(i=0; i<dim1; i++){
    for(j=0; j<dim2; j++){
      gsl_matrix_set(reParam, i, j, gsl_ran_gaussian(r, sd));
    }
  }
}

/* initialize alpha or beta based on random effects */
void initClineParameter(gsl_matrix * clineParam, gsl_matrix * locusEffect, gsl_matrix * popEffect, 
			int nloci, int npop){
  int i, j;
  double le;
  
  for(i=0; i<nloci; i++){
    le = gsl_matrix_get(locusEffect, i, 0);
    for(j=0; j<npop; j++){
      gsl_matrix_set(clineParam, i, j, le + gsl_matrix_get(popEffect, i, j));
    }
  }
}


/* copies current parameter values to last parameter values */
void copyCur2Last(paramcont * param, datadimcont * datadim, gsl_vector * locivec, 
		  gsl_vector * indvec, int linkMod){
  int i;

  gsl_matrix_memcpy(param->pi0Last, param->pi0Cur);
  gsl_matrix_memcpy(param->pi1Last, param->pi1Cur);
  
  for(i=0; i<datadim->nloci; i++){
    gsl_matrix_int_memcpy(param->zLast[i].mat, param->zCur[i].mat);
  }
  
  gsl_matrix_get_col(indvec, param->hi, 0);
  gsl_matrix_set_col(param->hi, 1, indvec);
  
  gsl_matrix_get_col(locivec, param->gamma, 0);
  gsl_matrix_set_col(param->gamma, 1, locivec);
  
  gsl_matrix_memcpy(param->etaLast, param->etaCur);

  gsl_matrix_get_col(locivec, param->zeta, 0);
  gsl_matrix_set_col(param->zeta, 1, locivec);
 
  gsl_matrix_memcpy(param->kappaLast, param->kappaCur);

  param->tauAlphaLast = param->tauAlphaCur;
  param->tauBetaLast = param->tauBetaCur;

  gsl_matrix_get_col(locivec, param->nu, 0);
  gsl_matrix_set_col(param->nu, 1, locivec);

  gsl_matrix_get_col(locivec, param->omega, 0);
  gsl_matrix_set_col(param->omega, 1, locivec);

  gsl_matrix_memcpy(param->alphaLast, param->alphaCur);
  gsl_matrix_memcpy(param->betaLast, param->betaCur);

  if (linkMod == 1){
    param->rhoLast = param->rhoCur;
  }

}
