#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_poly.h>
#include <float.h>
#include <omp.h>

#include <math.h>

#include "bgc.h"
#include "mvrandist.h"

using namespace std;

/* -------------------------------------------------------------- */
/*                          mcmc functions                        */
/* -------------------------------------------------------------- */

/* control function for mcmc updates */
void mcmcUpdate(paramcont * param, datacont * data, datadimcont * datadim, auxcont * auxvar,
		int allDiploid,	int onePrecision, int sumZ, int linkMod, int ngs, 
		double trunc){

  int i, j;

  // update parental allele frequencies, pop 0 then pop 1
  if (ngs == 1){ // update p0cnt and p1cnt
    updatePcnt(data,auxvar,datadim,data->p0cnt,data->ngsp0cnt,datadim->nindP0,param->pi0Cur);
    updatePcnt(data,auxvar,datadim,data->p1cnt,data->ngsp1cnt,datadim->nindP1,param->pi1Cur);
  }

  updatePi(data->p0cnt,param->pi0Cur,datadim->nallele,datadim->nloci,auxvar->vecAllelen1,auxvar->vecAllelen2);
  updatePi(data->p1cnt,param->pi1Cur,datadim->nallele,datadim->nloci,auxvar->vecAllelen1,auxvar->vecAllelen2);

  // update z, ancestry, uses gibbs sampling
  if (ngs == 0){
    updateZ(param,data,datadim,allDiploid);
  }
  else if (ngs == 1){ // ngs update of ancestry
    updateZngs(param,data,datadim,auxvar);
  }

  // update hi, hybrid index, uses random walk metropolis
  updateHi(param,datadim,auxvar,allDiploid);

  /* CLEAN UP */

  // update fis, inbreeding coefficient for ancestry, uses random walk metropolis
  //   updateFis(h,zCur,alphaCur,betaCur,zero2one,zero2one1,zero2one2,zero2one3,M,pops,
  // 	    tuneFis,nind,nloci,ploidyVec,allDiploid,fis);
  
  // set fis = 0
  for (i=0; i<datadim->nind; i++){
    gsl_vector_set(param->fis, i, 0);
  }

  // update gamma, locus effect on genomic cline paramter alpha
  if (linkMod == 0){
    updateGamma(param,datadim,auxvar,allDiploid,sumZ);
  }
  else if (linkMod == 1){
    setSigma(param->sigmaAll,data->proxMat,param->rhoCur,datadim->nloci);
    updateGammaCAR(param,datadim,auxvar,allDiploid);
  }

  // update zeta, locus effect on genomic cline paramter beta
  if (linkMod == 0){
    updateZeta(param,datadim,auxvar,allDiploid,sumZ);
  }
  else if (linkMod == 1){
    updateZetaCAR(param,datadim,auxvar,allDiploid);
  }

  // update nested population effects for cline center and cline width
  if (datadim->npop == 1){ // set to zero if npop = 1
    for (i=0; i<datadim->nloci; i++){
      gsl_matrix_set(param->etaCur,i,0,0); // set to zero
      gsl_matrix_set(param->kappaCur,i,0,0); // set to zero
    }
  }
  else { // update if npop != 1, uses MVN
    for (i=0; i<datadim->nloci; i++){
      updateEta(param,datadim,auxvar,i,allDiploid,onePrecision);
      updateKappa(param,datadim,auxvar,i,allDiploid,onePrecision);
    }
  }

  // update precision parameters, tauAlpha, tauBeta
  if (linkMod == 0){
    updateTau(param->gamma,&param->tauAlphaCur,datadim->nloci,trunc);
    updateTau(param->zeta,&param->tauBetaCur,datadim->nloci,trunc);
  }
  else if (linkMod == 1){
    updateTauCAR(param->gamma,&param->tauAlphaCur,param->tauAlphaLast,param->sigmaAll,param->rhoCur,datadim->nloci,param->tuneTau);
    updateTauCAR(param->zeta,&param->tauBetaCur,param->tauBetaLast,param->sigmaAll,param->rhoCur,datadim->nloci,param->tuneTau);
    // update rho if linkage model is used
    updateRhoCAR(param,data,datadim);
  }

  // update precision parameters nu and omega
  if (datadim->npop == 1){ // no need to update
    ;
  }
  else { // update if npop != 1
    updateNuOmega(param->etaCur,param->nu,datadim->nloci,datadim->npop,onePrecision);
    updateNuOmega(param->kappaCur,param->omega,datadim->nloci,datadim->npop,onePrecision);
  }

  // calculate alpha and beta based on current mcmc samples
  for(i=0; i<datadim->nloci; i++){
    for(j=0; j<datadim->npop; j++){
      gsl_matrix_set(param->alphaCur,i,j, gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j));
      gsl_matrix_set(param->betaCur,i,j, gsl_matrix_get(param->zeta,i,0) + gsl_matrix_get(param->kappaCur,i,j));
    }
  }
  // calculate phi based on mcmc samples for this iteration
  calcPhi(param,datadim,auxvar);

}

/* function to update pi, parental allele frequencies, samples directly from posterior */
void updatePi(gsl_matrix_int * pcnt, gsl_matrix * piCur, gsl_vector_int * nallele,
	      int nloci, gsl_vector * alpha, gsl_vector * theta){
  int i, k, nk;

  for(i=0; i<nloci; i++){
    nk = gsl_vector_int_get(nallele, i);
    for(k=0; k<nk; k++){
      gsl_vector_set(alpha, k, gsl_matrix_int_get(pcnt, i, k) + 1); // 1 added for dirichlet prior
    }
    gsl_ran_dirichlet(r, nk, alpha->data, theta->data);
    for(k=0; k<nk; k++){
      gsl_matrix_set(piCur, i, k, gsl_vector_get(theta, k));
    }
  }

}

/* function to update z, ancestry of hybrids, uses gibbs sampling similar to structure */
void updateZ(paramcont * param, datacont * data, datadimcont * datadim, int allDiploid){
  int i, n, k, nk;
  double p00 = 0, p01 = 0, p10 = 0, p11 = 0, ptot;
  double pr[4]; // array of probabilities
  unsigned int sam[4]; // sampled ancestry for a locus
  int obs[2]; // obs alleles
  int acopy; // allele copy 0 or 1
  int pl; // ploidy

  // only have diploid data
  if (allDiploid == 1){
    // loop through loci
    for(i=0; i<datadim->nloci; i++){
      nk = gsl_vector_int_get(datadim->nallele, i);
      // loop through individuals
      for(n=0; n<datadim->nind; n++){
	acopy = 0;
	p00 = 0;
	p01 = 0;
	p11 = 0;
	// missing data, set probabilities based on phi
	if(gsl_matrix_int_get(data->phcnt[i].mat, n, 0) == -9){
	  p00 = gsl_pow_2(1 - gsl_matrix_get(param->phi, i, n)); // (1-phi)^2
	  p01 = (1 - gsl_matrix_get(param->phi, i, n)) * gsl_matrix_get(param->phi, i, n); // (1-phi)phi
	  p10 = p01; // both hets ancestry have equal prob.
	  p11 = gsl_pow_2(gsl_matrix_get(param->phi, i, n)); // phi^2	  
	}
	// loop through alleles
	else{ // individual has data for this locus
	  for(k=0; k<nk; k++){
	    if(gsl_matrix_int_get(data->phcnt[i].mat, n, k) == 1){
	      obs[acopy] = k;
	      acopy++;
	    }
	    else if (gsl_matrix_int_get(data->phcnt[i].mat, n, k) == 2){
	      obs[acopy] = k;
	      acopy++;
	      obs[acopy] = k;
	      acopy++;	      
	    }
	  }
	  // calculate p(z|x,pi,phi) = p(x|z,pi) p(z|phi)
	  if (acopy==2){
	    if (obs[0] == obs[1]){ // homozyogous genotype
	      // p00 = p_0^2 * (1-phi)^2
	      p00 = gsl_pow_2(gsl_matrix_get(param->pi0Cur, i, obs[0])) * 
		gsl_pow_2(1 - gsl_matrix_get(param->phi, i, n));
	      // p01 = p_0 * p_1 (1-phi)phi
	      p01 = gsl_matrix_get(param->pi0Cur, i, obs[0]) * gsl_matrix_get(param->pi1Cur, i, obs[0]) * 
		(1 - gsl_matrix_get(param->phi, i, n)) * gsl_matrix_get(param->phi, i, n);
	      p10 = p01; // probs of het. ancestry are equal
	      // p1 = p_1^2 *  phi^2
	      p11 = gsl_pow_2(gsl_matrix_get(param->pi1Cur, i, obs[0])) * 
		gsl_pow_2(gsl_matrix_get(param->phi, i, n));
	    }
	    else { // heterozygous genotype
	      // p00 = 2 p_0 q_0 * (1-phi)^2
	      p00 = 2 * gsl_matrix_get(param->pi0Cur, i, obs[0]) * gsl_matrix_get(param->pi0Cur, i, obs[1]) * 
		gsl_pow_2(1 - gsl_matrix_get(param->phi, i, n));
	      // p01 = ((p_0 q_1) + (p_1 q_0)) p_1 (1-phi)phi
	      p01 = (gsl_matrix_get(param->pi0Cur, i, obs[0]) * gsl_matrix_get(param->pi1Cur, i, obs[1]) +
		     gsl_matrix_get(param->pi1Cur, i, obs[0]) * gsl_matrix_get(param->pi0Cur, i, obs[1])) *
		(1 - gsl_matrix_get(param->phi, i, n)) * gsl_matrix_get(param->phi, i, n);
	      p10 = p01;
	      // p11 = 2 p_1 q_1 * phi^2
	      p11 = 2 * gsl_matrix_get(param->pi1Cur, i, obs[0]) * gsl_matrix_get(param->pi1Cur, i, obs[1]) * 
			   gsl_pow_2(gsl_matrix_get(param->phi, i, n));
	      
	    }
	  }
	  else {
	    cerr << "Indiviudal " << n << " is missing informaiton for locus " << i << endl;
	    exit(0);
	  }
	  ptot = p00 + p01 + p10 + p11;
	  //cerr << p00 << " " << p01 << " " << p11 << " : " << ptot << endl;
	  pr[0] = p00 / ptot;
	  pr[1] = p01 / ptot;
	  pr[2] = p10 / ptot;
	  pr[3] = p11 / ptot;
	  // sample ancestry pair based on probabilities
	  gsl_ran_multinomial(r, 4, 1, pr, sam);
	  if(sam[0] == 1){ // homo pop 0
	    gsl_matrix_int_set(param->zCur[i].mat, n, 0, 0);
	    gsl_matrix_int_set(param->zCur[i].mat, n, 1, 0);
	  }
	  else if(sam[1] == 1){ // het
	    gsl_matrix_int_set(param->zCur[i].mat, n, 0, 0);
	    gsl_matrix_int_set(param->zCur[i].mat, n, 1, 1);
	  }
	  else if(sam[2] == 1){ // het
	    gsl_matrix_int_set(param->zCur[i].mat, n, 0, 1);
	    gsl_matrix_int_set(param->zCur[i].mat, n, 1, 0);
	  }
	  else{ // homo pop 1
	    gsl_matrix_int_set(param->zCur[i].mat, n, 0, 1);
	    gsl_matrix_int_set(param->zCur[i].mat, n, 1, 1);
	  }
	}
      }
    }
  }
  else { // some loci could be haploid
    for(i=0; i<datadim->nloci; i++){
      nk = gsl_vector_int_get(datadim->nallele, i);
      pl = gsl_vector_int_get(datadim->ploidyVec, i);
      if (pl==2){ // diploid
	// loop through individuals
	for(n=0; n<datadim->nind; n++){
	  acopy = 0;
	  p00 = 0;
	  p01 = 0;
	  p11 = 0;
	  // missing data, set probabilities based on phi and fis
	  if(gsl_matrix_int_get(data->phcnt[i].mat, n, 0) == -9){
	    p00 = gsl_pow_2(1 - gsl_matrix_get(param->phi, i, n)) * (1 - gsl_vector_get(param->fis, n)) + 
	      (1 - gsl_matrix_get(param->phi, i, n)) * gsl_vector_get(param->fis, n); // (1-phi)^2 (1-f) + f(1-phi)
	    p01 = (1 - gsl_matrix_get(param->phi, i, n)) * gsl_matrix_get(param->phi, i, n) * 
	      (1 - gsl_vector_get(param->fis, n)); // (1-phi)phi (1-f)
	    p10 = p01; // both hets ancestry have equal prob.
	    p11 = gsl_pow_2(gsl_matrix_get(param->phi, i, n)) * (1 - gsl_vector_get(param->fis, n)) + 
	      gsl_matrix_get(param->phi, i, n) * gsl_vector_get(param->fis, n); // phi^2 (1-f) + f(phi)	  
	  }
	  // loop through alleles
	  else{ // individual has data for this locus
	    for(k=0; k<nk; k++){
	      if(gsl_matrix_int_get(data->phcnt[i].mat, n, k) == 1){
		obs[acopy] = k;
		acopy++;
	      }
	      else if (gsl_matrix_int_get(data->phcnt[i].mat, n, k) == 2){
		obs[acopy] = k;
		acopy++;
		obs[acopy] = k;
		acopy++;	      
	      }
	    }
	    // calculate p(z|x,pi,phi,f) = p(x|z,pi) p(z|phi,f)
	    if (acopy==2){
	      if (obs[0] == obs[1]){ // homozyogous genotype
		// p00 = p_0^2 * ((1-phi)^2 (1-f) + f(1-phi))
		p00 = gsl_pow_2(gsl_matrix_get(param->pi0Cur, i, obs[0])) * 
		  (gsl_pow_2(1 - gsl_matrix_get(param->phi, i, n)) * (1 - gsl_vector_get(param->fis, n)) + 
		   (1 - gsl_matrix_get(param->phi, i, n)) * gsl_vector_get(param->fis, n));
		// p01 = p_0 * p_1 (1-phi)phi (1-f)
		p01 = gsl_matrix_get(param->pi0Cur, i, obs[0]) * gsl_matrix_get(param->pi1Cur, i, obs[0]) * 
		  ((1 - gsl_matrix_get(param->phi, i, n)) * gsl_matrix_get(param->phi, i, n) * 
		   (1 - gsl_vector_get(param->fis, n)));
		p10 = p01; // probs of het. ancestry are equal
		// p1 = p_1^2 *  (phi^2 (1-f) + f(phi))
		p11 = gsl_pow_2(gsl_matrix_get(param->pi1Cur, i, obs[0])) * 
		  (gsl_pow_2(gsl_matrix_get(param->phi, i, n)) * (1 - gsl_vector_get(param->fis, n)) + 
		   gsl_matrix_get(param->phi, i, n) * gsl_vector_get(param->fis, n));
	      }
	      else { // heterozygous genotype
		// p00 = 2 p_0 q_0 * ((1-phi)^2 (1-f) + f(1-phi))
		p00 = 2 * gsl_matrix_get(param->pi0Cur, i, obs[0]) * gsl_matrix_get(param->pi0Cur, i, obs[1]) * 
		  (gsl_pow_2(1 - gsl_matrix_get(param->phi, i, n)) * (1 - gsl_vector_get(param->fis, n)) + 
		   (1 - gsl_matrix_get(param->phi, i, n)) * gsl_vector_get(param->fis, n));
		// p01 = ((p_0 q_1) + (p_1 q_0)) p_1 (1-phi)phi (1-f)
		p01 = (gsl_matrix_get(param->pi0Cur, i, obs[0]) * gsl_matrix_get(param->pi1Cur, i, obs[1]) +
		       gsl_matrix_get(param->pi1Cur, i, obs[0]) * gsl_matrix_get(param->pi0Cur, i, obs[1])) *
		  ((1 - gsl_matrix_get(param->phi, i, n)) * gsl_matrix_get(param->phi, i, n) * 
		   (1 - gsl_vector_get(param->fis, n)));
		p10 = p01;
		// p11 = 2 p_1 q_1 * (phi^2 (1-f) + f(phi))
		p11 = 2 * gsl_matrix_get(param->pi1Cur, i, obs[0]) * gsl_matrix_get(param->pi1Cur, i, obs[1]) * 
		  (gsl_pow_2(gsl_matrix_get(param->phi, i, n)) * (1 - gsl_vector_get(param->fis, n)) + 
		   gsl_matrix_get(param->phi, i, n) * gsl_vector_get(param->fis, n));
	      
	      }
	    }
	    else {
	      cerr << "Indiviudal " << n << " is missing informaiton for locus " << i << endl;
	      exit(0);
	    }
	    ptot = p00 + p01 + p10 + p11;
	    //cerr << ptot << endl;
	    pr[0] = p00 / ptot;
	    pr[1] = p01 / ptot;
	    pr[2] = p10 / ptot;
	    pr[3] = p11 / ptot;
	    // sample ancestry pair based on probabilities
	    gsl_ran_multinomial(r, 4, 1, pr, sam);
	    if(sam[0] == 1){ // homo pop 0
	      gsl_matrix_int_set(param->zCur[i].mat, n, 0, 0);
	      gsl_matrix_int_set(param->zCur[i].mat, n, 1, 0);
	    }
	    else if(sam[1] == 1){ // het
	      gsl_matrix_int_set(param->zCur[i].mat, n, 0, 0);
	      gsl_matrix_int_set(param->zCur[i].mat, n, 1, 1);
	    }
	    else if(sam[2] == 1){ // het
	      gsl_matrix_int_set(param->zCur[i].mat, n, 0, 1);
	      gsl_matrix_int_set(param->zCur[i].mat, n, 1, 0);
	    }
	    else{ // homo pop 1
	      gsl_matrix_int_set(param->zCur[i].mat, n, 0, 1);
	      gsl_matrix_int_set(param->zCur[i].mat, n, 1, 1);
	    }
	  }
	}
      }
      else if (pl == 1){ // haploid
	for(n=0; n<datadim->nind; n++){
	  acopy = 0;
	  p00 = 0;
	  p11 = 0;
	  // missing data, set probabilities based on phi and fis
	  if(gsl_matrix_int_get(data->phcnt[i].mat, n, 0) == -9){
	    p00 = (1 - gsl_matrix_get(param->phi, i, n)); // (1-phi)
	    p11 = gsl_matrix_get(param->phi, i, n); // phi	  
	  }
	  // loop through alleles
	  else{ // individual has data for this locus
	    for(k=0; k<nk; k++){
	      if(gsl_matrix_int_get(data->phcnt[i].mat, n, k) == 1){
		obs[acopy] = k;
		acopy++;
	      } 
	    }
	    if (acopy==1){
	      // p00 = p_0 * (1-phi)
	      p00 = gsl_matrix_get(param->pi0Cur, i, obs[0]) * (1 - gsl_matrix_get(param->phi, i, n));
	      // p11 = p_1 phi
	      p11 = gsl_matrix_get(param->pi1Cur, i, obs[0]) * gsl_matrix_get(param->phi, i, n);
	    }
	    else {
	      cerr << "Missing or too much data for haploid locus " << i << " individual " << n << endl;
	      exit(0);
	    }
	  }
	  ptot = p00 + p01;
	  pr[0] = p00 / ptot;
	  pr[1] = p11 / ptot;
	  pr[2] = 0;
	  pr[3] = 0;
	  gsl_ran_multinomial(r, 2, 1, pr, sam);
	  if(sam[0] == 1){ // pop 0 hap, but set both anyway
	    gsl_matrix_int_set(param->zCur[i].mat, n, 0, 0);
	    gsl_matrix_int_set(param->zCur[i].mat, n, 1, 0);
	  }
	  else{ // homo pop 1
	    gsl_matrix_int_set(param->zCur[i].mat, n, 0, 1);
	    gsl_matrix_int_set(param->zCur[i].mat, n, 1, 1);
	  }
	}
      }	    
    }
  }
}
    

/* function to update hi, uses random-walk metropolis */
void updateHi(paramcont * param, datadimcont * datadim, auxcont * auxvar, int allDiploid){
  int n;
  double value, u, logMratio, dev;
  
  // loop through individuals
  for(n=0; n<datadim->nind; n++){
    // propose hi = old value + random deviate from uniform on - hiprop to + hiprop
    dev = gsl_ran_flat(r, (-1 * param->hiprop), param->hiprop);
    value = gsl_matrix_get(param->hi, n, 0) + dev;
    if (0 <= value && value <= 1){ // otherwise posterior = 0, so need to set to current
      gsl_matrix_set(param->hi, n, 0, value);
      logMratio = calcLogMratioHi(param,datadim,auxvar,gsl_vector_int_get(datadim->pops,n),n,allDiploid);
      u = log(gsl_ran_flat(r, 0, 1));
      if (logMratio < u){ // not keeping proposed value
	gsl_matrix_set(param->hi, n, 0, gsl_matrix_get(param->hi, n, 1));
      }
    }
    else {
      gsl_matrix_set(param->hi, n, 0,  gsl_matrix_get(param->hi, n, 1));
    }
  }
}

/* calculates metropolis ratio for hybrid index */
double calcLogMratioHi(paramcont * param, datadimcont * datadim, auxcont * auxvar, 
		       int j, int n, int allDiploid){
  double pNew = 0;
  double pOld = 0;
  int i;
  double onePhi, y, mMax, mMin, logMratio;
  int splice;
  double prob;
  double f;

  // loop through loci and calc P(z|phi), for proposed
  f = gsl_vector_get(param->fis, n); // fis for individual n
  for(i=0; i<datadim->nloci; i++){
    gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
    gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
    gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
    splice = needSplice(gsl_matrix_get(param->alphaCur,i,j),gsl_matrix_get(param->betaCur,i,j),auxvar->zero2one1,
			auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
    onePhi = calcOnePhi(gsl_matrix_get(param->alphaCur,i,j),gsl_matrix_get(param->betaCur,i,j),
			gsl_matrix_get(param->hi, n, 0),y,mMin,mMax,splice);
     
    if (allDiploid == 1){    // include both allele copies
      // individual has homozygous ancestry for pop 0
      if ((gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1)) && 
	  (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0)){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if ((gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1)) &&
	       (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1)){
	prob = gsl_pow_2(onePhi) * (1 - f) +  onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
    else if (gsl_vector_int_get(datadim->ploidyVec,i) == 2){
      // individual has homozygous ancestry for pop 0
      if ((gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1)) 
	  && (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0)){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	       gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
    else {
      prob = gsl_ran_bernoulli_pdf(gsl_matrix_int_get(param->zCur[i].mat, n, 0), onePhi);
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
  }
  // prior is uniform no need to add

  // loop through loci and calc P(z|phi), for old
  for(i=0; i<datadim->nloci; i++){
    onePhi = gsl_matrix_get(param->phi, i, n);
 
    if (allDiploid == 1){    // include both allele copies
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	       gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
    else if (gsl_vector_int_get(datadim->ploidyVec,i) == 2){
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	 gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) +  (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	       gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
    else {
      prob = gsl_ran_bernoulli_pdf(gsl_matrix_int_get(param->zCur[i].mat, n, 0), onePhi);
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
  }
  // prior is uniform no need to add
  logMratio = pNew - pOld;
  //cerr << "new: " << pNew << ", old: " << pOld << ", ratio: " << logMratio << endl;

  return (logMratio);
}
/* function to update fis, uses random-walk metropolis */
void updateFis(gsl_matrix * h, mat_container_int * zCur, gsl_matrix * alphaCur,
	       gsl_matrix * betaCur, gsl_vector * zero2one, gsl_vector * zero2one1,
	       gsl_vector * zero2one2, gsl_vector * zero2one3, gsl_vector * M, 
	       gsl_vector_int * pops, double tuneFis, int nind, 
	       int nloci, gsl_vector_int * ploidyVec, int allDiploid, 
	       gsl_vector * fis){
  int n;
  double value, u, logMratio, dev;
  
  // loop through individuals
  for(n=0; n<nind; n++){
    // propose fis = old value + random deviate from uniform on - tuneFis to + tuenFis
    dev = gsl_ran_flat(r, (-1 * tuneFis), tuneFis);
    value = gsl_vector_get(fis, n) + dev;
    if (-1 <= value && value <= 1){ // otherwise posterior = 0, so need to set to current
      logMratio = calcLogMratioFis(h,zCur,alphaCur,betaCur,zero2one,zero2one1,zero2one2,
				   zero2one3,M,gsl_vector_int_get(pops,n),n,nloci,
				   ploidyVec, allDiploid, fis, value);
      u = log(gsl_ran_flat(r, 0, 1));
      if (logMratio >= u){ // assign proposed value to fis
	gsl_vector_set(fis, n, value);
      }
    }
  }
}

/* calculates metropolis ratio for fis */
double calcLogMratioFis(gsl_matrix * h, mat_container_int * zCur, gsl_matrix * alphaCur,
			gsl_matrix * betaCur, gsl_vector * zero2one, gsl_vector * zero2one1,
			gsl_vector * zero2one2, gsl_vector * zero2one3, gsl_vector * M,
		        int j, int n, int nloci, gsl_vector_int * ploidyVec, 
			int allDiploid, gsl_vector * fis, double newf){
  double pNew = 0;
  double pOld = 0;
  int i;
  double onePhi, y, mMax, mMin, logMratio;
  int splice;
  double prob;

  double oldf = gsl_vector_get(fis, n); // old fis for individual n

  // loop through loci and calc P(z|phi), for proposed
  for(i=0; i<nloci; i++){
    gsl_vector_memcpy(zero2one1,zero2one);
    gsl_vector_memcpy(zero2one2,zero2one);
    gsl_vector_memcpy(zero2one3,zero2one);
    splice = needSplice(gsl_matrix_get(alphaCur,i,j),gsl_matrix_get(betaCur,i,j),zero2one1,
			zero2one2,zero2one3,&y,&mMax,&mMin,M);
    onePhi = calcOnePhi(gsl_matrix_get(alphaCur,i,j),gsl_matrix_get(betaCur,i,j),
			gsl_matrix_get(h, n, 0),y,mMin,mMax,splice);
     
    if (allDiploid == 1){    // include both allele copies
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(zCur[i].mat, n, 0) == gsl_matrix_int_get(zCur[i].mat, n, 1) && 
	  gsl_matrix_int_get(zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - newf) + (1 - onePhi) * newf;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(zCur[i].mat, n, 0) == gsl_matrix_int_get(zCur[i].mat, n, 1) &&
	       gsl_matrix_int_get(zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - newf) + onePhi * newf;
     }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - newf);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
    else if (gsl_vector_int_get(ploidyVec,i) == 2){
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(zCur[i].mat, n, 0) == gsl_matrix_int_get(zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - newf) + (1 - onePhi) * newf;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(zCur[i].mat, n, 0) == gsl_matrix_int_get(zCur[i].mat, n, 1) 
	       && gsl_matrix_int_get(zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - newf) + onePhi * newf;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - newf);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
    else {
      prob = gsl_ran_bernoulli_pdf(gsl_matrix_int_get(zCur[i].mat, n, 0), onePhi);
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
  }
  // prior is uniform no need to add

  // loop through loci and calc P(z|phi), for old
  for(i=0; i<nloci; i++){
    if (allDiploid == 1){    // include both allele copies
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(zCur[i].mat, n, 0) == gsl_matrix_int_get(zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - oldf) + (1 - onePhi) * oldf;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(zCur[i].mat, n, 0) == gsl_matrix_int_get(zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - oldf) + onePhi * oldf;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - oldf);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
    else if (gsl_vector_int_get(ploidyVec,i) == 2){
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(zCur[i].mat, n, 0) == gsl_matrix_int_get(zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - oldf) + (1 - onePhi) * oldf;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(zCur[i].mat, n, 0) == gsl_matrix_int_get(zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - oldf) + onePhi * oldf;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - oldf);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
    else {
      prob = gsl_ran_bernoulli_pdf(gsl_matrix_int_get(zCur[i].mat, n, 0), onePhi);
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
  }
  // prior is uniform no need to add
  logMratio = pNew - pOld;
  //cerr << "new: " << pNew << ", old: " << pOld << ", ratio: " << logMratio << endl;
  return (logMratio);
}

/* function to determine whether splice is necessary for the cline function */
int needSplice(double a, double b, gsl_vector * zero2one1, gsl_vector * zero2one2,
	       gsl_vector * zero2one3, double * y, double * mMax, double * mMin,
	       gsl_vector * M){
  double dh0, dh1, sumConst, dmin;
  double coA, coB, coC, coD, d0one, d0two;
  double m1, m2, m3;

  if (b < 0){
    dh0 = 2 * a - 2 * b + 1;
    dh1 = 2 * a - 4 * a - 2 * b + 1;
    sumConst = 2 * a - 2 * b + 1;
    gsl_vector_scale(zero2one1, (-4 * a));
    gsl_vector_scale(zero2one2, 12 * b);
    gsl_vector_mul(zero2one3, zero2one3);
    gsl_vector_scale(zero2one3, -12 * b);
    gsl_vector_add(zero2one1, zero2one2);
    gsl_vector_add(zero2one1, zero2one3);
    gsl_vector_add_constant(zero2one1, sumConst);
    dmin = gsl_vector_min(zero2one1);
    if ((dmin < 0) && (dh0 > 0) && (dh1 > 0)){ // conditions for local extrema
      coA = -12 * b;
      coB = -4 * a + 12 * b;
      coC = 1 + 2 * a - 2 * b;
      gsl_poly_solve_quadratic(coA, coB, coC, &d0one, &d0two);
      // calculate theta at d0one and d0two and take mean
      d0one = 2 * d0one * (1 - d0one) * (a + (2 * d0one -1) * b);
      d0two = 2 * d0two * (1 - d0two) * (a + (2 * d0two -1) * b);
      *y = (d0one + d0two) / 2; // value of theta for splice
      coD = -4 * b;
      coA = (-2 * a + 6 * b) / coD;
      coB = (2 * a - 2 * b + 1) / coD;
      coC = (-1 * *y) / coD;
      gsl_poly_solve_cubic(coA, coB, coC, &m1, &m2, &m3);
      gsl_vector_set(M, 0, m1);
      gsl_vector_set(M, 1, m2);
      gsl_vector_set(M, 2, m3);
      gsl_vector_minmax(M,mMin,mMax);
      return(1);
    }
  }
  return(0); // do not need to splice
}

/* calculates phi for a single combination of hybrid index alpha and beta */
double calcOnePhi(double a, double b, double h, double y, double mMin, double mMax,
		  int splice){
  double theta = 0;
  
  if (splice == 0) {
    theta = h + 2 * h * (1 - h) * (a + b * (2 * h - 1));
  }
  else if (splice == 1){
    if ((mMin <= h) && (h <= mMax)){
      theta = y;
    }
    else {
      theta = h + 2 * h * (1 - h) * (a + b * (2 * h - 1));     
    }
  }
  if ((0 <= theta)  && (theta <= 1)){
    return(theta);
  }
  else if (theta < 0){
    return(0);
  }
  else {
    return(1);
  }
}

/* function to update gamma, the locus effect on cline center, uses metropolis random walk */
void updateGamma(paramcont * param, datadimcont * datadim, auxcont * auxvar, int allDiploid, 
		 int sumZ){
  int i;
  double prop, logMratio, u;
  double paramsum = 0;

  for(i=0; i<datadim->nloci; i++){
    // gaussian proposal for gamma centered on current value
    prop = gsl_ran_gaussian(r, param->tuneGamma) + gsl_matrix_get(param->gamma, i, 1);
    gsl_matrix_set(param->gamma, i, 0, prop);
    
    logMratio = calcLogMratioGamma(param,datadim,auxvar,i,allDiploid);
    u = log(gsl_ran_flat(r, 0, 1));
    if (logMratio < u){ // not keeping proposed value
      gsl_matrix_set(param->gamma, i, 0, gsl_matrix_get(param->gamma, i, 1));
    }
    paramsum += gsl_matrix_get(param->gamma,i,0);
  }
  paramsum = paramsum / datadim->nloci;
  if (sumZ == 1){
    for (i=0; i<datadim->nloci; i++){
      gsl_matrix_set(param->gamma,i,0, gsl_matrix_get(param->gamma,i,0) - paramsum);
    }
  }
}

/* calculates metropolis ratio for gamma */
double calcLogMratioGamma(paramcont * param, datadimcont * datadim, auxcont * auxvar,
			  int i, int allDiploid){
  double pNew = 0;
  double pOld = 0;
  int n, j, splice;
  double onePhi, y, mMax, mMin, logMratio;
  double oneAlpha;
  double prob;
  double f;

  gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
  
  j = gsl_vector_int_get(datadim->pops, 0);
  oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
  splice = needSplice(oneAlpha,gsl_matrix_get(param->betaCur,i,j),auxvar->zero2one1,
		      auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
  // for proposed value; loop through individuals
  for(n=0; n<datadim->nind; n++){
    // only need to recalculate alpha and splice if pop changes
    if (j != gsl_vector_int_get(datadim->pops, n)){
      j = gsl_vector_int_get(datadim->pops, n);
      gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
      oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
      splice = needSplice(oneAlpha,gsl_matrix_get(param->betaCur,i,j),auxvar->zero2one1,
			  auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
    }
    onePhi = calcOnePhi(oneAlpha,gsl_matrix_get(param->betaCur,i,j),
			gsl_matrix_get(param->hi, n, 0),y,mMin,mMax,splice);
    f = gsl_vector_get(param->fis, n); // fis for individual n

    if (allDiploid == 1){    // include both allele copies
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) +  onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
    else if (gsl_vector_int_get(datadim->ploidyVec,i) == 2){
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
    else {
      prob = gsl_ran_bernoulli_pdf(gsl_matrix_int_get(param->zCur[i].mat, n, 0), onePhi);
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
  }
  // conditional prior
  pNew += log(gsl_ran_gaussian_pdf(gsl_matrix_get(param->gamma,i,0),sqrt(1/param->tauAlphaCur)));
 
  // for old value; loop through individuals

  gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
  
  j = gsl_vector_int_get(datadim->pops, 0);
  splice = needSplice(gsl_matrix_get(param->alphaCur,i,j),gsl_matrix_get(param->betaCur,i,j),
		      auxvar->zero2one1,auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,
		      auxvar->M);
  for(n=0; n<datadim->nind; n++){
    if (j != gsl_vector_int_get(datadim->pops, n)){ // only need to recalculate alpha and splice if pop changes
      j = gsl_vector_int_get(datadim->pops, n);
      gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
      splice = needSplice(gsl_matrix_get(param->alphaCur,i,j),gsl_matrix_get(param->betaCur,i,j),
			  auxvar->zero2one1,auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,
			  auxvar->M);
    }
    onePhi = calcOnePhi(gsl_matrix_get(param->alphaCur,i,j),gsl_matrix_get(param->betaCur,i,j),
			gsl_matrix_get(param->hi, n, 0),y,mMin,mMax,splice);
    f = gsl_vector_get(param->fis, n); // fis for individual n

    if (allDiploid == 1){    // include both allele copies
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) +  (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
    else if (gsl_vector_int_get(datadim->ploidyVec,i) == 2){
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
    else {
      prob = gsl_ran_bernoulli_pdf(gsl_matrix_int_get(param->zCur[i].mat, n, 0), onePhi);
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
  }
  // conditional prior
  pOld += log(gsl_ran_gaussian_pdf(gsl_matrix_get(param->gamma,i,1),sqrt(1/param->tauAlphaCur)));

  logMratio = pNew - pOld;

  //cerr << "gamma: " << pNew << "," << pOld << "," << logMratio << endl;
  return (logMratio);
}

/* function to update zeta, the locus effect on cline width, uses metropolis random walk */
void updateZeta(paramcont * param, datadimcont * datadim, auxcont * auxvar, int allDiploid, int sumZ){
  int i;
  double prop, logMratio, u;
  double paramsum = 0;

  for(i=0; i<datadim->nloci; i++){
    // gaussian proposal for zeta centered on current value
    prop = gsl_ran_gaussian(r, param->tuneZeta) + gsl_matrix_get(param->zeta, i, 1);
    gsl_matrix_set(param->zeta, i, 0, prop);
    //cerr << gsl_matrix_get(zeta,i,0) << ","  << gsl_matrix_get(zeta,i,1) << endl;

    logMratio = calcLogMratioZeta(param,datadim,auxvar,i,allDiploid);
    u = log(gsl_ran_flat(r, 0, 1));
    if (logMratio < u){ // not keeping proposed value
      gsl_matrix_set(param->zeta, i, 0, gsl_matrix_get(param->zeta, i, 1));
    }
    paramsum += gsl_matrix_get(param->zeta,i,0);
  }
  paramsum = paramsum / datadim->nloci;
  if (sumZ == 1){
    for (i=0; i<datadim->nloci; i++){
      gsl_matrix_set(param->zeta,i,0, gsl_matrix_get(param->zeta,i,0) - paramsum);
    }
  }
}

/* calculates metropolis ratio for zeta */
double calcLogMratioZeta(paramcont * param, datadimcont * datadim, auxcont * auxvar,
			 int i, int allDiploid){
  double pNew = 0;
  double pOld = 0;
  int n, j, splice;
  double onePhi, y, mMax, mMin, logMratio;
  double oneBeta, oneAlpha;
  double prob;
  double f;

  gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
  
  j = gsl_vector_int_get(datadim->pops, 0);
  oneBeta = gsl_matrix_get(param->zeta,i,0) + gsl_matrix_get(param->kappaCur,i,j);
  oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
  splice = needSplice(oneAlpha,oneBeta,auxvar->zero2one1,
		      auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);

  // for proposed value; loop through individuals
  for(n=0; n<datadim->nind; n++){
    // only need to recalculate alpha, beta  and splice if pop changes
    if (j != gsl_vector_int_get(datadim->pops, n)){ 
      j = gsl_vector_int_get(datadim->pops, n);
      gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
      oneBeta = gsl_matrix_get(param->zeta,i,0) + gsl_matrix_get(param->kappaCur,i,j);
      oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
      splice = needSplice(oneAlpha,oneBeta,auxvar->zero2one1,
			  auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
    }
    onePhi = calcOnePhi(oneAlpha,oneBeta,
			gsl_matrix_get(param->hi, n, 0),y,mMin,mMax,splice);
    f = gsl_vector_get(param->fis, n); // fis for individual n
    if (allDiploid == 1){    // include both allele copies
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) +  onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
    else if (gsl_vector_int_get(datadim->ploidyVec,i) == 2){
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
    else {
      prob = gsl_ran_bernoulli_pdf(gsl_matrix_int_get(param->zCur[i].mat, n, 0), onePhi);
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
  }
  // conditional prior
  pNew += log(gsl_ran_gaussian_pdf(gsl_matrix_get(param->zeta,i,0),sqrt(1/param->tauBetaCur)));
 
  gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
  
  j = gsl_vector_int_get(datadim->pops, 0);
  oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
  splice = needSplice(oneAlpha,gsl_matrix_get(param->betaCur,i,j),auxvar->zero2one1,
		      auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);

  // for old value; loop through individuals
  for(n=0; n<datadim->nind; n++){
    // only need to recalculate alpha, beta  and splice if pop changes
    if (j != gsl_vector_int_get(datadim->pops, n)){
      j = gsl_vector_int_get(datadim->pops, n);
      gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
      oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
      splice = needSplice(oneAlpha,gsl_matrix_get(param->betaCur,i,j),auxvar->zero2one1,
			  auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
    }
    onePhi = calcOnePhi(oneAlpha,gsl_matrix_get(param->betaCur,i,j),
			gsl_matrix_get(param->hi, n, 0),y,mMin,mMax,splice);
    f = gsl_vector_get(param->fis, n); // fis for individual n
    if (allDiploid == 1){    // include both allele copies
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
    else if (gsl_vector_int_get(datadim->ploidyVec,i) == 2){
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
    else {
      prob = gsl_ran_bernoulli_pdf(gsl_matrix_int_get(param->zCur[i].mat, n, 0), onePhi);
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }    
  }

  // conditional prior
  pOld += log(gsl_ran_gaussian_pdf(gsl_matrix_get(param->zeta,i,1),sqrt(1/param->tauBetaCur)));

  logMratio = pNew - pOld;
  //cerr << "zeta: " << pNew << "," << pOld << "," << logMratio << endl;
  return (logMratio);
}

/* update eta, zero centered nested population effects, uses metropolis with MVN */
void updateEta(paramcont * param, datadimcont * datadim, auxcont * auxvar,
	       int i, int allDiploid, int onePrecision){
  int j;
  double logMratio, u;
  double sumdev = 0;
  // proposes nested population effects, proposes vector of effects from MVN
  rmvnorm(r,auxvar->vecPopn1,param->etaPropMat,auxvar->vecPopn2); 
  for(j=0; j<(datadim->npop-1); j++){
    gsl_matrix_set(param->etaCur, i, j, gsl_matrix_get(param->etaLast, i, j) + 
		   gsl_vector_get(auxvar->vecPopn2, j));
  }
  for(j=0; j<(datadim->npop-1); j++){ // not currently doing anything with sumdev, could be used for s2z
    sumdev += gsl_matrix_get(param->etaCur, i, j);
  }

  gsl_matrix_set(param->etaCur, i, j, -1 * sumdev); // set last, reparameterization by sweeping

  logMratio = calcLogMratioEta(param,datadim,auxvar,i,allDiploid,onePrecision);
  u = log(gsl_ran_flat(r, 0, 1));
  if (logMratio < u){ // not keeping proposed value
    for(j=0; j<datadim->npop; j++){
      gsl_matrix_set(param->etaCur, i, j, gsl_matrix_get(param->etaLast, i, j));
    }
  }
}

/* calculates metropolis ratio for eta */
double calcLogMratioEta(paramcont * param, datadimcont * datadim, auxcont * auxvar, int i, 
			int allDiploid,int onePrecision){
  double pNew = 0;
  double pOld = 0;
  int n, j, w, splice;
  double onePhi, y, mMax, mMin, logMratio;
  double oneAlpha, oneBeta;
  double prob;
  double f;

  gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
  
  j = gsl_vector_int_get(datadim->pops, 0);
  oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
  oneBeta = gsl_matrix_get(param->zeta,i,0) + gsl_matrix_get(param->kappaCur,i,j);
  splice = needSplice(oneAlpha,oneBeta,auxvar->zero2one1,
		      auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
  // for proposed value; loop through individuals
  for(n=0; n<datadim->nind; n++){
    if (j != gsl_vector_int_get(datadim->pops, n)){ // only need to recalculate alpha and splice if pop changes
      j = gsl_vector_int_get(datadim->pops, n);
      gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
      oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
      oneBeta = gsl_matrix_get(param->zeta,i,0) + gsl_matrix_get(param->kappaCur,i,j);
      splice = needSplice(oneAlpha,oneBeta,auxvar->zero2one1,
			  auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
    }
    onePhi = calcOnePhi(oneAlpha,oneBeta,
			gsl_matrix_get(param->hi, n, 0),y,mMin,mMax,splice);
    f = gsl_vector_get(param->fis, n); // fis for individual n
    if (allDiploid == 1){    // include both allele copies
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) +  onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
    else if (gsl_vector_int_get(datadim->ploidyVec,i) == 2){
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
    else {
      prob = gsl_ran_bernoulli_pdf(gsl_matrix_int_get(param->zCur[i].mat, n, 0), onePhi);
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
  }

  // conditional prior, from MVN
  for(j=0; j<(datadim->npop-1); j++){
    gsl_vector_set(auxvar->vecPopn2, j, gsl_matrix_get(param->etaCur, i, j));
  }		   
  for(j=0; j<(datadim->npop-1); j++){
    for(w=0; w<(datadim->npop-1); w++){
      if(w==j){
	if (onePrecision == 1){
	  // just use first nu, one precision across loci
	  gsl_matrix_set(param->popCoVaMat, j, w, 
			 (1 / gsl_matrix_get(param->nu, 0, 0)) * (1 - 1 / datadim->npop));
	}
	else{
	  gsl_matrix_set(param->popCoVaMat, j, w, 
			 (1 / gsl_matrix_get(param->nu, i, 0)) * (1 - 1 / datadim->npop));
	}
      }
      else{
	gsl_matrix_set(param->popCoVaMat, j, w, (-1 / datadim->npop));
      }
    }
  }
  pNew += log(dmvnorm(auxvar->vecPopn2, auxvar->vecPopn1, param->popCoVaMat));


  // for old value; loop through individuals
  gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
  
  j = gsl_vector_int_get(datadim->pops, 0);
  oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaLast,i,j);
  oneBeta = gsl_matrix_get(param->zeta,i,0) + gsl_matrix_get(param->kappaCur,i,j);
  splice = needSplice(oneAlpha,oneBeta,auxvar->zero2one1,
		      auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
  for(n=0; n<datadim->nind; n++){
    if (j != gsl_vector_int_get(datadim->pops, n)){ // only need to recalculate alpha and splice if pop changes
      j = gsl_vector_int_get(datadim->pops, n);
      gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
      oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaLast,i,j);
      oneBeta = gsl_matrix_get(param->zeta,i,0) + gsl_matrix_get(param->kappaCur,i,j);
      splice = needSplice(oneAlpha,oneBeta,auxvar->zero2one1,
			  auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
    }

    onePhi = calcOnePhi(oneAlpha,oneBeta,
			gsl_matrix_get(param->hi, n, 0),y,mMin,mMax,splice);
    f = gsl_vector_get(param->fis, n); // fis for individual n
    if (allDiploid == 1){    // include both allele copies
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) +  onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
    else if (gsl_vector_int_get(datadim->ploidyVec,i) == 2){
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
    else {
      prob = gsl_ran_bernoulli_pdf(gsl_matrix_int_get(param->zCur[i].mat, n, 0), onePhi);
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
  }
  // conditional prior, from MVN
  for(j=0; j<(datadim->npop-1); j++){
    gsl_vector_set(auxvar->vecPopn2, j, gsl_matrix_get(param->etaLast, i, j));
  }	
  //gsl_matrix_get_row(vecPopn2, etaLast, i);
  for(j=0; j<(datadim->npop-1); j++){
    for(w=0; w<(datadim->npop-1); w++){
      if(w==j){
	if (onePrecision == 1){
	  // share precision across loci, use nu[0,0]
	  gsl_matrix_set(param->popCoVaMat, j, w, 
			 (1 / gsl_matrix_get(param->nu, 0, 0))  * (1 - 1 / datadim->npop));
	}
	else{
	  gsl_matrix_set(param->popCoVaMat, j, w, 
			 (1 / gsl_matrix_get(param->nu, i, 0))  * (1 - 1 / datadim->npop));
	}
      }
      else{
	gsl_matrix_set(param->popCoVaMat, j, w, (-1 / datadim->npop));
      }
    }
  }
  pOld += log(dmvnorm(auxvar->vecPopn2, auxvar->vecPopn1, param->popCoVaMat));
  logMratio = pNew - pOld;
  //cerr << "eta:  " << pNew << "," << pOld << endl;
  return (logMratio);
}
/* update kappa, zero centered nested population effects, uses metropolis with MVN */
void updateKappa(paramcont * param, datadimcont * datadim, auxcont * auxvar,
		 int i, int allDiploid, int onePrecision){
  int j;
  double logMratio, u;
  double sumdev = 0;
  // proposes nested population effects, proposes vector of effects from MVN
  rmvnorm(r,auxvar->vecPopn1,param->etaPropMat,auxvar->vecPopn2); 
  for(j=0; j<(datadim->npop-1); j++){
    gsl_matrix_set(param->kappaCur, i, j, gsl_matrix_get(param->kappaLast, i, j) + 
		   gsl_vector_get(auxvar->vecPopn2, j));
  }
  for(j=0; j<(datadim->npop-1); j++){
    sumdev += gsl_matrix_get(param->kappaCur, i, j);
  }
  gsl_matrix_set(param->kappaCur, i, j, -1 * sumdev); // set last, reparameterization by sweepingv

  logMratio = calcLogMratioKappa(param,datadim,auxvar,i,allDiploid,onePrecision);
  u = log(gsl_ran_flat(r, 0, 1));
  if (logMratio < u){ // not keeping proposed value
    for(j=0; j<datadim->npop; j++){
      gsl_matrix_set(param->kappaCur, i, j, gsl_matrix_get(param->kappaLast, i, j));
    }
  }
}

/* calculates metropolis ratio for kappa */
double calcLogMratioKappa(paramcont * param, datadimcont * datadim, auxcont * auxvar,
			  int i, int allDiploid,int onePrecision){
  double pNew = 0;
  double pOld = 0;
  int n, j, w, splice;
  double onePhi, y, mMax, mMin, logMratio;
  double oneAlpha, oneBeta;
  double prob;
  double f;

  gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
  
  j = gsl_vector_int_get(datadim->pops, 0);
  oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
  oneBeta = gsl_matrix_get(param->zeta,i,0) + gsl_matrix_get(param->kappaCur,i,j);
  splice = needSplice(oneAlpha,oneBeta,auxvar->zero2one1,
		      auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
  // for proposed value; loop through individuals
  for(n=0; n<datadim->nind; n++){
    if (j != gsl_vector_int_get(datadim->pops, n)){ // only need to recalculate alpha and splice if pop changes
      j = gsl_vector_int_get(datadim->pops, n);
      gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
      oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
      oneBeta = gsl_matrix_get(param->zeta,i,0) + gsl_matrix_get(param->kappaCur,i,j);
      splice = needSplice(oneAlpha,oneBeta,auxvar->zero2one1,
			  auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
    }
    onePhi = calcOnePhi(oneAlpha,oneBeta,
			gsl_matrix_get(param->hi, n, 0),y,mMin,mMax,splice);
    f = gsl_vector_get(param->fis, n); // fis for individual n
    if (allDiploid == 1){    // include both allele copies
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) +  onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
    else if (gsl_vector_int_get(datadim->ploidyVec,i) == 2){
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
    else {
      prob = gsl_ran_bernoulli_pdf(gsl_matrix_int_get(param->zCur[i].mat, n, 0), onePhi);
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
    }
  }
  // conditional prior, from MVN
  for(j=0; j<(datadim->npop-1); j++){
    gsl_vector_set(auxvar->vecPopn2, j, gsl_matrix_get(param->kappaCur, i, j));
  }	
  //gsl_matrix_get_row(vecPopn2, kappaCur, i);
  for(j=0; j<(datadim->npop-1); j++){
    for(w=0; w<(datadim->npop-1); w++){
      if(w==j){
	if(onePrecision == 1){
	  // use one precision omega across loci
	  gsl_matrix_set(param->popCoVaMat, j, w, 
			 (1 / gsl_matrix_get(param->omega, 0, 0)) * (1 - 1 / datadim->npop));
	}
	else{
	  gsl_matrix_set(param->popCoVaMat, j, w, 
			 (1 / gsl_matrix_get(param->omega, i, 0)) * (1 - 1 / datadim->npop));
	}
      }
      else{
	gsl_matrix_set(param->popCoVaMat, j, w, (-1 / datadim->npop));
      }
    }
  }
  pNew += log(dmvnorm(auxvar->vecPopn2, auxvar->vecPopn1, param->popCoVaMat));
 
  // for old value; loop through individuals
  gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
  
  j = gsl_vector_int_get(datadim->pops, 0);
  oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
  oneBeta = gsl_matrix_get(param->zeta,i,0) + gsl_matrix_get(param->kappaLast,i,j);
  splice = needSplice(oneAlpha,oneBeta,auxvar->zero2one1,
		      auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);

  for(n=0; n<datadim->nind; n++){
    if (j != gsl_vector_int_get(datadim->pops, n)){ // only need to recalculate alpha and splice if pop changes
      j = gsl_vector_int_get(datadim->pops, n);
      gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
      oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
      oneBeta = gsl_matrix_get(param->zeta,i,0) + gsl_matrix_get(param->kappaLast,i,j);
      splice = needSplice(oneAlpha,oneBeta,auxvar->zero2one1,
			  auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
    }

    onePhi = calcOnePhi(oneAlpha,oneBeta,
			gsl_matrix_get(param->hi, n, 0),y,mMin,mMax,splice);
    f = gsl_vector_get(param->fis, n); // fis for individual n
    if (allDiploid == 1){    // include both allele copies
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) +  (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
    else if (gsl_vector_int_get(datadim->ploidyVec,i) == 2){
      // individual has homozygous ancestry for pop 0
      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	prob = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
      }
      // individual has homozygous ancestry for pop 1
      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	prob = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
      }
      // individual has heterozygous ancestry
      else { 
	prob = 2 * (1 - onePhi) * onePhi * (1 - f);
      }
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }
    else {
      prob = gsl_ran_bernoulli_pdf(gsl_matrix_int_get(param->zCur[i].mat, n, 0), onePhi);
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }    
  }
  // conditional prior, from MVN
  for(j=0; j<(datadim->npop-1); j++){
    gsl_vector_set(auxvar->vecPopn2, j, gsl_matrix_get(param->kappaLast, i, j));
  }	
  //gsl_matrix_get_row(vecPopn2, kappaLast, i);
  for(j=0; j<(datadim->npop-1); j++){
    for(w=0; w<(datadim->npop-1); w++){
      if(w==j){
	if(onePrecision == 1){
	  // use one precision parameter omega across loci
	  gsl_matrix_set(param->popCoVaMat, j, w,
			 (1 / gsl_matrix_get(param->omega, 0, 0)) * (1  - 1 / datadim->npop));
	}
	else{
	  gsl_matrix_set(param->popCoVaMat, j, w, 
			 (1 / gsl_matrix_get(param->omega, i, 0)) * (1  - 1 / datadim->npop));
	}
      }
      else{
	gsl_matrix_set(param->popCoVaMat, j, w, (-1 / datadim->npop));
      }
    }
  }
  pOld += log(dmvnorm(auxvar->vecPopn2, auxvar->vecPopn1, param->popCoVaMat));
  logMratio = pNew - pOld;
  //cerr <<  "kappa:  " << pNew << "," << pOld << endl;
  return (logMratio);
}

/* update zero centered precision parameters for locus effects (i.e., tauAlpha and tauBeta), uses gibbs sampling */
/* trunc > 0 for truncated gamma */
void updateTau(gsl_matrix * parameters, double * precision, int nloci, double trunc){
  int i; 
  double v = 0;
  double a, b;
  double gampriora = 0.001;  // shape parameter
  double gampriorb = 0.001;  // rate paramerer
  
  *precision = -9; // set to -9 to enter loop

  for(i=0; i<nloci; i++){
    v += gsl_pow_2(gsl_matrix_get(parameters, i, 0));    
  }

  a = nloci * 0.5 + gampriora;
  b = 0.5 * v + gampriorb; //
  while (*precision <= trunc){ // loop until precision is within the bounds of the gamma prior
    *precision = gsl_ran_gamma(r, a, 1/b); // 1/b makes it a scale parameter for gsl
  }
}

/* update zero centered precision parameters for the nested population effects, uses gibbs sampling */
void updateNuOmega(gsl_matrix * parameters, gsl_matrix * precision, int nloci, int npop,
		   int onePrecision){
  int i, j; 
  double v = 0;
  double a, b;
  double gampriora = 0.001; // shape parameter
  double gampriorb = 0.001; // rate paramerer
 
  if (onePrecision == 0){ // seperate precision parameter for each locus
    for(i=0; i<nloci; i++){
      v = 0;
      for(j=0; j<npop; j++){
	v += gsl_pow_2(gsl_matrix_get(parameters, i, j));    
      }

      a = (npop / 2) + gampriora;
      b = 0.5 * v + gampriorb;
      gsl_matrix_set(precision, i, 0, gsl_ran_gamma(r, a, 1/b)); // 1/b makes it a scale parameter for gsl
    }
  }

  else{ // one precision parameter for all loci
    v = 0;
    for(i=0; i<nloci; i++){
      for(j=0; j<npop; j++){
	v += gsl_pow_2(gsl_matrix_get(parameters, i, j));    
      }
    }
    a = ((npop * nloci) / 2) + gampriora;
    b = 0.5 * v + gampriorb;
    gsl_matrix_set(precision, 0, 0, gsl_ran_gamma(r, a, 1/b)); // 1/b makes it a scale parameter for gsl
  }
}

/* calculate phi based on current mcmc samples */
void calcPhi(paramcont * param, datadimcont * datadim, auxcont * auxvar){
  int i, j, n;
  double onePhi, y, mMax, mMin;
  int splice;

  for(i=0; i<datadim->nloci; i++){
    gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
    gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
    gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
    j = gsl_vector_int_get(datadim->pops, 0);
    splice = needSplice(gsl_matrix_get(param->alphaCur,i,j),gsl_matrix_get(param->betaCur,i,j),auxvar->zero2one1,
			auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
    for(n=0; n<datadim->nind; n++){
      if (j != gsl_vector_int_get(datadim->pops,n)){
	gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
	gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
	gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
	j = gsl_vector_int_get(datadim->pops, n);
	splice = needSplice(gsl_matrix_get(param->alphaCur,i,j),gsl_matrix_get(param->betaCur,i,j),auxvar->zero2one1,
			    auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
      }
      onePhi = calcOnePhi(gsl_matrix_get(param->alphaCur,i,j),gsl_matrix_get(param->betaCur,i,j),
			  gsl_matrix_get(param->hi, n, 0),y,mMin,mMax,splice);
      gsl_matrix_set(param->phi,i,n,onePhi);
    }
  }

}

/* calculate Ln Likelihood, which includes all model terms except hyperpriors and parental allele frequencies */
double calcLnL(paramcont * param, datacont * data, datadimcont * datadim, int allDiploid){
  int i, n, k, nk;
  double prob = 0;
  int acopy, pl;
  int obs[2]; // obs alleles
  double ptemp = 0;
  double f, onePhi;

 
  // p(x|pi,z)
  // only have diploid data
  if (allDiploid == 1){
    for(i=0; i<datadim->nloci; i++){
      nk = gsl_vector_int_get(datadim->nallele, i);
      for(n=0; n<datadim->nind; n++){
	acopy = 0;
	f = gsl_vector_get(param->fis, n); // fis for individual n
	onePhi = gsl_matrix_get(param->phi, i, n);
	for(k=0; k<nk; k++){
	  if(gsl_matrix_int_get(data->phcnt[i].mat, n, k) == 1){
	    obs[acopy] = k;
	    acopy++;
	  }
	  else if (gsl_matrix_int_get(data->phcnt[i].mat, n, k) == 2){
	    obs[acopy] = k;
	    acopy++;
	    obs[acopy] = k;
	    acopy++;	      
	  }
	}
	// individual has homozygous ancestry for pop 0
	if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	  ptemp = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
	  if (acopy == 2){
	    if (obs[0] == obs[1]){ // homozy. allele
	      ptemp *= gsl_pow_2(gsl_matrix_get(param->pi0Cur, i, obs[0]));
	    }
	    else {
	      ptemp *= 2 * gsl_matrix_get(param->pi0Cur, i, obs[0]) * gsl_matrix_get(param->pi0Cur, i, obs[1]);
	    }
	  }
	}
	// individual has homozygous ancestry for pop 1
	else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	  ptemp = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
	  if (acopy == 2){
	    if (obs[0] == obs[1]){ // homozy. allele
	      ptemp *= gsl_pow_2(gsl_matrix_get(param->pi1Cur, i, obs[0]));
	    }
	    else {
	      ptemp *= 2 * gsl_matrix_get(param->pi1Cur, i, obs[0]) * gsl_matrix_get(param->pi1Cur, i, obs[1]);
	    }
	  }
	}
	// individual has heterozygous ancestry
	else { 
	  ptemp = (1 - onePhi) * onePhi * ((1 - f) * 0.5);
	  if (acopy == 2){
	    if (obs[0] == obs[1]){ // homozy. allele
	      ptemp *= 2 * gsl_matrix_get(param->pi0Cur, i, obs[0]) * gsl_matrix_get(param->pi0Cur, i, obs[1]);
	    }
	    else {
	      ptemp *= (gsl_matrix_get(param->pi0Cur, i, obs[0]) * gsl_matrix_get(param->pi1Cur, i, obs[1]) + 
			gsl_matrix_get(param->pi1Cur, i, obs[0]) * gsl_matrix_get(param->pi0Cur, i, obs[1]));
	    }
	  }
	}
	if (gsl_isnan(ptemp) == 1){
	  ptemp = DBL_MIN;
	}
	if (ptemp <= 0){
	  ptemp = DBL_MIN;
	}
	prob += log(ptemp);
      }
    }
  }

  // diploid and haploid
  else{
    for(i=0; i<datadim->nloci; i++){
      nk = gsl_vector_int_get(datadim->nallele, i);
      pl = gsl_vector_int_get(datadim->ploidyVec, i);
      if (pl==2){ // diploid
	for(n=0; n<datadim->nind; n++){
	  acopy = 0;
	  f = gsl_vector_get(param->fis, n); // fis for individual n
	  onePhi = gsl_matrix_get(param->phi, i, n);
	  for(k=0; k<nk; k++){
	    if(gsl_matrix_int_get(data->phcnt[i].mat, n, k) == 1){
	      obs[acopy] = k;
	      acopy++;
	    }
	    else if (gsl_matrix_int_get(data->phcnt[i].mat, n, k) == 2){
	      obs[acopy] = k;
	      acopy++;
	      obs[acopy] = k;
	      acopy++;	      
	    }
	  }
	  // individual has homozygous ancestry for pop 0
	  if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	    ptemp = gsl_pow_2(1 - onePhi) * (1 - f) + (1 - onePhi) * f;
	    if (acopy == 2){
	      if (obs[0] == obs[1]){ // homozy. allele
		ptemp *= gsl_pow_2(gsl_matrix_get(param->pi0Cur, i, obs[0]));
	      }
	      else {
		ptemp *= 2 * gsl_matrix_get(param->pi0Cur, i, obs[0]) * gsl_matrix_get(param->pi0Cur, i, obs[1]);
	      }
	    }
	  }
	  // individual has homozygous ancestry for pop 1
	  else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
	  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	    ptemp = gsl_pow_2(onePhi) * (1 - f) + onePhi * f;
	    if (acopy == 2){
	      if (obs[0] == obs[1]){ // homozy. allele
		ptemp *= gsl_pow_2(gsl_matrix_get(param->pi1Cur, i, obs[0]));
	      }
	      else {
		ptemp *= 2 * gsl_matrix_get(param->pi1Cur, i, obs[0]) * gsl_matrix_get(param->pi1Cur, i, obs[1]);
	      }
	    }
	  }
	  // individual has heterozygous ancestry
	  else { 
	    ptemp = (1 - onePhi) * onePhi * ((1 - f) * 0.5);
	    if (acopy == 2){
	      if (obs[0] == obs[1]){ // homozy. allele
		ptemp *= 2 * gsl_matrix_get(param->pi0Cur, i, obs[0]) * gsl_matrix_get(param->pi1Cur, i, obs[1]);
	      }
	      else {
		ptemp *= (gsl_matrix_get(param->pi0Cur, i, obs[0]) * gsl_matrix_get(param->pi1Cur, i, obs[1]) + 
			  gsl_matrix_get(param->pi1Cur, i, obs[0]) * gsl_matrix_get(param->pi0Cur, i, obs[1]));
	      }
	    }
	  }
	}
	if (gsl_isnan(ptemp) == 1){
	  ptemp = DBL_MIN;
	}
	if (ptemp <= 0){
	  ptemp = DBL_MIN;
	}
	prob += log(ptemp);
      }
      else if (pl==1){
    	for(n=0; n<datadim->nind; n++){
	  acopy = 0;
	  onePhi = gsl_matrix_get(param->phi, i, n);
	  for(k=0; k<nk; k++){
	    if(gsl_matrix_int_get(data->phcnt[i].mat, n, k) == 1){
	      obs[acopy] = k;
	      acopy++;
	    }
	  }
	  // p0 ancestry
	  if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
	    ptemp =  (1 - onePhi);
	    if (acopy == 1){
	      ptemp *=  gsl_matrix_get(param->pi0Cur, i, obs[0]);
	    }
	  }
	  else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
	    ptemp = onePhi;
	    if (acopy == 1){
	      ptemp *=  gsl_matrix_get(param->pi1Cur, i, obs[0]);
	    }
	  }
	}
	if (gsl_isnan(ptemp) == 1){
	  ptemp = DBL_MIN;
	}
	if (ptemp <= 0){
	  ptemp = DBL_MIN;
	}
	prob += log(ptemp);
      }
    }
  }
  return(prob);
}

/* calculates the quantile of each random effect based on the conditional prior and writes to file */
void calcQuants(paramcont * param, datadimcont * datadim, FILE * fpQgamma, FILE * fpQzeta, 
		FILE * fpQeta, FILE * fpQkappa){
  
  int i, j;
  double q, oneNu, oneOmega;

  // quants for gamma and zeta
  // first time no comma
  q = gsl_cdf_gaussian_P(gsl_matrix_get(param->gamma,0,0),sqrt(1/param->tauAlphaCur));
  fprintf(fpQgamma, "%.5f", q);
  q = gsl_cdf_gaussian_P(gsl_matrix_get(param->zeta,0,0),sqrt(1/param->tauBetaCur));
  fprintf(fpQzeta, "%.5f", q);
  for(i=1; i<datadim->nloci; i++){
    // cdf for gamma  and zeta based on the lower tail
    q = gsl_cdf_gaussian_P(gsl_matrix_get(param->gamma,i,0),sqrt(1/param->tauAlphaCur));
    fprintf(fpQgamma, ",%.5f", q);
    q = gsl_cdf_gaussian_P(gsl_matrix_get(param->zeta,i,0),sqrt(1/param->tauBetaCur));
    fprintf(fpQzeta, ",%.5f", q);
  }
  fprintf(fpQgamma, "\n");
  fprintf(fpQzeta, "\n");
  
  // quants for eta and kappa
  if(datadim->npop > 1){ // only if there is more than one popualtion
    for(i=0; i<datadim->nloci; i++){
      fprintf(fpQeta, "%i",i);
      fprintf(fpQkappa, "%i",i);
      oneNu = gsl_matrix_get(param->nu,i,0);
      oneOmega = gsl_matrix_get(param->omega,i,0);
      for(j=0; j<datadim->npop; j++){
	q = gsl_cdf_gaussian_P(gsl_matrix_get(param->etaCur,i,j),sqrt(1/oneNu));
	fprintf(fpQeta, ",%.5f",q);
	q = gsl_cdf_gaussian_P(gsl_matrix_get(param->kappaCur,i,j),sqrt(1/oneOmega));
	fprintf(fpQkappa, ",%.5f",q);
      }
      fprintf(fpQeta, "\n");
      fprintf(fpQkappa, "\n");
    }
  }
}

/* calcualtes and prints interspecific heterozygosity based on ancestry (z) */
void calcHet(mat_container_int * zCur, int nloci, int nind, gsl_vector_int * ploidyVec,
	     int allDiploid, FILE * fpHet, int step){
  int i, n;
  double tot, het;
  double phet;
  

  if (allDiploid == 1){
    n = 0;
    tot = 0;
    het = 0;
    for(i=0; i<nloci; i++){
      if (gsl_matrix_int_get(zCur[i].mat,n,0) != gsl_matrix_int_get(zCur[i].mat,n,1)){
	het++;
      }
      tot++;
    }
    phet = het / tot;
    fprintf(fpHet, "%.5f", phet);
    for(n=1; n<nind; n++){
      tot = 0;
      het = 0;
      for(i=0; i<nloci; i++){
	if (gsl_matrix_int_get(zCur[i].mat,n,0) != gsl_matrix_int_get(zCur[i].mat,n,1)){
	  het++;
	}
	tot++;
      }
      phet = het / tot;
      fprintf(fpHet, ",%.5f", phet);
    }
    fprintf(fpHet, "\n");
  }
  else {
    n = 0;
    tot = 0;
    het = 0;
    for(i=0; i<nloci; i++){
      if (gsl_vector_int_get(ploidyVec,i) == 2){
	if (gsl_matrix_int_get(zCur[i].mat,n,0) != gsl_matrix_int_get(zCur[i].mat,n,1)){
	  het++;
	}
	tot++;
      }
    }
    phet = het / tot;
    fprintf(fpHet, "%.5f", phet);
    for(n=1; n<nind; n++){
      tot = 0;
      het = 0;
      for(i=0; i<nloci; i++){
	if (gsl_vector_int_get(ploidyVec,i) == 2){
	  if (gsl_matrix_int_get(zCur[i].mat,n,0) != gsl_matrix_int_get(zCur[i].mat,n,1)){
	    het++;
	  }
	  tot++;
	}
      }
      phet = het / tot;
      fprintf(fpHet, ",%.5f", phet);
    }
    fprintf(fpHet, "\n");
  }
}
