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

#include <math.h>

#include "bgc.h"
#include "mvrandist.h"

using namespace std;

/* -------------------------------------------------------------- */
/*                   linkage model functions                      */
/* -------------------------------------------------------------- */

/* gets the number of chromosomes */
void getNchrom(string filename, int * nchrom){
  ifstream infile;
  int wordCtr;
  string line, oneword;
  // open file
  infile.open(filename.c_str());
  if (!infile){
    cerr << "Cannot open file " << filename << endl;
    exit(1);
  }

  *nchrom = 0; // set n chrom to zero

  // get data from file
  while (getline(infile, line)){ // read a line of input into line
      wordCtr = 0;
      istringstream stream(line);
      while (stream >> oneword){ // read a word at a time
	if (wordCtr == 1){ // chromosome number
	  if (atof(oneword.c_str()) > *nchrom){
	    *nchrom = atoi(oneword.c_str());
	  }
	}
	wordCtr++;
      }
  }
  infile.close();


}

/* reads in map locations and calculates a pair-wise distance matrix */
void getMap(string filename, gsl_matrix * proxMat, gsl_matrix * mapData, int nloci, double maxDist,
	    gsl_vector_int * ncloci, int nchrom){
  int i, j;
  int wordCtr, lineCtr;
  string line, oneword;
  ifstream infile;
  double prox;
  double curChrom = 1; 
  double chromNloci;

  // open file
  infile.open(filename.c_str());
  if (!infile){
    cerr << "Cannot open file " << filename << endl;
    exit(1);
  }

  chromNloci = 0;
  lineCtr = 0;

  // get data from file
  while (getline(infile, line)){ // read a line of input into line
      wordCtr = 0;
      istringstream stream(line);
      while (stream >> oneword){ // read a word at a time
	if (wordCtr == 1){ // chromosome number
	  gsl_matrix_set(mapData,lineCtr,0,atof(oneword.c_str()));
	  if (atof(oneword.c_str()) > nchrom){
	    nchrom = atoi(oneword.c_str());
	  }
	  if (atof(oneword.c_str()) != curChrom){
	    gsl_vector_int_set(ncloci, ((int) curChrom - 1), (int) chromNloci);
	    curChrom = atof(oneword.c_str());
	    chromNloci = 0;
	  }
	  chromNloci++;
	}
	else if (wordCtr == 2){ // map location
	  gsl_matrix_set(mapData,lineCtr,1,atof(oneword.c_str()));
	}
	wordCtr++;
      }	    
      lineCtr++;
  }
  gsl_vector_int_set(ncloci, ((int) curChrom - 1), (int) chromNloci);
  infile.close();

  // calculate proximity matrix
  for (i=0; i<nloci; i++){
    for (j=0; j<nloci; j++){
      if (gsl_matrix_get(mapData,i,0) != gsl_matrix_get(mapData,j,0)){ // on different chrom, set prox to 0
	gsl_matrix_set(proxMat,i,j,0.0);
      }
      else { // non-zero proximity
	prox = 1 - (fabs(gsl_matrix_get(mapData,i,1) - gsl_matrix_get(mapData,j,1)) / maxDist);
	if (prox > 1){
	  prox = 1;
	}
	else if (prox < 0){
	  prox = 0;
	}
	gsl_matrix_set(proxMat,i,j,prox);
      }
    }
  }

}

/* separates the proximity matrix into one matrix per chromosome, also saves first and last locus per chrom. */
void sepMap(gsl_matrix * proxMat, gsl_matrix * mapData, int nloci, gsl_matrix_int * firstLast, 
	    mat_container * chromProxMat){
  int i, j;
  int curChrom = 1;
  int rowCtr = 0;
  int colCtr = 0;

  gsl_matrix_int_set(firstLast,curChrom-1,0,0);

  for (i=0; i<nloci; i++){
    if (curChrom != (int) gsl_matrix_get(mapData,i,0)){ // new chromosome
      gsl_matrix_int_set(firstLast,curChrom-1,1,i-1);
      curChrom = (int) gsl_matrix_get(mapData,i,0);
      gsl_matrix_int_set(firstLast,curChrom-1,0,i);
      rowCtr = 0;
    }
    colCtr = 0;
    for (j=0; j<nloci; j++){
      if (gsl_matrix_get(mapData,i,0) == gsl_matrix_get(mapData,j,0)){ // on the same chrom
	gsl_matrix_set(chromProxMat[curChrom-1].mat,rowCtr,colCtr, gsl_matrix_get(proxMat,i,j));
	colCtr++;
      }
    }
    rowCtr++;
  }
  gsl_matrix_int_set(firstLast,curChrom-1,1,i-1);
}

/* function to write rho mcmc samples */
void printRho(double rhoCur, FILE * fpRho){
  
  fprintf(fpRho, "%.5f\n", rhoCur);
}

/* sets the weight matrix, sigmaAll, based on rho and proxMat */
void setSigma(gsl_matrix * sigmaAll, gsl_matrix * proxMat, double rho, int nloci){
  int i, j; // locus
  double weight;
  double c3 = 20; // constants to make wi+ ~ 1 set to 20, 0.6, and 0.25, adjust for locus number
  double c2 = 0.6;
  double c1 = 0.25 / nloci;

  for (i=0; i<nloci; i++){
    for (j=0; j<nloci; j++){
      if (i == j){
	gsl_matrix_set(sigmaAll,i,j,0); // w_ii = 0
      }
      weight = c1 + c2 * exp(-c3 * (1 - gsl_matrix_get(proxMat,i,j)));
      gsl_matrix_set(sigmaAll,i,j,weight);
    }
  }
}


/* function to update gamma for linkage model , the locus effect on cline center, uses metropolis random walk */
void updateGammaCAR(paramcont * param, datadimcont * datadim, auxcont * auxvar, int allDiploid){
  double logMratio, u;
  int i;
  double prop;
  double paramsum = 0;

  for (i=0; i<datadim->nloci; i++){
    prop = gsl_ran_gaussian(r, param->tuneGamma) + gsl_matrix_get(param->gamma, i, 1);
    gsl_matrix_set(param->gamma, i, 0, prop);
      
    logMratio = calcLogMratioGammaCAR(param,datadim,auxvar,i,allDiploid);
    u = log(gsl_ran_flat(r, 0, 1));
    if (logMratio < u){ // not keeping proposed vector value
      gsl_matrix_set(param->gamma,i,0, gsl_matrix_get(param->gamma,i,1));
    }
    paramsum += gsl_matrix_get(param->gamma,i,0);
  }
  paramsum = paramsum / datadim->nloci;
  for (i=0; i<datadim->nloci; i++){
    gsl_matrix_set(param->gamma,i,0, gsl_matrix_get(param->gamma,i,0) - paramsum);
  }
}

/* calculates metropolis ratio for gamma */
double calcLogMratioGammaCAR(paramcont * param, datadimcont * datadim, auxcont * auxvar,
			     int i, int allDiploid){
  double pNew = 0;
  double pOld = 0;
  int n, j, splice, l;
  double onePhi, y, mMax, mMin, logMratio;
  double oneAlpha;
  double mu = 0;
  double s_weight = 0;
  double prob;
  double f;

  // calculate s_weight
  for (l=0; l<datadim->nloci; l++){ // 
    if (i != l){
      s_weight +=  gsl_matrix_get(param->sigmaAll,i,l);
    }
  }
  //cerr << "gamma: " << i << ", " << s_weight << endl;

  gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
  
  j = gsl_vector_int_get(datadim->pops, 0);
  oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
  splice = needSplice(oneAlpha,gsl_matrix_get(param->betaCur,i,j),auxvar->zero2one1,
		      auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
  // for proposed value; loop through individuals
  for(n=0; n<datadim->nind; n++){
    if (j != gsl_vector_int_get(datadim->pops, n)){ // only need to recalculate alpha and splice if pop changes
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

    if (allDiploid == 1){    
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

  // conditional autoregressive prior
  for (l=0; l<datadim->nloci; l++){ // loop through loci to calculate weighted mu
    if (i != l){
      mu += gsl_matrix_get(param->gamma,l,0) * (gsl_matrix_get(param->sigmaAll,i,l) / s_weight);
    }
  }
  mu = mu * param->rhoCur;
  pNew += log(gsl_ran_gaussian_pdf(gsl_matrix_get(param->gamma,i,0) - mu, 
				   sqrt((1 / param->tauAlphaCur) / s_weight)));
  mu = 0;

  // for old value
  gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
  
  j = gsl_vector_int_get(datadim->pops, 0);
  splice = needSplice(gsl_matrix_get(param->alphaCur,i,j),gsl_matrix_get(param->betaCur,i,j),auxvar->zero2one1,
		      auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
  for(n=0; n<datadim->nind; n++){
    if (j != gsl_vector_int_get(datadim->pops, n)){ // only need to recalculate alpha and splice if pop changes
      j = gsl_vector_int_get(datadim->pops, n);
      gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
      gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
      splice = needSplice(gsl_matrix_get(param->alphaCur,i,j),gsl_matrix_get(param->betaCur,i,j),auxvar->zero2one1,
			  auxvar->zero2one2,auxvar->zero2one3,&y,&mMax,&mMin,auxvar->M);
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

  // conditional autoregressive prior
  for (l=0; l<datadim->nloci; l++){ // loop through loci to calculate weighted mu
    if (i != l){
      mu += gsl_matrix_get(param->gamma,l,0) * (gsl_matrix_get(param->sigmaAll,i,l) / s_weight);
    }
  }
  mu = mu * param->rhoCur;
  pOld += log(gsl_ran_gaussian_pdf(gsl_matrix_get(param->gamma,i,1) - mu, 
				   sqrt((1 / param->tauAlphaCur) / s_weight)));
  logMratio = pNew - pOld;

//   cerr << "gamma: " << gsl_matrix_get(gamma,i,0) << ", " << gsl_matrix_get(gamma,i,1) << "," << mu << endl;
//   cerr << "dif: " << logMratio << ", new: " << pNew << ", old: " << pOld << endl;


  return (logMratio);
}

/* function to update zeta, the locus effect on cline width, uses metropolis random walk */
void updateZetaCAR(paramcont * param, datadimcont * datadim, auxcont * auxvar, int allDiploid){
  double logMratio, u;
  int i;
  double prop;
  double paramsum = 0;

  for (i=0; i<datadim->nloci; i++){
    prop = gsl_ran_gaussian(r, param->tuneZeta) + gsl_matrix_get(param->zeta, i, 1);
    gsl_matrix_set(param->zeta, i, 0, prop);
      
    logMratio = calcLogMratioZetaCAR(param,datadim,auxvar,i,allDiploid);
    u = log(gsl_ran_flat(r, 0, 1));
      
    if (logMratio < u){ // not keeping proposed vector value
      gsl_matrix_set(param->zeta,i,0, gsl_matrix_get(param->zeta,i,1));
    }
    paramsum += gsl_matrix_get(param->zeta,i,0);
  }
  paramsum = paramsum / datadim->nloci;
  for (i=0; i<datadim->nloci; i++){
    gsl_matrix_set(param->zeta,i,0, gsl_matrix_get(param->zeta,i,0) - paramsum);
  }
}

/* calculates metropolis ratio for zeta */
double calcLogMratioZetaCAR(paramcont * param, datadimcont * datadim, auxcont * auxvar,
			    int i, int allDiploid){
  double pNew = 0;
  double pOld = 0;
  int n, j, splice, l;
  double onePhi, y, mMax, mMin, logMratio;
  double oneBeta, oneAlpha;
  double mu = 0;
  double s_weight = 0;
  double prob;
  double f;

  // calculate s_weight
  for (l=0; l<datadim->nloci; l++){ // 
    if (i != l){
      s_weight +=  gsl_matrix_get(param->sigmaAll,i,l);
    }
  }

  //cerr << "zeta: " << i << ", " << s_weight << endl;

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

  // conditional autoregressive prior
  for (l=0; l<datadim->nloci; l++){ // loop through loci to calculate weighted mu
    if (i != l){
      mu += gsl_matrix_get(param->zeta,l,0) * (gsl_matrix_get(param->sigmaAll,i,l) / s_weight);
    }
  }
  mu = mu * param->rhoCur;
  pNew += log(gsl_ran_gaussian_pdf(gsl_matrix_get(param->zeta,i,0) - mu, sqrt((1 / param->tauBetaCur) / s_weight)));
  mu = 0;

  // for old value
  gsl_vector_memcpy(auxvar->zero2one1,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one2,auxvar->zero2one);
  gsl_vector_memcpy(auxvar->zero2one3,auxvar->zero2one);
  
  j = gsl_vector_int_get(datadim->pops, 0);
  oneAlpha = gsl_matrix_get(param->gamma,i,0) + gsl_matrix_get(param->etaCur,i,j);
  oneBeta = gsl_matrix_get(param->zeta,i,1) + gsl_matrix_get(param->kappaCur,i,j);
  splice = needSplice(oneAlpha,oneBeta,auxvar->zero2one1,
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
      oneBeta = gsl_matrix_get(param->zeta,i,1) + gsl_matrix_get(param->kappaCur,i,j);
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

  // conditional autoregressive prior
  for (l=0; l<datadim->nloci; l++){ // loop through loci to calculate weighted mu
    if (i != l){
      mu += gsl_matrix_get(param->zeta,l,0) * (gsl_matrix_get(param->sigmaAll,i,l) / s_weight);
    }
  }
  mu = mu * param->rhoCur;
  pOld += log(gsl_ran_gaussian_pdf(gsl_matrix_get(param->zeta,i,1) - mu, 
				   sqrt((1 / param->tauBetaCur) / s_weight)));
  logMratio = pNew - pOld;
//   cerr << "zeta: " << gsl_matrix_get(zeta,i,0) << "," << gsl_matrix_get(zeta,i,1) << "," << mu << endl;
//   cerr << "new: " << pNew << ", old: " << pOld << ", dif: " << logMratio << endl;

  return (logMratio);
}

/* update tau precision parameters with CAR linkage model, uses metropolis random-walk */
void updateTauCAR(gsl_matrix * parameters, double * precisionCur, double precisionLast, gsl_matrix * sigmaAll, 
		  double rho, int nloci, double tuneTau){
  double dev, u;
  double pNew = 0;
  double pOld = 0;
  //double gampriora = 0.001;  // shape parameter
  //double gampriorb = 0.001;  // rate paramerer
  double unipriorlb = 0.01;  // unif prior lb
  double unipriorub = 10000;  // unif prior ub
  double logMratio;
  int i,l;
  double mu = 0;
  double s_weight = 0;
  double prob;

  dev = gsl_ran_gaussian(r, tuneTau);
  *precisionCur = precisionLast + dev;

  //if (*precisionCur > 0){
  if (*precisionCur < unipriorub && *precisionCur > unipriorlb){
    // calculate pOld and pNew
    for (i=0; i<nloci; i++){
      s_weight = 0;
      mu = 0;
      // calculate s_weight
      for (l=0; l<nloci; l++){ //
	if (i != l){
	  s_weight +=  gsl_matrix_get(sigmaAll,i,l);
	}
      }
      // loop through loci to calculate weighted mu
      for (l=0; l<nloci; l++){ 
	if (i != l){
	  mu += gsl_matrix_get(parameters,l,0) * (gsl_matrix_get(sigmaAll,i,l) / s_weight);
	}
      }
      mu = mu * rho;
      prob = gsl_ran_gaussian_pdf(gsl_matrix_get(parameters,i,0) - mu, sqrt((1 / *precisionCur) / s_weight));
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pNew += log(prob);
      prob = gsl_ran_gaussian_pdf(gsl_matrix_get(parameters,i,0) - mu, sqrt((1 / precisionLast) / s_weight));
      if (gsl_isnan(prob) == 1){
	prob = DBL_MIN;
      }
      if (prob <= 0){
	prob = DBL_MIN;
      }
      pOld += log(prob);
    }

    // gamma prior for precision
    // pNew += log(gsl_ran_gamma_pdf(*precisionCur, gampriora, 1/gampriorb));
    // pOld += log(gsl_ran_gamma_pdf(precisionLast, gampriora, 1/gampriorb));
    // uniform prior
    pNew += log(gsl_ran_flat_pdf(*precisionCur, unipriorlb, unipriorub));
    pOld += log(gsl_ran_flat_pdf(precisionLast, unipriorlb, unipriorub));

    logMratio = pNew - pOld;
  
//     cerr << *precisionCur << ", " << precisionLast << endl;
//     cerr << "dif: " << logMratio << ", new: " << pNew << ", old: " << pOld << endl;

    u = log(gsl_ran_flat(r, 0, 1));
  
    if (logMratio < u){ // not keeping proposed vector value
      *precisionCur = precisionLast;
    }
  }
  else { // p(tau < 0) = 0
    *precisionCur = precisionLast;
  }
}

/* update correlation parameter rho, uses metropolis random-walk */
void updateRhoCAR(paramcont * param, datacont * data, datadimcont * datadim){
  double dev, u;
  double pNew = 0;
  double pOld = 0;
  double unipriorlb = 0; // uniform prior [0,1] 
  double unipriorub = 1; 
  double logMratio;
  int i,l;
  double mu_g = 0;
  double mu_z = 0;
  double s_weight = 0;

  dev = gsl_ran_gaussian(r, param->tuneRho);
  param->rhoCur = param->rhoLast + dev;
  if (param->rhoCur > unipriorub || param->rhoCur < unipriorlb){
    param->rhoCur = param->rhoLast;
  }
  else { // p(*rhoCur) > 0, calc M-ratio
    // pOld first, it uses old sigmaAll
    for (i=0; i<datadim->nloci; i++){
      s_weight = 0;
      mu_g = 0;
      mu_z = 0;
      // calculate s_weight
      for (l=0; l<datadim->nloci; l++){ //
	if (i != l){
	  s_weight +=  gsl_matrix_get(param->sigmaAll,i,l);
	}
      }
      // loop through loci to calculate weighted mu
      for (l=0; l<datadim->nloci; l++){ 
	if (i != l){
	  mu_g += gsl_matrix_get(param->gamma,l,0) * (gsl_matrix_get(param->sigmaAll,i,l) / s_weight);
	  mu_z += gsl_matrix_get(param->zeta,l,0) * (gsl_matrix_get(param->sigmaAll,i,l) / s_weight);
	}
      }
      mu_g = mu_g * param->rhoLast;
      mu_z = mu_z * param->rhoLast;
      pOld += log(gsl_ran_gaussian_pdf(gsl_matrix_get(param->gamma,i,0) - mu_g, 
				       sqrt((1 / param->tauAlphaCur) / s_weight)));
      pOld += log(gsl_ran_gaussian_pdf(gsl_matrix_get(param->zeta,i,0) - mu_z, 
				       sqrt((1 / param->tauBetaCur) / s_weight)));
    }

    // next recalculate sigmaAll based on new rho NO LONGER NECESSARY
    //setSigma(sigmaAll,proxMat,*rhoCur,nloci);
    // now calculate pNew
    for (i=0; i<datadim->nloci; i++){
      s_weight = 0;
      mu_g = 0;
      mu_z = 0;
      // calculate s_weight
      for (l=0; l<datadim->nloci; l++){ //
	if (i != l){
	  s_weight +=  gsl_matrix_get(param->sigmaAll,i,l);
	}
      }
      // loop through loci to calculate weighted mu
      for (l=0; l<datadim->nloci; l++){ 
	if (i != l){
	  mu_g += gsl_matrix_get(param->gamma,l,0) * (gsl_matrix_get(param->sigmaAll,i,l) / s_weight);
	  mu_z += gsl_matrix_get(param->zeta,l,0) * (gsl_matrix_get(param->sigmaAll,i,l) / s_weight);
	}
      }
      mu_g = mu_g * param->rhoCur;
      mu_z = mu_z * param->rhoCur;
      pNew += log(gsl_ran_gaussian_pdf(gsl_matrix_get(param->gamma,i,0) - mu_g,
				       sqrt((1 / param->tauAlphaCur) / s_weight)));
      pNew += log(gsl_ran_gaussian_pdf(gsl_matrix_get(param->zeta,i,0) - mu_z, 
				       sqrt((1 / param->tauBetaCur) / s_weight)));
    }

    // uniform prior
    pNew += log(gsl_ran_flat_pdf(param->rhoCur, unipriorlb, unipriorub));
    pOld += log(gsl_ran_flat_pdf(param->rhoLast, unipriorlb, unipriorub));
    if (gsl_isinf(pNew)==1){
      cerr << "pNew for rho is inf" << endl;
      exit(1);
    }

    logMratio = pNew - pOld;
    u = log(gsl_ran_flat(r, 0, 1));
 
    //cerr << "new: " << pNew << "  old: " << pOld << "  dif: " << logMratio << endl;

    if (logMratio < u){ // not keeping proposed vector value
      param->rhoCur = param->rhoLast;
      setSigma(param->sigmaAll,data->proxMat,param->rhoCur,datadim->nloci);// reset sigmaAll
    }
  }
}

/* calculates the quantile of each random effect based on the conditional prior and writes to file */
void calcQuantsCAR(paramcont * param, datadimcont * datadim, FILE * fpQgamma, FILE * fpQzeta){
  
  int i, l;
  // global and local outlier
  double qg, ql;
  double mu = 0;
  double s_weight = 0;

  // quants for gamma and zeta, first is global, second is local
  for (i=0; i<datadim->nloci; i++){
    for (l=0; l<datadim->nloci; l++){ //
      if (i != l){
	s_weight +=  gsl_matrix_get(param->sigmaAll,i,l);
      }
    }
    // cdf for gamma based on the lower tail
    for (l=0; l<datadim->nloci; l++){ 
      if (i != l){
	mu += gsl_matrix_get(param->gamma,l,0) * (gsl_matrix_get(param->sigmaAll,i,l) / s_weight);
      }
    }
    qg = gsl_cdf_gaussian_P(gsl_matrix_get(param->gamma,i,0), sqrt((1/param->tauAlphaCur) / s_weight));
    ql = gsl_cdf_gaussian_P(gsl_matrix_get(param->gamma,i,0) - mu, sqrt((1/param->tauAlphaCur) / s_weight));
    fprintf(fpQgamma, "%d,%.5f,%.5f\n", i,qg,ql);
    mu = 0;
    // cdf for zeta based on the lower tail
    for (l=0; l<datadim->nloci; l++){ 
      if (i != l){
	mu += gsl_matrix_get(param->zeta,l,0) * (gsl_matrix_get(param->sigmaAll,i,l) / s_weight);
      }
    }
    qg = gsl_cdf_gaussian_P(gsl_matrix_get(param->zeta,i,0), sqrt((1/param->tauBetaCur) / s_weight));
    ql = gsl_cdf_gaussian_P(gsl_matrix_get(param->zeta,i,0) - mu, sqrt((1/param->tauBetaCur) / s_weight));
    fprintf(fpQzeta, "%d,%.5f,%.5f\n", i,qg,ql);
    mu = 0;
    s_weight = 0;
  }
}


void setH(int nind, string hfile, gsl_matrix * hi){
  int n = 0;
  string line, oneword;
  ifstream infile;

  infile.open(hfile.c_str());
  while (getline(infile, line)){ // read a line of input into line
    istringstream stream(line);
    while (stream >> oneword){
      gsl_matrix_set(hi,n,0,atof(oneword.c_str()));
      n++;
    }
  }
  infile.close();
  for (n=0; n<nind; n++){
    cerr << gsl_matrix_get(hi,n,0) << endl;
  }
}



