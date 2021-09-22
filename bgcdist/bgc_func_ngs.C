#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <float.h>

#include <math.h>

#include "bgc.h"

using namespace std;

/* -------------------------------------------------------------- */
/*                     functions for ngs data                     */
/* -------------------------------------------------------------- */

/* read in the count data for the admixed and parental populations */
void getDatangs(string hybridfile, string p0file, string p1file, datacont * data,
		datadimcont * datadim){
  int i = -1;
  int j = 0, n = 0, k = 0;
  string line, oneword;
  ifstream infile;

  // open file for hybrid count data
  infile.open(hybridfile.c_str());
  if (!infile){
    cerr << "Cannot open file " << hybridfile << endl;
    exit(1);
  }

  // read in data for hybrids
  while (getline(infile, line)){ // read a line of input into line
    if (MAR == line[0]){ // MAR is the character designating a locus
      i++; // increment locus
      if (data->error >= 1){ // locus-specific error information is included in this fiel, which follows locus number
	istringstream stream(line);
	stream >> oneword;
	stream >> oneword; // throw away "locus" and locus number"
	stream >> oneword;
	gsl_vector_set(data->errVec,i,atof(oneword.c_str())); 
      }
      else { // set all errors probs. equal to data->error
	gsl_vector_set(data->errVec,i,data->error); 
      }	


      j = -1;
      n = 0;
    }
    else if (POPMAR == line[0]){
      j++; // increment population
    }
    else { // this is a line with count data
      k = 0;
      gsl_vector_int_set(datadim->pops, n, j);// set population for this ind x locus
      istringstream stream(line);
      while (stream >> oneword){ // read a word at a time
	gsl_matrix_int_set(data->phcnt[i].mat, n, k, atoi(oneword.c_str())); // count data for allele x ind x locus
	k++;
      } 
      n++; // increment individual
    }
  }
  infile.close();

  // open file for parental population 0  count data
  infile.open(p0file.c_str());
  if (!infile){
    cerr << "Cannot open file " << p0file << endl;
    exit(1);
  }

  i = -1;

  // read in data for parental 0
  while (getline(infile, line)){ // read a line of input into line
    if (MAR == line[0]){ // MAR is the character designating a locus
      i++; // increment locus
      n = 0;
    }
    else { // this is a line with count data
      k = 0;
      istringstream stream(line);
      while (stream >> oneword){ // read a word at a time
	gsl_matrix_int_set(data->ngsp0cnt[i].mat, n, k, atoi(oneword.c_str())); // count data 
	k++;
      } 
      n++;
    }
  }
  infile.close();

  // open file for parental population 1 count data
  infile.open(p1file.c_str());
  if (!infile){
    cerr << "Cannot open file " << p1file << endl;
    exit(1);
  }

  i = -1;

  // read in data for parental 1
  while (getline(infile, line)){ // read a line of input into line
    if (MAR == line[0]){ // MAR is the character designating a locus
      i++; // increment locus
      n = 0;
    }
    else { // this is a line with count data
      k = 0;
      istringstream stream(line);
      while (stream >> oneword){ // read a word at a time
	gsl_matrix_int_set(data->ngsp1cnt[i].mat, n, k, atoi(oneword.c_str())); // count data 
	k++;
      } 
      n++;
    }
  }
  infile.close();
}

/* function to update z, ancestry of hybrids, uses gibbs sampling similar to structure, but with ngs data */
/* we analytically sum across the unkown genotype */
void updateZngs(paramcont * param, datacont * data, datadimcont * datadim, auxcont * auxvar){
  int i, n, k, l, x, nk;
  double p00 = 0, p01 = 0, p10 = 0, p11 = 0, ptot = 0;
  double pr[4]; // array of probabilities
  unsigned int sam[4]; // sampled ancestry for a locus
  int acopy; // allele copy 0 or 1
  int finite;
  double err; // error rate to specific nucleotides
  double probtemp;

  // loop through loci
  for(i=0; i<datadim->nloci; i++){
    data->error = gsl_vector_get(data->errVec, i);
    nk = gsl_vector_int_get(datadim->nallele, i);
    err = (double) data->error / (nk - 1);
    // loop through individuals
    for(n=0; n<datadim->nind; n++){
      acopy = 0;
      for(k=0; k<nk; k++){
	// the individual has this allele
	if(gsl_matrix_int_get(data->phcnt[i].mat, n, k) > 0){
	  acopy++;
	}
      }
      // no data for this locus, set prob. based on phi
      if(acopy==0){
	p00 = gsl_pow_2(1 - gsl_matrix_get(param->phi, i, n)); // (1-phi)^2
	p01 = (1 - gsl_matrix_get(param->phi, i, n)) * gsl_matrix_get(param->phi, i, n); // (1-phi)phi
	p11 = gsl_pow_2(gsl_matrix_get(param->phi, i, n)); // phi^2	  
      }
      // there are reads for this locus
      else {
	p00 = 0;
	p01 = 0;
	p11 = 0;
	// copy count data to an auxillary vector
	for(k=0; k<nk; k++){
	  gsl_vector_uint_set(auxvar->vecAllele_int, k,
			      gsl_matrix_int_get(data->phcnt[i].mat, n, k));
	}
	// loop through all allele combinations
	for(k=0; k<nk; k++){
	  for(l=k; l<nk; l++){
	    // this is a homozygous genotype
	    if (k == l){
	      // set multinomial probs. for homozygous 
	      for (x=0; x<nk; x++){
		if (k == x){
		  gsl_vector_set(auxvar->vecAllelen1, x, (1.0 - data->error));
		}
		else {
		  gsl_vector_set(auxvar->vecAllelen1, x, err);
		}
	      } 
	      // probability of data given the genotype, i.e., p(d | g)
	      probtemp = gsl_ran_multinomial_pdf(nk, auxvar->vecAllelen1->data, auxvar->vecAllele_int->data);
	      // multiply by probability of genotype given ancestry
	      // p(d | g) * p(g | z), note k = l
	      p00 += probtemp * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l);
	      p01 += probtemp * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
	      p11 += probtemp * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
	    }
	    // this is a heterozygous genotype
	    else { 
	      // set multinomial probs. for homozygous 
	      for (x=0; x<nk; x++){
		if (x == k || x == l){
		  gsl_vector_set(auxvar->vecAllelen1, x, (0.5 - data->error + 0.5 * err));
		}
		else {
		  gsl_vector_set(auxvar->vecAllelen1, x, err);
		}
	      } 
	      // probability of data given the genotype, i.e., p(d | g)
	      probtemp = gsl_ran_multinomial_pdf(nk, auxvar->vecAllelen1->data, auxvar->vecAllele_int->data);
	      // multiply by probability of genotype given ancestry
	      // p(d | g) * p(g | z, pi), note k = l
	      p00 += probtemp * 2 * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l);
	      p01 += probtemp * (gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l) +
				 gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l));
	      p11 += probtemp * 2 * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
	    }
	  }
	}
      	// calculate p(z|phi) and mulitply by p(x|g) p(g|z,pi)
	p00 *= gsl_pow_2(1 - gsl_matrix_get(param->phi, i, n));
	p01 *= (1 - gsl_matrix_get(param->phi, i, n)) * gsl_matrix_get(param->phi, i, n);
	p11 *= gsl_pow_2(gsl_matrix_get(param->phi, i, n));
      }
      //cerr << "ancestry: " << i << " " << n << " " << p00 << "," << p01 << "," << p11 << endl;
      finite = gsl_isinf(p00);
      if (finite == 1){
	p00 = DBL_MAX;
      }
      else if (finite == -1){
	p00 = DBL_MIN;
      }      
      finite = gsl_isinf(p01);
      if (finite == 1){
	p01 = DBL_MAX;
      }
      else if (finite == -1){
	p01 = DBL_MIN;
      }      
      finite = gsl_isinf(p11);
      if (finite == 1){
	p11 = DBL_MAX;
      }
      else if (finite == -1){
	p11 = DBL_MIN;
      }

      p10 = p01;
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

/* initialize ancestry, this only works for ngs */
void initZngs(mat_container_int * zCur, mat_container_int * phcnt, gsl_matrix * pi0Cur,
	      gsl_matrix * pi1Cur, gsl_vector_int * nallele, int nloci, int nind){
  int i, n, k, nk, count = 0;
  double p00 = 0, p01 = 0, p10 = 0, p11 = 0, ptot = 0;
  double pr[4]; // array of probabilities
  unsigned int sam[4]; // sampled ancestry for a locus
  int obs[2]; // obs alleles
  int acopy; // allele copy 0 or 1
  int finite;
  int sample = 0; // bool, need to sample z

  // loop through loci
  for(i=0; i<nloci; i++){
    nk = gsl_vector_int_get(nallele, i);
    // loop through individuals
    for(n=0; n<nind; n++){
      acopy = 0;
      count = 0;
      for(k=0; k<nk; k++){
	// the individual has this allele
	if(gsl_matrix_int_get(phcnt[i].mat, n, k) > 0){
	  obs[acopy] = k;
	  // number of times obs, only relevant if one allele obs.
	  count = gsl_matrix_int_get(phcnt[i].mat, n, k); 
	  acopy++;
	}
      }
      // no data for this locus, sample with p = 0.5
      if(acopy==0){
	gsl_matrix_int_set(zCur[i].mat, n, acopy, gsl_ran_bernoulli(r, 0.5));
	acopy++;
	gsl_matrix_int_set(zCur[i].mat, n, acopy, gsl_ran_bernoulli(r, 0.5));
	acopy++;
	sample = 0;
      }
      // two alleles were observed the genotype is known
      else if(acopy==2){
	// calclate p(x|z,pi) * p(z)
	p00 =  2 * gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi0Cur, i, obs[0])) * gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi0Cur, i, obs[1]));
	p11 = 2 * gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi1Cur, i, obs[0])) * gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi1Cur, i, obs[1]));
	p01 = (gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi0Cur, i, obs[0])) * gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi1Cur, i, obs[1])) + 
	       gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi1Cur, i, obs[0])) * gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi0Cur, i, obs[1])));
	sample = 1;
      }

      // one allele was observed, genotype is not known
      // loop through possible genotypes
      else{
	p00 = 0;
	p01 = 0;
	p11 = 0;
	// calculate sum_g p(x|g) p(g|z,pi)
	for(k=0; k<nk; k++){
	  // this is a potential homozygous genetopye
	  if(k == obs[0]){
	    // calclate p(x|g) = 1 * p(g|pi,z)
	    p00 += gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi0Cur, i, obs[0])) * gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi0Cur, i, k));
	    p01 += gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi0Cur, i, obs[0])) * gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi1Cur, i, k));
	    p11 += gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi1Cur, i, obs[0])) * gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi1Cur, i, k));
	  }
	  // this is a heterozygous genotype
	  else{
	    // calclate calculate p(x|g) p(g|z,pi)
	    p00 += gsl_pow_int(0.5,count) * 2 * gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi0Cur, i, obs[0])) *
	      gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi0Cur, i, k));
	    p11 += gsl_pow_int(0.5,count) * 2 * gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi1Cur, i, obs[0])) * 
	      gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi1Cur, i, k));
	    p01 += gsl_pow_int(0.5,count) * (gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi0Cur, i, obs[0])) *
					     gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi1Cur, i, k)) + 
					     gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi1Cur, i, obs[0])) *
					     gsl_ran_bernoulli_pdf(1, gsl_matrix_get(pi0Cur, i, k)));

	  }
	}
	sample = 1;
      }
      //cerr << "ancestry: " << i << " " << n << " " << p00 << "," << p01 << "," << p11 << endl;
      if (sample == 1){ // still need to sample z
	finite = gsl_isinf(p00);
	if (finite == 1){
	  p00 = DBL_MAX;
	}
	else if (finite == -1){
	  p00 = DBL_MIN;
	}      
	finite = gsl_isinf(p01);
	if (finite == 1){
	  p01 = DBL_MAX;
	}
	else if (finite == -1){
	  p01 = DBL_MIN;
	}      
	finite = gsl_isinf(p11);
	if (finite == 1){
	  p11 = DBL_MAX;
	}
	else if (finite == -1){
	  p11 = DBL_MIN;
	}

	p10 = p01;
	ptot = p00 + p01 + p10 + p11;
	pr[0] = p00 / ptot;
	pr[1] = p01 / ptot;
	pr[2] = p10 / ptot;
	pr[3] = p11 / ptot;
	// sample ancestry pair based on probabilities
	gsl_ran_multinomial(r, 4, 1, pr, sam);
	if(sam[0] == 1){ // homo pop 0
	  gsl_matrix_int_set(zCur[i].mat, n, 0, 0);
	  gsl_matrix_int_set(zCur[i].mat, n, 1, 0);
	}
	else if(sam[1] == 1){ // het
	  gsl_matrix_int_set(zCur[i].mat, n, 0, 0);
	  gsl_matrix_int_set(zCur[i].mat, n, 1, 1);
	}
	else if(sam[2] == 1){ // het
	  gsl_matrix_int_set(zCur[i].mat, n, 0, 1);
	  gsl_matrix_int_set(zCur[i].mat, n, 1, 0);
	}
	else{ // homo pop 1
	  gsl_matrix_int_set(zCur[i].mat, n, 0, 1);
	  gsl_matrix_int_set(zCur[i].mat, n, 1, 1);
	}
      }
    }
  }
}

/* for ngs, p0cnt and p1cnt are the unobserved genotypes, which we initilize based on the observed data */
void initPcnts(datacont * data, datadimcont * datadim, gsl_vector * pvec, 
	       gsl_vector_uint * sam){
  int i, k, n, nk;
  int nobs = 0;
  double pr;
  int count = 0;
  int obs = 0;

  for(i=0; i<datadim->nloci; i++){
    nk = gsl_vector_int_get(datadim->nallele, i); // number of alleles for this locus
    // parental population 0
    for(n=0; n<datadim->nindP0; n++){
      nobs = 0;
      for(k=0; k<nk; k++){
	// must have this allele, increment by 1
	if(gsl_matrix_int_get(data->ngsp0cnt[i].mat, n, k) > 0){ 
	  gsl_matrix_int_set(data->p0cnt, i, k, (gsl_matrix_int_get(data->p0cnt, i, k) + 1));
	  nobs++;
	  obs = k;
	  count = gsl_matrix_int_get(data->ngsp0cnt[i].mat, n, k);
	}
      }
      // if no alleles were observed
      if (nobs == 0){
	pr = double (1 / double(nk));
	for (k=0; k<nk; k++){
	  gsl_vector_set(pvec, k, pr);
	}
	gsl_ran_multinomial(r, nk, 2, pvec->data, sam->data); // sample 2 alleles
	for (k=0; k<nk; k++){
	  if (gsl_vector_uint_get(sam, k) == 1){
	    gsl_matrix_int_set(data->p0cnt, i, k, (gsl_matrix_int_get(data->p0cnt, i, k) + 1));
	  }
	}
      }
      else if (nobs == 1){
	for (k=0; k<nk; k++){
	  if (k == obs){ // potential homozygote
	    pr = 1;
	  }
	  else {
	    pr = gsl_pow_int(0.5, count);
	  }
	  gsl_vector_set(pvec, k, pr);
	}
	gsl_ran_multinomial(r, nk, 1, pvec->data, sam->data);
	for (k=0; k<nk; k++){
	  if (gsl_vector_uint_get(sam, k) == 1){
	    gsl_matrix_int_set(data->p0cnt, i, k, (gsl_matrix_int_get(data->p0cnt, i, k) + 1));
	  }
	}
      }
    }
    // parental population 1
    for(n=0; n<datadim->nindP1; n++){
      nobs = 0;
      for(k=0; k<nk; k++){
	// must have this allele, increment by 1
	if(gsl_matrix_int_get(data->ngsp1cnt[i].mat, n, k) > 0){ 
	  gsl_matrix_int_set(data->p1cnt, i, k, (gsl_matrix_int_get(data->p1cnt, i, k) + 1));
	  nobs++;
	  obs = k;
	  count = gsl_matrix_int_get(data->ngsp1cnt[i].mat, n, k);
	}	
      }
      // if no alleles were observed
      if (nobs == 0){
	pr = double (1 / double(nk));
	for (k=0; k<nk; k++){
	  gsl_vector_set(pvec, k, pr);
	}
	gsl_ran_multinomial(r, nk, 2, pvec->data, sam->data); // sample 2 alleles
	for (k=0; k<nk; k++){
	  if (gsl_vector_uint_get(sam, k) == 1){
	    gsl_matrix_int_set(data->p1cnt, i, k, (gsl_matrix_int_get(data->p1cnt, i, k) + 1));
	  }
	}
     }
      else if (nobs == 1){
	for (k=0; k<nk; k++){
	  if (k == obs){ // potential homozygote
	    pr = 1;
	  }
	  else {
	    pr = gsl_pow_int(0.5, count);
	  }
	  gsl_vector_set(pvec, k, pr);
	}
	gsl_ran_multinomial(r, nk, 1, pvec->data, sam->data);
	for (k=0; k<nk; k++){
	  if (gsl_vector_uint_get(sam, k) == 1){
	    gsl_matrix_int_set(data->p1cnt, i, k, (gsl_matrix_int_get(data->p1cnt, i, k) + 1));
	  }
	}
      }
    }
  }
}

/* gibbs sample update of parental population counts */
void updatePcnt(datacont * data, auxcont * auxvar, datadimcont * datadim, gsl_matrix_int * pcnt, 
		mat_container_int * ngsPcnt, int nind, gsl_matrix * pi){
  int i, k, l, x, n, nk;
  int nobs = 0;
  unsigned int dip = 2; // diploid
  double err, probtemp;

  // set pcnt to zero
  gsl_matrix_int_set_zero(pcnt);

  for(i=0; i<datadim->nloci; i++){
    data->error = gsl_vector_get(data->errVec, i);
    nk = gsl_vector_int_get(datadim->nallele,i); // number of alleles for this locus
    err = (double) data->error / (nk - 1);
    for(n=0; n<nind; n++){
       nobs = 0;
      // determine wether there is data for this individual
      for(k=0; k<nk; k++){
	if(gsl_matrix_int_get(ngsPcnt[i].mat, n, k) > 0){ // data present
	  nobs++;
	}	
       }
      // no data for this individual
      if(nobs == 0){
	// update based on population allele frequences
	for(k=0; k<nk; k++){
	  gsl_vector_set(auxvar->vecAllelen1, k, gsl_matrix_get(pi, i, k));
	}
	// sample
	gsl_ran_multinomial(r, nk, dip, auxvar->vecAllelen1->data, auxvar->vecAllele_int->data);
	// record values
	for(k=0; k<nk; k++){
	  gsl_matrix_int_set(pcnt, i, k, (gsl_matrix_int_get(pcnt, i, k) + 
					  gsl_vector_uint_get(auxvar->vecAllele_int, k)));
	}
      }
      // data present for this individual
      else {
	// loop through all allele combinations
	for(k=0; k<nk; k++){
	  gsl_vector_uint_set(auxvar->vecAllele_int, k, gsl_matrix_int_get(ngsPcnt[i].mat, n, k));
	}
	for(k=0; k<nk; k++){
	  for(l=0; l<nk; l++){
	    // this is a homozygous genotype
	    if (k == l){
	      // set multinomial probs. for homozygous 
	      for (x=0; x<nk; x++){
		if (k == x){
		  gsl_vector_set(auxvar->vecAllelen1, x, (1.0 - data->error));
		}
		else {
		  gsl_vector_set(auxvar->vecAllelen1, x, err);
		}
	      } 	
      	      // p(x | g) * p(g | pi)
	      probtemp = gsl_ran_multinomial_pdf(nk, auxvar->vecAllelen1->data, auxvar->vecAllele_int->data) * 
		gsl_pow_2(gsl_matrix_get(pi, i, k));
	      gsl_vector_set(auxvar->vecnkbynk, (nk * k + l), probtemp);
	    }
	    // this is a heterozygous genotype
	    else { 
	      // set multinomial probs. for homozygous 
	      for (x=0; x<nk; x++){
		if (x == k || x == l){
		  gsl_vector_set(auxvar->vecAllelen1, x, (0.5 - data->error + 0.5 * err));
		}
		else {
		  gsl_vector_set(auxvar->vecAllelen1, x, err);
		}
	      }
     	      // p(x | g) * p(g | pi), just pq not 2pq, bucause we are storing both hets.
	      probtemp = gsl_ran_multinomial_pdf(nk, auxvar->vecAllelen1->data, auxvar->vecAllele_int->data) * 
		gsl_matrix_get(pi, i, k) * gsl_matrix_get(pi, i, l);
	      gsl_vector_set(auxvar->vecnkbynk, (nk * k + l), probtemp);
	    }
	  }
	}
	gsl_vector_uint_set_zero(auxvar->vecnkbynk_int);
	gsl_ran_multinomial(r, (nk*nk), 1, auxvar->vecnkbynk->data, auxvar->vecnkbynk_int->data);
	for (k=0; k<nk; k++){
	  for (l=0; l<nk; l++){
	    // this is the sampled genotype
	    if (gsl_vector_uint_get(auxvar->vecnkbynk_int, (nk * k + l)) == 1){
	      gsl_matrix_int_set(pcnt, i, k, (gsl_matrix_int_get(pcnt, i, k) + 1));
	      gsl_matrix_int_set(pcnt, i, l, (gsl_matrix_int_get(pcnt, i, l) + 1));
	    }
	  }
	}
      }
    }
  }
}


/* calculate Ln Likelihood of the data from the admixed populations = p(x | z, pi) = \prod_i,n sum_g p(x | g) p(g | pi, z) */
double calcLnLngs(paramcont * param, datacont * data, datadimcont * datadim, auxcont * auxvar){

  int i, n, k, l, x, nk;
  int acopy = 0;
  double probsum, probtemp, prob = 0;
  double err;
  
  // loop through loci
  for(i=0; i<datadim->nloci; i++){
    data->error = gsl_vector_get(data->errVec, i);
    nk = gsl_vector_int_get(datadim->nallele, i);
    err = (double) data->error / (nk - 1);
    // loop through individuals
    for(n=0; n<datadim->nind; n++){
      probsum = 0;
      acopy = 0;
      for(k=0; k<nk; k++){
	gsl_vector_uint_set(auxvar->vecAllele_int, k,
			    gsl_matrix_int_get(data->phcnt[i].mat, n, k));
	// the individual has this allele
	if(gsl_matrix_int_get(data->phcnt[i].mat, n, k) > 0){
	  acopy++;
	}
      }
      // if acopy == 0 no data for this locus, do not include in likelihood calculation
      if (acopy > 0){
	for(k=0; k<nk; k++){
	  for(l=k; l<nk; l++){
	    // this is a homozygous genotype
	    if (k == l){
	      // set multinomial probs. for homozygous 
	      for (x=0; x<nk; x++){
		if (k == x){
		  gsl_vector_set(auxvar->vecAllelen1, x, (1.0 - data->error));
		}
		else {
		  gsl_vector_set(auxvar->vecAllelen1, x, err);
		}
	      } 
	      // probability of data given the genotype, i.e., p(x | g)
	      probtemp = gsl_ran_multinomial_pdf(nk, auxvar->vecAllelen1->data, auxvar->vecAllele_int->data);
	      // probability of genotype | ancestry
	      // ancestry only pop. 0
	      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
		  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
		probtemp = probtemp * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l);
	      }
	      // ancestry only pop. 1	      
	      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
		       gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
		probtemp = probtemp * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
	      }
	      // ancestry from each population 
	      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) != gsl_matrix_int_get(param->zCur[i].mat, n, 1)){	      
		probtemp = probtemp * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
	      }
	      if (probtemp == 0){
		probtemp = DBL_MIN;
	      }
	      probsum += probtemp;
	    }
	    // this is a heterozygous genotype
	    else { 
	      // set multinomial probs. for homozygous 
	      for (x=0; x<nk; x++){
		if (x == k || x == l){
		  gsl_vector_set(auxvar->vecAllelen1, x, (0.5 - data->error + 0.5 * err));
		}
		else {
		  gsl_vector_set(auxvar->vecAllelen1, x, err);
		}
	      } 
	      // probability of data given the genotype, i.e., p(x | g)
	      probtemp = gsl_ran_multinomial_pdf(nk, auxvar->vecAllelen1->data, auxvar->vecAllele_int->data);
	      // probability of genotype | ancestry
	      // ancestry only pop. 0
	      if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
		  gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 0){
		probtemp = probtemp * 2 * gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l);
	      }
	      // ancestry only pop. 1	      
	      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) == gsl_matrix_int_get(param->zCur[i].mat, n, 1) &&
		       gsl_matrix_int_get(param->zCur[i].mat, n, 0) == 1){
		probtemp = probtemp * 2 * gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l);
	      }
	      // ancestry from each population 
	      else if (gsl_matrix_int_get(param->zCur[i].mat, n, 0) != gsl_matrix_int_get(param->zCur[i].mat, n, 1)){	      
		probtemp =  probtemp * (gsl_matrix_get(param->pi0Cur, i, k) * gsl_matrix_get(param->pi1Cur, i, l) +
					gsl_matrix_get(param->pi1Cur, i, k) * gsl_matrix_get(param->pi0Cur, i, l));
	      }
	      if (probtemp == 0){
		probtemp = DBL_MIN;
	      }
	      probsum += probtemp;
	    }
	  }
	}
	prob += log(probsum);
      }
    }
  }
  return(prob);
}
