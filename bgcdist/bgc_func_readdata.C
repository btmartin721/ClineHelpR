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

/* prints out the proper usage of the program */
void usage(char * name){
  fprintf(stderr,"\n%s version %s\n\n", name, VERSION); 
  fprintf(stderr, "Usage:   bgc -a p0infile -b p1infile -h admixinfile [options]\n");
  fprintf(stderr, "-a     Infile with genetic data for parental population 0\n");
  fprintf(stderr, "-b     Infile with genetic data for parental population 1\n");
  fprintf(stderr, "-h     Infile with genetic data for admixed population(s)\n");
  fprintf(stderr, "-M     Infile with genetic map, only used for ICARrho model\n");
  fprintf(stderr, "-F     Prefix added to all outfiles\n");
  fprintf(stderr, "-O     Format to write MCMC samples [default = 0]\n");
  fprintf(stderr, "       0 = HDF5\n");
  fprintf(stderr, "       1 = ascii text\n");
  fprintf(stderr, "       2 = HDF5 and ascii text\n");
  fprintf(stderr, "-x     Number of MCMC steps for the analysis [default = 1000]\n");
  fprintf(stderr, "-n     Discard the first n MCMC samples as a burn-in [default = 0]\n");
  fprintf(stderr, "-t     Thin MCMC samples by recording every nth value [default = 1]\n");
  fprintf(stderr, "-p     Species which parameter samples to print [default = 0]\n");
  fprintf(stderr, "       0 = print log likelihood, alpha, beta, and hybrid index\n");
  fprintf(stderr, "       1 = also print precision parameters\n");
  fprintf(stderr, "       2 = also print eta and kappa\n");
  fprintf(stderr, "-q     Boolean, calculate and print cline parameter quantiles [default = 0]\n");
  fprintf(stderr, "-i     Boolean, calculate and print interspecific-heterozygosity [default = 0]\n");
  fprintf(stderr, "-N     Boolean, use genotype-uncertainty model [default = 0]\n");
  fprintf(stderr, "-E     Boolean, use sequence error model, only valid in conjunction with the\n");
  fprintf(stderr, "       genotype-uncertainty model [default = 0]\n");
  fprintf(stderr, "-m     Boolean, use ICARrho model for linked loci [default = 0]\n");
  fprintf(stderr, "-d     Boolean, all loci are diploid [default = 1]\n");
  fprintf(stderr, "-s     Boolean, sum-to-zero constraint on locus cline parameters [default = 1]\n");
  fprintf(stderr, "-o     Boolean, assume a constant population-level cline parameter variance\n");
  fprintf(stderr, "       for all loci [default = 0]\n");
  fprintf(stderr, "-I     Select algorithm for to initialize MCMC [default = 1]\n");
  fprintf(stderr, "-D     Maximum distance between loci, free recombination [default = 0.5]\n");
  fprintf(stderr, "-T     If non-zero, use a truncated gamma prior for tau with this upper bound [default = 0]\n");
  fprintf(stderr, "-u     MCMC tunning paramter, maximum deviate from uniform for proposed hybrid index\n");
  fprintf(stderr, "       hybrid index [default = 0.1]\n");
  fprintf(stderr, "-g     MCMC tunning paramter, standard deviation for Gaussian proposal of cline\n");
  fprintf(stderr, "       parameter gamma [default = 0.05]\n");
  fprintf(stderr, "-z     MCMC tunning paramter, standard deviation for Gaussian proposal of cline\n");
  fprintf(stderr, "       parameter zeta [default = 0.05]\n");
  fprintf(stderr, "-e     MCMC tunning paramter, standard deviation for Gaussian proposal of cline\n");
  fprintf(stderr, "       parameters eta and kappa [default = 0.02]\n");
  fprintf(stderr, "-v     Display software version\n");
  exit(1);
}

/*------------------------------------------------------*/
/*         data reading functions                       */
/*------------------------------------------------------*/

/* determine the number of loci based on hybridfile */
int getNloci(string filename){
  string line;
  int lociCtr = 0;
  ifstream infile;

  // open file
  infile.open(filename.c_str());
  if (!infile){
    cerr << "Cannot open file " << filename << endl;
    exit(1);
  }
  // get number of loci
  while (getline(infile, line)){ // read a line of input into line
    if ( MAR == line[0]){ // MAR is the character designating a locus
      lociCtr++;
    }
  }
  infile.close();
  return(lociCtr);
}

/* determine and record the number of alleles per locus, based on hybridfile */
int getNallele(string filename, gsl_vector_int * nallele){
  string line, oneword;
  int lociCtr = -1;
  int wordCtr = 0;
  int max = 0;
  int first = 1;
  ifstream infile;

  // open file
  infile.open(filename.c_str());
  if (!infile){
    cerr << "Cannot open file " << filename << endl;
    exit(1);
  }

  // determine number of alleles per locus
  while (getline(infile, line)){ // read a line of input into line
    if (MAR == line[0]){ // MAR is the character designating a locus
      lociCtr++;
      first = 1;
    }
    else if (POPMAR == line[0]){
      ;
    }
    else if (first == 1){ // this is a line with count data
      first = 0;
      wordCtr = 0;
      istringstream stream(line);
      while (stream >> oneword){ // read a word at a time
	  wordCtr++;
      }
      gsl_vector_int_set(nallele, lociCtr, wordCtr);
      if (wordCtr > max){
	max = wordCtr;
      }
    }
  }
  infile.close();
  return(max);
}

/* determine and record the number of populations, based on hybridfile */
int getNpop(string filename){
  string line;
  int popCtr = 0;
  int first = 1;
  ifstream infile;

  // open file
  infile.open(filename.c_str());
  if (!infile){
    cerr << "Cannot open file " << filename << endl;
    exit(1);
  }

  // determine number of alleles per locus
  while (getline(infile, line)){ // read a line of input into line
    if (MAR == line[0]){ // MAR is the character designating a locus
      if (first == 0){
	break;
      }
      first = 0;
    }
    else if (POPMAR == line[0]){
      popCtr++;
    }
  }
  infile.close();
  return(popCtr);
}

/* determine and record the number of individuals per locus, based on hybridfile */
int getNind(string filename){
  string line;
  int lociCtr = -1;
  int indCtr = 0;
  ifstream infile;

  // open file
  infile.open(filename.c_str());
  if (!infile){
    cerr << "Cannot open file " << filename << endl;
    exit(1);
  }

  // determine number of alleles per locus
  while (getline(infile, line)){ // read a line of input into line
    if (MAR == line[0]){ // MAR is the character designating a locus
      if (lociCtr > -1){
	break;
      }
      indCtr = 0;
      lociCtr++;
    }
    else if (POPMAR == line[0]){
      ;
    }
    else { // this is a line with count data
      indCtr++;
    }
  }
  infile.close();
  return(indCtr);
}

/* read in the count data for the admixed and parental populations */
void getData(string hybridfile,string p0file,string p1file, datacont *data,
	     datadimcont * datadim, int allDiploid){
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
      if (allDiploid != 1){ // need to get ploidy for locus, which follows locus number
	istringstream stream(line);
	stream >> oneword;
	stream >> oneword; // throw away "locus" and locus number"
	stream >> oneword;
	gsl_vector_int_set(datadim->ploidyVec,i,atoi(oneword.c_str())); // set ploidy number
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
    cerr << "Cannot open file " << hybridfile << endl;
    exit(1);
  }

  i = -1;

  // read in data for parental 0
  while (getline(infile, line)){ // read a line of input into line
    if (MAR == line[0]){ // MAR is the character designating a locus
      i++; // increment locus
    }
    else { // this is a line with count data
      k = 0;
      istringstream stream(line);
      while (stream >> oneword){ // read a word at a time
	gsl_matrix_int_set(data->p0cnt, i, k, atoi(oneword.c_str())); // count data 
	k++;
      } 
    }
  }
  infile.close();

  // open file for parental population 1 count data
  infile.open(p1file.c_str());
  if (!infile){
    cerr << "Cannot open file " << hybridfile << endl;
    exit(1);
  }

  i = -1;

  // read in data for parental 1
  while (getline(infile, line)){ // read a line of input into line
    if (MAR == line[0]){ // MAR is the character designating a locus
      i++; // increment locus
    }
    else { // this is a line with count data
      k = 0;
      istringstream stream(line);
      while (stream >> oneword){ // read a word at a time
	gsl_matrix_int_set(data->p1cnt, i, k, atoi(oneword.c_str())); // count data 
	k++;
      } 
    }
  }
  infile.close();
}
