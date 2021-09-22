#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "hdf5.h"

// macro definitions
#define VERSION "1.03 -- 03 April 2014"
#define MAR 'l' //character that identifies a locus in the data
#define POPMAR 'p' //character that identifies population data 
#define PLOIDY 2 // number of allele copies per locus, may not work when != 2
#define INITSD 0.01 // standard deviation for inintialization of locus random effects
#define INITSDP 0.01 // standard deviation for inintialization of population random effects
#define INITRHO 3.0 // shape for gamma proposal of rho
#define HPROP 0.1 // bounds of deviate from unifrom for hi proposal, DEFAULT only
#define TGAMMA 0.05 // sd of gamma proposals, DEFAULT only
#define TZETA 0.05 // sd of zeta proposals, DEFAULT only
#define TETA 0.02 // sd of eta proposals, also used for kappa, DEFAULT only
#define TTAU 2 // sd of tau proposals, DEFAULT only
#define TRHO 0.1 // sd of tau proposals, DEFAULT only
#define TFIS 0.02 // bounds of deviate from uniform for fis proposal, DEFAULT ONLY
extern gsl_rng * r;

// struct definitions
struct mat_vec_container_int { 
  gsl_matrix_int * mat; 
  gsl_vector_int * vec; 
}; 

struct mat_container_int { 
  gsl_matrix_int * mat; 
}; 

struct mat_container { 
  gsl_matrix * mat; 
}; 

struct vec_container { 
  gsl_vector * vec; 
}; 

struct datacont{
  // count variables
  gsl_matrix_int * p0cnt;
  gsl_matrix_int * p1cnt;
  mat_container_int * phcnt;

  // variables for ngs parental populations, input file is different
  mat_container_int * ngsp0cnt;
  mat_container_int * ngsp1cnt;

  // map location variables
  gsl_matrix * proxMat;
  mat_container * chromProxMat;
  gsl_matrix * mapData; // column 1 is chrom. number, column 2 is location
  gsl_matrix_int * firstLast;

  // sequence error
  double error; // fixed error rate, set to 9 for variable error rate
  gsl_vector * errVec; // locus error rate
};

struct paramcont{

  // variables to store MCMC output
  gsl_matrix * pi0Cur; // parental allele frequences
  gsl_matrix * pi1Cur; 
  gsl_matrix * pi0Last; // parental allele frequences
  gsl_matrix * pi1Last; 
  mat_container_int * zCur; // ancestry
  mat_container_int * zLast; // ancestry
  gsl_matrix * hi; // hybrid index, current and last
  gsl_matrix * gamma; //  locus effect for cline center (alpha), current and last  
  gsl_matrix * etaCur; //  nested population effect for cline center (alpha)
  gsl_matrix * etaLast; //  nested population effect for cline center (alpha)  
  gsl_matrix * zeta; // locus effect for cline width (beta), current and last
  gsl_matrix * kappaCur; // nested population effect for cline width (beta)
  gsl_matrix * kappaLast; // nested population effect for cline width (beta)
  double tauAlphaCur; // precision for gamma random effects
  double tauAlphaLast; // precision for gamma random effects
  double tauBetaCur; // precision for zeta random effects
  double tauBetaLast; // precision for zeta random effects
  gsl_matrix * nu; //  precision for eta random effects
  gsl_matrix * omega; // precision for kappa random effects
  gsl_matrix * alphaCur; // cline center paramer
  gsl_matrix * alphaLast; // cline center parameter
  gsl_matrix * betaCur; // cline width paramer
  gsl_matrix * betaLast; // cline width parameter

  gsl_matrix * phi; // ancestry probability

  gsl_vector * fis;  // inbreeding coefficient

  gsl_matrix * popCoVaMat;

  // variables for linkage model
  double rhoCur; // correlation coefficient for linked markers
  double rhoLast;
  mat_container * sigmaMat; // per chrom
  gsl_matrix * sigmaAll; // on matrix for all loci

  // proposal variables
  gsl_matrix * etaPropMat;
  gsl_matrix * gamzetPropMat;
  double hiprop;
  double tuneGamma;
  double tuneZeta;
  double tuneEta;
  double tuneTau;
  double tuneRho;
  double tuneFis;


};

struct datadimcont{
  int nloci;
  int nind;
  int maxNallele;
  int npop;
  int nchrom;

  gsl_vector_int * nallele;
  gsl_vector_int * pops;
  gsl_vector_int * ncloci;
  gsl_vector_int * ploidyVec;

  // variables for ngs parental populations
  int nindP0;
  int nindP1;

};

struct auxcont{
  // auxillarly vectors
  gsl_vector * vecAllelen1;
  gsl_vector * vecAllelen2;
  gsl_vector_uint * vecAllele_int;
  gsl_vector * vecLocin;
  gsl_vector * vecLocin2;
  gsl_vector * vecIndn;
  gsl_vector * vecPopn1;
  gsl_vector * vecPopn2;
  gsl_vector * vecUnit;

  gsl_vector * vecnkbynk; // number of alleles by number of alleles
  gsl_vector_uint * vecnkbynk_int; 

  gsl_vector * zero2one;
  gsl_vector * zero2one1;
  gsl_vector * zero2one2;
  gsl_vector * zero2one3;

  gsl_vector * M;

  vec_container * chromVecLocin;
  vec_container * chromVecLocin2;
};

struct hdf5cont{
  hid_t file;
  hid_t datatype, dataspacelocus, dataspaceindividual, dataspacemcmc, dataspacelocuspop;
  hid_t datasetAlpha, datasetBeta, datasetHi, datasetLnL;
  hid_t datasetTauA, datasetTauB, datasetEta, datasetKappa, datasetRho;
  hid_t datasetQgamma, datasetQzeta, datasetQeta, datasetQkappa;
  hid_t datasetQgammaCAR, datasetQzetaCAR, datasetHet;
  hid_t locusvector, indvector, unitvector;
  hsize_t dimslocus[2], dimsindividual[2], dims[1], dimspop[3], dimsnull[1]; /* dataset dimension */ 
  herr_t status; /* not sure what this does, but it appears to be important?? */
  hsize_t startr3[3]; /* Start of hyperslab */
  hsize_t startr2[2];
  hsize_t startr1[1];
  hsize_t countr3[3]; /* Block count */
  hsize_t countr2[2];
  hsize_t countr1[1]; 
  hsize_t blockr3[3]; /* Block size */
  hsize_t blockr2[2]; 
  hsize_t blockr1[1];
};

// function declarations
using namespace std;

/* ------- functions from bgc_func_readdata.C -------- */
void usage(char * name);

int getNloci(string filename);
int getNallele(string filename, gsl_vector_int * nallele);
int getNpop(string filename);
int getNind(string filename);
void getData(string hybridfile,string p0file,string p1file, datacont * data,
	     datadimcont * datadim, int allDiploid);
/* ------- functions from bgc_func_initialize.C -------- */
void initParams(paramcont * param, datacont * data, datadimcont * datadim,
		auxcont * auxvar, int ngs, int linkMod);
void initParams2(paramcont * param, datacont * data, datadimcont * datadim,
		auxcont * auxvar, int ngs, int linkMod);
void initRandomEffect(gsl_matrix * reParam, int dim1, int dim2, double sd);
void initClineParameter(gsl_matrix * clineParam, gsl_matrix * locusEffect, gsl_matrix * popEffect, 
			int nloci, int npop);
void copyCur2Last(paramcont * param, datadimcont * datadim, gsl_vector * locivec, 
		  gsl_vector * indvec, int linkMod);

/* ------- functions from bgc_func_mcmc.C -------- */
void mcmcUpdate(paramcont * param, datacont * data, datadimcont * datadim, auxcont * auxvar,
		int allDiploid,	int onePrecision, int sumZ, int linkMod, int ngs, 
		double trunc);
void updatePi(gsl_matrix_int * pcnt, gsl_matrix * piCur, gsl_vector_int * nallele,
	      int nloci, gsl_vector * alpha, gsl_vector * theta);
void updateZ(paramcont * param, datacont * data, datadimcont * datadim, int allDiploid);
void updateHi(paramcont * param, datadimcont * datadim, auxcont * auxvar, int allDiploid);
double calcLogMratioHi(paramcont * param, datadimcont * datadim, auxcont * auxvar, 
		       int j, int n, int allDiploid);
void updateFis(gsl_matrix * h, mat_container_int * zCur, gsl_matrix * alphaCur,
	       gsl_matrix * betaCur, gsl_vector * zero2one, gsl_vector * zero2one1,
	       gsl_vector * zero2one2, gsl_vector * zero2one3, gsl_vector * M, 
	       gsl_vector_int * pops, double tuneFis, int nind, 
	       int nloci, gsl_vector_int * ploidyVec, int allDiploid, 
	       gsl_vector * fis);
double calcLogMratioFis(gsl_matrix * h, mat_container_int * zCur, gsl_matrix * alphaCur,
			gsl_matrix * betaCur, gsl_vector * zero2one, gsl_vector * zero2one1,
			gsl_vector * zero2one2, gsl_vector * zero2one3, gsl_vector * M,
		        int j, int n, int nloci, gsl_vector_int * ploidyVec, 
			int allDiploid, gsl_vector * fis, double newf);
int needSplice(double a, double b, gsl_vector * zero2one1, gsl_vector * zero2one2,
	       gsl_vector * zero2one3, double * y, double * mMax, double * mMin,
	       gsl_vector * M);
double calcOnePhi(double a, double b, double h, double y, double mMin, double mMax,
		  int splice);

void updateGamma(paramcont * param, datadimcont * datadim, auxcont * auxvar,  
		 int allDiploid, int sumZ);
double calcLogMratioGamma(paramcont * param, datadimcont * datadim, auxcont * auxvar, 
			  int i,  int allDiploid);
void updateZeta(paramcont * param, datadimcont * datadim, auxcont * auxvar,
		int allDiploid, int sumZ);
double calcLogMratioZeta(paramcont * param, datadimcont * datadim, auxcont * auxvar,
			 int i, int allDiploid);
void updateEta(paramcont * param, datadimcont * datadim, auxcont * auxvar,
	       int i, int allDiploid, int onePrecision);
double calcLogMratioEta(paramcont * param, datadimcont * datadim, auxcont * auxvar,
			int i, int npop, int allDiploid);
void updateKappa(paramcont * param, datadimcont * datadim, auxcont * auxvar,
		 int i, int allDiploid, int onePrecision);
double calcLogMratioKappa(paramcont * param, datadimcont * datadim, auxcont * auxvar,
			  int i, int allDiploid,int onePrecision);
void updateTau(gsl_matrix * parameters, double * precision, int nloci, double trunc);
void updateNuOmega(gsl_matrix * parameters, gsl_matrix * precision, int nloci, int npop,
		   int onePrecision);

void calcPhi(paramcont * param, datadimcont * datadim, auxcont * auxvar);
double calcLnL(paramcont * param, datacont * data, datadimcont * datadim, 
	       int allDiploid);
void calcQuants(paramcont * param, datadimcont * datadim, FILE * fpQgamma, FILE * fpQzeta, 
		FILE * fpQeta, FILE * fpQkappa);
void calcHet(mat_container_int * zCur, int nloci, int nind, gsl_vector_int * ploidyVec,
	     int allDiploid, FILE * fpHet, int step);

/* ------- functions from bgc_func_write.C -------- */
void printSamples(paramcont * param, datadimcont * datadim, double LnL,
		  FILE * fpHi, FILE * fpAlpha, FILE * fpBeta, FILE * fpLnL);
void printTaus(double tauAlphaCur, double tauBetaCur, FILE * fpTau);
void printPop(gsl_matrix * etaCur, gsl_matrix * kappaCur, int npop, int nloci, 
	      FILE * fpEta, FILE * fpKappa);

/* ------- functions from bgc_func_linkage.C -------- */
void getNchrom(string filename, int * nchrom);
void getMap(string filename, gsl_matrix * proxMat, gsl_matrix * mapData, int nloci, double maxDist,
	    gsl_vector_int * ncloci, int nchrom);
void sepMap(gsl_matrix * proxMat, gsl_matrix * mapData, int nloci, gsl_matrix_int * firstLast, 
	    mat_container * chromProxMat);
			    
void printRho(double rhoCur, FILE * fpRho);

void setSigma(gsl_matrix * sigmaAll, gsl_matrix * proxMat, double rho, int nloci);
void updateGammaCAR(paramcont * param, datadimcont * datadim, auxcont * auxvar, 
		    int allDiploid);
double calcLogMratioGammaCAR(paramcont * param, datadimcont * datadim, auxcont * auxvar,
			     int i, int allDiploid);
void updateZetaCAR(paramcont * param, datadimcont * datadim, auxcont * auxvar, 
		   int allDiploid);
double calcLogMratioZetaCAR(paramcont * param, datadimcont * datadim, auxcont * auxvar,
			    int i, int allDiploid);
void updateTauCAR(gsl_matrix * parameters, double * precisionCur, double precisionLast, 
		  gsl_matrix * sigmaAll, double rho, int nloci, double tuneTau);
void updateRhoCAR(paramcont * param, datacont * data, datadimcont * datadim);
void calcQuantsCAR(paramcont * param, datadimcont * datadim, FILE * fpQgamma, FILE * fpQzeta);

void setH(int nind, string hfile, gsl_matrix * hi);

/* ------- functions from bgc_func_ngs.C -------- */
void getDatangs(string hybridfile,string p0file,string p1file, datacont * data,
		datadimcont * datadim);
void updateZngs(paramcont * param, datacont * data, datadimcont * datadim, auxcont * auxvar);
void initZngs(mat_container_int * zCur, mat_container_int * phcnt, gsl_matrix * pi0Cur,
	      gsl_matrix * pi1Cur, gsl_vector_int * nallele, int nloci, 
	      int nind);
void initPcnts(datacont * data, datadimcont * datadim, gsl_vector * pvec, 
	       gsl_vector_uint * sam);
void updatePcnt(datacont * data, auxcont * auxvar, datadimcont * datadim, gsl_matrix_int * pcnt, 
		mat_container_int * ngsPcnt, int nind, gsl_matrix * pi);
double calcLnLngs(paramcont * param, datacont * data, datadimcont * datadim, auxcont * auxvar);

/* -------- functions for writing to hdf5 format ------------ */
void writeHDF5basic(paramcont * param, datadimcont * datadim, auxcont * auxvar, double LnL, hdf5cont * hdf5, 
		    int step, int burn, int thin);
void writeHDF5tau(paramcont * param, auxcont * auxvar, double talpha, double tbeta, hdf5cont * hdf5, 
		  int step, int burn, int thin);
void writeHDF5pop(paramcont * param, datadimcont * datadim, auxcont * auxvar, hdf5cont * hdf5, 
		  int step, int burn, int thin);
void writeHDF5rho(auxcont * auxvar, double rho, hdf5cont * hdf5, int step, int burn, int thin);
void writeHDF5quants(paramcont * param, datadimcont * datadim, auxcont * auxvar, hdf5cont * hdf5, 
		     int step, int burn, int thin);
void writeHDF5quantsCAR(paramcont * param, datadimcont * datadim, auxcont * auxvar, hdf5cont * hdf5, 
		     int step, int burn, int thin);
void writeHDF5het(paramcont * param, datadimcont * datadim, auxcont * auxvar, hdf5cont * hdf5, 
		  int step, int burn, int thin, int allDiploid);
