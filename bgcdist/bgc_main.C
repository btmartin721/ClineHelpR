// file: bgc_main.C
// Time-stamp: <Thursday, 03 April 2014, 11:12 MDT -- zgompert>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <time.h>
#include <omp.h>
#include <getopt.h>
#include "bgc.h"
#include "mvrandist.h"
#include "hdf5.h"

/* To compile on OSX, requires GSL, no longer uses -framework Accelerate? */
/* h5c++ -Wall -O2 -o bgc bgc_main.C bgc_func_readdata.C bgc_func_initialize.C bgc_func_mcmc.C bgc_func_write.C bgc_func_linkage.C bgc_func_ngs.C bgc_func_hdf5.C mvrandist.c -lgsl -lm  */

/* To compile on Ubuntu Linux, i.e., my laptop. Remove -Wall to turn warnings off */
/* h5c++ -Wall -O2 -o bgc bgc_main.C bgc_func_readdata.C bgc_func_initialize.C bgc_func_mcmc.C bgc_func_write.C bgc_func_linkage.C bgc_func_ngs.C bgc_func_hdf5.C mvrandist.c -lgsl -lgslcblas */

/* h5cc -Wall -O3 -o flat2hdf5 flat2hdf5.c -lgsl -lgslcblas */

/* To compile on Linux cluster Seismic against atlas, requires GSL */
/* g++ -O2 -I/usr/local/gsl-1.14/include -I/usr/local/atlas-3.9.28/include  -L/usr/local/gsl-1.14/lib -L/usr/local/atlas-3.9.28/lib -Wall -o bgc bgc_main.C bgc_func_readdata.C bgc_func_initialize.C bgc_func_mcmc.C bgc_func_write.C bgc_func_linkage.C bgc_func_ngs.C mvrandist.c  -lgsl -lgslcblas -latlas -lm */


using namespace std;

gsl_rng * r;  /* global state variable for random number generator */

/* ----------------- */
/* beginning of main */
/* ----------------- */

int main(int argc, char *argv[]) {

  time_t start = time(NULL);
  time_t end;
  int rng_seed = 0, ch = 0, i = 0, j = 0, w = 0, step = 0;
  int mcmcL = 1000;
  int format = 0; // 0 = HDF5, 1 = text, 2 = HDF5 and text
  int linkMod = 0; // default, do not use the linkage model
  int thin = 1; // frequency to print mcmc samples
  int printLevel = 0; // 0 = LnL, alpha, beta, hi; 1 adds tau, 2 adds eta and kappa
  int printQ = 0; // 1 turns on calculation and printing of quantiles for cline parameters
  int printHet = 0; // 1 turns on calculation and printing of interspecific heterozygosity
  int burn = 0;
  double trunc = 0; // minimum bounds of the gamma prior on tau, > 0 = truncated gamma prior
  int allDiploid = 1; // set to zero with command line and provide ploidy info.
  int onePrecision = 0; // share precision prameter for nu and omega across loci
  int sumZ = 1; // use sum to zero constraints on locus effects, default is 0 = false
  int ngs = 0; // boolean, set to 1 to use ngs model
  double maxDist = 0.5; // default maximum distance, edit with -D
  int simpleInit = 1; // boolean, set to 0 for initialization based on data
  double avalue = 0;
  double LnL = 0;

  // main variables, these are structures that include many variables
  datacont data;
  paramcont param;
  datadimcont datadim;
  auxcont auxvar;
  hdf5cont hdf5;

  // set variables to 0
  datadim.nloci = 0;
  datadim.nind = 0;
  datadim.maxNallele = 0;
  datadim.npop = 0;
  datadim.nchrom = 0;
  datadim.nindP0 = 0;
  datadim.nindP1 = 0;
  param.tauAlphaCur = 0; 
  param.tauAlphaLast = 0;
  param.tauBetaCur = 0;
  param.tauBetaLast = 0;

  param.rhoCur = 0; // correlation coefficient for linked markers
  param.rhoLast = 0;

  // tunning parameters for proposal distribution
  param.hiprop = HPROP; // adjust with -u
  param.tuneGamma = TGAMMA; // adjust with -g
  param.tuneZeta = TZETA; // adjust with -z
  param.tuneEta = TETA; // adjust with -e
  param.tuneTau = TTAU; // adjust ?
  param.tuneRho = TRHO; // adjust ?
  param.tuneFis = TFIS; // adjust -f

  // fixed error rate
  data.error = 0.0; // set to 9 for locus-specific error rates

  // infiles
  string hybridfile = "undefined"; // inputfiles
  string p0file = "undefined";
  string p1file = "undefined";
  string mapfile = "undefined";
  string fileinfo = "OF";
  string onefile = "undefined";

  // files for MCMC output
  FILE * fpHi = NULL;
  FILE * fpAlpha = NULL;
  FILE * fpBeta = NULL;
  FILE * fpLnL = NULL;
  FILE * fpTau = NULL;
  FILE * fpEta = NULL;
  FILE * fpKappa = NULL; 
  FILE * fpQgamma = NULL;
  FILE * fpQzeta = NULL;
  FILE * fpQeta = NULL;
  FILE * fpQkappa = NULL;
  FILE * fpHet = NULL;
  FILE * fpRho = NULL;
  //FILE * fpFis = NULL;

  // variables for getopt_long
  static struct option long_options[] = {
    {"version", no_argument,       0, 'v'},
    //     {"append",  no_argument,       0, 'b'},
    //     {"delete",  required_argument, 0, 'd'},
    //     {"create",  required_argument, 0, 'c'},
    //     {"file",    required_argument, 0, 'f'},
    {0, 0, 0, 0}
  };
  int option_index = 0;

  // get command line arguments
  if (argc < 2) {
    usage(argv[0]);
  }

  /* Create a new file using H5F_ACC_TRUNC access, default file
     creation properties, and default file access properties. */
  char * hdf5outfile = (char*) "mcmcout.hdf5";

  while ((ch = getopt_long(argc, argv, "h:a:b:x:u:g:z:e:t:p:n:d:q:o:s:v:i:m:M:D:F:H:N:I:f:T:E:O:",
			   long_options, &option_index)) != -1){
    switch(ch){
    case 'h':
      hybridfile = optarg;
      break;
    case 'a':
      p0file = optarg;
      break;
    case 'b':
      p1file = optarg;
      break;
    case 'x':
      mcmcL = atoi(optarg);
      break;
    case 'u':
      param.hiprop = atof(optarg);
      break;
    case 'g':
      param.tuneGamma = atof(optarg);
      break;
     case 'z':
      param.tuneZeta = atof(optarg);
      break;
     case 'e':
      param.tuneEta = atof(optarg);
      break;
     case 't':
      thin = atoi(optarg);
      break;
    case 'p':
      printLevel = atoi(optarg);
      break;
    case 'n':
      burn = atoi(optarg);
      break;
    case 'd':
      allDiploid = atoi(optarg);
      break;
    case 'q':
      printQ = atoi(optarg);
      break;
    case 'o':
      onePrecision = atoi(optarg);
      break;
    case 's':
      sumZ = atoi(optarg);
      break;
    case 'v':
      printf("bgc (Bayesian genomic clines) version %s\n", VERSION); 
      // VERSION is a macro that is definied in bgc.h
      exit(0); // note program will exit if this option is specified
    case 'i':
      printHet = atoi(optarg);
      break;
    case 'm':
      linkMod = atoi(optarg);
      break;
    case 'M':
      mapfile = optarg;
      break;
    case 'D':
      maxDist = atof(optarg);
      break;
    case 'F':
      fileinfo = optarg;
      hdf5outfile = strcat(optarg, ".hdf5");
      break;
    case 'N':
      ngs = atoi(optarg);
      break;
    case 'I':
      simpleInit = atoi(optarg);
      break;
    case 'f':
      param.tuneFis = atof(optarg);
      break;
    case 'T':
      trunc = atof(optarg);
      break;
    case 'E':
      data.error = atof(optarg);
      break;
    case 'O':
      format = atoi(optarg);
      break;
    case '?':
    default:
      usage(argv[0]);
    }
  }


  if (format != 1){
    hdf5.file = H5Fcreate(hdf5outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }


  if (allDiploid == 0 && ngs == 1){
    cerr << "WARNING: the ngs model assumes diploid data\n" << endl;
  }
  
  // set up gsl random number generation 
  gsl_rng_env_setup();
  r = gsl_rng_alloc (gsl_rng_default);
  srand(time(NULL));
  rng_seed = rand();
  gsl_rng_set(r, rng_seed); /* seed gsl_rng with output of rand, which
			       was seeded with result of time(NULL) */
  
  cerr << "Reading input files" << endl;

  // determine number of loci and number of alleles per locus
  datadim.nloci = getNloci(hybridfile);
  datadim.nallele = gsl_vector_int_calloc(datadim.nloci);
  datadim.maxNallele = getNallele(hybridfile, datadim.nallele);
  datadim.npop = getNpop(hybridfile);

  // determine number of individuals per locus
  datadim.nind = getNind(hybridfile);

  // for ngs, need to get number for each parent too
  if (ngs == 1){
    datadim.nindP0 = getNind(p0file);
    datadim.nindP1 = getNind(p1file);
  }

  cerr << "Number of loci: " << datadim.nloci << endl;
  cerr << "Number of admixed populations: " << datadim.npop << endl;
  cerr << "Number of individuals: " << datadim.nind << endl;
  if (linkMod == 1){
    cerr << "Using the linkage model for locus effects" << endl;
  }

  if (ngs == 1){
    cerr << "Allowing for uncertainty in allele counts" << endl;
  }

  cerr << "Allocating memory" << endl;

  // allocate memory for data
  data.p0cnt = gsl_matrix_int_calloc(datadim.nloci,datadim.maxNallele);
  data.p1cnt = gsl_matrix_int_calloc(datadim.nloci,datadim.maxNallele);
  if (ngs == 1){
    data.ngsp0cnt = new mat_container_int[datadim.nloci];
    data.ngsp1cnt = new mat_container_int[datadim.nloci];
    for (i=0; i<datadim.nloci; i++){
      data.ngsp0cnt[i].mat = gsl_matrix_int_calloc(datadim.nindP0,
					      gsl_vector_int_get(datadim.nallele,i));
      data.ngsp1cnt[i].mat = gsl_matrix_int_calloc(datadim.nindP1,
					      gsl_vector_int_get(datadim.nallele,i));
    }
  }
  data.phcnt = new mat_container_int[datadim.nloci];
  for (i=0; i<datadim.nloci; i++){
    data.phcnt[i].mat = gsl_matrix_int_calloc(datadim.nind,
					 gsl_vector_int_get(datadim.nallele,i));
  }
  datadim.pops = gsl_vector_int_calloc(datadim.nind);
  datadim.ploidyVec = gsl_vector_int_calloc(datadim.nloci);

  // allocate memory for MVN matrixes used for multiple population analysis
  if (datadim.npop > 1){
    param.popCoVaMat = gsl_matrix_calloc(datadim.npop-1,datadim.npop-1);
    param.etaPropMat = gsl_matrix_calloc(datadim.npop-1,datadim.npop-1);
    for(j=0; j<(datadim.npop-1); j++){
      for(w=0; w<(datadim.npop-1); w++){
	if(w==j){
	  gsl_matrix_set(param.etaPropMat, j, w, gsl_pow_2(param.tuneEta));
	}
	else{
	  gsl_matrix_set(param.etaPropMat, j, w, -1/(datadim.npop - 1));
	}
      }
    }
  }

  // allocate memory for auxillarly variables
  auxvar.vecAllelen1 = gsl_vector_calloc(datadim.maxNallele);
  auxvar.vecAllelen2 = gsl_vector_calloc(datadim.maxNallele);
  auxvar.vecAllele_int = gsl_vector_uint_calloc(datadim.maxNallele);
  auxvar.vecIndn = gsl_vector_calloc(datadim.nind);
  auxvar.vecLocin = gsl_vector_calloc(datadim.nloci);
  if (datadim.npop > 1){
    auxvar.vecPopn1 = gsl_vector_calloc(datadim.npop-1);
    auxvar.vecPopn2 = gsl_vector_calloc(datadim.npop-1);
  }
  if (linkMod == 1){
    auxvar.vecLocin2 = gsl_vector_calloc(datadim.nloci);
  }
  auxvar.vecUnit = gsl_vector_calloc(1);
  auxvar.vecnkbynk = gsl_vector_calloc(datadim.maxNallele * datadim.maxNallele);
  auxvar.vecnkbynk_int = gsl_vector_uint_calloc(datadim.maxNallele * datadim.maxNallele);


  // 100 divides 0..1 into 100 intervals for determinig whether the
  // splicing fuction is needed
  auxvar.zero2one = gsl_vector_calloc(100);
  auxvar.zero2one1 = gsl_vector_calloc(100);
  auxvar.zero2one2 = gsl_vector_calloc(100);
  auxvar.zero2one3 = gsl_vector_calloc(100);  

  avalue = 0;
  for(i=0; i<100; i++){
    gsl_vector_set(auxvar.zero2one, i, avalue);
    avalue += 0.01;
  }

  auxvar.M = gsl_vector_calloc(3);

  // allocate memory for sequence error rates
  data.errVec = gsl_vector_calloc(datadim.nloci);

  // allocate memory for MCMC output, only storing current (0) and last (1) iteration
  param.pi0Cur = gsl_matrix_calloc(datadim.nloci,datadim.maxNallele);
  param.pi1Cur = gsl_matrix_calloc(datadim.nloci,datadim.maxNallele);
  param.pi0Last = gsl_matrix_calloc(datadim.nloci,datadim.maxNallele);
  param.pi1Last = gsl_matrix_calloc(datadim.nloci,datadim.maxNallele);
  param.zCur = new mat_container_int[datadim.nloci];
  param.zLast = new mat_container_int[datadim.nloci];
  for (i=0; i<datadim.nloci; i++){
    param.zCur[i].mat = gsl_matrix_int_calloc(datadim.nind,PLOIDY);
    param.zLast[i].mat = gsl_matrix_int_calloc(datadim.nind,PLOIDY);
  }
  param.hi = gsl_matrix_calloc(datadim.nind,2); // current (0) and last (1)
  param.gamma = gsl_matrix_calloc(datadim.nloci,2); // current (0) and last (1)
  param.etaCur = gsl_matrix_calloc(datadim.nloci,datadim.npop);
  param.etaLast = gsl_matrix_calloc(datadim.nloci,datadim.npop);
  param.zeta = gsl_matrix_calloc(datadim.nloci,2); // current (0) and last (1)
  param.kappaCur = gsl_matrix_calloc(datadim.nloci,datadim.npop);
  param.kappaLast = gsl_matrix_calloc(datadim.nloci,datadim.npop);
  param.nu = gsl_matrix_calloc(datadim.nloci,2); // current (0) and last (1)
  param.omega = gsl_matrix_calloc(datadim.nloci,2); // current (0) and last (1)
  param.alphaCur = gsl_matrix_calloc(datadim.nloci,datadim.npop);
  param.alphaLast = gsl_matrix_calloc(datadim.nloci,datadim.npop);
  param.betaCur = gsl_matrix_calloc(datadim.nloci,datadim.npop);
  param.betaLast = gsl_matrix_calloc(datadim.nloci,datadim.npop);

  param.phi = gsl_matrix_calloc(datadim.nloci,datadim.nind);
  param.fis = gsl_vector_calloc(datadim.nind);

  if (format != 1){ // HDF5 output
    /* HDF5, create datasets */
    /*
     * Describe the size of the array and create the data space for
     * fixed size dataset.  I reuse dataspace and datatype, because they
     * are the same for each of the parameters.
     */
    // set dimensions
    hdf5.dimslocus[0] = datadim.nloci; // rows
    hdf5.dimslocus[1] = (mcmcL - burn) / thin; // cols
    hdf5.dimsindividual[0] = datadim.nind;
    hdf5.dimsindividual[1] = (mcmcL - burn) / thin;
    hdf5.dims[0] = (mcmcL - burn) / thin;
    hdf5.dimsnull[0] = 1;
    hdf5.dimspop[0] = datadim.nloci;
    hdf5.dimspop[1] = (mcmcL - burn) / thin;
    hdf5.dimspop[2] = datadim.npop;

    // create dataspace
    hdf5.dataspacelocus = H5Screate_simple(2,  hdf5.dimslocus, NULL); /* RANK is 2 */
    hdf5.dataspacelocuspop = H5Screate_simple(3,  hdf5.dimspop, NULL); /* RANK is 3 */
    hdf5.dataspaceindividual = H5Screate_simple(2,  hdf5.dimsindividual, NULL); /* RANK is 2 */
    hdf5.dataspacemcmc = H5Screate_simple(1,  hdf5.dims, NULL); /* RANK is 1 */

    // define datatype by copying an existing datatype, little endian double
    hdf5.datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    hdf5.status = H5Tset_order(hdf5.datatype, H5T_ORDER_LE);

    // create datasets, for now just alpha, beta, and hi
    hdf5.datasetAlpha = H5Dcreate2(hdf5.file, "alpha", hdf5.datatype, hdf5.dataspacelocuspop, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetBeta = H5Dcreate2(hdf5.file, "beta", hdf5.datatype, hdf5.dataspacelocuspop, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetHi = H5Dcreate2(hdf5.file, "hi", hdf5.datatype, hdf5.dataspaceindividual, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetLnL = H5Dcreate2(hdf5.file, "LnL", hdf5.datatype, hdf5.dataspacemcmc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetTauA = H5Dcreate2(hdf5.file, "tau-alpha", hdf5.datatype, hdf5.dataspacemcmc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetTauB = H5Dcreate2(hdf5.file, "tau-beta", hdf5.datatype, hdf5.dataspacemcmc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetEta = H5Dcreate2(hdf5.file, "eta", hdf5.datatype, hdf5.dataspacelocuspop, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetKappa = H5Dcreate2(hdf5.file, "kappa", hdf5.datatype, hdf5.dataspacelocuspop, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetRho = H5Dcreate2(hdf5.file, "rho", hdf5.datatype, hdf5.dataspacemcmc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetQgamma = H5Dcreate2(hdf5.file, "gamma-quantile", hdf5.datatype, hdf5.dataspacelocus, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetQzeta = H5Dcreate2(hdf5.file, "zeta-quantile", hdf5.datatype, hdf5.dataspacelocus, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetQeta = H5Dcreate2(hdf5.file, "eta-quantile", hdf5.datatype, hdf5.dataspacelocuspop, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetQkappa = H5Dcreate2(hdf5.file, "kappa-quantile", hdf5.datatype, hdf5.dataspacelocuspop, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetQgammaCAR = H5Dcreate2(hdf5.file, "gamma-quantile-local", hdf5.datatype, hdf5.dataspacelocus, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetQzetaCAR = H5Dcreate2(hdf5.file, "zeta-quantile-local", hdf5.datatype, hdf5.dataspacelocus, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hdf5.datasetHet = H5Dcreate2(hdf5.file, "interspecific-het", hdf5.datatype, hdf5.dataspaceindividual, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* need a buffer for storing a vector for each parameter
       for use in the hdf5write */
    hdf5.locusvector = H5Screate_simple(1, &hdf5.dimslocus[0], NULL); 
    hdf5.indvector = H5Screate_simple(1, &hdf5.dimsindividual[0], NULL); 
    hdf5.unitvector = H5Screate_simple(1, &hdf5.dimsnull[0], NULL); 
  }

  // read in data
  if (ngs == 0){
    getData(hybridfile,p0file,p1file,&data,&datadim,allDiploid);
  }
  else if (ngs == 1){
    getDatangs(hybridfile,p0file,p1file,&data,&datadim);
  }

  if (linkMod == 1){ // use linkage model, so get map data and allocate memory
    data.proxMat = gsl_matrix_calloc(datadim.nloci,datadim.nloci);
    data.mapData = gsl_matrix_calloc(datadim.nloci,2);

    getNchrom(mapfile,&datadim.nchrom);

    datadim.ncloci = gsl_vector_int_calloc(datadim.nchrom);
    data.firstLast = gsl_matrix_int_calloc(datadim.nchrom,2);

    getMap(mapfile,data.proxMat,data.mapData,datadim.nloci,maxDist,datadim.ncloci,datadim.nchrom);

    auxvar.chromVecLocin = new vec_container[datadim.nchrom];
    auxvar.chromVecLocin2 = new vec_container[datadim.nchrom];
    for (i=0; i<datadim.nchrom; i++){
      auxvar.chromVecLocin[i].vec = gsl_vector_calloc(gsl_vector_int_get(datadim.ncloci,i));
      auxvar.chromVecLocin2[i].vec = gsl_vector_calloc(gsl_vector_int_get(datadim.ncloci,i));
    }

    data.chromProxMat = new mat_container[datadim.nchrom]; // one matrix per chrom.
    for (i=0; i<datadim.nchrom; i++){
      data.chromProxMat[i].mat = gsl_matrix_calloc(gsl_vector_int_get(datadim.ncloci,i),
						   gsl_vector_int_get(datadim.ncloci,i));
    }

    sepMap(data.proxMat,data.mapData,datadim.nloci,data.firstLast,data.chromProxMat);

    param.gamzetPropMat = gsl_matrix_calloc(datadim.nloci,datadim.nloci);
    for (i=0; i<datadim.nloci; i++){
      for (j=0; j<datadim.nloci; j++){
	if (i == j){
	  gsl_matrix_set(param.gamzetPropMat,i,j,0.03 * 0.03);
	}
	else {
	  gsl_matrix_set(param.gamzetPropMat,i,j,0);
	}
      }
    }
    param.sigmaMat = new mat_container[datadim.nchrom]; // one matrix per chrom.
    param.sigmaAll = gsl_matrix_calloc(datadim.nloci,datadim.nloci);
    for (i=0; i<datadim.nchrom; i++){
      param.sigmaMat[i].mat = gsl_matrix_calloc(gsl_vector_int_get(datadim.ncloci,i),
						gsl_vector_int_get(datadim.ncloci,i));
    }
  }

  if (format != 0){ // text file ouput
    /* open output files */
    onefile = fileinfo + "_hi_output.txt";
    fpHi = fopen(onefile.c_str(), "w");  // "w" write rather than "a" append
    onefile = fileinfo + "_alpha_output.txt";
    fpAlpha = fopen(onefile.c_str(), "w");
    onefile = fileinfo + "_beta_output.txt";
    fpBeta = fopen(onefile.c_str(), "w");
    onefile = fileinfo + "_LnL_output.txt";
    fpLnL = fopen(onefile.c_str(), "w");
    //onefile = fileinfo + "_fis_output.txt";
    //fpFis = fopen(onefile.c_str(), "w");
    if (printLevel > 0){
      onefile = fileinfo + "_tau_output.txt";
      fpTau = fopen(onefile.c_str(), "w");
    }
    if (printLevel > 1){
      onefile = fileinfo + "_eta_output.txt";
      fpEta = fopen(onefile.c_str(), "w");
      onefile = fileinfo + "_kappa_output.txt";
      fpKappa = fopen(onefile.c_str(), "w");
    }
    if (printQ == 1){
      onefile = fileinfo + "_q_gamma_output.txt";
      fpQgamma = fopen(onefile.c_str(), "w");
      onefile = fileinfo + "_q_zeta_output.txt";
      fpQzeta = fopen(onefile.c_str(), "w");
      onefile = fileinfo + "_q_etq_output.txt";
      fpQeta = fopen(onefile.c_str(), "w");
      onefile = fileinfo + "_q_kappa_output.txt";
      fpQkappa = fopen(onefile.c_str(), "w");
    }
    if (printHet == 1){
      onefile = fileinfo + "_heterozygosity_output.txt";
      fpHet = fopen(onefile.c_str(), "w");
    }
    if (linkMod == 1){
      onefile = fileinfo + "_rho_output.txt";
      fpRho = fopen(onefile.c_str(), "w");
    }
  }

  // initialize parameters
  cerr << "Initializing MCMC chain" << endl;

  if (simpleInit == 1){
    initParams(&param,&data,&datadim,&auxvar,ngs,linkMod);
  }
  
  else{ // uses data to initialize hybrid index and ancestry
    initParams2(&param,&data,&datadim,&auxvar,ngs,linkMod);
  }

  copyCur2Last(&param,&datadim,auxvar.vecLocin,auxvar.vecIndn,linkMod);

  calcPhi(&param,&datadim,&auxvar);

  cerr << "Initialization complete" << endl;

  // mcmc loop
  for(step=0; step<mcmcL; step++){
    // one mcmc step for each parameter
    mcmcUpdate(&param,&data,&datadim,&auxvar,allDiploid,onePrecision,sumZ,linkMod,ngs,trunc);

    if (step % thin == 0 && step > (burn -1)){
      if (ngs == 0){
	LnL = calcLnL(&param,&data,&datadim,allDiploid);
      }
      else if (ngs == 1){
	LnL = calcLnLngs(&param,&data,&datadim,&auxvar);
      }
      // always print hi,alpha,beta,LnL
      if (format != 0){
	printSamples(&param,&datadim,LnL,fpHi,fpAlpha,fpBeta,fpLnL);
      }
      if (format != 1){
	writeHDF5basic(&param,&datadim,&auxvar,LnL,&hdf5,step,burn,thin);
      }
      if (printLevel > 0){ // also print taus
	if (format != 0){
	  printTaus(param.tauAlphaCur,param.tauBetaCur,fpTau);
	}
	if (format != 1){
	  writeHDF5tau(&param,&auxvar,param.tauAlphaCur,param.tauBetaCur,&hdf5,step,burn,thin);
	}
      }
      if (printLevel > 1){ // also print pop effects
	if (format != 0){
	  printPop(param.etaCur,param.kappaCur,datadim.npop,datadim.nloci,fpEta,fpKappa);
	}
	if (format != 1){
	  writeHDF5pop(&param,&datadim,&auxvar,&hdf5,step,burn,thin);
	}
      }
      if (linkMod == 1){ // also print rho
	if (format != 0){
	  printRho(param.rhoCur,fpRho);
	}
	if (format != 1){
	  writeHDF5rho(&auxvar,param.rhoCur,&hdf5,step,burn,thin);
	}
      }
      // calculate and print quantiles for gamma, zeta, eta, and kappa cline parameters
      if (printQ == 1){ 
	if (linkMod == 0){
	  if (format != 0){
	    calcQuants(&param,&datadim,fpQgamma,fpQzeta,fpQeta,fpQkappa);
	  }
	  if (format != 1){
	    writeHDF5quants(&param,&datadim,&auxvar,&hdf5,step,burn,thin);
	  }
	}
	else if (linkMod == 1){
	  if (format != 0){
	    calcQuantsCAR(&param,&datadim,fpQgamma,fpQzeta);
	  }
	  if (format != 1){
	    writeHDF5quantsCAR(&param,&datadim,&auxvar,&hdf5,step,burn,thin);
	  }
	}
      }
      if (printHet == 1){
	if (format != 0){
	  calcHet(param.zCur,datadim.nloci,datadim.nind,datadim.ploidyVec,allDiploid,fpHet,step);
	}
	if (format != 1){
	  writeHDF5het(&param,&datadim,&auxvar,&hdf5,step,burn,thin,allDiploid);
	}
      }
    }

    if (step % 1000 == 0){
      if (ngs == 0){
	LnL = calcLnL(&param,&data,&datadim,allDiploid);
      }
      else if (ngs == 1){
	LnL = calcLnLngs(&param,&data,&datadim,&auxvar);
      }
      cerr << "mcmc iteration: " << step << " ........ LnL: " << LnL << endl;
    }
    // copy parameters over 
    copyCur2Last(&param,&datadim,auxvar.vecLocin,auxvar.vecIndn,linkMod);
  }

  cerr << "Finished mcmc" << endl;

  // close text output files
  if (format != 0){
    fclose(fpHi);
    fclose(fpBeta);
    fclose(fpAlpha);
    fclose(fpLnL);
    //fclose(fpFis);
    if (printLevel > 0){
      fclose(fpTau);
    }

    if (printLevel > 1){
      fclose(fpEta);
      fclose(fpKappa);
    }

    if (printQ == 1){
      fclose(fpQgamma);
      fclose(fpQzeta);
      fclose(fpQeta);
      fclose(fpQkappa);

    }

    if (printHet == 1){
      fclose(fpHet);
    }
    if (linkMod == 1){
      fclose(fpRho);
    }
  }
  // free memory
  gsl_matrix_int_free(data.p0cnt);
  gsl_matrix_int_free(data.p1cnt);
  for (i=0; i<datadim.nloci; i++){
    gsl_matrix_int_free(data.phcnt[i].mat);
  }
  delete data.phcnt;
  if (linkMod == 1){
    gsl_matrix_free(data.proxMat);
    gsl_matrix_free(param.gamzetPropMat);
    gsl_matrix_free(param.sigmaAll);
    for (i=0; i<datadim.nchrom; i++){
      gsl_matrix_free(data.chromProxMat[i].mat);
      gsl_matrix_free(param.sigmaMat[i].mat);
      gsl_vector_free(auxvar.chromVecLocin[i].vec);
      gsl_vector_free(auxvar.chromVecLocin2[i].vec);
    }
    delete data.chromProxMat;
    delete param.sigmaMat;
    delete auxvar.chromVecLocin;
    delete auxvar.chromVecLocin2;
  }
  gsl_vector_int_free(datadim.pops);
  gsl_vector_int_free(datadim.ploidyVec);

  gsl_vector_free(auxvar.vecAllelen1);
  gsl_vector_free(auxvar.vecAllelen2);
  gsl_vector_uint_free(auxvar.vecAllele_int);
  gsl_vector_free(auxvar.vecIndn);
  gsl_vector_free(auxvar.vecLocin);

  if (datadim.npop > 1){
    gsl_vector_free(auxvar.vecPopn1);
    gsl_vector_free(auxvar.vecPopn2);
    gsl_matrix_free(param.popCoVaMat);
    gsl_matrix_free(param.etaPropMat);
  }

  gsl_vector_free(auxvar.zero2one);
  gsl_vector_free(auxvar.zero2one1);
  gsl_vector_free(auxvar.zero2one2);
  gsl_vector_free(auxvar.zero2one3);
  gsl_vector_free(auxvar.M);

  gsl_matrix_free(param.pi0Cur);
  gsl_matrix_free(param.pi0Last);
  gsl_matrix_free(param.pi1Cur);
  gsl_matrix_free(param.pi1Last);
  for (i=0; i<datadim.nloci; i++){
    gsl_matrix_int_free(param.zCur[i].mat);
    gsl_matrix_int_free(param.zLast[i].mat);
  }
  delete param.zCur;
  delete param.zLast;
  gsl_matrix_free(param.hi);
  gsl_matrix_free(param.gamma);
  gsl_matrix_free(param.etaCur);
  gsl_matrix_free(param.etaLast);
  gsl_matrix_free(param.zeta);
  gsl_matrix_free(param.kappaCur);
  gsl_matrix_free(param.kappaLast);
  gsl_matrix_free(param.nu);
  gsl_matrix_free(param.omega);
  gsl_matrix_free(param.alphaCur);
  gsl_matrix_free(param.alphaLast);
  gsl_matrix_free(param.betaCur);
  gsl_matrix_free(param.betaLast);
  gsl_matrix_free(param.phi);

  // close hdf5 output file
  if (format != 1){
    H5Fclose(hdf5.file);
  }
  // prints run time
  end = time(NULL);
  cerr << "Runtime: " << (end-start)/3600 << " hr " << (end-start)%3600/60 << " min ";
  cerr << (end-start)%60 << " sec" << endl;
  return 0;
}
