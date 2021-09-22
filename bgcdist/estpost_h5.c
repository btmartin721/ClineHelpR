/* read hdf5 and write data summary */
/* Time-stamp: <Tuesday, 21 February 2012, 14:45 CST -- zgompert> */

/* Compilation for Linux */
/* h5cc -Wall -O3 -o estpost estpost_h5.c -lgsl -lgslcblas */
/* Compilation for OSX */
/* h5cc -Wall -O3 -o estpost estpost_h5.c -lgsl -lm */


#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_histogram.h>


#include "hdf5.h"
#include <getopt.h>

/* header, function declerations */

#define VERSION "1.0 - 19 February 2012"

void usage(char * name);
void estpost(hid_t file, const char * param, double credint, int burn, int nbins, int sumtype, int wids);
void calclocuspop(const char * param, hid_t dataspace, hid_t dataset, gsl_vector * sample, gsl_histogram * hist,
		  double credint, int nloci, int nsamples, int npop, int burn, int nbins, int sumtype, int wids);
void calcparam(const char * param, hid_t dataspace, hid_t dataset, gsl_vector * sample, gsl_histogram * hist,
	       double credint, int nparam, int nsamples, int burn, int nbins, int sumtype, int wids);
void calcsimple(const char * param, hid_t dataspace, hid_t dataset, gsl_vector * sample, gsl_histogram * hist, 
		double credint, int nsamples, int burn, int nbins, int sumtype, int wids);
void calcci(gsl_vector * sample, double credint, int nsamples);
void calchist(gsl_vector * sample,  gsl_histogram * hist, int nsamples, int nbins);
void writetext(gsl_vector * sample,  int nsamples);

FILE * outfp; 

/* beginning of main */

int main (int argc, char **argv) {
  int sumtype = 0; /* summary to perform
		      0 = estimate and ci, 1 = histogram, 2 = convert to text */
  int ch = 0;
  int burn = 0; /* discard the first burn samples as a burn-in */
  int nbins = 20; /* number of bins for histogram */
  int wids = 1; /* boolean, write ids in first column */

  char * infile = "undefined"; /* filename */
  char * outfile = "postout.txt";
  char * param = "undefined"; /* paramter to summarize */

  double credint = 0.95; /* default = 95% credible interval */

  /* variables for getopt_long */
  static struct option long_options[] = {
    {"version", no_argument, 0, 'v'},
    {0, 0, 0, 0} 
  };
  int option_index = 0;

  /* variables for hdf5 */
  hid_t       file;         /* file handle */
  
  /*  get command line arguments */
  if (argc < 2) {
    usage(argv[0]);
  }
  
  while ((ch = getopt_long(argc, argv, "i:o:v:p:c:b:h:s:w:", 
			   long_options, &option_index)) != -1){
    switch(ch){
    case 'i':
      infile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'p':
      param = optarg;
      break;
    case 'c':
      credint = atof(optarg);
      break;
    case 'b':
      burn = atoi(optarg);
      break;
    case 'h':
      nbins = atoi(optarg);
      break;
    case 's':
      sumtype = atoi(optarg);
      break;
    case 'w':
      wids = atoi(optarg);
      break;
    case 'v':
      printf("%s version %s\n", argv[0], VERSION); 
      /* VERSION is a macro  */
      exit(0); /* note program will exit if this option is specified */
    case '?':
    default:
      usage(argv[0]);
    }
  }
  /* open the h5 file, read only */
  file = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);

 /* open the outfile */
  outfp = fopen(outfile, "w");
  if ( !outfp ){
    fprintf(stderr, "Can't open %s for writing!\n", outfile);
    exit(1);
  }
  if ((sumtype == 0) && (wids == 1)){
    fprintf(outfp,"param,mean,median,%.3f_CI_LB,%.3f_CI_UB\n", credint, credint);
  }
  
  /* main function */
  estpost(file, param, credint, burn, nbins, sumtype, wids);
  
  fclose(outfp);
  H5Fclose(file);

  return 0;

}

/* ---------------- Functions ------------------ */

/* Prints usage */
void usage(char * name){
  fprintf(stderr,"\n%s version %s\n\n", name, VERSION); 
  fprintf(stderr, "Usage:   estpost -i inputfile -p parameter [options]\n");
  fprintf(stderr, "-i     Infile, MCMC results from bgc in HDF5 format\n");
  fprintf(stderr, "-o     Outfile [default = postout]\n");
  fprintf(stderr, "-p     Name of parameter to summarize, e.g., 'alpha'\n");
  fprintf(stderr, "-c     Credible interval to calculate [default = 0.95]\n");
  fprintf(stderr, "-b     Number of additinal MCMC samples to discard for burn-in [default = 0]\n");
  fprintf(stderr, "-h     Number of bins for posterior sample histogram [default = 20]\n");
  fprintf(stderr, "-s     Which summary to perform: 0 = posterior estimates and credible intervals\n");
  fprintf(stderr, "                                 1 = histogram of posterior samples\n");
  fprintf(stderr, "                                 2 = convert to plain text\n");
  fprintf(stderr, "-w     Write parameter identification to file, boolean [default = 1]\n");
  fprintf(stderr, "-v     Display software version\n");
  exit(1);
}

/* Function estimates the mean, median and credible interval 
   of the posterior distribution */
void estpost(hid_t file, const char * param, double credint, int burn, int nbins, int sumtype, int wids){
  gsl_vector * sample;
  gsl_histogram * hist;
  hid_t dataset, dataspace;
  hsize_t dims[3]; /* note 3 dimensions is maximum */
  int rank = 0, status_n = 0;
  /* dimensions for mcmc samples */
  int nloci = 0, nind = 0, npop = 0, nsamples = 0;

  hist = gsl_histogram_alloc((size_t) nbins);

  dataset = H5Dopen2(file, param, H5P_DEFAULT);
  dataspace = H5Dget_space(dataset); 
  rank = H5Sget_simple_extent_ndims(dataspace);
  status_n = H5Sget_simple_extent_dims(dataspace, dims, NULL);

  /* interpret dimensions based on parameter, then calculate desired quantities */
  if ((strcmp(param, (const char *) "alpha") == 0) || (strcmp(param, (const char *) "beta") == 0) || 
      (strcmp(param, (const char *) "eta") == 0) || (strcmp(param, (const char *) "kappa") == 0) || 
      (strcmp(param, (const char *) "eta-quantile") == 0) || 
      (strcmp(param, (const char *) "kappa-quantile") == 0)){
    nloci = dims[0];
    nsamples = dims[1];
    if (burn >= nsamples){
      printf("Burnin exceeds number of samples\n");
      exit(1);
    }
    npop = dims[2];
    sample = gsl_vector_calloc(nsamples - burn);
    printf("parameter dimensions for %s: loci = %d, samples = %d, pouplations = %d\n", 
	   param, nloci, nsamples, npop);
    calclocuspop(param, dataspace, dataset, sample, hist, credint, nloci, nsamples, npop, burn, nbins, sumtype,
		 wids);
    
  }
  else if ((strcmp(param, (const char *) "gamma-quantile") == 0) || 
	   (strcmp(param, (const char *) "zeta-quantile") == 0) || 
	   (strcmp(param, (const char *) "gamma-quantile-local") == 0) || 
	   (strcmp(param, (const char *) "zeta-quantile-local") == 0)){
    nloci = dims[0];
    nsamples = dims[1];
    if (burn >= nsamples){
      printf("Burnin exceeds number of samples\n");
      exit(1);
    }
    sample = gsl_vector_calloc(nsamples - burn);
    printf("parameter dimensions for %s: loci = %d, samples = %d\n", 
	   param, nloci, nsamples);
    calcparam(param, dataspace, dataset, sample, hist, credint, nloci, nsamples, burn, nbins, sumtype, wids);
  }
  else if ((strcmp(param, (const char *) "hi") == 0) || (strcmp(param, (const char *) "interspecific-het") == 0)){
    nind = dims[0];
    nsamples = dims[1];
    if (burn >= nsamples){
      printf("Burnin exceeds number of samples\n");
      exit(1);
    }
    sample = gsl_vector_calloc(nsamples - burn);
    printf("parameter dimensions for %s: individuals = %d, samples = %d\n", 
	   param, nind, nsamples);
    calcparam(param, dataspace, dataset, sample, hist, credint, nind, nsamples, burn, nbins, sumtype, wids);
  }
  else if ((strcmp(param, (const char *) "LnL") == 0) || (strcmp(param, (const char *) "tau-alpha") == 0) || 
      (strcmp(param, (const char *) "tau-beta") == 0) || (strcmp(param, (const char *) "rho") == 0)){
    nsamples = dims[0];
    if (burn >= nsamples){
      printf("Burnin exceeds number of samples\n");
      exit(1);
    }
    sample = gsl_vector_calloc(nsamples - burn);
    printf("parameter dimensions for %s: samples = %d\n", 
	   param, nsamples);
    calcsimple(param, dataspace, dataset, sample, hist, credint, nsamples, burn, nbins, sumtype, wids);
  }
  else {
    printf("The parameter you specified does not exist\n");
    exit(1);
  }
}

/* calculate summaries for parameters indexed by locus and
   population */
void calclocuspop(const char * param, hid_t dataspace, hid_t dataset, gsl_vector * sample, gsl_histogram * hist,
		  double credint, int nloci, int nsamples, int npop, int burn, int nbins, int sumtype, int wids){
  int i, j;
  hid_t mvector;
  hsize_t start[3];  /* Start of hyperslab */
  hsize_t count[3] = {1,1,1};  /* Block count */
  hsize_t block[3];
  herr_t ret;
  hsize_t samdim[1];

  /* create vector for buffer */
  samdim[0] = nsamples - burn;
  mvector = H5Screate_simple(1, samdim, NULL);
  /* loop through population, then locus */
  for (j=0; j<npop; j++){
    for (i=0; i<nloci; i++){
      if (wids == 1){
	fprintf(outfp,"%s_loc_%d_pop_%d,",param,i,j);
      }
      start[0] = i; start[1] = burn; start[2] = j;
      block[0] = 1; block[1] = nsamples - burn; block[2] = 1;
      ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, block);
      ret = H5Dread(dataset, H5T_NATIVE_DOUBLE, mvector, dataspace, H5P_DEFAULT, sample->data);
      /* estimate and credible intervals */
      if (sumtype == 0){
	calcci(sample, credint, (nsamples - burn));
      }
      /* generate histogram */
      else if (sumtype == 1){
	calchist(sample, hist, (nsamples - burn), nbins);
      }
      /* write to text file */
      else if (sumtype == 2){
	writetext(sample, (nsamples - burn));
      }
    }
  } 
}

/* calculate summaries for parameters indexed by locus or individual */
void calcparam(const char * param, hid_t dataspace, hid_t dataset, gsl_vector * sample, gsl_histogram * hist,
	       double credint, int nparam, int nsamples, int burn, int nbins, int sumtype, int wids){
  int i;
  hid_t mvector;
  hsize_t start[2];  /* Start of hyperslab */
  hsize_t count[2] = {1,1};  /* Block count */
  hsize_t block[2];
  herr_t ret;
  hsize_t samdim[1];

  /* create vector for buffer */
  samdim[0] = nsamples - burn;
  mvector = H5Screate_simple(1, samdim, NULL);
  /* loop through parameters */
 
  for (i=0; i<nparam; i++){
    if (wids == 1){
      fprintf(outfp,"%s_param_%d,",param,i);
    }
    start[0] = i; start[1] = burn;
    block[0] = 1; block[1] = (nsamples - burn);
    ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, block);
    ret = H5Dread(dataset, H5T_NATIVE_DOUBLE, mvector, dataspace, H5P_DEFAULT, sample->data);
    /* estimate and credible intervals */
    if (sumtype == 0){
      calcci(sample, credint, nsamples - burn);
    }
    /* generate histogram */
    else if (sumtype == 1){
      calchist(sample, hist, (nsamples - burn), nbins);
    }
    /* write to text file */
    else if (sumtype == 2){
      writetext(sample, (nsamples - burn));
    }
  }
}

/* calculate summaries for single parameters */
void calcsimple(const char * param, hid_t dataspace, hid_t dataset, gsl_vector * sample, gsl_histogram * hist,
		double credint, int nsamples, int burn, int nbins, int sumtype, int wids){
  hid_t mvector;
  hsize_t start[1];  /* Start of hyperslab */
  hsize_t count[1] = {1};  /* Block count */
  hsize_t block[1];
  herr_t ret;
  hsize_t samdim[1];

  /* create vector for buffer */
  samdim[0] = nsamples - burn;
  mvector = H5Screate_simple(1, samdim, NULL);
  /* loop through parameters */
 
  if (wids == 1){
    fprintf(outfp,"%s,",param);
  }
  start[0] = burn;
  block[0] = (nsamples - burn);
  ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, block);
  ret = H5Dread(dataset, H5T_NATIVE_DOUBLE, mvector, dataspace, H5P_DEFAULT, sample->data);
  /* estimate and credible intervals */
  if (sumtype == 0){
    calcci(sample, credint, (nsamples - burn));
  }
  /* generate histogram */
  else if (sumtype == 1){
    calchist(sample, hist, (nsamples - burn), nbins);
  }
  /* write to text file */
  else if (sumtype == 2){
    writetext(sample, (nsamples - burn));
  }
}

/* Function calculates and prints mean, median, and ci */
void calcci(gsl_vector * sample, double credint, int nsamples){
  double lb = (1.0 - credint)/2.0;
  double ub = 1.0 - lb;

  /* sort samples in ascending order */
  gsl_sort_vector(sample);
  fprintf(outfp,"%.6f,", gsl_stats_mean(sample->data, 1, nsamples));
  fprintf(outfp,"%.6f,", gsl_stats_median_from_sorted_data(sample->data, 1, nsamples));
  fprintf(outfp,"%.6f,", gsl_stats_quantile_from_sorted_data(sample->data, 1, nsamples,lb));
  fprintf(outfp,"%.6f\n", gsl_stats_quantile_from_sorted_data(sample->data, 1, nsamples,ub));
}

/* Function generates and prints a histogram of the posterior */
void calchist(gsl_vector * sample,  gsl_histogram * hist, int nsamples, int nbins){
  int i;
  double min = 0, max = 0;
  double upper = 0, lower = 0, midpoint = 0;
  
  min = gsl_vector_min(sample);
  max = gsl_vector_max(sample);
  max += max * 0.0001; /* add a small bit to max, as upper bound is exclusive */
  (void) gsl_histogram_set_ranges_uniform(hist, min, max);
  
  /* generate histogram */
  for (i=0; i<nsamples; i++){
    (void) gsl_histogram_increment(hist, gsl_vector_get(sample, i));
  }
  
  /* write histogram */
  for (i=0; i<nbins; i++){
    gsl_histogram_get_range(hist, i, &lower, &upper);
    midpoint = (lower + upper) / 2.0;
    if (i == 0){
      fprintf(outfp,"%.6f,%d", midpoint, (int) gsl_histogram_get(hist, i));
    }
    else{
      fprintf(outfp,",%.6f,%d", midpoint, (int) gsl_histogram_get(hist, i));
    }
  }
  fprintf(outfp, "\n");

}

/* Function prints the ordered samples from posterior as plain text */
void writetext(gsl_vector * sample,  int nsamples){
  int i;

  fprintf(outfp,"%.5f", gsl_vector_get(sample, i));
   for (i=1; i<nsamples; i++){
    fprintf(outfp,",%.5f", gsl_vector_get(sample, i));
  }    
  fprintf(outfp,"\n");
}
