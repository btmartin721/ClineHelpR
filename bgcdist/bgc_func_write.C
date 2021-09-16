#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <float.h>

#include "bgc.h"

using namespace std;

/*------------------------------------------------------*/
/*             mcmc sample writing functions            */
/*------------------------------------------------------*/

/* function to write main mcmc results */
void printSamples(paramcont * param, datadimcont * datadim, double LnL, 
		  FILE * fpHi, FILE * fpAlpha, FILE * fpBeta, FILE * fpLnL){
  int i, j, n;
  
  fprintf(fpHi, "%.5f", gsl_matrix_get(param->hi,0,0));
  //fprintf(fpFis, "%.5f", gsl_vector_get(param->fis,0));
  for(n=1; n<datadim->nind; n++){
    fprintf(fpHi, ",%.5f", gsl_matrix_get(param->hi,n,0));
    //fprintf(fpFis, ",%.5f", gsl_vector_get(param->fis,n));
  }
  fprintf(fpHi, "\n");
  //fprintf(fpFis, "\n");

  for(i=0; i<datadim->nloci; i++){
    fprintf(fpAlpha, "%i", i);
    for(j=0; j<datadim->npop; j++){
      fprintf(fpAlpha, ",%.5f", gsl_matrix_get(param->alphaCur,i,j));
    }
    fprintf(fpAlpha, "\n");
  }

  for(i=0; i<datadim->nloci; i++){
    fprintf(fpBeta, "%i", i);
    for(j=0; j<datadim->npop; j++){
      fprintf(fpBeta, ",%.5f", gsl_matrix_get(param->betaCur,i,j));
    }
    fprintf(fpBeta, "\n");
  }
  
  fprintf(fpLnL, "%.5f\n", LnL);

}

/* function to write tau mcmc samples */
void printTaus(double tauAlphaCur, double tauBetaCur, FILE * fpTau){
  
  fprintf(fpTau, "%.5f,%.5f\n", tauAlphaCur,tauBetaCur);
}

/* function to write nested pop eta and kappa mcmc samples */
void printPop(gsl_matrix * etaCur, gsl_matrix * kappaCur, int npop, int nloci, 
	       FILE * fpEta, FILE * fpKappa){
  int i, j;
  for(i=0; i<nloci; i++){
    fprintf(fpEta, "%i", i);
    fprintf(fpKappa, "%i", i);
    for(j=0; j<npop; j++){
      fprintf(fpEta, ",%.5f", gsl_matrix_get(etaCur,i,j));
      fprintf(fpKappa, ",%.5f", gsl_matrix_get(kappaCur,i,j));
    }
    fprintf(fpEta, "\n");
    fprintf(fpKappa, "\n");
  }
}
