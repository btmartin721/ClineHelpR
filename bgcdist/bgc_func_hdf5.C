#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cdf.h>
#include <float.h>

#include "hdf5.h"
#include "bgc.h"

using namespace std;

/*------------------------------------------------------*/
/*                  hdf5 writing functions              */
/*------------------------------------------------------*/

/* write hybrid index, alpha, beta, and log likelihood in hdf5 format */
void writeHDF5basic(paramcont * param, datadimcont * datadim, auxcont * auxvar, double LnL, hdf5cont * hdf5, 
		    int step, int burn, int thin){
  int j;
  // select hyperlab,  H5S_SELECT_SET overwrites old selection, stride is NULL
  hdf5->startr2[0] = 0; 
  hdf5->startr2[1] = (step - burn) / thin;
  hdf5->countr2[0] = 1;
  hdf5->countr2[1] = 1;
  hdf5->blockr2[0] = datadim->nloci;
  hdf5->blockr2[1] = 1;

  hdf5->startr3[0] = 0; 
  hdf5->startr3[1] = (step - burn) / thin;
  hdf5->countr3[0] = 1;
  hdf5->countr3[1] = 1;
  hdf5->countr3[2] = 1;
  hdf5->blockr3[0] = datadim->nloci;
  hdf5->blockr3[1] = 1;
  hdf5->blockr3[2] = 1;

  // write alpha and beta
  for (j=0; j<datadim->npop; j++){
    hdf5->startr3[2] = j;
    hdf5->status = H5Sselect_hyperslab(hdf5->dataspacelocuspop, H5S_SELECT_SET, hdf5->startr3, NULL, hdf5->countr3, hdf5->blockr3);

    gsl_matrix_get_col(auxvar->vecLocin, param->alphaCur, j);
    hdf5->status = H5Dwrite(hdf5->datasetAlpha, H5T_NATIVE_DOUBLE, hdf5->locusvector, hdf5->dataspacelocuspop, H5P_DEFAULT,  
			    auxvar->vecLocin->data);

    gsl_matrix_get_col(auxvar->vecLocin, param->betaCur, 0);
    hdf5->status = H5Dwrite(hdf5->datasetBeta, H5T_NATIVE_DOUBLE, hdf5->locusvector, hdf5->dataspacelocuspop, H5P_DEFAULT,  
			    auxvar->vecLocin->data);
  }
    
  // write hybrid index
  hdf5->blockr2[0] = datadim->nind;
  hdf5->status = H5Sselect_hyperslab(hdf5->dataspaceindividual, H5S_SELECT_SET, hdf5->startr2, NULL, 
				     hdf5->countr2, hdf5->blockr2);
  gsl_matrix_get_col(auxvar->vecIndn, param->hi, 0);
  hdf5->status = H5Dwrite(hdf5->datasetHi, H5T_NATIVE_DOUBLE, hdf5->indvector, hdf5->dataspaceindividual, 
			  H5P_DEFAULT, auxvar->vecIndn->data);

  // write lnL
  hdf5->startr1[0] = (step - burn) / thin;
  hdf5->countr1[0] = 1;
  hdf5->blockr1[0] = 1;
  hdf5->status = H5Sselect_hyperslab(hdf5->dataspacemcmc, H5S_SELECT_SET, hdf5->startr1, NULL, 
				     hdf5->countr1, hdf5->blockr1);
  gsl_vector_set(auxvar->vecUnit, 0, LnL);
  hdf5->status = H5Dwrite(hdf5->datasetLnL, H5T_NATIVE_DOUBLE, hdf5->unitvector, hdf5->dataspacemcmc, 
			  H5P_DEFAULT, auxvar->vecUnit->data);
}

/* write tau parameters in hdf5 format */
void writeHDF5tau(paramcont * param, auxcont * auxvar, double talpha, double tbeta, hdf5cont * hdf5, 
		  int step, int burn, int thin){
      
  // write tau alpha and tau beta
  hdf5->startr1[0] = (step - burn) / thin;
  hdf5->countr1[0] = 1;
  hdf5->blockr1[0] = 1;
  hdf5->status = H5Sselect_hyperslab(hdf5->dataspacemcmc, H5S_SELECT_SET, hdf5->startr1, NULL, 
				     hdf5->countr1, hdf5->blockr1);
  gsl_vector_set(auxvar->vecUnit, 0, talpha);
  hdf5->status = H5Dwrite(hdf5->datasetTauA, H5T_NATIVE_DOUBLE, hdf5->unitvector, hdf5->dataspacemcmc, 
			  H5P_DEFAULT, auxvar->vecUnit->data);
  gsl_vector_set(auxvar->vecUnit, 0, tbeta);
  hdf5->status = H5Dwrite(hdf5->datasetTauB, H5T_NATIVE_DOUBLE, hdf5->unitvector, hdf5->dataspacemcmc, 
			  H5P_DEFAULT, auxvar->vecUnit->data);
}


/* write cline parameter population effects, eta and kappa, in hdf5 format */
void writeHDF5pop(paramcont * param, datadimcont * datadim, auxcont * auxvar, hdf5cont * hdf5, 
		  int step, int burn, int thin){
  int j;
  hdf5->startr3[0] = 0; 
  hdf5->startr3[1] = (step - burn) / thin;
  hdf5->countr3[0] = 1;
  hdf5->countr3[1] = 1;
  hdf5->countr3[2] = 1;
  hdf5->blockr3[0] = datadim->nloci;
  hdf5->blockr3[1] = 1;
  hdf5->blockr3[2] = 1;

  // write alpha and beta
  for (j=0; j<datadim->npop; j++){
    hdf5->startr3[2] = j;
    hdf5->status = H5Sselect_hyperslab(hdf5->dataspacelocuspop, H5S_SELECT_SET, hdf5->startr3, NULL, hdf5->countr3, hdf5->blockr3);

    gsl_matrix_get_col(auxvar->vecLocin, param->etaCur, j);
    hdf5->status = H5Dwrite(hdf5->datasetEta, H5T_NATIVE_DOUBLE, hdf5->locusvector, hdf5->dataspacelocuspop, H5P_DEFAULT,  
			    auxvar->vecLocin->data);

    gsl_matrix_get_col(auxvar->vecLocin, param->kappaCur, 0);
    hdf5->status = H5Dwrite(hdf5->datasetKappa, H5T_NATIVE_DOUBLE, hdf5->locusvector, hdf5->dataspacelocuspop, H5P_DEFAULT,  
			    auxvar->vecLocin->data);
  }
}

/* write rho parameter in hdf5 format */
void writeHDF5rho(auxcont * auxvar, double rho, hdf5cont * hdf5, int step, int burn, int thin){
  // write rho
  hdf5->startr1[0] = (step - burn) / thin;
  hdf5->countr1[0] = 1;
  hdf5->blockr1[0] = 1;
  hdf5->status = H5Sselect_hyperslab(hdf5->dataspacemcmc, H5S_SELECT_SET, hdf5->startr1, NULL, 
				     hdf5->countr1, hdf5->blockr1);
  gsl_vector_set(auxvar->vecUnit, 0, rho);
  hdf5->status = H5Dwrite(hdf5->datasetRho, H5T_NATIVE_DOUBLE, hdf5->unitvector, hdf5->dataspacemcmc, 
			  H5P_DEFAULT, auxvar->vecUnit->data);
}


/* calculates the quantile of each random effect based on the conditional prior and writes to a hdf5 file */
void writeHDF5quants(paramcont * param, datadimcont * datadim, auxcont * auxvar, hdf5cont * hdf5, 
		     int step, int burn, int thin){
  
  int i, j;
  double q, oneNu, oneOmega;

  hdf5->startr2[0] = 0; 
  hdf5->startr2[1] = (step - burn) / thin;
  hdf5->countr2[0] = 1;
  hdf5->countr2[1] = 1;
  hdf5->blockr2[0] = datadim->nloci;
  hdf5->blockr2[1] = 1;

  hdf5->startr3[0] = 0; 
  hdf5->startr3[1] = (step - burn) / thin;
  hdf5->countr3[0] = 1;
  hdf5->countr3[1] = 1;
  hdf5->countr3[2] = 1;
  hdf5->blockr3[0] = datadim->nloci;
  hdf5->blockr3[1] = 1;
  hdf5->blockr3[2] = 1;

  // quants for gamma 
  for(i=0; i<datadim->nloci; i++){
     q = gsl_cdf_gaussian_P(gsl_matrix_get(param->gamma,i,0),sqrt(1/param->tauAlphaCur));
     gsl_vector_set(auxvar->vecLocin, i, q);
  }
  hdf5->status = H5Sselect_hyperslab(hdf5->dataspacelocus, H5S_SELECT_SET, hdf5->startr2, NULL, 
				     hdf5->countr2, hdf5->blockr2);
  hdf5->status = H5Dwrite(hdf5->datasetQgamma, H5T_NATIVE_DOUBLE, hdf5->locusvector, hdf5->dataspacelocus, 
			  H5P_DEFAULT, auxvar->vecLocin->data);
  
  // quants for zeta 
  for(i=0; i<datadim->nloci; i++){
    q = gsl_cdf_gaussian_P(gsl_matrix_get(param->zeta,i,0),sqrt(1/param->tauBetaCur));
    gsl_vector_set(auxvar->vecLocin, i, q);
  }
  hdf5->status = H5Dwrite(hdf5->datasetQzeta, H5T_NATIVE_DOUBLE, hdf5->locusvector, hdf5->dataspacelocus, 
			  H5P_DEFAULT, auxvar->vecLocin->data);

  // quantiles for eta and kappa if multiple populations are included
  if(datadim->npop > 1){
    for (j=0; j<datadim->npop; j++){
      hdf5->startr3[2] = j;
      hdf5->status = H5Sselect_hyperslab(hdf5->dataspacelocuspop, H5S_SELECT_SET, hdf5->startr3, NULL, hdf5->countr3, hdf5->blockr3);
      // quantiles for eta
      for(i=0; i<datadim->nloci; i++){
	oneNu = gsl_matrix_get(param->nu,i,0);
	q = gsl_cdf_gaussian_P(gsl_matrix_get(param->etaCur,i,j),sqrt(1/oneNu));
	gsl_vector_set(auxvar->vecLocin, i, q);
      }
      hdf5->status = H5Dwrite(hdf5->datasetQeta, H5T_NATIVE_DOUBLE, hdf5->locusvector, hdf5->dataspacelocuspop, H5P_DEFAULT,  
			      auxvar->vecLocin->data);
      // quantiles for kappa
      for(i=0; i<datadim->nloci; i++){
	oneOmega = gsl_matrix_get(param->omega,i,0);
	q = gsl_cdf_gaussian_P(gsl_matrix_get(param->kappaCur,i,j),sqrt(1/oneOmega));
	gsl_vector_set(auxvar->vecLocin, i, q);
      }
      hdf5->status = H5Dwrite(hdf5->datasetQkappa, H5T_NATIVE_DOUBLE, hdf5->locusvector, hdf5->dataspacelocuspop, H5P_DEFAULT,  
			      auxvar->vecLocin->data);
    }
  }
}

/* calculates the quantile of each random effect based on the conditional prior and writes to a hdf5 file */
void writeHDF5quantsCAR(paramcont * param, datadimcont * datadim, auxcont * auxvar, hdf5cont * hdf5, 
			int step, int burn, int thin){

  int i, l;
  // global and local outlier
  double qg, ql;
  double mu = 0;
  double s_weight = 0;
  
  hdf5->startr2[0] = 0; 
  hdf5->startr2[1] = (step - burn) / thin;
  hdf5->countr2[0] = 1;
  hdf5->countr2[1] = 1;
  hdf5->blockr2[0] = datadim->nloci;
  hdf5->blockr2[1] = 1;

  hdf5->status = H5Sselect_hyperslab(hdf5->dataspacelocus, H5S_SELECT_SET, hdf5->startr2, NULL, 
				     hdf5->countr2, hdf5->blockr2);


  // quantiles for gamma, global than local
  for (i=0; i<datadim->nloci; i++){
    // calculate weights
    mu = 0;
    s_weight = 0;
    for (l=0; l<datadim->nloci; l++){ 
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
    gsl_vector_set(auxvar->vecLocin, i, qg);
    gsl_vector_set(auxvar->vecLocin2, i, ql);
  }
  hdf5->status = H5Dwrite(hdf5->datasetQgamma, H5T_NATIVE_DOUBLE, hdf5->locusvector, hdf5->dataspacelocus, 
			  H5P_DEFAULT, auxvar->vecLocin->data);
  hdf5->status = H5Dwrite(hdf5->datasetQgammaCAR, H5T_NATIVE_DOUBLE, hdf5->locusvector, hdf5->dataspacelocus, 
			  H5P_DEFAULT, auxvar->vecLocin2->data);
  

  // quantiles for zeta, global than local
  for (i=0; i<datadim->nloci; i++){
    // calculate weights
    mu = 0;
    s_weight = 0;
    for (l=0; l<datadim->nloci; l++){ 
      if (i != l){
	s_weight +=  gsl_matrix_get(param->sigmaAll,i,l);
      }
    }
    // cdf for zeta based on the lower tail
    for (l=0; l<datadim->nloci; l++){ 
      if (i != l){
	mu += gsl_matrix_get(param->zeta,l,0) * (gsl_matrix_get(param->sigmaAll,i,l) / s_weight);
      }
    }
    qg = gsl_cdf_gaussian_P(gsl_matrix_get(param->zeta,i,0), sqrt((1/param->tauBetaCur) / s_weight));
    ql = gsl_cdf_gaussian_P(gsl_matrix_get(param->zeta,i,0) - mu, sqrt((1/param->tauBetaCur) / s_weight));
    gsl_vector_set(auxvar->vecLocin, i, qg);
    gsl_vector_set(auxvar->vecLocin2, i, ql);
  }
  hdf5->status = H5Dwrite(hdf5->datasetQzeta, H5T_NATIVE_DOUBLE, hdf5->locusvector, hdf5->dataspacelocus, 
			  H5P_DEFAULT, auxvar->vecLocin->data);
  hdf5->status = H5Dwrite(hdf5->datasetQzetaCAR, H5T_NATIVE_DOUBLE, hdf5->locusvector, hdf5->dataspacelocus, 
			  H5P_DEFAULT, auxvar->vecLocin2->data);
  
}

/* calcualtes and prints interspecific heterozygosity based on ancestry (z) and writes to a hdf5 file */
void writeHDF5het(paramcont * param, datadimcont * datadim, auxcont * auxvar, hdf5cont * hdf5, 
		  int step, int burn, int thin, int allDiploid){

  int i, n;
  double tot, het;
  double phet;
  
  hdf5->startr2[0] = 0; 
  hdf5->startr2[1] = (step - burn) / thin;
  hdf5->countr2[0] = 1;
  hdf5->countr2[1] = 1;
  hdf5->blockr2[0] = datadim->nind;
  hdf5->blockr2[1] = 1;

  hdf5->status = H5Sselect_hyperslab(hdf5->dataspaceindividual, H5S_SELECT_SET, hdf5->startr2, NULL, 
				     hdf5->countr2, hdf5->blockr2);

  if (allDiploid == 1){
    for(n=0; n<datadim->nind; n++){
      tot = 0;
      het = 0;
      for(i=0; i<datadim->nloci; i++){
	if (gsl_matrix_int_get(param->zCur[i].mat,n,0) != gsl_matrix_int_get(param->zCur[i].mat,n,1)){
	  het++;
	}
	tot++;
      }
      phet = het / tot;
      gsl_vector_set(auxvar->vecIndn, n, phet);
    }
  }
  else {
    for(n=0; n<datadim->nind; n++){
      tot = 0;
      het = 0;
      for(i=0; i<datadim->nloci; i++){
	if (gsl_vector_int_get(datadim->ploidyVec,i) == 2){
	  if (gsl_matrix_int_get(param->zCur[i].mat,n,0) != gsl_matrix_int_get(param->zCur[i].mat,n,1)){
	    het++;
	  }
	  tot++;
	}
      }
      phet = het / tot;
      gsl_vector_set(auxvar->vecIndn, n, phet);
    }
  }
  hdf5->status = H5Dwrite(hdf5->datasetHet, H5T_NATIVE_DOUBLE, hdf5->indvector, hdf5->dataspaceindividual, 
			  H5P_DEFAULT, auxvar->vecIndn->data);
}
