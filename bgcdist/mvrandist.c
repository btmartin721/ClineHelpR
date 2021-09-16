/***************************************************************************************
 *  Multivariate Gaussian probability density function
 *  Multivariate Gaussian random variate
 *  Multivariate Student-t probability density function
 *  Multivariate Student-t random variate
 *  Wishart random variate
 *  Using the GSL - GNU Scientific Library. Visit www.gnu.org/software/gsl
 *
 *  Copyright (C) 2007  Ralph dos Santos Silva
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  AUTHOR 
 *     Ralph dos Santos Silva
 *     September, 2007
***************************************************************************************/
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
/*****************************************************************************************************************/
int rmvnorm(const gsl_rng *r, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *output){
/* This function returns a random variate from the multivariate Gaussian distribution.  
*	 mean   = the mean vector, size n
*	 var    = variance-covariance matrix, dimension n x n
*	 output = vector of sampled values, size n
*/
unsigned int k;
gsl_matrix *work;

if( (mean->size != var->size1) || (var->size1 != var->size2) || (mean->size != output->size) ){
  GSL_ERROR("Incompatible dimensions in rmvnorm", GSL_EINVAL);
}

work = gsl_matrix_alloc(mean->size,mean->size);

gsl_matrix_memcpy(work,var);
gsl_linalg_cholesky_decomp(work);

for(k=0; k<mean->size; k++)
	gsl_vector_set( output, k, gsl_ran_ugaussian(r) );

gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, output);
gsl_vector_add(output,mean);

gsl_matrix_free(work);

return GSL_SUCCESS;
}
/*****************************************************************************************************************/
double dmvnorm(const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var){
/* This function computes the probability density p(x) at x for a multivariate Gaussian distribution.
*	 mean   = the mean vector, size n
*	 var    = variance-covariance matrix, dimension n x n
*/
int s;
double ax,ay;
gsl_vector *ym, *xm;
gsl_matrix *work,*winv;
gsl_permutation *p;

if( (x->size != mean->size) || (mean->size != var->size1) || (var->size1 != var->size2)){
  GSL_ERROR("Incompatible dimensions in dmvnorm", GSL_EINVAL);
}

work = gsl_matrix_alloc(mean->size,mean->size);
winv = gsl_matrix_alloc(mean->size,mean->size);
p = gsl_permutation_alloc(mean->size);

gsl_matrix_memcpy( work, var );
gsl_linalg_LU_decomp( work, p, &s );
gsl_linalg_LU_invert( work, p, winv );
ax = gsl_linalg_LU_det( work, s );
gsl_matrix_free( work );
gsl_permutation_free( p );

xm = gsl_vector_alloc(mean->size);
gsl_vector_memcpy( xm, x);
gsl_vector_sub( xm, mean );
ym = gsl_vector_alloc(mean->size);
gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
gsl_matrix_free( winv );
gsl_blas_ddot( xm, ym, &ay);
gsl_vector_free(xm);
gsl_vector_free(ym);
ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),mean->size)*ax );

return ay;
}
/*****************************************************************************************************************/
int rmvt(const gsl_rng *r, const gsl_vector *location, const gsl_matrix *scale, const double dof, gsl_vector *output){
/* This function returns a random variate from the multivariate Student-t distribution.
*	 location  = the location vector, size n
*	 scale     = scale matrix, dimension n x n
*	 dof	     = degrees of freedom
*	 output    = vector of sampled values, size n
*/
unsigned int k;
double ax = 0.5*dof; 
gsl_matrix *work;

if( (location->size != scale->size1) || (scale->size1 != scale->size2) || (location->size != output->size) ){
  GSL_ERROR("Incompatible dimensions in rmvt", GSL_EINVAL);
}

work = gsl_matrix_alloc(location->size,location->size);

ax = gsl_ran_gamma(r,ax,(1/ax));     /* gamma distribution */

gsl_matrix_memcpy(work,scale);
gsl_matrix_scale(work,(1/ax));       /* scaling the matrix */
gsl_linalg_cholesky_decomp(work);

for(k=0; k<location->size; k++)
	gsl_vector_set( output, k, gsl_ran_ugaussian(r) );

gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, output);
gsl_vector_add(output, location);

gsl_matrix_free(work);

return GSL_SUCCESS;
}
/*****************************************************************************************************************/
double dmvt(const gsl_vector *x, const gsl_vector *location, const gsl_matrix *scale, const double dof){
/* This function computes the probability density p(x) at x for a multivariate Student-t distribution.
*	 location  = the location vector, size n
*	 scale     = scale matrix, dimension n x n
*	 dof	     = degrees of freedom
*/
int s;
double ax,ay,az;
gsl_vector *ym, *xm;
gsl_matrix *work,*winv;
gsl_permutation *p;

if( (x->size != location->size) || (location->size != scale->size1) || (scale->size1 != scale->size2)){
  GSL_ERROR("Incompatible dimensions in dmvt", GSL_EINVAL);
}

az=0.5*(dof + location->size);

work = gsl_matrix_alloc(location->size,location->size);
winv = gsl_matrix_alloc(location->size,location->size);
p = gsl_permutation_alloc(location->size);

gsl_matrix_memcpy( work, scale );
gsl_linalg_LU_decomp( work, p, &s );
gsl_linalg_LU_invert( work, p, winv );
ax = gsl_linalg_LU_det( work, s );
gsl_matrix_free( work );
gsl_permutation_free( p );

xm = gsl_vector_alloc(location->size);
gsl_vector_memcpy( xm, x);
gsl_vector_sub( xm, location );
ym = gsl_vector_alloc(location->size);
gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
gsl_matrix_free( winv );
gsl_blas_ddot( xm, ym, &ay);
gsl_vector_free(xm);
gsl_vector_free(ym);

ay = pow((1+ay/dof),-az)*gsl_sf_gamma(az)/(gsl_sf_gamma(0.5*dof)*sqrt( pow((dof*M_PI),location->size)*ax ));

return ay;
}
/*****************************************************************************************************************/
int rwishart(const gsl_rng *r, const int dof, const gsl_matrix *scale, gsl_matrix *output){
/* This function returns a random variate from the (matrix) Wishart distribution.
*	 dof	     = degrees of freedom
*	 scale     = scale matrix, dimension n x n
*	 output    = matrix of sampled values, dimension n x n
*/
unsigned int k,l;
gsl_matrix *work;

if( (scale->size1 != scale->size2) || (output->size1 != output->size2) || (scale->size1 != output->size1) ){
  GSL_ERROR("Incompatible dimensions in rwishart", GSL_EINVAL);
}

work = gsl_matrix_calloc(scale->size1,scale->size1);

for(k=0; k<scale->size1; k++){
	gsl_matrix_set( work, k, k, sqrt( gsl_ran_chisq( r, (dof-k) ) ) );
	for(l=0; l<k; l++){
		gsl_matrix_set( work, k, l, gsl_ran_ugaussian(r) );
	}
}
gsl_matrix_memcpy(output,scale);
gsl_linalg_cholesky_decomp(output);
gsl_blas_dtrmm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1.0,output,work);
gsl_blas_dsyrk(CblasUpper,CblasNoTrans,1.0,work,0.0,output);

return GSL_SUCCESS;
}

