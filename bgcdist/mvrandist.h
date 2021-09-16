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
int rmvnorm(const gsl_rng *r, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *output);
double dmvnorm(const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var);
int rmvt(const gsl_rng *r, const gsl_vector *location, const gsl_matrix *scale, const double dof, gsl_vector *output);
double dmvt(const gsl_vector *x, const gsl_vector *location, const gsl_matrix *scale, const double dof);
int rwishart(const gsl_rng *r, const int dof, const gsl_matrix *scale, gsl_matrix *output);
