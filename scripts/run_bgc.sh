#!/bin/bash

# ------------------------------------------------------------------------------
# --- Run bgc
# ------------------------------------------------------------------------------

mcmc=$3
burnin=$4
thin=$5
gt_uncert=$6
error=$7
linkage=$8

# Arguments set
# -O (0, 1, or 2)     : Format to write MCMC samples (default=0)
#                 :     0 -> use HDF5 output (requires program estpost); 
#                 :     1 -> ascii text
#                 :     2 -> HDF5 and ascii text
# -x (integer)    : Number of MCMC generations
# -n (integer)    : discard n steps as burn-in
# -t (integer)    : sample MCMC chain every t generations
# -N (0 or 1)     : Whether to use genotype uncertainty model (recommended) (default=0)
# -E (float)      : Use sequence error model and genotype uncertainty model (default=0)
#                 :     Error rate can be specified with a float (e.g. 0.001 = .1% error rate)
# -q (0 or 1)     : Whether to calculate cline parameter quantiles (default=0)
# -I (0 or 1)     : Whether to initialize hybrid index/ancestry based on data (default=1); 
#                 :     0 -> initialize
#                 :     1 -> do not initialize
# -p (0 or 1)     : Whether to print default parameters + precision parameters (default=0)
# -m (0 or 1)     : Whether to use linkage model (default=0)
# -M (string)     : Filename for linkage map input file, only used for ICARrho model
# -d (0 or 1)     : Whether all loci are diploid (default=1)
#                 :     0 -> Not diploid; only valid for the known genotype model
#                 :     1 -> all loci are diploid
# -o (0 or 1)     : Whether to assume a constant population-level cline parameter variance for all loci (default=0)
# -i (0 or 1)     : Whether to calculate and print interspecific heterozygosity (default=0)
# -s (0 or 1)     : Whether to apply sum-to-zero constraint on locus cline parameters (default=1)
# -T (0 or float) : If non-zero, use a truncated gamme prior for tau with this upper bound (default=0)
# -u (float)      : MCMC tuning parameter, maximum deviate from uniform for proposed hybrid index  (default=0.1)
# -g (float)      : standard deviation for Gaussian proposal of cline parameter gamma (default=0.05)
# -z (float)      : standard deviation for Gaussian proposal of cline parameter zeta (default=0.05)
# -e (float)      : standard deviation for Gaussian proposal of cline parameters eta and kappa (default=0.02)
# -D (float)      : Maximum distance between loci, free recombination. Only used for linkage model (default=0.5)

# Not sure how BGC manages seeding. So in case it's random by clock time, I staggered them here.
sleepytime=`shuf -i 10-240 -n 1`
sleep "$sleepytime"

if [ $linkage = 0 ]; then
	$HOME/local/bin/bgc -a data/${2}_bgc_p0in.txt -b data/${2}_bgc_p1in.txt -h data/${2}_bgc_admixedin.txt -O 0 -x $mcmc -n $burnin -t $thin -N $gt_uncert -E $error -m $linkage -q 1 -I 1 -p 1 -F results/${2}_bgc_mcmcout_${1};

else
	$HOME/local/bin/bgc -a data/${2}_bgc_p0in.txt -b data/${2}_bgc_p1in.txt -h data/${2}_bgc_admixedin.txt -O 0 -x $mcmc -n $burnin -t $thin -N $gt_uncert -E $error -m $linkage -M data/${2}_bgc_map.txt -q 1 -I 1 -p 1 -F results/${2}_bgc_mcmcout_${1};

fi

# ------------------------------------------------------------------------------
# --- Run estpost
# ------------------------------------------------------------------------------

# Arguments set
# -p (string)      : Name of parameter to summarize (required)
#                  :      Possibilities: 
#                  :          LnL, alpha, beta, eta, eta-quantile,
#                  :          gamma-quantile, gamma-quantile-local, hi, interspecific-het, kappa,
#                  :          kappa-quantile, rho, tau-alpha, tau-beta, zeta-quantile, zeta-quantile-local
# -s (0, 1, or 2)  : Which summary to perform (required)
#                  :      0 -> posterior estimates and credibility intervals; 
#                  :      1 -> histogram of posterior samples; 
#                  :      2 -> plain text
# -w (0 or 1)      : Whether to write parameter identification and headers to file
#                  :      0 -> Do not write to file
#                  :      1 -> Write to file (default)

# -c (float)       : Credibility interval (default=0.95)
# -b (integer)     : Discard b generations as burn-in (based on number of thinned samples) (default=0)
# -h (integer)     : Number of bins for posterior sample histogram (default=20)

# Get log-likelihood
$HOME/local/bin/estpost -i results/${2}_bgc_mcmcout_${1}.hdf5 -p LnL -o results/${2}_bgc_stat_ln1_${1} -s 2 -w 0

# Get alpha
$HOME/local/bin/estpost -i results/${2}_bgc_mcmcout_${1}.hdf5 -p alpha -o results/${2}_bgc_stat_a0_${1} -s 2 -w 0

# Get beta
$HOME/local/bin/estpost -i results/${2}_bgc_mcmcout_${1}.hdf5 -p beta -o results/${2}_bgc_stat_b0_${1} -s 2 -w 0

# Get hybrid index
$HOME/local/bin/estpost -i results/${2}_bgc_mcmcout_${1}.hdf5 -p hi -o results/${2}_bgc_stat_hi_${1} -s 2 -w 0

# Get gamma-quantile (qa)
$HOME/local/bin/estpost -i results/${2}_bgc_mcmcout_${1}.hdf5 -p gamma-quantile -o results/${2}_bgc_stat_qa_${1} -s 2 -w 0

# Get zeta-quantile (qb)
$HOME/local/bin/estpost -i results/${2}_bgc_mcmcout_${1}.hdf5 -p zeta-quantile -o results/${2}_bgc_stat_qb_${1} -s 2 -w 0

exit

