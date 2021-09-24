#!/bin/bash

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "Run bgc and estpost."
   echo
   echo "Usage: run_bgc.sh -s <PATH_TO_BGC_SETTINGS_FILE>"
   echo
   echo "Options:"
   echo "-h     Print this help menu."
   echo "-s     Path to BGC settings file."
   echo
}

# Get the options
while getopts "hs:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      s) # BGC settings file
	 BGC_SETTINGS=${OPTARG}
	 ;;
     \?) # incorrect option
	 echo "Error: Invalid option"
	 Help
	 exit;;
   esac
done

# ------------------------------------------------------------------------------
# --- Run bgc
# ------------------------------------------------------------------------------

if [ -z "$BGC_SETTINGS" ]; then
        echo "You must supply the path to the bgc_settings.txt file"
        echo "Usage: run_bgc.sh -s <PATH_TO_BGC_SETTINGS_FILE>"
	Help
        exit 1;
fi

if [ ! -f "$BGC_SETTINGS" ]; then
        echo "Couldn't find specified BGC settings file!"
	echo "Usage: run_bgc.sh -s <PATH_TO_BGC_SETTINGS_FILE>"
	Help
        exit 1;
fi

# Source bgc settings into environment from file
source $BGC_SETTINGS
export $(sed '/^#/d' $BGC_SETTINGS | cut -d= -f1)

# Arguments set
# -O (0, 1, or 2)     : Format to write MCMC samples (default=0)
#                 :     0 -> use HDF5 output (requires program estpost); 
#                 :     1 -> ascii text
#                 :     2 -> HDF5 and ascii text
# -x (integer)    : Number of MCMC generations
# -n (integer)    : discard n steps as burn-in
# -t (integer)    : sample MCMC chain every t generations
# -N (0 or 1)     : Whether to use genotype uncertainty model (recommended for NGS data) (default=0)
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

p0in=${INPUT_DATA_DIR}/${prefix}_p0in.txt
p1in=${INPUT_DATA_DIR}/${prefix}_p1in.txt
admixedin=${INPUT_DATA_DIR}/${prefix}_admixedin.txt

if [ ! -f "$p0in" ]; then
       	echo "Could not find parent0 input file!"
	echo "Check the prefix in the bgc settings file,"
	echo " make sure it is in the location you specified in the bgc_settings file,"
	echo " and make sure it has the filename suffix *_p0in.txt"
	exit 2;
fi

if [ ! -f "$p1in" ]; then
	echo "Could not find parent1 input file!"
        echo "Check the prefix in the bgc settings file,"
        echo " make sure it is in the location you specified in the bgc_settings file,"
        echo " and make sure it has the filename suffix *_p1in.txt"
        exit 2;
fi

if [ ! -f "$admixedin" ]; then
	echo "Could not find admixed input file!"
        echo "Check the prefix in the bgc settings file,"
        echo "make sure it is in the location you specified in the bgc_settings file,"
        echo " and make sure it has the filename suffix *_admixedin.txt"
        exit 2;
fi

# Not sure how BGC manages seeding. So in case it's random by clock time, I staggered them here.
sleepytime=`shuf -i 10-240 -n 1`
sleep "$sleepytime"

echo "Running BGC. Results will be located in the ${RESULTS_DIR}/ directory..."

if [ $linkage -eq 0 ]; then
	${BGC_PATH}/bgc -a ${INPUT_DATA_DIR}/${prefix}_p0in.txt \
		-b ${INPUT_DATA_DIR}/${prefix}_p1in.txt \
		-h ${INPUT_DATA_DIR}/${prefix}_admixedin.txt \
		-O 0 -x $mcmc -n $burnin -t $thin -N $gt_uncert \
		-E $error -m $linkage -q 1 -I 1 -p 1 \
		-F ${RESULTS_DIR}/${prefix}_mcmcout_${run} \
		-g $stddev_gamma -z $stddev_zeta \
		-e $stddev_eta_kappa -u $mcmc_tuning;

else
	linkagefile=${INPUT_DATA_DIR}/${prefix}_map.txt
	if [ ! -f "$linkagefile" ]; then
		echo "Could not find linkage map input file!"
        	echo "Check the prefix in the bgc settings file,"
        	echo " make sure it is in the location you specified in the bgc_settings file,"
        	echo " and make sure it has the filename suffix *_bgc_map.txt"
        	exit 2;
	fi
	${BGC_PATH}/bgc -a ${INPUT_DATA_DIR}/${prefix}_p0in.txt \
		-b ${INPUT_DATA_DIR}/${prefix}_p1in.txt \
		-h ${INPUT_DATA_DIR}/${prefix}_admixedin.txt \
		-O 0 -x $mcmc -n $burnin -t $thin -N $gt_uncert \
		-E $error -m $linkage -M ${INPUT_DATA_DIR}/${prefix}_map.txt \
		-D $linkage_dist -q 1 -I 1 -p 1 \
		-F ${RESULTS_DIR}/${prefix}_mcmcout_${run} \
		-g $stddev_gamma -z $stddev_zeta \
		-e $stddev_eta_kappa -u $mcmc_tuning;
fi

echo "Done running BGC!"

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

echo "Running estpost; Results will be located in the ${RESULTS_DIR}/ directory"

# Get log-likelihood
${ESTPOST_PATH}/estpost -i ${RESULTS_DIR}/${prefix}_mcmcout_${run}.hdf5 \
	-p LnL -o ${RESULTS_DIR}/${prefix}_stat_ln1_${run} -s 2 -w 0

# Get alpha
${ESTPOST_PATH}/estpost -i ${RESULTS_DIR}/${prefix}_mcmcout_${run}.hdf5 \
	-p alpha -o ${RESULTS_DIR}/${prefix}_stat_a0_${run} -s 2 -w 0

# Get beta
${ESTPOST_PATH}/estpost -i ${RESULTS_DIR}/${prefix}_mcmcout_${run}.hdf5 \
	-p beta -o ${RESULTS_DIR}/${prefix}_stat_b0_${run} -s 2 -w 0

# Get hybrid index
${ESTPOST_PATH}/estpost -i ${RESULTS_DIR}/${prefix}_mcmcout_${run}.hdf5 \
	-p hi -o ${RESULTS_DIR}/${prefix}_stat_hi_${run} -s 2 -w 0

# Get gamma-quantile (qa)
$ESTPOST_PATH/estpost -i ${RESULTS_DIR}/${prefix}_mcmcout_${run}.hdf5 \
	-p gamma-quantile -o ${RESULTS_DIR}/${prefix}_stat_qa_${run} -s 2 -w 0

# Get zeta-quantile (qb)
$ESTPOST_PATH/estpost -i ${RESULTS_DIR}/${prefix}_mcmcout_${run}.hdf5 \
	-p zeta-quantile -o ${RESULTS_DIR}/${prefix}_stat_qb_${run} -s 2 -w 0

echo "Done with estpost!"

exit 0

