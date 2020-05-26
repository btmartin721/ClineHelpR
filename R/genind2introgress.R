#!/usr/local/bin/R
#' Aggregate BGC Output
#'
#' This function aggregates multiple BGC runs into one file. Runs are joined
#' by columns and must all contain the same number of post-burnin MCMC
#' iterations. It assumes that input files are contained in a single directory
#' and have specific file suffixes. See the README.md file for details on
#' file suffixes. Files can have population prefixes.
#'
#' @param results.dir Path to a directory containing all the BGC results
#' @param prefix Prefix to the input files; e.g., "population1"
#' @param thin Thin to every n MCMC samples. E.g. thin=2 cuts samples in half
#' @param discard Discard first N samples from each BGC run
#' @return A list of data.frames (each data.frame is a BGC parameter)
#' @export
#' @examples
#' aggregate.results <- combine_bgc_output(results.dir = "./results",
#'                                         prefix = "population1")
#'
#' aggregate.results <- combine_bgc_output(results.dir = "./results",
#'                                         prefix = "population2",
#'                                         thin = 2)
#'
#' aggregate.results <- combine_bgc_output(results.dir = "./results",
#'                                         prefix = "pop3",
#'                                         discard = 2000)
genind2introgress <- function(
	genind,
	popA,
	popB,
	admix,
	missingPerLoc=0.5,
	missingPerInd=0.5,
	subset=NULL,
	drop=TRUE
){
	
	ret<-list()
	
	all_pops<-c(popA, popB, admix)
	
	#drop individuals not in selected pops
	sub<-genind[pop(genind) %in% all_pops,]
	
	#drop loci having more than missingPerLoc missing data
	sub<-filter_missingByLoc(sub, prop=missingPerLoc)
	
	#drop individuals having more than missingPerInd missing data
	#NOTE: Also dropping loci that become monomorphic
	sub<-filter_missingByInd(sub, prop=missingPerInd, drop=drop)
	
	#if necessary, randomly subset loci to the number requested
	if (!is.null(subset)){
		sub<-filter_randomSubsetLoci(sub, sample=subset)
	}
	
	#begin parsing remaining data to create inputs for introgress 
	#set up loci information table ()
	loci.data<-data.frame(cbind(c(names(genind$loc.n.all)), genind$type), stringsAsFactors = FALSE)
	colnames(loci.data)<-c("locus", "type")
	loci.data[loci.data["type"]=="codom", "type"]<-"C"
	loci.data[loci.data["type"]=="dom", "type"]<-"D"
	#print(loci.data)
}

###############################################
#Utility functions
filter_missingByLoc<-function(genind, prop=0.5){
	missing<-propTyped(genind, by="loc")
	#print(genind[loc=c(missing>prop)])
	return(genind[loc=c(missing>prop)])
}

filter_missingByInd<-function(genind, prop=0.5, drop=TRUE){
	missing<-propTyped(genind, by="ind")
	#print(missing>prop)
	return(genind[c(missing>prop), drop=drop])
}

filter_randomSubsetLoci<-function(genind, sample){
	if (sample < nLoc(genind)){
		locs <- sort(sample(1:nLoc(genind), sample, replace=FALSE))
		return(genind[loc=c(locs)])
	}
}

library(adegenet)
data(nancycats)
dat<-genind2introgress(nancycats, popA=c("P01", "P02", "P03"), 
                      popB=c("P15", "P14", "P13"), 
                      admix=c("P10", "P11", "P09"),
                      missingPerInd=0.5, missingPerLoc=0.5,
                      subset=5)
