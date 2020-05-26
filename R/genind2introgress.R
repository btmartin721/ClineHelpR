#!/usr/local/bin/R
#'Generate inputs for Introgress
#'
#' This function parses a genind object to create input data frames for 
#' running genomic cline analysis in the Introgress R package
#'
#' @param gen Genind object containing data for parental and admixed populations
#' @param p1 Character or vector containing population names for the p1 group
#' @param p2 Character or vector containing population names for the p2 group
#' @param admix Character or vector containing population names for the admixed individuals
#' @param missingPerInd Proportion of missing genotypes allowed to retain an individual 
#' @param missingPerLoc Proportion of missing genotypes allowed to retain a locus
#' @param subset An optional parameter defining a number of loci to randomly retain
#' @param drop Boolean value defining whether or not to remove monomorphic loci
#' @return A list of data.frames (each data.frame is an introgress input)
#' @export
#' @examples
#' introgress_input <- genind2introgress(genind, p1=c("PopA", "PopB"),
#'                                         p2="PopC", admix=c("PopD1", "PopD2"))
#'
#' introgress_input <- genind2introgress(genind, p1=c("PopA", "PopB"),
#'                                         p2="PopC", admix=c("PopD1", "PopD2"),
#'                                         missingPerInd=0.75, missingPerLoc=0.25,
#'                                         drop=TRUE)
#'
#' introgress_input <- genind2introgress(genind, p1=c("PopA", "PopB"),
#'                                         p2="PopC", admix=c("PopD1", "PopD2"),
#'                                         subset=500)
genind2introgress <- function(
	gen,
	p1,
	p2,
	admix,
	missingPerLoc=0.5,
	missingPerInd=0.5,
	subset=NULL,
	drop=TRUE
){
	
	ret<-list()
	
	all_pops<-c(p1, p2, admix)
	
	#drop individuals not in selected pops
	sub<-gen[pop(gen) %in% all_pops,]
	
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
	loci.data<-data.frame(cbind(c(names(sub$loc.n.all)), sub$type), stringsAsFactors = FALSE)
	colnames(loci.data)<-c("locus", "type")
	loci.data[loci.data["type"]=="codom", "type"]<-"C"
	loci.data[loci.data["type"]=="dom", "type"]<-"D"
	
	#genotype table for p1
	p1_genind<-sub[pop(sub) %in% p1,]
	p1.data<-adegenet::genind2df(p1_genind, sep="/")
	p1.data<-subset(p1.data, select=-c(pop))
	p1.data[is.na(p1.data)] <- "NA/NA"
	
	#genotype table for p2
	p2_genind<-sub[pop(sub) %in% p2,]
	p2.data<-adegenet::genind2df(p2_genind, sep="/")
	p2.data<-subset(p2.data, select=-c(pop))
	p2.data[is.na(p2.data)] <- "NA/NA"
	
	#genotype table for admix population
}

###############################################
#Utility functions
filter_missingByLoc<-function(gen, prop=0.5){
	missing<-propTyped(gen, by="loc")
	#print(genind[loc=c(missing>prop)])
	return(gen[loc=c(missing>prop)])
}

filter_missingByInd<-function(gen, prop=0.5, drop=TRUE){
	missing<-propTyped(gen, by="ind")
	#print(missing>prop)
	return(genind[c(missing>prop), drop=drop])
}

filter_randomSubsetLoci<-function(gen, sample){
	if (sample < nLoc(gen)){
		locs <- sort(sample(1:nLoc(gen), sample, replace=FALSE))
		return(gen[loc=c(locs)])
	}
}

library(adegenet)
data(nancycats)
dat<-genind2introgress(nancycats, p1=c("P01", "P02", "P03"), 
                      p2=c("P15", "P14", "P13"), 
                      admix=c("P10", "P11", "P09"),
                      missingPerInd=0.5, missingPerLoc=0.5,
                      subset=5)
