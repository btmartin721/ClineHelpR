#!/usr/bin/env python

"""
Script to convert VCF files to BGC (Bayesian Genomic Cline) format, using the genotype uncertainties.

Uses the read depth for each allele and each individual
which is output in the ipyrad VCF file.

BGC format looks like this:

Parental population files:

locus_1
22 4
1 21
4 55
33 0
locus_2
33 5
22 3
0 1
0 0

Admixed population files:
locus_1
pop_0
0 0
0 43
33 5
0 33
locus_2
pop_0
22 3
33 5
33 0
0 0


So: the loci contain 4 individuals, and each column represents the read depth
for each allele. Must be bi-allelic data. The admixed file requires a population ID line.

The ipyrad VCF output file contains a column that includes:
GT (Genotype):DP (total reads):CATG (# reads per allele) delimited by a colon.

Dependencies:
	PyVCF (I used version 0.6.8)
"""

# Load necessary modules.
import argparse
import sys
import vcf # PyVCF module


def main():

	# Parse command-line arguments.
	args = Get_Arguments()

	# Read popmap into dictionary. {'indID': 'popID'}
	popmap = read_popmap(args.popmap)
	popsamples = get_samples_by_pop(popmap, args.admixed, args.p1, args.p2)

	# Read VCF using PyVCF module.
	vcf_reader = vcf.Reader(open(args.vcf), "r")

	admix_file = "{}_admixedin.txt".format(args.outprefix)
	p1_file = "{}_p0in.txt".format(args.outprefix)
	p2_file = "{}_p1in.txt".format(args.outprefix)
	loci_order_file = "{}_loci.txt".format(args.outprefix)
	linkage_file = "{}_map.txt".format(args.outprefix)

	# Open file handles.
	admix = open(admix_file, "w")
	p1 = open(p1_file, "w")
	p2 = open(p2_file, "w")
	loci = open(loci_order_file, "w")
	loci.write("#CHROM POS\n")
	linkage_fh = open(linkage_file, "w")

	# Initialize variables used in for loop.
	previous_chrom = None
	chrom_number = 1
	pos_list = list()
	locus_list = list()

	# Count number of loci in VCF file.
	nrecords = 0
	for record in vcf_reader:
		nrecords += 1
	print("Processing {} records in VCF file...\n".format(nrecords))

	vcf_reader = vcf.Reader(open(args.vcf), "r")

	# For each locus:
	for j, record in enumerate(vcf_reader, start = 1):
		locus = "locus_{}".format(j)

		# For linkage map. Check if CHROM == previous CHROM and if so, add to pos_list.
		# If CHROM is different than previous, add the previous
		# chromosome as a chunk. If last record, add current at end.
		if previous_chrom is None:
			pos_list.append(record.POS)
			locus_list.append(locus)
			diff_counter = 0
			
		elif previous_chrom and str(record.CHROM) == str(previous_chrom):
			pos_list.append(record.POS)
			diff_counter = 0
			locus_list.append(locus)
			# If last record
			if j == nrecords:
				pos_min = min(pos_list) # Get min/max of pos_list.
				pos_max = max(pos_list)
				if args.linkage:
					normalize_linkagemap(pos_list, pos_min, pos_max, chrom_number, linkage_fh, locus_list)
		
		# else current CHROM different than previous_chrom:
		else:
			diff_counter += 1
			pos_min = min(pos_list) # Get min and max of previous chrom.
			pos_max = max(pos_list)
			if j == nrecords: # If last record and different:
				if diff_counter > 1:
					if args.linkage:
						normalize_linkagemap(pos_list, pos_min, pos_max, chrom_number, linkage_fh, locus_list, 1.0)

					locus_list.clear() # Adding last locus.
					locus_list.append(locus)

					# Now add the last record.
					if args.linkage:
						normalize_linkagemap(pos_list, pos_min, pos_max, chrom_number+1, linkage_fh, locus_list, 1.0)
				else:
					if args.linkage:
						normalize_linkagemap(pos_list, pos_min, pos_max, chrom_number, linkage_fh, locus_list)

					locus_list.clear()
					locus_list.append(locus)

					if args.linkage:
						normalize_linkagemap(pos_list, pos_min, pos_max, chrom_number+1, linkage_fh, locus_list, 1.0)

			else: # If different than previous_chrom but not last record:
				if diff_counter > 1 and args.linkage: # If multiple different in a row.
					normalize_linkagemap(pos_list, pos_min, pos_max, chrom_number, linkage_fh, locus_list, 1.0)
				else: # Normalize 0 -> 1 and write to file.
					if args.linkage:
						normalize_linkagemap(pos_list, pos_min, pos_max, chrom_number, linkage_fh, locus_list)

			# Reset pos_list and then append current chrom value.
			chrom_number += 1 # Increase chromosome count.
			pos_list.clear()
			locus_list.clear()
			locus_list.append(locus)
			pos_list.append(record.POS)

		# Get ref and alt alleles.
		ref = record.REF
		alt = record.ALT

		# Make sure all sites are bi-allelic. Required for BGC.
		if len(alt) > 1:
			raise ValueError("All SNPs must be bi-allelic. >2 alleles detected for locus {} {} ...Terminating execution.\n".format(record.CHROM, record.POS))

		# Get alternate allele. pyvcf has it as a list.
		alt = str(alt[0])

		write_output(record, popsamples, ref, alt, locus, args.outprefix, admix, p1, p2)

		loci.write("{} {}\n".format(record.CHROM, record.POS))

		# For next iteration.
		previous_chrom = record.CHROM
		previous_pos = record.POS

	admix.close()
	p1.close()
	p2.close()
	loci.close()
	linkage_fh.close()

	print("DONE!\n\n")

def normalize_linkagemap(mylist, nmin, nmax, number, fh, loc, one_pos=None):
	"""
	Normalize chromosome positions from 0 to 1 and write to file.
	Input:
		mylist: list of positions for each chromosome.
		nmin: minimum value in mylist.
		nmax: maximum value in mylist.
		number: chromosome number.
		fh: filehandle to write output to.
		loc: list of strings identifying the locus: "locus_N".
		one_pos: 1.0 if chromosome only has one locus; else==None.
	"""
	if one_pos:
		fh.write("{} {} {}\n".format(loc[0], number, one_pos))
	else:
		for i, val in enumerate(mylist):
			mylist[i] = (val - nmin) / (nmax - nmin)
			fh.write("{} {} {:.20f}\n".format(loc[i], number, mylist[i]))

def get_allele_depth(record, pop, ref, alt, sampledict):
	"""
	Get read depths for each allele of each sample. Should be within for loop.
	Input:
		record: pyvcf reader object.
		pop: Population ID that is dictionary key in sampledict.
		ref: reference allele
		alt: alternate allele.
		sampledict: dict(list) {popID: [inds]}
	"""
	possible = ["C", "A", "T", "G"]

	# For each sample:
	result = list()
	for call in record.samples:

		# CallData from PyVCF has depth counts as comma-delimited.
		alleles = call.data[2].split(",")

		# cast depth counts to integers.
		alleles = [int(x) for x in alleles]

		# Make into dict object: {'C': DepthCount, 'A': DC, 'T': DC, 'G': DC}
		allele_depth = dict(zip(possible, alleles))

		try:
			idx = sampledict[pop].index(call.sample)
			result.append("{}  {}".format(allele_depth[ref], allele_depth[alt]))
		except:
			continue

	return result


def write_output(record, sampledict, ref, alt, locus, prefix, admix, p1, p2):
	"""
	Write to the three output files: admixed, p1, and p2. This function should be within a for loop.
	Input:
		record: pyvcf reader object.
		sampledict: dict(list) {population: [inds]}
		ref: reference allele (string)
		alt: alternate allele (string)
		locus: locus name (string)
		prefix: prefix for output files. Specified by user at command-line
		admix: admixed file handle.
		p1: p1 file handle.
		p2: p2 file handle.
	Returns:
		Writes to three output files.

	"""
	# Get possible alleles in order of CATG.
	#print(sampledict["P2"])
	# For each sample: get depth counts.
	admix_output = get_allele_depth(record, "Admixed", ref, alt, sampledict)
	p1_output = get_allele_depth(record, "P1", ref, alt, sampledict)
	p2_output = get_allele_depth(record, "P2", ref, alt, sampledict)

	admix.write("{}\npop_0\n".format(locus))
	p1.write("{}\n".format(locus))
	p2.write("{}\n".format(locus))
	#print("{} {}\n".format(allele_depth[ref], allele_depth[alt]))
	for ind in admix_output:
		admix.write("{}\n".format(ind))

	for ind in p1_output:
		p1.write("{}\n".format(ind))

	for ind in p2_output:
		p2.write("{}\n".format(ind))

def get_samples_by_pop(d, admix, p1, p2):
	"""
	Finds samples associated with each of the three user-specified populations:
	Input:
		d: popmap dictionary {indids: popids}.
		admix: admixed population string specified by user at command-line.
		p1: p1 population string specified by user at command-line.
		p2: p2 population string specified by user at command-line.
	Returns:
		popdict: dictionary {'Admixed': [inds], 'P1': [inds], 'P2': [inds]}
	"""
	popdict = dict()
	admix_inds = list()
	p1_inds = list()
	p2_inds = list()
	for k,v in d.items():
		if str(v) == str(admix):
			admix_inds.append(k)
		elif str(v) == str(p1):
			p1_inds.append(k)
		elif str(v) == str(p2):
			p2_inds.append(k)
		else:
			raise ValueError("\n\nPopulation {} is different than P1, P2 and Admixed! Aborting.")

	# Error handling if wrong population specified.
	if not admix_inds:
		raise ValueError("\n\nPopulation {} was not found in the popmap file! Aborting\n\n".format(admix))
	if not p1_inds:
		raise ValueError("\n\nPopulation {} was not found in the popmap file! Aborting\n\n".format(p1))
	if not p2_inds:
		raise ValueError("\n\nPopulation {} was not found in the popmap file! Aborting\n\n".format(p2))

	popdict["Admixed"] = admix_inds
	popdict["P1"] = p1_inds
	popdict["P2"] = p2_inds

	# Sanity checks. Print number of individuals in each population.
	print("\n\nP1 population has {} individuals...\n".format(len(popdict["P1"])))

	print("P2 population has {} individuals...\n".format(len(popdict["P2"])))

	print("Admixed populalation has {} individuals...\n\n".format(len(popdict["Admixed"])))

	return popdict

def read_popmap(filename):
	"""
	Read a population map file.
	Input:
		filename: file should be two-columns, tab-separated. IndIDs\tPopIDs. No header.
	Returns:
		results: dict of {indIDs: popIDs}
	"""
	results = dict()
	with open(filename, "r") as popfile:
		for line in popfile:
			line = line.strip()
			inds, pops = line.split()
			results[inds] = pops
	return results

def Get_Arguments():
	"""
	Parse command-line arguments. Imported with argparse.
	Returns: object of command-line arguments.
	"""
	parser = argparse.ArgumentParser(description="Convert VCF file to BGC format (with genotype uncertainties). Currently only handles three populations maximum (P1, P2, and Admixed).", add_help=False)

	required_args = parser.add_argument_group("Required Arguments")
	optional_args = parser.add_argument_group("Optional Arguments")

	## Required Arguments
	required_args.add_argument("-v", "--vcf",
								type=str,
								required=True,
								help="Input VCF file")
	required_args.add_argument("-m", "--popmap",
								type=str,
								required=True,
								help="Two-column tab-separated population map file: inds\tpops. No header line.")
	required_args.add_argument("--p1",
								type=str,
								required=True,
								help="Parental population 1")
	required_args.add_argument("--p2",
								type=str,
								required=True,
								help="Parental population 2")
	required_args.add_argument("--admixed",
								type=str,
								required=True,
								help="Admixed population (limit=1 population)")
	optional_args.add_argument("-o", "--outprefix",
								type=str,
								required=False,
								default="bgc",
								help="Specify output prefix for BGC files.")
	optional_args.add_argument("-l", "--linkage",
								default = True, action="store_false")
	optional_args.add_argument("-h", "--help",
								action="help",
								help="Displays this help menu")

	if len(sys.argv)==1:
		print("\nExiting because no command-line options were called.\n")
		parser.print_help(sys.stderr)
		sys.exit(1)

	args = parser.parse_args()
	return args

if __name__ == "__main__":
	main()
