#!/usr/bin/python
# 
# 
# This program simulates the F2 offspring of a cross between two different isogenic lines and the
# selection of individuals with a particular trait. The original purpose was to create mapping
# populations. It takes two chromosomes from Arabidopsis thaliana and creates recombinant versions.
# The recombination frequency distribution is based in Salome et al 2012, Heredity (2012) 108,
# 447-455; doi:10.1038/hdy.2011.95. (The program can be modified to be used with other species
# mainly by modifying the list(s) 'chr_xo_freq'). The positions of the crossovers are random.
# The program then selects the recombinant chromosomes based on whether they carry or not a given
# mutation selected by the user (the phenotype-causing mutation).
# 
# PARAMETERS EXPLANATION:
# 
# -outdir: Name of directory relative to current location where to place the recombinant chromosomes
# 		created.
# 
# -atc (integer) [1, 2, 3, 4, 5]: Chromosome number from Arabidopsis thaliana. Used to apply the
# 		recombination frequency distribution corresponding to that particular chromosome. It must match
# 		the chromosomes provided in 'parmut' and 'parpol'.
# 
# -parmut: Parental A. In modes 'r', 'd', and 'di', this parental must be the mutant parental. It
# 		must be a fasta-formatted file with a single contig. The length of the contig must be equal to the
# 		length of -parpol.
# 
# -parpol: Parental B. Polymorphic to parental A. In modes 'r', 'd', and 'di', this parental must
# 		be the non-mutant parental. It must be a fasta-formatted file with a single contig. The length of
# 		the contig must be equal to the length of -parmut.
# 
# -pcm1 (integer): Position of the causal mutation in -parmut.
# 
# -pcm2 (integer): Position of the causal mutation in -parpol. Only required if mode is 'dr'.
# 
# -mod [re, dm, dw, tr]: Selection mode. 're': Recessive mutation and selection of the mutant phenotype
# 		(recessive phenotypic class); 'dm': Dominant mutation and selection of the mutant phenotype (dominant
# 		phenotypic class); 'dw': Dominant mutation and selection of the wild type phenotype (recessive
# 		phenotypic class); 'tr': Two recessive mutations and selection of the double mutant phenotype.
# 		In the first three modes, the causal mutation is provided in parental -parmut. In the last mode,
# 		each causal mutation is provided in a different parental.
# 
# -nrc (integer): Number of recombinant chromosomes to create. This number corresponds to the number
# 		of chromosomes after the selection has been performed.



import argparse, os, shutil
from random import randint


# Parse command arguments
parser = argparse.ArgumentParser()
parser.add_argument('-outdir', action="store", dest='out_dir', required=True)
parser.add_argument('-atc', action="store", dest='at_chr', type=int, required=True, choices=set((1,2,3,4,5))) # Choose between 1, 2, 3, 4, 5
parser.add_argument('-parmut', action="store", dest='parental_a_mutant', required=True) #mutated genome
parser.add_argument('-parpol', action="store", dest='parental_b_polymorphic', required=True) #polymorphic genome
parser.add_argument('-pcm1', action="store", dest='mut_a_pos', type=int, required=True)
parser.add_argument('-pcm2', action="store", dest='mut_b_pos', type=int)
parser.add_argument('-mod', action="store", dest='selection_mode', required=True, choices=set(('re','dm','dw','tr'))) #Choose between...
parser.add_argument('-nrc', action="store", dest='nbr_rec_chrs', type=int, required=True)
args = parser.parse_args()

out_dir = args.out_dir
at_chr = int(args.at_chr)
parental_a_mutant = args.parental_a_mutant
parental_b_polymorphic = args.parental_b_polymorphic
mut_a_pos = args.mut_a_pos
mut_b_pos = args.mut_b_pos
selection_mode = args.selection_mode # "re" = recessive, "dm" = dominant mt-phe, "dw" = dominant wt-phe, "tr" = double recessive
nbr_rec_chrs = args.nbr_rec_chrs

if selection_mode == 'tr' and mut_b_pos is None:
	quit('Quit. Selected mode is "double recessive" but no position for mutation b was provided. See program description.')		


# Function to parse fasta files
def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith('>'):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))


# Function to create a recombinant chromosome from two parental chromosomes and
# the list 'crossover_positions' with XO positions.
# I decide to group this code in function to be used in the different selection modes
def create_rec_seq():
	rec_chr = ''
	for key,val in enumerate(crossover_positions):	
		if starting_parental == 0:
			if key % 2 == 0:
				try: rec_chr += seq_parental_a[crossover_positions[key]:crossover_positions[key+1]]
				except IndexError: pass
			if key % 2 != 0:
				try: rec_chr += seq_parental_b[crossover_positions[key]:crossover_positions[key+1]]
				except IndexError: pass
		if starting_parental == 1:
			if key % 2 == 0:
				try: rec_chr += seq_parental_b[crossover_positions[key]:crossover_positions[key+1]]
				except IndexError: pass
			if key % 2 != 0:
				try: rec_chr += seq_parental_a[crossover_positions[key]:crossover_positions[key+1]]
				except IndexError: pass					
	return rec_chr

# Function to divide a long string ('data') into chunks of defined length ('batch_size')
def batch_gen(data, batch_size):
	for i in range(0, len(data), batch_size):
		yield data[i:i+batch_size]

# Function to create Fasta file and write the sequence of a recombinant chromosome to it.
def create_rec_chr_file():
	output_file = open(out_dir + '/rec_chr_' + str(iter1 + 1) + '.fa', 'w')
	output_file.write('>rec_chr_' + str(iter1 + 1))
	for chunk in batch_gen(rec_chr, 80): # Write to file a small chunk of the contig sequence in each line
		output_file.write('\n' + chunk)
	output_file.close()


# Create a subdirectory to place the recombinant chromosomes. If it already exists, remove it first
if os.path.exists(out_dir):
	shutil.rmtree(out_dir)
os.makedirs(out_dir)


# Read fasta files of parentals
with open(parental_a_mutant) as fp:
	for name_parental_a, seq_parental_a in read_fasta(fp):
		#print(name_parental_a, seq_parental_a)
		pass

with open(parental_b_polymorphic) as fp:
	for name_parental_b, seq_parental_b in read_fasta(fp):
		#print(name_parental_b, seq_parental_b)
		pass

# Check that both parentals have identical length. If not, quit.
if len(seq_parental_a) != len(seq_parental_b):
	quit('Quit. The lengths of the two parental sequences provided are not identical')

# Determine length of provided chromosome
chr_len = len(seq_parental_a)

# Define recombination frequencies of chromosomes in Arabidopsis.
# Data from Salome et al 2012, Heredity (2012) 108, 447-455; doi:10.1038/hdy.2011.95
# In each sublist [x,y,z]:
#	y - x = Percentage of chromosomes
#	z = Number of crossover events
#	x and y: Used to compare with a random number [1, 100] to simulate chromosomes
#				following the frequencies of crossover events published by
#				Salome et al 2012.
#	x = left interval (open)
#	y = right interval (closed)

if at_chr == 1:
	chr_xo_freq = [[0,14,0],[14,45,1],[45,78,2],[78,93,3],[93,98,4],[98,100,5]]
if at_chr == 2:
	chr_xo_freq = [[0,26,0],[26,68,1],[68,93,2],[93,98,3],[98,99,4],[99,100,5]]
if at_chr == 3:
	chr_xo_freq = [[0,20,0],[20,59,1],[59,87,2],[87,97,3],[97,99,4],[99,100,5]]
if at_chr == 4:
	chr_xo_freq = [[0,24,0],[24,67,1],[67,92,2],[92,98,3],[98,99,4],[99,100,5]]
if at_chr == 5:
	chr_xo_freq = [[0,16,0],[16,50,1],[50,81,2],[81,95,3],[95,99,4],[99,100,5]]


# Create a defined number of recombinant chromosomes 
iter1 = 0
while iter1 < nbr_rec_chrs:
	
	# Reset variables to 0 ('does not contain the mutation') at the beggining of each loop
	chr_carries_mutation_a = False
	chr_carries_mutation_b = False

	# Randomly choose which of the two parentals the program will choose as "seed".
	# The probabilities are 0.5/0.5. Use later.
	starting_parental = randint(0, 1) # 0: start with mutated (parmut) genome, 1: start with polymorphic (parpol) genome
	
	# Create a radom number and compare it to the XO frequency table of the user-chosen chr 
	# The XO frequency table represents the frequency of each number of XOs in a chromosome
	# By comparing a series of random numbers with this table, it can be created a series of XOs
	# numbers that follow exactly the real frequencies
	rand_nbr = randint(1, 100) # 100 different options (values 0 and 100 are included)
	for i in chr_xo_freq:
		if rand_nbr > i[0] and rand_nbr <= i[1]:
			nbr_crossovers = i[2] # This is the number of XOs the current rec chr will have 
	
	# Randomly create the genomic position of each XO event
	iter2 = 0
	crossover_positions = list()
	while iter2 < nbr_crossovers:
		crossover_positions.append(randint(1, chr_len))
		iter2 +=1
	
	# Add the beggining and end coordinates of the chromosome to the list 'crossover_positions'
	crossover_positions.append(0)
	crossover_positions.append(chr_len)
	crossover_positions.sort()
	
	# Determine if the resulting recombinant chr carries the primary (a) causal mutation
	for key,val in enumerate(crossover_positions):
		
		if starting_parental == 0 and key % 2 == 0:
			if mut_a_pos > crossover_positions[key] and mut_a_pos <= crossover_positions[key+1]:
				chr_carries_mutation_a = True
		
		if starting_parental == 1 and key % 2 != 0:
			if mut_a_pos > crossover_positions[key] and mut_a_pos <= crossover_positions[key+1]:
				chr_carries_mutation_a = True
	
	# Determine if the resulting recombinant chr carries the secondary (b) causal mutation
	for key,val in enumerate(crossover_positions):
		
		if starting_parental == 0 and key % 2 != 0:
			if mut_b_pos > crossover_positions[key] and mut_b_pos <= crossover_positions[key+1]:
				chr_carries_mutation_b = True

		if starting_parental == 1 and key % 2 == 0:
			if mut_b_pos > crossover_positions[key] and mut_b_pos <= crossover_positions[key+1]:
				chr_carries_mutation_b = True

	
	# Chromosome selection
	# If recombinant chr contains the desired mutation(s), execute 'create_rec_seq()'
	# and 'create_rec_chr_file()' functions.
	# The first function creates a Fasta file using the info in 'crossover_positions', and
	# the second writes the output to a file.
	
	# Select all chromosomes that carry mutation A (all phenotypically mutant plants)
	if selection_mode == 're':
		if chr_carries_mutation_a == True:
			rec_chr = create_rec_seq()
			create_rec_chr_file()
			iter1 +=1
	
	# Select all chromosomes that carry mutation A and also some that do not, so the final
	# proportion of chrs that carry the mutation is 0.67 (all phenotyoically mutant plants).
	# This happens when selecting mutants with dominant mutation based on phenotype.
	if selection_mode == 'dm':
		if chr_carries_mutation_a == True:
			rec_chr = create_rec_seq()
			create_rec_chr_file()
			iter1 +=1
		else:
			chr_filtering_threshold = randint(1,100)
			if chr_filtering_threshold > 50:
				rec_chr = create_rec_seq()
				create_rec_chr_file()
				iter1 +=1
	
	# Select all chromosomes that do not carry mutation A (all phenotypically wild type plants)
	if selection_mode == 'dw':
		if chr_carries_mutation_a == False:
			rec_chr = create_rec_seq()
			create_rec_chr_file()
			iter1 +=1
	
	# Select all chromosomes that carry mutations A and B (all phenotipically double recessive mutants)
	if selection_mode == 'tr':
		if chr_carries_mutation_a == True and chr_carries_mutation_b == True:
			rec_chr = create_rec_seq()
			create_rec_chr_file()
			iter1 +=1
	

print 'Done!'
