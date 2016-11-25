#!/usr/bin/python
#
# This program takes two chromosomes from Arabidopsis thaliana that have the exact
# same length and are in Fasta format and creates recombinanat versions.
# It simulates selfing of an heterozygous line that comes from crossing two isogenic lines.
# The recombination frequency distribution is based in Salome et al 2012, Heredity (2012)
# 108, 447-455; doi:10.1038/hdy.2011.95. The position of the crossovers is random.
# The program then selects the recombinanat chromosomes based on whether they carry ir not
# a given mutation selected by the user (the phenotype-causing mutation).
#
# The program provides several modes to select recombinanat chromosomes
# (a) Recesive mutation. Selection of the mutant phenotype (recessive phenotypic class).
# (b) Dominant mutation. Selection of the mutant phenotype (dominant phenotypic class).
# (c) Dominant mutation. Selection of the wild type phenotype (recessive phenotypic class).
# (d) Two recessive mutations. Selection of the double mutant phenotype.
# In the first 3 modes, the causal mutation is provided in parental "parmut". In the last mode,
# each causal mutation is provided in a different parental.
#
# The user has to provide:
# (a) Arabidopsis thaliana chromosome to model (1, 2, 3, 4, 5). It must match the chromosome.
# provided in the next two arguments.
# (b) Parental A (mutant) sequence in fasta format.
# (c) Parental B (polymorphic) sequence in fasta format.
# (d) Position of the causal mutation in parental A.
# (d) [Only required in 'Two recessive mutations' mode] Position of the causal mutation in parental B.
# (e) Selection mode (described above).
# (f) Number of recombinant chromosome to create (after selection).
#
# Due to the pourpose of this program, it only supports fasta files with a single
# contig.
#
# 2016 - David Wilson - dws1985@hotmail.com
#


import argparse
import random


#Parse command arguments
parser = argparse.ArgumentParser()
parser.add_argument('-achr', action="store", dest='at_chr', 
required=True, type=int, choices=set((1,2,3,4,5))) # Choose between 1, 2, 3, 4, 5
parser.add_argument('-parmut', action="store", dest='parental_a_mutant', required=True) #mutated genome
parser.add_argument('-parpol', action="store", dest='parental_b_polymorphic', required=True) #polymorphic genome
parser.add_argument('-mutapos', action="store", dest='mut_a_pos', required=True)
parser.add_argument('-mutbpos', action="store", dest='mut_b_pos', required=True)
parser.add_argument('-smod', action="store", dest='selection_mode',
required=True, choices=set(('r','d','di','dr'))) #Choose between...
parser.add_argument('-nrec', action="store", dest='nbr_rec_chrs', required=True)
args = parser.parse_args()

at_chr = int(args.at_chr)
parental_a_mutant = args.parental_a_mutant
parental_b_polymorphic = args.parental_b_polymorphic
mut_a_pos = int(args.mut_a_pos)
mut_b_pos = int(args.mut_b_pos)
selection_mode = args.selection_mode #"r" = recessive, "d" = dominant mt-phe, "di" = dominant wt-phe, "dr" = double recessive
nbr_rec_chrs = int(args.nbr_rec_chrs)


#Define recombination frequencies of chromosomes in Arabidopsis.
#Data from Salome et al 2012, Heredity (2012) 108, 447-455; doi:10.1038/hdy.2011.95
#In each sublist [x,y,z]:
#	y - x = Percentage of chromosomes
#	z = Number of crossover events
#	x and y: Used to compare with a random number [1, 100] to simulate chromosomes
#				following the frequencies of crossover events published by
#				Salome et al 2012.
#	x = left interval (open)
#	y = right interval (closed)

if at_chr == 1:
	chr_len = 50
	chr_xo_freq = [[0,14,0],[14,45,1],[45,78,2],[78,93,3],[93,98,4],[98,100,5]]
if at_chr == 2:
	chr_len = 50
	chr_xo_freq = [[0,26,0],[26,68,1],[68,93,2],[93,98,3],[98,99,4],[99,100,5]]
if at_chr == 3:
	chr_len = 50
	chr_xo_freq = [[0,20,0],[20,59,1],[59,87,2],[87,97,3],[97,99,4],[99,100,5]]
if at_chr == 4:
	chr_len = 50
	chr_xo_freq = [[0,24,0],[24,67,1],[67,92,2],[92,98,3],[98,99,4],[99,100,5]]
if at_chr == 5:
	chr_len = 50
	chr_xo_freq = [[0,16,0],[16,50,1],[50,81,2],[81,95,3],[95,99,4],[99,100,5]]


#Function to parse fasta files
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


#Function to create a recombinant chromosome from two parental chromosomes and
#the list 'crossover_positions' with XO positions.
#I decide to group this code in function to be used in the different selection modes
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


#Function to create Fasta file and write the sequence of a recombinant chromosome to it.
def create_rec_chr_file():
	output_file = open('output/rec_chr_' + str(iter1 + 1) + '.fa', 'w')
	output_file.write('rec_chr_' + str(iter1 + 1) + '\n' + rec_chr)
	output_file.close()
			

#read fasta files of parentals
with open(parental_a_mutant) as fp:
	for name_parental_a, seq_parental_a in read_fasta(fp):
		#print(name_parental_a, seq_parental_a)
		print 'Sequence of parental A parsed.'

with open(parental_b_polymorphic) as fp:
	for name_parental_b, seq_parental_b in read_fasta(fp):
		#print(name_parental_b, seq_parental_b)
		print 'Sequence of parental B parsed.'


#Check that both parentals have identical length. If not, quit.
if len(seq_parental_a) != len(seq_parental_b):
	quit('Quit. The lengths of the two parental sequences provided are not identical')


#Create a defined number of recombinant chromosomes 
iter1 = 0
while iter1 < nbr_rec_chrs:
	
	#Reset variables to 0 ('does not contain the mutation') at the beggining of each loop
	chr_carries_mutation_a = False
	chr_carries_mutation_b = False

	#Randomly choose which of the two parentals the program will choose as "seed".
	#The probabilities are 0.5/0.5. Use later.
	starting_parental = random.randint(0, 1) #0: start with mutated genome, 1: start with polymorphic genome
	
	#Create a radom number and compare it to the XO frequency table of the user-chosen chr 
	#The XO frequency table represents the frequency of each number of XOs in a chromosome
	#By comparing a series of random numbers with this table, it can be created a series of XOs
	#numbers that follow exactly the real frequencies
	rand_nbr = random.randint(1, 100) #100 different options (values 0 and 100 are included)
	for i in chr_xo_freq:
		if rand_nbr > i[0] and rand_nbr <= i[1]:
			nbr_crossovers = i[2] #This is the number of XOs the current rec chr will have 
	
	#Randomly create the genomic position of each XO event
	iter2 = 0
	crossover_positions = list()
	while iter2 < nbr_crossovers:
		crossover_positions.append(random.randint(1, chr_len))
		iter2 +=1
	
	#add the beggining and end coordinates of the chromosome to the list 'crossover_positions'
	crossover_positions.append(0)
	crossover_positions.append(chr_len)
	crossover_positions.sort()
		
	#Determine if the resulting recombinant chr carries the primary causal mutation
	for key,val in enumerate(crossover_positions):
		
		if starting_parental == 0 and key % 2 != 0:
			if mut_a_pos > crossover_positions[key] and mut_a_pos <= crossover_positions[key+1]:
				chr_carries_mutation_a = True

		if starting_parental == 1 and key % 2 == 0:
			if mut_a_pos > crossover_positions[key] and mut_a_pos <= crossover_positions[key+1]:
				chr_carries_mutation_a = True

	#Determine if the resulting recombinant chr carries the secondary causal mutation
	for key,val in enumerate(crossover_positions):
		
		if starting_parental == 0 and key % 2 == 0:
			if mut_b_pos > crossover_positions[key] and mut_b_pos <= crossover_positions[key+1]:
				chr_carries_mutation_b = True

		if starting_parental == 1 and key % 2 != 0:
			if mut_b_pos > crossover_positions[key] and mut_b_pos <= crossover_positions[key+1]:
				chr_carries_mutation_b = True

	
	#Chromosome selection
	#If recombinant chr contains the desired mutation(s), execute 'create_rec_seq()'
	#and 'create_rec_chr_file()' functions.
	#The first function creates a Fasta file using the info in 'crossover_positions', and
	#the second writes the output to a file.
	
	#Select all chromosomes that carry mutation A (all phenotypically mutant plants)
	if selection_mode == 'r':
		if chr_carries_mutation_a == True:
			rec_chr = create_rec_seq()
			create_rec_chr_file()
			iter1 +=1
			print 'rec_chr: ' + rec_chr
	
	#Select all chromosomes that carry mutation A and also some that do not, so the final
	#proportion of chrs that carry the mutation is 0.67 (all phenotyoically mutant plants).
	#This happens when selecting mutants with dominant mutation based on phenotype.
	if selection_mode == 'd':
		if chr_carries_mutation_a == True:
			rec_chr = create_rec_seq()
			create_rec_chr_file()
			iter1 +=1
			print 'rec_chr: ' + rec_chr
		else:
			chr_filtering_threshold = random.randint(1,100)
			if chr_filtering_threshold > 50:
				rec_chr = create_rec_seq()
				create_rec_chr_file()
				iter1 +=1
				print 'rec_chr: ' + rec_chr
	
	#Select all chromosomes that do not carry mutation A (all phenotypically wild type plants)
	if selection_mode == 'di':
		if chr_carries_mutation_a == False:
			rec_chr = create_rec_seq()
			create_rec_chr_file()
			iter1 +=1
			print 'rec_chr: ' + rec_chr
	
	#Select all chromosomes that carry mutations A and B (all phenotipically double recessive mutants)
	if selection_mode == 'dr':
		if chr_carries_mutation_a == True and chr_carries_mutation_b == True:
			rec_chr = create_rec_seq()
			create_rec_chr_file()
			iter1 +=1
			print 'rec_chr: ' + rec_chr
	

print 'Done!'
