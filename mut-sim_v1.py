#!/usr/bin/python
#
# This program takes the sequence of a DNA contig in Fasta format, duplicates it,
# and introduceds mutations. It offers three modes: 
# (a) Drift SNPs: Any base can be mutated and any other base can be the variant, like in genetic
# drift.
# (b) EMS SNPs: Only CG > AT transitions, as caused by ethyl methanesulfonate.
# (c) Large insertions: A sequence provided by the user is inserted at random positions
#
# The user has to provide:
# (a) Contig sequence in fasta format
# (b) [Only required in mode 'large insertions'] Insertion sequence in fasta format
# (c) Number of mutations to introduce
# (c) Mode chosen
#
# Due to the pourpose of this program, it only supports fasta files with a single
# contig.
#
# It does not check wether the number of mutations asked to introduce is bigger
# than the number of possible bases to mutate. In that case, the program runs
# forever. This was not implemented because that situation is extremely unlikely if 
# simulating real life experiments.
#
# 2016 - David Wilson - dws1985@hotmail.com
#


import argparse
import random


#Parse command arguments
parser = argparse.ArgumentParser()
parser.add_argument('-nbr', action="store", dest='total_nbr_mutations', required=True)
parser.add_argument('-mod', action="store", dest='mutator_mode',
choices=set(('d','e','li')), required=True) #Choose between d, e, li (d = genetic drift SNP mutations, e = EMS SNP mutations GC>AT, li = long insertion (t-dna, transposon))
parser.add_argument('-con', action="store", dest='contig_source', required=True)
parser.add_argument('-ins', action="store", dest='insertion_source')
args = parser.parse_args()

total_nbr_mutations = int(args.total_nbr_mutations)
mutator_mode = args.mutator_mode   
contig_source = args.contig_source
insertion_source = args.insertion_source

#function to parse fasta file (based on one of the Biopython IOs)
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

#read contig fasta file
with open(contig_source) as fp:
	for name_contig, seq_contig in read_fasta(fp):
		#print(name_contig, seq_contig)
		pass

#calculate contig length
contig_length = len(seq_contig)
seq_contig_lowercase = seq_contig.lower()

#Randomly choose mutation positions, add them to a list, and finally sort the list
all_mut_pos = list() #will store the list of mutated positions

iter = 0
while iter < total_nbr_mutations:
	mut_pos = random.randint(1, contig_length - 2) #I deliberately exclude first and last nucleotides
	if mut_pos not in all_mut_pos:
		if mutator_mode == 'd' or mutator_mode == 'li': #if mode is drift, any base can be mutated
			iter +=1
			all_mut_pos.append(mut_pos)
		if mutator_mode == 'e': #if mode is EMS, ony G and A can be mutated
			wt_base = seq_contig_lowercase[mut_pos:mut_pos+1]
			if wt_base == 'g' or wt_base == 'c':
				iter +=1
				all_mut_pos.append(mut_pos)

all_mut_pos.sort()


#If mode is SNPs (drift or EMS), use list 'all_mut_pos' and determine the mutant allele at each position
if mutator_mode == 'd' or mutator_mode == 'e':
	all_mut_info = list() #will store mutated positions, wt base, and mut base

	#Iterate over list, determine the wt base for each position in list, and randomly choose a mt base
	for i in all_mut_pos:
		wt_base = seq_contig_lowercase[i:i+1]

		if mutator_mode == 'd': #if mode is drift, the mutation can be any base
			options_drift = {'a': ['t', 'c', 'g'], 't': ['c', 'g', 'a'], 'c': ['g', 'a', 't'], 'g': ['a', 't', 'c']}
			try:
				possible_mt_bases = options_drift[wt_base] #lookup the alterantive bases of "wt_base" in "options_drift"
				mt_base = possible_mt_bases[random.randint(0, 2)] #randomly pick between the possible mut bases
			except KeyError:
				mt_base = 'n'

		if mutator_mode == 'e': #if mode is ems, only GC>AT possible
			options_ems = {'g': 'a', 'c': 't'}
			try:
				mt_base = options_ems[wt_base]
			except KeyError:
				mt_base = 'n'

		mut_info = (i, wt_base, mt_base) 
		all_mut_info.append(mut_info)
	
	#iterate over all_mut_info list, mutate the wt sequence, and write the mutant version to file	
	mutant_seq_contig = list(seq_contig_lowercase) #I use a list because python strings do not support item assignment
		
	for mut_info in all_mut_info:
		#write to info file
		mutant_seq_contig[mut_info[0]] = mut_info[2] #for each mutation in list "all_mut_info", appropriately substitute in list "mutant_seq"
	
	output2 = open('outputs/mutant_sequence.fa', 'w')
	output2.write(name_contig + '__mutated' + '\n' + ''.join(mutant_seq_contig))
	output2.close()
		

if mutator_mode == 'li':
	#how insertion positions are be handled:
	#  ATCG     .ATCG     A.TCG     AT.CG     ATC.G     ATCG.
	#  01234    0          1          2          3          4

	mutant_seq_contig = ''

	#read insertion fasta file
	with open(insertion_source) as fp:
		for name_ins, seq_ins in read_fasta(fp):
			#print (name_ins, seq_ins)
			pass
	
	#create a new list that contains the insertion points and the beggining and end coordinates 
	all_mut_pos_plus_ends = list(all_mut_pos)
	
	all_mut_pos_plus_ends.insert(0, 0) #add the value "0" at position 0 of the list
	all_mut_pos_plus_ends.append(len(seq_contig)) #add the last base of the contig to the end of the list 
	
	#Iterate over list "all_mut_pos_plus_ends", take the coord values and substring the contig sequence, and build the mutant sequence
	iter_ins = 0
	while iter_ins < len(all_mut_pos_plus_ends):
		try:
			fragment = seq_contig[all_mut_pos_plus_ends[iter_ins]:all_mut_pos_plus_ends[iter_ins + 1]]
			#print fragment_len
			mutant_seq_contig += fragment
			if iter_ins < len(all_mut_pos_plus_ends) - 2:
				mutant_seq_contig += seq_ins
		except IndexError:
			pass
		
		iter_ins +=1
	
	output2 = open('outputs/mutant_sequence.fa', 'w')
	output2.write(name_contig + '__mutated' + '\n' + mutant_seq_contig)
	output2.close()	


#Iterate over list "all_mut_info" and create file with info of mutations
output1 = open('outputs/info_all_mutations.txt', 'w')

#If mode is SNPs, write header with 3 columns and use list "all_mut_info".
if mutator_mode == 'd' or mutator_mode == 'e':
	output1.write('position\twt_base\tmt_base\n')
	for mut_info in all_mut_info:
		pos = int(mut_info[0]) + 1
		output1.write(str(pos) + '\t' + str(mut_info[1]) + '\t' + str(mut_info[2]) + '\n')

#If mode is large insertions, write name of inserted sequence, single column header, and just use list "all_mut_pos".
if mutator_mode == 'li':
	output1.write('name of inserted sequence: ' + name_ins[1:] + '\n\npositions of insertions\n')
	for mut_pos in all_mut_pos:
		output1.write(str(mut_pos) + '\n')
output1.close()


print 'Done!'
