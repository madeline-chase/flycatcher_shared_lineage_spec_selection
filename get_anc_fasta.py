#! /usr/bin/env python3

### This script is used to print the ancestral state fasta file. It takes a fasta sequence and a bed formatted file with
### the ancestral alleles, and replaces positions in the original fasta with the ancestral allele. The format of the ancestral
### allele file should be: Scaffold, Start, Stop, Ancestral_allele (without a header). The script then prints the new fasta file.
### Run the script: python3 get_anc_fasta.py <orig_fasta.fasta> <ancestral_state.bed>  >  ancestral_state.fasta
### Written 9/11/19 by Madeline Chase


from sys import argv

## Open the orignal fasta file.
orig_fasta = open(argv[1],"r")

##Save the sequences in the original file to a dictionary
fasta_dict = {}
l = list()
header = None
for line in orig_fasta:
	line = line[:-1] ## Get rid of new line character
	if line.startswith('>'): # a new record
        # save the previous record to the dict
		if header:
			fasta_dict[header] = ''.join(l)
			l=[]    # empty the list
		header = line[1:]
	else:
		l.append(line)

# save the last record
fasta_dict[header] = ''.join(l)
orig_fasta.close()


## open the ancestral allele file
anc_state = open(argv[2],"r")

## save the ancestral alleles to a list
for line in anc_state:
	scaff = line.strip().split('\t')[0]
	position = int(line.strip().split('\t')[1])
	anc_allele = line.strip().split('\t')[3]

	fasta_dict[scaff] = fasta_dict[scaff][:position] + anc_allele + fasta_dict[scaff][position+1:]

## For each allele replace the nucleotide in the original fasta with the ancestral allele
#for entry in range(0, len(anc_list)):
#	scaff = anc_list[entry][0]
#	position = int(anc_list[entry][1])
#	anc_allele = anc_list[entry][3]

	#print(fasta_dict[scaff][position])
#	fasta_dict[scaff] = fasta_dict[scaff][:position] + anc_allele + fasta_dict[scaff][position+1:]

## print the header and sequence for each entry in the updated fasta dictionary
for x in fasta_dict:
	print("> " + x, end = '\n')
	print(fasta_dict[x], end = '\n')
