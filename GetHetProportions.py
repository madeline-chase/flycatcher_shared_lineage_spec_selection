## This script will take in a genotype file, and output the proportion of heterzygous genotypes for each SNP, divided into provided groups.
## To run the script, type: python3 GetHetProportions2.py -i <input.geno> -g Par_fem Sample_32,Sample_4,Sample_53,Sample_54,Sample_55,Sample_77,Sample_92 -g Par_male Sample_31,Sample_48,Sample_49,Sample_51,Sample_52,Sample_56,Sample_82,Sample_83 > <het_proportions.txt>
## The -i flag specifies the input file. The -g flag can be provided multiple times, and provides both the (arbitrary) group name, and the Sample names belonging to that group (must match what is in header of file).
## Written by Madeline Chase Oct 24 2019

#! /usr/bin/env python3

from sys import argv
import argparse

## Add argument parser to get groups with sample names and input file
parser = argparse.ArgumentParser()
parser.add_argument('-g', dest = 'groups', action='append', nargs="*", help = "Name of group to calculate heterzygous proportions for, followed by comma separated list of samples in group. Must be provided at least once, can be provided multiple times. Sample names must match file header.")
parser.add_argument('-i', '--input', dest = 'input', help = "Input file in .geno format")
args = parser.parse_args()

group_list = args.groups
input_file = args.input


## Populate two lists with the groups to calculate heterzygous proportions for, and the sample names in each group.
groupNames = []
sampNames = []

for p in args.groups:
    groupNames.append(p[0])
    if len(p) > 1: sampNames.append(p[1].split(','))

## Iterate through .geno file
count = 0
for x in open(input_file):
    line = x.strip().split('\t')
    ## On the first pass through file, get the indexes for each sample from the header
    if count == 0:
        group_index_list = []
        for group in range(0,len(sampNames)):
            group_index = []
            for sample in range(0,len(sampNames[group])):
                group_index.append(line.index(sampNames[group][sample]))
            group_index_list.append(group_index)
        ## Print the header of the output file
        print("Scaff","Position",'\t'.join(groupNames), sep = '\t')
    else:
        ## For each SNP, calculate the proportion of heterozygous genotypes in each group
        print(line[0], line[1], end = '\t',sep='\t')
        for group in range(0,len(group_index_list)):
            
            het_len = 0
            total_len = 0
            for sample in range(0, len(group_index_list[group])):
                samp_index = group_index_list[group][sample]
                A1 = line[samp_index].split('/')[0]
                A2 = line[samp_index].split('/')[1]

                if A1 != A2:
                    het_len += 1
                if A1 and A2 != "N":
                    total_len += 1
            if total_len == 0:
                print('NA', end = '\t')
            else:
                print(het_len/total_len, sep = '\t', end = '\t')

    print('\n', end = '')
    count+=1
