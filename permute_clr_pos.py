import numpy as np
from sys import argv
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = 'input', help = "Input file containing CLR data")
parser.add_argument('-b', '--bed', dest = 'bed_in', help = "File containing bed file with regions to find significant clr in")
parser.add_argument('-o', '--out', dest = 'out_name', help = "Prefix to append to outfile names")
args = parser.parse_args()

input = args.input
bed_file = args.bed_in
out_name = args.out_name

clr_input = input
clr_data = pd.read_csv(clr_input, sep = '\t', names = ['Chr', 'Start', 'Stop', 'Clr','Alpha'])

clr_sub_out_name = '{}_shuffled_clr.bed'.format(out_name)
bed_intersect_out = '{}_bed_intersect_out.bed'.format(out_name)
final_out_name = '{}_permuted_clr_olaps.txt'.format(out_name)

for iter in range(0,1000):
    clr_data['shuff_clr'] = np.random.permutation(clr_data.Clr)

    clr_data_subset = clr_data[clr_data.shuff_clr>46.25]

    clr_data_subset.to_csv(clr_sub_out_name, sep = '\t', header = False, index = False)

    bedtools_str = "bedtools intersect -a {0} -b {1} -c > {2}_bed_intersect_out_{3}.bed".format(bed_file, clr_sub_out_name, out_name, iter)
    #cat_str = "cat {0} {1} > tmp".format(final_out_name, bed_intersect_out)
    #rename_str = "mv tmp {}".format(final_out_name)

    os.system(bedtools_str)
    #os.system(cat_str)
    #os.system(rename_str)

    print('Finished with iteration .. {}'.format(iter))
