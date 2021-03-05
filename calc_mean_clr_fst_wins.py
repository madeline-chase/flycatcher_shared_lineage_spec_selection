import csv
import numpy as np
import allel
from sys import argv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest = 'SF_out')
parser.add_argument('-c', '--chroms', dest = 'chr_input', help = "File containing list of chromosomes to iterate over")
parser.add_argument('-w', '--wins', dest = 'win_input', help = "File containing window coordinates. Format is Chr Start Stop, with no header")
args = parser.parse_args()

SF_input = args.SF_out
chr_input = args.chr_input
win_input = args.win_input

chr_list = []
for line in open(chr_input):
    chr_list.append(line.strip())

with open(win_input) as winput:
    reader1 = csv.reader(winput, delimiter = '\t')
    win_data = list(zip(*reader1))

win_chr = list(win_data[0])
start = win_data[1]
stop = win_data[2]

win_coords = []

for site in range(0,len(start)):
     site_start = int(start[site])
     site_stop = int(stop[site])
     win_list = [site_start,site_stop]
     win_coords.append(win_list)


with open(SF_input) as inf:
    reader = csv.reader(inf, delimiter='\t')
    data = list(zip(*reader))

chr = list(data[0])


clr = []
pos = []

for j in range(0,len(chr)):
    clr.append(float(data[2][j]))
    pos.append(int(data[1][j]))

idx = allel.ChromPosIndex(chr, pos)
win_idx = allel.ChromPosIndex(win_chr,start)



print('Chr','Start','Stop','Sites','Mean_CLR', sep = '\t')
for x in range(0, len(chr_list)):
    current_chr = chr_list[x]

    if current_chr in chr:
        chr_pos = pos[idx.locate_range(current_chr)]
        chr_clr = clr[idx.locate_range(current_chr)]
    else:
        chr_pos = []
        chr_freq_counts = []

    current_win_coords = win_coords[win_idx.locate_range(current_chr)]

    mean_clr, wins, counts = allel.windowed_statistic(chr_pos,chr_clr,statistic=np.mean,windows=current_win_coords, fill = "NA")

    for win in range(0,len(mean_clr)):
        print(current_chr,wins[win][0],wins[win][1],counts[win],mean_clr[win], sep = '\t')

