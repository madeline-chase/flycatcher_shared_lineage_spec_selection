import csv
import numpy as np
import allel
from sys import argv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--pi', dest = 'pi_input', help = "Input file containing snp based pi")
parser.add_argument('-c', '--chroms', dest = 'chr_input', help = "File containing list of chromosomes to iterate over")
parser.add_argument('-w', '--win', dest = 'win_input', help = "File containing window coordinates and total callable sites. Format is Chr Start Stop Sites. Coordinates should be 1-based and inclusive.")
parser.add_argument('-s', '--stat', dest = 'calc_stat', choices = ['pi','dxy'], help = "Specify whether input is snp based pi or dxy. Affects only format of output header." )
args = parser.parse_args()

pi_input = args.pi_input
chr_input = args.chr_input
win_input = args.win_input
calc_stat = args.calc_stat

chr_list = []

for line in open(chr_input):
    chr_list.append(line.strip())

with open(win_input) as winput:
    reader1 = csv.reader(winput, delimiter = ' ')
    win_data = list(zip(*reader1))

win_chr = list(win_data[0])
start = list(win_data[1])
stop = list(win_data[2])
total_sites = list(win_data[3])

win_coords = []
win_sites = []

for site in range(0, len(win_chr)):
    site_bases = int(total_sites[site])
    site_start = int(start[site])
    site_stop = int(stop[site])
    win_list = [site_start,site_stop]
    win_coords.append(win_list)
    win_sites.append(site_bases)

#print(win_coords)

with open(pi_input) as inf:
    reader = csv.reader(inf, delimiter = ' ')
    data = list(zip(*reader))

chr = list(data[0])

pos = []
freq = []

for j in range(0, len(chr)):
    pos.append(int(data[1][j]))
    freq.append(float(data[2][j]))

idx = allel.ChromPosIndex(chr, pos)
win_idx = allel.ChromPosIndex(win_chr,start)
if calc_stat == "pi":
    print("Scaff", "Start", "Stop", "Seg_sites","Pi_sum", "Total_sites", "Pi", sep = '\t')
elif calc_stat == "dxy":
    print("Scaff", "Start", "Stop", "Seg_sites","dxy_sum", "Total_sites", "dxy", sep = '\t')

for x in range(0, len(chr_list)):
    current_chr = chr_list[x]
    if current_chr in chr:
        chr_pos = pos[idx.locate_range(current_chr)]
        chr_freq = freq[idx.locate_range(current_chr)]
    else:
        chr_pos = []
        chr_freq = []
    current_win_coords = win_coords[win_idx.locate_range(current_chr)]
    current_win_sites = win_sites[win_idx.locate_range(current_chr)]
    #nnz, windows, counts = allel.windowed_statistic(chr_pos,chr_freq,statistic=np.sum,size = 200000)
    nnz, wins, counts = allel.windowed_statistic(chr_pos,chr_freq,statistic=np.sum,windows = current_win_coords)

    for y in range(0, len(wins)):
        if nnz[y] is np.nan:
            print(current_chr, wins[y][0], wins[y][1], counts[y], 0, current_win_sites[y], 0, sep = '\t')
        else:
            print(current_chr, wins[y][0], wins[y][1], counts[y], nnz[y], current_win_sites[y], nnz[y]/current_win_sites[y],sep = '\t')
