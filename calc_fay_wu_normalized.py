import csv
import numpy as np
import allel
from sys import argv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--counts', dest = 'counts', help = "Input file containing allele counts with format: Chr, Pos, N_allele, N_chr, Anc_allele, Anc_count, Derived_allele, Derived_count")
parser.add_argument('-c', '--chroms', dest = 'chr_input', help = "File containing list of chromosomes to iterate over")
parser.add_argument('-w', '--win', dest = 'win_input', help = "File containing window coordinates and total callable sites. Format is Chr Start Stop")
parser.add_argument('-n', '--n_chr', dest = 'n_chr', help = "Integer specifying the number of sampled chromosomes", type=int)
args = parser.parse_args()

count_input = args.counts
chr_input = args.chr_input
win_input = args.win_input
n_chr = (args.n_chr)

chr_list = []
for line in open(chr_input):
    chr_list.append(line.strip())

with open(win_input) as winput:
    reader1 = csv.reader(winput, delimiter = '\t')
    win_data = list(zip(*reader1))

win_chr = list(win_data[0])
start = list(win_data[1])
stop = list(win_data[2])

win_coords = []
win_sites = []

for site in range(0, len(win_chr)):
    site_start = int(start[site])
    site_stop = int(stop[site])
    win_list = [site_start,site_stop]
    win_coords.append(win_list)

with open(count_input) as inf:
    reader = csv.reader(inf, delimiter = '\t')
    data = list(zip(*reader))

chr = list(data[0])

pos = []
freq_counts = []
for j in range(0,len(chr)):
    pos.append(int(data[1][j]))
    freq_counts.append(int(data[-1][j]))

def win_sfs(x):
    return(allel.sfs(x,n=n_chr))


idx = allel.ChromPosIndex(chr, pos)
win_idx = allel.ChromPosIndex(win_chr,start)

def a_n(n):
    an_sum = 0
    for x in range(1,n-1):
        an_sum += 1/x
    return(an_sum)

def b_n(n):
    bn_sum = 0
    for x in range(1,n-1):
        bn_sum += 1/(x*x)
    return(bn_sum)

for x in range(0, len(chr_list)):
    current_chr = chr_list[x]

    if current_chr in chr:
        chr_pos = pos[idx.locate_range(current_chr)]
        chr_freq_counts = freq_counts[idx.locate_range(current_chr)]
    else:
        chr_pos = []
        chr_freq_counts = []

    current_win_coords = win_coords[win_idx.locate_range(current_chr)]

    der_sfs, wins, counts = allel.windowed_statistic(chr_pos, chr_freq_counts, statistic=win_sfs, windows = current_win_coords)
    for win in range(0, len(der_sfs)):
        #print(der_sfs[win])
        theta_h_sum = 0
        theta_pi_sum = 0
        theta_l_sum = 0
        s_sum = 0
        if not isinstance(der_sfs[win], float):
            for i in range(1,len(der_sfs[win])-1):
                theta_h_site = (2*der_sfs[win][i]*i*i)/(n_chr*(n_chr-1))
                theta_h_sum+=theta_h_site

                theta_pi_site = (2*der_sfs[win][i]*i*(n_chr - i))/(n_chr*(n_chr-1))
                theta_pi_sum += theta_pi_site

                s_sum += der_sfs[win][i]

                theta_l_site = i * der_sfs[win][i]
                theta_l_sum += theta_l_site

            theta_l_sum = theta_l_sum * 1/(n_chr-1)

            theta_w = s_sum/a_n(n_chr)


            var_pi_l = ((n_chr -2)/(6*(n_chr-1)))*theta_w + ((18*n_chr*n_chr*(3*n_chr+2)*b_n(n_chr+1) - (88*n_chr*n_chr*n_chr+9*n_chr*n_chr-13*n_chr+6))/(9*n_chr*((n_chr-1)*(n_chr-1))))*((s_sum*(s_sum-1))/((a_n(n_chr)*a_n(n_chr))+b_n(n_chr)))
            #h_stat = theta_pi_sum - theta_h_sum
            h_stat = (theta_pi_sum - theta_l_sum)/(np.sqrt(var_pi_l))
            win_sites = counts[win]
            print(current_chr, wins[win][0], wins[win][1], win_sites, s_sum, theta_w, theta_pi_sum, theta_h_sum, theta_l_sum, var_pi_l, h_stat, sep = '\t')
            #print(der_sfs[win])
        #print(current_chr, wins[win][0], wins[win][1], win_sites,der_sfs[win])
