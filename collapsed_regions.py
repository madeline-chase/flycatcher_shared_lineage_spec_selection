#! /usr/bin/env python3

from sys import argv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = 'input', help = "Input file of heterozygote proportions")
parser.add_argument('--excess_snps', '-e', dest = 'excess', type = int, help = "Number of SNPs that must show excess heterozygosity")
parser.add_argument('--nonexcess_snps', '-n', dest = 'nonexcess', type = int, help = "Number of SNPs without excess heterozygosity that can be in between two SNPs with high heterozygosity.")
parser.add_argument('--cutoff', '-c', dest = 'cutoff', type = float, help = "A SNP must have greater than this proportion of heterozygotes to be considered an excess.")
parser.add_argument('--distance', '-d', dest = 'dist', type = int, help = "Maximum distance between two SNPs with excess heterozygosity")
parser.add_argument('--extend', '-x', dest = 'extend', type = int, help = "Extend the collapsed regions identified by this many bps. Can help to deal with the gradual rise in excess heterozygosity at the ends of collapsed regions.")
args = parser.parse_args()


excess = args.excess
y = args.nonexcess
step = y + 2
input = args.input
cutoff = args.cutoff
dist = args.dist
data = []
extend = args.extend

for position in open(input):
    data.append(position.strip().split('\t'))

collapsed = []
for position in range(0, len(data)):
    if len(collapsed) == 0:
        if float(data[position][2]) > cutoff:
            extra_excess = 1
            for z in range(position+1, position+step):
                if z < len(data):
                #print(data[z])
                    if float(data[z][2]) > cutoff and data[z][0] == data[position][0] and int(data[z][1]) - int(data[position][1]) < dist:
                        extra_excess += 1
            if extra_excess >= 2:
                collapsed.append(data[position][0:2])
        #print(collapsed)
    elif len(collapsed) > 0:
        if float(data[position][2]) > cutoff:
            extra_excess = 1
            for z in range(position+1, position+step):
                if z < len(data):
                    if float(data[z][2]) > cutoff and data[z][0] == data[position][0] and int(data[z][1]) - int(data[position][1]) < dist:
                        extra_excess += 1
            if extra_excess >= 2:
                collapsed.append(data[position][0:2])
            elif extra_excess < 2:
                collapsed.append(data[position][0:2])
                if len(collapsed) >= excess:
                    if int(collapsed[0][1])-extend-1 < 0:
                        start = 0
                    else:
                        start = int(collapsed[0][1])-extend-1
                    print(collapsed[0][0], start, int(data[position][1])+extend, sep = '\t')
                collapsed = []
