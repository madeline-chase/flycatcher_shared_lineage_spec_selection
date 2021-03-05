#! /usr/bin/env python3

from sys import argv

input = argv[1]

seqs = []


for line in open(input):
     seqs.append(line.strip())

#print(len(seqs))

print("Scaff","Start","Stop","Sites", sep = '\t')

for window in range(0,len(seqs)-1,2):
    coordinate = seqs[window]
    fasta = seqs[window+1]

    scaff = coordinate.split(":")[0][1:]
    #print(scaff)
    start = coordinate.split(":")[1].split("-")[0]
    stop = coordinate.split(":")[1].split("-")[1]

    fasta_N_rmvd = fasta.replace("N","")
    num_sites = len(fasta_N_rmvd)

    print(scaff, start, stop, num_sites, sep = '\t')


    #print(window)
