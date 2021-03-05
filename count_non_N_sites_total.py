from sys import argv

fasta = argv[1]

fasta_list = []
for line in open(fasta):
    fasta_list.append(line.strip())

total_length = 0
for seq in range(0,len(fasta_list)):
    if not fasta_list[seq][0].startswith(">"):
        orig_seq = fasta_list[seq]
        seq_n_replaced = orig_seq.replace("N", "")
        total_length+=len(seq_n_replaced)

print(total_length)

