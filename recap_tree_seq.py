import msprime, pyslim
import numpy as np
from sys import argv

# Get the slim simulation output and simulation number from the command line
infile = argv[1]
outnum = str(argv[2])

# Load the tree sequence from slim
ts = pyslim.load(infile)

# Define function to remove a fraction of the mutations added in tskit
# This is used because I have variable genomic regions with different forms of mutations
# Coding regions and CNEs have deleterious mutations that were added in SLiM but I need to
# add neutral mutations into noncoding/non-CNE regions. This function can be used after
# adding neutral mutations to the tree sequence to remove a proportion of mutations from
# specified regions.

def remove_mutations(ts, start, end, proportion):
    '''
    This function will return a new tree sequence the same as the input,
    but after removing each non-SLiM mutation within regions specified in lists
    start and end with probability `proportion`, independently. So then, if we
    want to add neutral mutations with rate 1.0e-8 within the regions and 0.7e-8
    outside the regions, we could do
      ts = pyslim.load("my.trees")
      first_mut_ts = msprime.mutate(ts, rate=1e-8)
      mut_ts = remove_mutations(first_mut_ts, start, end, 0.3)

    :param float proportion: The proportion of mutations to remove.
    '''
    new_tables = ts.dump_tables()
    new_tables.mutations.clear()
    mutation_map = [-1 for _ in range(ts.num_mutations)]
    for j, mut in enumerate(ts.mutations()):
        keep_mutation = True
        for i in range(len(start)):
            left = start[i]
            right = end[i]
            assert(left < right)
            if i > 0:
                assert(end[i - 1] <= left)
            if mut.position >= left and mut.position < right and len(mut.metadata) == 0:
                keep_mutation = (random.uniform(0, 1) > proportion)
        if keep_mutation:
            mutation_map[j] = new_tables.mutations.num_rows
            if mut.parent < 0:
                new_parent = -1
            else:
                new_parent = mutation_map[mut.parent]
            new_tables.mutations.add_row(site = mut.site, node = mut.node,
                    derived_state = mut.derived_state,
                    parent = new_parent,
                    metadata = mut.metadata)
    return new_tables.tree_sequence()

# Read in recombination map used for simulations

positions = []
rates = []

with open("Rec.200kb.5kGap.chr11_scaffs.chrompos.bed") as file:
    for line in file:
        components = line.split('\t')
        positions.append(float(components[2]))
        # Need to convert from cM/Mb, in the simulations I multiplied the actual
        # recombination rate by 10 to scale for the downscaled population size
        # so I multiply rates by 1e-7 (rather than 1e-8)
        rates.append(1e-7*float(components[3]))

# These additions to the positions and rates lists are necessary for the msprime format
positions.insert(0,0)
rates.append(0.0)

# Make recombination map
recomb_map = msprime.RecombinationMap(positions,rates)

# Read in gene + CNE locations
gene_starts = []
gene_stops = []
with open("fAlb15.introns_exons_intergenic.CNEs_included.ensemble_e73.chrom_pos.chrom11.sorted.bed") as gene_file:
    for line in gene_file:
        components = line.split('\t')
        if components[-1] == 'CNE' or components[-1] == 'exonic':
            gene_starts.append(int(components[1]))
            gene_stops.append(int(components[2]))

# Recapitate tree sequences using the recombinaton map specified
recap_ts = ts.recapitate(recombination_map = recomb_map, Ne = 10000)

# Add neutral mutations to the tree with mu = 4.6e-8. This is 10x the estimated
# mu for collared flycatcher, again to scale for the reduction in pop size in
# simulations
ts = pyslim.SlimTreeSequence(msprime.mutate(recap_ts, rate = 4.6e-8, keep = True))

# Now, given the list of coordinates for coding regions and CNEs, I remove 90% of
# the neutral mutations that were added. Within the simulations I scaled the mutation
# rate of deleterious mutations such that 10% of mutations would be neutral
for region in range(0,len(gene_starts)):
    ts = remove_mutations(ts, gene_starts[region], gene_stops[region], 0.9)


# Get a list of individuals alive at the end of the simulation
alive = ts.individuals_alive_at(0)

# Choose a random sample of 15 individuals from alive and get nodes for these individuals
sample_inds = [np.random.choice(alive, 15, replace = False)]
keep_inds = np.concatenate(sample_inds)
keep_nodes = []

for i in keep_inds:
    keep_nodes.extend(ts.individual(i).nodes)

# Simplify the treesequence with the keep_nodes
ts_simp = ts.simplify(keep_nodes)

sfs = ts_simp.allele_frequency_spectrum(polarised=True)

# Write an output file in the format required for SweepFinder2
outfile = open('slim_sim_out_for_SF2_{}.txt'.format(outnum), "w")
header = "position\tx\tn\tfolded\n"
outfile.write(header)
for variant in ts_simp.variants():
    # Some sites will have multiple mutations (i.e. alleles != 0 or 1), I want to
    # remove these sites
    if ((variant.genotypes <= 1).sum() == variant.genotypes.size).astype(np.int) > 0:
        var_string = str(variant.site.position) + '\t' + str(variant.genotypes.sum()) + '\t' + str(variant.genotypes.size) + '\t0\n'
        outfile.write(var_string)
outfile.close()
