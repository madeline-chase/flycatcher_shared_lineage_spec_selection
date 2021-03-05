# Simulations and permutation tests for sweep signatures

To determine significance thresholds for the sweep measures, Fay and Wu's H and CLR, I ran simulations in SLiM. I simulated the action of background selection on an approximately 21Mbp chromosome, based on Chromosome 11 of the collared flycatcher reference assembly.

## Simulations in SLiM

I ran 1000 simulations in SLiM with the script `bgs_n10k_sfixed.slim` with the following command.

```Bash
for i in {1..1000}
do
  slim -seed $RANDOM -d out=$i bgs_n10k_sfixed.slim
done
```

## Recapitation and adding neutral mutations

After simulations completed, I used the python script `recap_tree_seq.py` to perform recapitation and add neutral mutations using [msprime](https://msprime.readthedocs.io/en/stable/index.html) and [tskit](https://tskit.readthedocs.io/en/latest/index.html). I then sampled 15 individuals from the final population and output allele count data in the format needed for SweepFinder2.

```Bash
for i in {1..1000}
do
  python recap_tree_seq.py bgs_test_N10k_gen100k_neut_mut_rmvd_fixed_sel_sim$i.trees sim$i
done
```



## Running SweepFinder2 on simulated data

I ran SweepFinder2 on the simulated data following the same procedure as the real data, where I first obtained the background SFS from the full dataset, and then defined a grid file with every position.

```Bash
for i in {1..1000}
do
  SweepFinder2 -f slim_sim_out_for_SF2_sim$i.txt slim_sim_out_for_SF2_sim$i.spec_file
  grep -v 'folded' slim_sim_out_for_SF2_sim$i.txt | cut -f1 > slim_sim_out_for_SF2_sim$i.grid_file.txt

  SweepFinder2 -lu slim_sim_out_for_SF2_sim$i.grid_file.txt slim_sim_out_for_SF2_sim$i.txt slim_sim_out_for_SF2_sim$i.spec_file slim_sim_out_for_SF2_sim$i.sf_out_file.txt
done
```

I then combined the output for all 1000 simulations, and used the python script `get_clr_sig_threshold.py` to get significance thresholds. We used a p-value < 0.001 in the paper.

```Bash
cat *.sf_out_file.txt > slim_sims.sf_out_file.combined.txt
python get_clr_sig_threshold.py slim_sims.sf_out_file.combined.txt > slim_sims.sf_out.clr_sig_thresholds.txt
```

## Running Fay and Wu's H on simulated data

I also ran Fay and Wu's H on all 1000 simulated datasets. This script is written to take the same input format as SweepFinder2. The `-n 30` option specifies that there are 30 chromosomes present in the population.

```Bash
python3 calc_fay_wu_normalized_from_sims.py -i slim_sim_out_for_SF2_sim$i.txt -n 30 > fay_wu_h_sim_out$i.txt
```

I then combined the output and used the script `get_H_sig_threshold.py` to get significance thresholds. We also used a p-value < 0.001 for Fay and Wu's H.

```Bash
cat fay_wu_h_sim_out* > fay_wu_h_sims_combined.txt
python3 get_H_sig_threshold.py fay_wu_h_sims_combined.txt > fay_wu_h_sig_thresholds.txt
```

## Permutation tests to define significance thresholds in genomic windows

To identify genomic windows showing a significant sweep signature, I performed a permutation test with the observed CLR values for each species to generate a null distribution of how many significant CLR values we would expect to find in a 50kb or 200kb window by chance. I used a p-value < 0.05 to determine the threshold number of sites.

```Bash
python3 permute_clr_pos.py -i coll.ug.SF_out_combined_scaffs_no_chrZ.chrompos_sorted.bed -b pop_gen_stats_peak_class_only.bed -o coll_fst_w200k
```

## Permutation tests to define significance thresholds in genes

I also used a permutation test with observed CLR values to determine the threshold number of sites that could be expected by chance in a gene to find genes that showed significant sweep signatures. Using the coordinates of genes from the collared flycatcher annotation, I performed 1000 permutations for each species and used a p-value < 0.05.

```Bash
python3 permute_clr_pos.py -i coll.ug.SF_out_combined_scaffs_no_chrZ.chrompos_sorted.bed -b full.genes.on.fAlb15.chromPos_strict_sorted.from_gtf.bed -o coll_genes
```
