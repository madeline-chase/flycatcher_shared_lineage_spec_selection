# SNP polarization and sweep analysis

We looked for evidence of selective sweeps across the genome using Fay and Wu's H and SweepFinder2. We first polarized SNPs to obtain the derived site frequency spectrum, using *Ficedula hyperthra* as one outgroup.

## Identifying ancestral alleles

To identify the ancestral state, we used *F. hyperythra* as one outgroup, and then combined collared and pied into one group and red-breasted and taiga into another. We alternated each of these groups as the second outgroup to polarize sites for both as ingroups respectively. We obtained allele counts for all SNPs for the three groups (*F. hyperythra*, collared+pied, and red-breasted+taiga) with VCFtools with the following command, which removes sites that have greater than 10 % missing data in any of the four species.

```Bash
vcftools --vcf par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.all_chr.scaff_sorted.vcf --counts --keep par_taig_inds.txt --out par_taig_freq.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.par_taig_coll_pied_max_miss_10_perc  --exclude-positions missing_sites.CM.fem_het_Z.par_taig_coll_pied_max_miss_10_perc.txt
```

I then combined the allele counts for all three groups, and then identified sites for which two of the three groups were fixed for the same allele.

```Bash
paste hyp_freq.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc_par_taig_coll.split.txt coll_pied_freq.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.par_taig_coll_max_miss_10_perc.split.txt | paste - par_taig_freq.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.par_taig_coll_max_miss_10_perc.split.txt | awk '{if($4!=0 && $8==0 && $16==0) print $1,$2,$5; else if($4!=0 && $6==0 && $14==0) print $1,$2,$7; else if($4!=0 && $8==0 && $24==0) print $1,$2,$5; else if($4!=0 && $6==0 && $22==0) print $1,$2,$7; else if($16==0 && $24==0) print $1,$2,$13; else if($14==0 && $22==0) print $1,$2,$15}' > hyp.coll_pied.par_taig.inferred_anc_state_unfilt.bed
```

To account for misidentification of the ancestral state due to multiple mutations, I identified sites where the ancestral site formed a CpG or CpG-prone site, as these are known to be hypermutable.

I first reconstructed the putative ancestral genome by modifying the collared flycatcher reference. I masked positions that were not callable in each if the four species, as well as positions where SNPs were called as the ancestral state is equivocal at these positions unless two of the three groups were fixed for the same alleles. I then replaced The reference base with the ancestral allele for all positions where ancestral states could be identified.

Masking sites was performed with the following command using BEDTools maskfasta.

```Bash
bedtools maskfasta -fi $REFERENCE -bed coll.pied.par.taig.NONCALLABLE_LOCI.merged_with_SNPs.bed -fo anc_fasta_noncall_snps_masked.fasta
```

Replacing reference nucleotide with ancestral state was performed with the python script `get_anc_fasta.py`

```Bash
get_anc_fasta.py anc_fasta_noncall_snps_masked.fasta hyp.coll_pied.par_taig.inferred_anc_state_unfilt.bed > anc_fasta_noncall_snps_masked_anc_allele.fasta
```

Once the ancestral genome has been reconstructed, I retrieve the trinucleotide sequences for all the polarized SNPs from the ancestral reference using BEDTools getfasta, and a bed formatted list of the trinucleotide positions for each SNP.

I got the trinucleotide positions with the following command. This works because the original file is in bed format.

```Bash
 awk -v OFS='\t' '{print $1,$2-1,$3+1}' hyp.coll_pied.par_taig.inferred_anc_state_unfilt.bed > hyp.coll_pied.par_taig.inferred_anc_state_unfilt.trinuc_pos.bed
```

and then obtained the trinucleotide sequences using BEDTools getfasta.

```Bash
bedtools getfasta -fi anc_fasta_noncall_snps_masked_anc_allele.fasta -bed hyp.coll_pied.par_taig.inferred_anc_state_unfilt.trinuc_pos.bed -fo hyp.coll_pied.par_taig.inferred_anc_state_unfilt.trinuc_seq.fasta
```

I then used the python script `CpG_filt_remove_sites.py` to remove polarized sites that form CpG and CpG-prone sites.

```Bash
python3 CpG_filt_remove_sites hyp.coll_pied.par_taig_join.inferred_anc_state_unfilt.trinuc_seq.fasta > hyp.coll_pied.par_taig.RM.CM.minDP.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.anc_state_cpg_filt.bed
```

This leaves me with my final set of polarized sites.

## Obtaining derived allele counts

I then obtained derived allele counts for each species with VCFtools. I used the following command to get allele counts for only the sites where it was possible to identify the ancestral state.

```Bash
vcftools --vcf par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.all_chr.scaff_sorted.vcf --counts --keep gotland.1993.2015.samples.txt --out coll_freq.minDP5.maxDP200.minGQ.polarized_sites --positions hyp.coll_pied.par_taig.RM.CM.minDP.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.anc_state_cpg_filt.txt
```

I then reordered the sites so the ancestral allele was always listed first with the following command.

```Bash
paste hyp.coll_pied.par_taig.RM.CM.minDP.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.anc_state_cpg_filt.txt coll_freq.minDP5.maxDP200.minGQ.polarized_sites_with_pied.split.txt | awk -v OFS='\t' '{if($3==$8) print $4,$5,$6,$7,$8,$9,$10,$11; else if($3==$10) print $4,$5,$6,$7,$10,$11,$8,$9}' | awk '!($8==0) {print}' > coll.freq.minDP5.maxDP200.minGQ.polar_sites_sorted.anc_fixed_rmvd.txt
```

These data were then used to calculate Fay and Wu's H and CLR with SweepFinder2.

## Calculating Fay and Wu's H

Fay and Wu's H was calculated in 200kb windows using a custom script that made use of the scikit-allel package. The formula for Fay and Wu's H was taken from equations 11 and 12 of [Zeng et al. 2006](https://www.genetics.org/content/174/3/1431), which improves upon the original estimation from Fay and Wu 2000 by normalizing the statistic.

The script was run with the following command, which provides the derived allele count with fixed sites removed for a single species, the scaffolds for which there are data, a file with window coordinates specified as well as the total number of callable sites within a window (see `callable_loci.md` for how called sites were estimated), and the number of total chromosomes present within the species.

```Bash
python3 calc_fay_wu_normalized.py -i coll.freq.minDP5.maxDP200.minGQ.polar_sites_sorted.fixed_rmvd.txt  -c coll.win_chrs_uniq.txt -w coll.callable_sites.w200k_s200k.depth_only.q10.minDP5.maxDP200.50_perc_called.txt -n 190 > coll.polar_sites.fay_wu_h_normalized.w200k_s200k.min_50_perc_win_called.fixed_rmvd.txt
```

The significance threshold for Fay and Wu's H was determined to be -1.60 based on simulations of background selection (see simulations.md). After aligning scaffolds to the collared flycatcher chromosomes, I identified windows with significantly negative Fay and Wu's H, and used BEDTools to intersect these with the Fst peak class data.

```Bash
awk '{if($11<-1.6) print}' coll.polar_sites.fay_wu_h_normalized.w200k_s200k.min_50_perc_win_called.fixed_rmvd.chrpos.bed | bedtools intersect -a pop_gen_wins_200k_s200k_fst_peak_class.bed -b stdin -c  
```

## Calculating CLR

SweepFinder2 was used to estimate CLR for each of the four species using the derived allele counts. I first removed sites that were fixed for the ancestral allele, as this can increase the frequency of false positives due to low mutation rates. Then I used the following command to run SweepFinder2 for each scaffold, using the autosomal SFS as the user defined background spectrum, and a grid file of all variants for the focal species.

```Bash
chrs=`cat coll_chrs_chrZ_rmvd.txt`

for chr in $chrs
do
  awk -v var="$chr" -v OFS='\t' '{if($1==var) {print $2,$8,$4,0}}' coll.freq.minDP5.maxDP200.minGQ.polar_sites_sorted.anc_fixed_rmvd.chrZ_rmvd.txt > coll.$chr.unfolded.SF_input.no_header.txt
  cat SF_header.txt coll.$chr.unfolded.SF_input.no_header.txt > coll.$chr.unfolded.SF_input.header.txt
  grep -v 'folded' coll.$chr.unfolded.SF_input.header.txt | cut -f1 > coll.$chr.unfolded.grid_file.txt
done

awk -v OFS='\t' '{print $2,$8,$4,0}' coll.freq.minDP5.maxDP200.minGQ.polar_sites_sorted.anc_fixed_rmvd.chrZ_rmvd.txt > coll.whole_genome.SF_input.no_header.txt
cat SF_header.txt coll.whole_genome.SF_input.no_header.txt > coll.whole_genome.SF_input.header.txt
SweepFinder2 -f coll.whole_genome.SF_input.header.txt coll.whole_genome.SF_spec

for chr in $chrs
do
    SweepFinder2 -lu coll.$chr.unfolded.grid_file.txt coll.$chr.unfolded.SF_input.header.txt coll.whole_genome.SF_spec coll.$chr.ug.SF_out.txt
done
```

The average CLR in 200kb windows was calculated using the python script `calc_mean_clr_fst_wins.py`

```Bash
python3 calc_mean_clr_fst_wins.py -i coll.ug.SF_out_combined_scaffs_no_chrZ.txt -c coll_chrs_chrZ_rmvd.txt -w pop_gen_wins_200k_s200k_fst_peak_class.txt > coll.ug.SF_out_combined_scaffs_no_chrZ.w200k_s200k_mean_clr.txt
```

I also identified the significance threshold for CLR based on the same simulations of background selection, which was 46.25. I additionally used a permutation test for each species to determine the expected number of significant sites overlapping with 200kb genomic windows by chance (see `simulations.md`). This was used to look for genomic windows showing a significant signature of a selective sweep.

I used BEDTools to get the counts of significant CLR sites in genomic windows with Fst peak class data, and then counted the number of windows more than threshold number of overlaps for a significant sweep signature for each Fst peak class.

```Bash
awk '{if($4>46.25) print}' coll.ug.SF_out_combined_scaffs_no_chrZ.chrompos_sorted.bed | bedtools intersect -a pop_gen_wins_200k_s200k_fst_peak_class.bed -b stdin -c | awk '{if($5>0) print}' > coll_fst_wins_w50k_s50k_sig_clr_olap_permuted_05_thres.bed
```

Additionally, I identified collared flycatcher genes that showed a significant sweep signature in each species and pulled out the ensembl IDs.

```Bash
awk -v OFS='\t' '{if($4>46.25) print}' coll.ug.SF_out_combined_scaffs_no_chrZ.chrompos_sorted.bed | bedtools intersect -a full.genes.on.fAlb15.chromPos_strict_sorted.from_gtf.bed -b stdin -c | awk '{if($7>0) print}' > ../results/sig_gene_out/coll_genes_sig_clr_overlaps_permuted_05_thres.bed
```

I then identified genes that showed a sweep signature in both of the independent species pairs (i.e. shared between either of collared and pied with either of red-breasted and taiga).

```Bash
cut -f4 coll_genes_sig_clr_overlaps_permuted_05_thres.bed | grep -Fwf - par_genes_sig_clr_overlaps_permuted_05_thres.bed > coll_par_sig_genes_overlap_05_thres.bed
cat coll_par_sig_genes_overlap_05_thres.bed coll_taig_sig_genes_overlap_05_thres.bed pied_par_sig_genes_overlap_05_thres.bed pied_taig_sig_genes_overlap_05_thres.bed | cut -f4 | sort | uniq > sig_genes_shared_across_sister_comps_05_thres.txt
```

Finally, I identified the genes that showed sweep signatures in only a single species.

```Bash
cat par_genes_sig_clr_overlaps_permuted_05_thres.bed pied_genes_sig_clr_overlaps_permuted_05_thres.bed taig_genes_sig_clr_overlaps_permuted_05_thres.bed | cut -f4 | sort | grep -Fwvf - coll_genes_sig_clr_overlaps_permuted_05_thres.bed > coll_unique_genes_05_thres.txt
```

## Correlation between Fay and Wu's H and CLR

Using the mean CLR estimated in 200kb windows, I calculated the correlation between Fay and Wu's H and CLR for each species in R.

```R

#### Read in CLR mean data in Fst windows ####

coll_mean_clr <- read.table('coll.ug.SF_out_combined_scaffs_no_chrZ.w200k_s200k_mean_clr.txt', header = T)
pied_mean_clr <- read.table('pied.ug.SF_out_combined_scaffs_no_chrZ.w200k_s200k_mean_clr.txt', header = T)
par_mean_clr <- read.table('par.ug.SF_out_combined_scaffs_no_chrZ.w200k_s200k_mean_clr.txt', header = T)
taig_mean_clr <- read.table('taig.ug.SF_out_combined_scaffs_no_chrZ.w200k_s200k_mean_clr.txt', header = T)

coll_mean_clr$Mid <- (coll_mean_clr$Start+coll_mean_clr$Stop)/2
pied_mean_clr$Mid <- (pied_mean_clr$Start+pied_mean_clr$Stop)/2
par_mean_clr$Mid <- (par_mean_clr$Start+par_mean_clr$Stop)/2
taig_mean_clr$Mid <- (taig_mean_clr$Start+taig_mean_clr$Stop)/2

colnames(coll_mean_clr) <- c("Chr", "Start","Stop", "Coll_clr_sites", "Coll_mean_CLR", "Mid")
colnames(pied_mean_clr) <- c("Chr", "Start","Stop", "pied_clr_sites", "pied_mean_CLR", "Mid")
colnames(par_mean_clr) <- c("Chr", "Start","Stop", "par_clr_sites", "par_mean_CLR", "Mid")
colnames(taig_mean_clr) <- c("Chr", "Start","Stop", "taig_clr_sites", "taig_mean_CLR", "Mid")

clr_joined <- join(taig_mean_clr, join(par_mean_clr, join(coll_mean_clr, pied_mean_clr, type = 'inner', by = c("Chr", "Mid")), type = 'inner', by = c("Chr", "Mid")), type = 'inner', by = c("Chr","Mid"))

## Combine CLR wins with other population genomic stats data

all_stats_joined_clr_wins <- join(all_stats_joined, clr_joined, by = c('Chr','Stop'), type = 'inner')


#### Read in Fay and Wu's H data for each species ####

# Parva
par_h_wins <- read.table('par.polar_sites.fay_wu_h_normalized.w200k_s200k.min_50_perc_win_called.fixed_rmvd.chrpos.txt')
colnames(par_h_wins) <- c('Chr','Start','Stop', 'Sites_par_h','S_sites_par', 'Theta_w_par','Theta_pi_par','Theta_h_par','Theta_l_par','Var_par', 'H_par')
par_h_wins$Mid <- (par_h_wins$Start+par_h_wins$Stop)/2
# Taiga
taig_h_wins <- read.table('taig.polar_sites.fay_wu_h_normalized.w200k_s200k.min_50_perc_win_called.fixed_rmvd.chrpos.txt')
colnames(taig_h_wins) <- c('Chr','Start','Stop', 'Sites_taig_h','S_sites_taig', 'Theta_w_taig','Theta_pi_taig','Theta_h_taig','Theta_l_taig','Var_taig', 'H_taig')
taig_h_wins$Mid <- (taig_h_wins$Start+taig_h_wins$Stop)/2

# Collared
coll_h_wins <- read.table('coll.polar_sites.fay_wu_h_normalized.w200k_s200k.min_50_perc_win_called.chrpos.fixed_rmvd.txt')
colnames(coll_h_wins) <- c('Chr','Start','Stop', 'Sites_coll_h','S_sites_coll', 'Theta_w_coll','Theta_pi_coll','Theta_h_coll','Theta_l_coll','Var_coll', 'H_coll')
coll_h_wins$Mid <- (coll_h_wins$Start+coll_h_wins$Stop)/2

# Pied
pied_h_wins <- read.table('pied.polar_sites.fay_wu_h_normalized.w200k_s200k.min_50_perc_win_called.fixed_rmvd.chrpos.txt')
colnames(pied_h_wins) <- c('Chr','Start','Stop', 'Sites_pied_h','S_sites_pied', 'Theta_w_pied','Theta_pi_pied','Theta_h_pied','Theta_l_pied','Var_pied', 'H_pied')
pied_h_wins$Mid <- (pied_h_wins$Start+pied_h_wins$Stop)/2


## For each species I want to remove the Z-chromosome and LGE22

par_h_wins <- filter(par_h_wins, Chr!='ChrLGE22' & Chr!='ChrZ')
taig_h_wins <- filter(taig_h_wins, Chr!='ChrLGE22' & Chr!='ChrZ')
coll_h_wins <- filter(coll_h_wins, Chr!='ChrLGE22' & Chr!='ChrZ')
pied_h_wins <- filter(pied_h_wins, Chr!='ChrLGE22' & Chr!='ChrZ')

h_stats_combined <- join(join(join(par_h_wins, taig_h_wins, by = c('Chr','Mid'), type = 'inner'), coll_h_wins, by = c('Chr', 'Mid'), type = 'inner'), pied_h_wins, by = c("Chr","Mid"), type = 'inner')
h_stats_combined <- h_stats_combined[,-c(2,3,13,14,23,24,33,34)]

## And combine h stat data with the other statistics
all_stats_joined_clr_h_stats <- join(all_stats_joined_clr_wins, h_stats_combined, by = c("Chr","Mid"), type = "inner")


#### Relationships between H and CLR ####

cor.test(all_stats_joined_clr_h_stats$H_coll, all_stats_joined_clr_h_stats$Coll_mean_CLR)
cor.test(all_stats_joined_clr_h_stats$H_pied, all_stats_joined_clr_h_stats$pied_mean_CLR)
cor.test(all_stats_joined_clr_h_stats$H_par, all_stats_joined_clr_h_stats$par_mean_CLR)
cor.test(all_stats_joined_clr_h_stats$H_taig, all_stats_joined_clr_h_stats$taig_mean_CLR)
```

I then used BEDTools to count how many significant Fay and Wu's H windows also showed a significant sweep signature with the CLR method.

```Bash
awk '{if($11<-1.60) print}' coll.polar_sites.fay_wu_h_normalized.w200k_s200k.min_50_perc_win_called.chrpos.fixed_rmvd.bed | bedtools intersect -a pop_gen_wins_200k_s200k_fst_peak_class.bed -b stdin | bedtools intersect -a stdin -b coll_fst_wins_w200k_s200k_sig_clr_overlaps_permuted_05_thres.bed
```
