# Genome-wide and window-based diversity and divergence statistics

I estimated Dxy, pi, and Fst across the whole genome, and in 50kb and 200kb genomic windows, to visualize at multiple scales.


## Calculating pi and Dxy

I followed the method to calculate pi and Dxy from [Irwin et al 2016](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.13792), using allele frequencies.

Allele counts were obtained using vcftools with the following command. This gets allele counts for a single species, and removes sites that were identified as missing within that species.

```Bash
vcftools --vcf $VCFDIR/par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.all_chr.scaff_sorted.vcf --counts --keep $SAMPDIR/gotland.1993.2015.samples.txt --out coll_freq.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.par_taig_coll_max_miss_10_perc --maxDP 200 --remove-filtered-geno-all --exclude-positions missing_sites.CM.fem_het_Z.par_taig_coll_max_miss_10_perc.txt
```

To calculate pi I first used the following command to estimate pi on a per SNP basis.

```Bash
awk '{p=$6/$4; q=$8/$4; print $1,$2,(2*p*q)}' coll_freq.minDP5.maxDP200.RM.CM.minGQ.max_miss_10_perc.fem_het_rmvd.hard_filt.SNPs_only.frq.txt > coll_freq.minDP5.maxDP200.RM.CM.minGQ.max_miss_10_perc.fem_het_rmvd.hard_filt.SNPs_only.snp_pi.txt
```

I then estimated average pi across the whole genome by taking the sum of pi for each SNP, and dividing by the total number of callable sites for each species with the following command. See `callable_loci.md` for how the number of callable sites was estimated.

```Bash
awk '{sum+=$3}' END '{print sum/844581092}' coll_freq.minDP5.maxDP200.RM.CM.minGQ.max_miss_10_perc.fem_het_rmvd.hard_filt.SNPs_only.snp_pi.txt
```

To get window based estimates of pi, I used the python script `calc_pi_dxy.py` that makes use of the [scikit-allel](https://scikit-allel.readthedocs.io/en/stable/) python library to estimate pi or dxy in windows. I provide a list of window coordinates in which to estimate pi, and the total number of callable sites estimated in each window.

```Bash
python3 calc_pi_dxy.py -p coll_freq.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.snp_pi.txt -c coll.uniq_scaffs.txt -w callable_sites.coll.w200kb_s200kb.min_50_perc_window_called.txt -s pi > coll_pi.snps_only.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.txt
```

A similar method was used to calculate Dxy. To begin, I combined the count data for two species and estimated Dxy on a per SNP basis with the following command (note: this requires both species have count data for the exact same set of SNPs, so I remove missing sites further downstream).

```Bash
paste coll_freq.minDP5.maxDP200.RM.CM.minGQ.max_miss_10_perc.fem_het_rmvd.hard_filt.split.txt taig_freq.minDP5.maxDP200.RM.CM.minGQ.max_miss_10_perc.fem_het_rmvd.hard_filt.split.txt | awk '{p1=$6/$4; q1=$8/$4; p2=$14/$12; q2=$16/$12; print $1,$2, (p1*q2)+(p2*q1)}' | awk '!($3==0) {print}' | grep -Fwvf missing_sites.CM.fem_het_Z.coll_taig_max_miss_10_perc.txt - > coll_taig.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.snp_dxy.txt
```

I then also took the sum of Dxy across all SNPs and then divided by the total number of callable sites in those two species.

```Bash
awk '{sum+=$3}' END print '{sum/826836105}'
```

Finally, I used the same script, `calc_pi_dxy.py`, to estimate Dxy in windows by giving the snp dxy input, the scaffolds for which there are snp data, and the coordinates for the windows to calculate dxy.

```Bash
python3 calc_pi_dxy.py -p coll_taig.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.snp_dxy.txt -c coll.taig.uniq_scaffs.txt -w  callable_sites.coll.taig.w200kb_s200kb.min_50_perc_window_called.txt -s dxy > coll_taig_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.txt
```

## Calculating pi and Dxy with exons and CNEs masked

We also estimated pi and Dxy in windows after masking sites that may be direct targets of selection (exons and conserved noncoding elements) to estimate the extent of linked selection. We used the same SNP based estimates of pi and Dxy, but removed sites that overlapped with exons or CNEs based on collared flycatcher annotations.

```Bash
grep -Fwvf par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.all_chr.scaff_sorted.exon.CNE.intersect.txt coll_freq.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.snp_pi.txt > coll_freq.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.snp_pi.auto.exons_cnes_masked.txt
```

I then used the same python script to estimate pi and Dxy in windows, this time with the total sites also having exons and CNEs excluded.

```Bash
python3 calc_pi_dxy.py -p coll_freq.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.snp_pi.exon_cnes_masked.txt -c coll.unique_scaffs.exons_cnes_masked.txt -w  callable_sites.coll.exons_cnes_masked.w200k_s200k.txt -s pi > coll_pi.snps_only.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.exon_cnes_masked.all_win.txt
```

## Estimating the proportion of shared SNPs

To estimate the proportion of shared SNPs in each species pair, I randomly downsampled each species to the smallest sample size (11) and retrieved the allele counts in VCFtools.

```Bash
vcftools --vcf par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.all_chr.scaff_sorted.vcf --counts --keep gotland.1993.2015.samples.txt --out coll_freq.RM.CM.minDP5.maxDP200.minGQ.down_samp_11.coll_max_miss_10_perc --max-indv 11 --exclude-positions  missing_sites.CM.fem_het_Z.coll_max_miss_10_perc.txt
```

## Calculating Fst

Fst was calculated using VCFtools.

The following command was used to estimate Weir and Cockerham's Fst in 200kb windows, and removes sites with more than 10% missing data in either of the two species being compared.

```Bash
vcftools --vcf par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.all_chr.scaff_sorted.vcf --weir-fst-pop gotland.1993.2015.samples.txt --weir-fst-pop pied_sweden.post_filt.txt --fst-window-size 200000 --fst-window-step 200000  --exclude-positions missing_sites.CM.fem_het_Z.coll_pied_max_miss_10_perc.txt --out coll_pied.minGQ.minDP5.maxDP200.biallelic.maxMiss10perc.femHetRmvd.RM.CM.fst_w200k_s200k_auto_scaffs_only
```

```Bash
vcftools --vcf par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.all_chr.scaff_sorted.vcf --weir-fst-pop par_inds.txt --weir-fst-pop alb_inds.postFilt.txt --fst-window-size 200000 --fst-window-step 200000  --exclude-positions missing_sites.CM.fem_het_Z.par_taig_max_miss_10_perc.txt --out par_taig.minGQ.minDP5.maxDP200.biallelic.maxMiss10perc.femHetRmvd.RM.CM.fst_w200k_s200k_auto_scaffs_only
```


## Identifying Fst peaks

I identified Fst peaks in 50kb and 200kb windows with a permutation test in R, after scaffolds were translated to chromosome positions based on the collared flycatcher assembly and linkage map.

```R
## Load some packages

library(signal)
library(plyr)
library(dplyr)


## Read the data and perform some filtering

par_taig_fst_auto <- read.table('par_taig.minGQ.minDP5.maxDP200.biallelic.maxMiss10perc.femHetRmvd.RM.CM.fst_w200k_s200k_auto_scaffs_only.windowed.weir.fst.chrompos')
colnames(par_taig_fst_auto) <- c('Chr', 'PT_Start', 'PT_Stop', 'PT_fst_Sites', 'PT_weight_fst', 'PT_fst')
par_taig_fst_auto$Mid <- (par_taig_fst_auto$PT_Start + par_taig_fst_auto$PT_Stop)/2
par_taig_fst_auto <- par_taig_fst_auto[order(par_taig_fst_auto$Chr, par_taig_fst_auto$Mid),]
par_taig_fst_auto <- filter(par_taig_fst_auto, Chr!='ChrLGE22') ## Removing Chroms with few windows
par_taig_fst_auto <- filter(par_taig_fst_auto, Chr!='Chr25')
par_taig_fst_auto <- filter(par_taig_fst_auto, PT_fst_Sites >= 200) ## Removing windows with few sites

pied_coll_fst_auto <- read.table('coll_pied.minGQ.minDP5.maxDP200.biallelic.maxMiss10perc.femHetRmvd.RM.CM.fst_w200k_s200k_auto_scaffs_only.windowed.weir.fst.chrompos')
colnames(pied_coll_fst_auto) <- c('Chr', 'PC_Start', 'PC_Stop', 'PC_fst_Sites', 'PC_weight_fst', 'PC_fst')
pied_coll_fst_auto$Mid <- (pied_coll_fst_auto$PC_Start+pied_coll_fst_auto$PC_Stop)/2
pied_coll_fst_auto <- pied_coll_fst_auto[order(pied_coll_fst_auto$Chr, pied_coll_fst_auto$Mid),]
pied_coll_fst_auto <- filter(pied_coll_fst_auto, Chr != 'ChrLGE22')
pied_coll_fst_auto <- filter(pied_coll_fst_auto, Chr != 'Chr25')
pied_coll_fst_auto <- filter(pied_coll_fst_auto, PC_fst_Sites >= 200)

##### Writing functions necessary to locate peaks

## Calc z_fst function
z_fst <- function(x){
  (x-mean(x))/sd(x)
}

## Function to take zfst and calculate across each chromosome
z_fst_chr <- function(y){
  chr_list <- unique(y[,1])
  chr <- y[,1]
  z_fst_vec <- c()
  for(i in 1:length(chr_list)){
    fst <- y[,5]
    current_chr <- chr_list[i]
    z_fst_chr_vals <- z_fst(fst[chr==current_chr])
    z_fst_vec <- append(z_fst_vec, z_fst_chr_vals)
  }
  return(z_fst_vec)
}

## Write a function to smooth windows by chromosome
## Takes two variables x = values of fst to smooth
## y = dataframe with chromosome values as first column
smooth_windows <- function(x,y){
  chr <- y[,1]
  zfst <- x
  chr_list <- unique(chr)

  zfst_smoothed <- c()
  for(m in 1:length(chr_list)){
    current_chr <- chr_list[m]
    zfst_smooth_chr <- sgolayfilt(zfst[chr==current_chr], p = smooth_p, n = smooth_n)
    zfst_smoothed <- append(zfst_smoothed, zfst_smooth_chr)
  }
  return(zfst_smoothed)
}

## Write a function to perform permutation test, applying smoothing function written above.
## x = vector of raw (unsmoothed) Fst values
## y = dataframe with fst wins output (chromosome is first column)
permute_zfst <- function(x,y){sapply(1:n_permutations, function(z){
  randomized_zfst <- sample(x, size = length(x), replace = F)
  smooth_windows(randomized_zfst, y)
})}

## Function to calculate p-values for observed Fst. x = matrix of permuted values across wins, y = vector of observed (smoothed) fst.
pval_zfst <- function(x,y){
  n_wins <- length(x[,1])
  pvals <- c()
  for(i in 1:n_wins){
    obs_val <- y[i]
    perm_vals <- x[i,]
    pval_win <- length(perm_vals[perm_vals>=obs_val])/n_permutations
    pvals <- append(pvals, pval_win)
  }
  return(pvals)
}

## Variables that need to be set to run the diff functions:

## Define smoothing function
# p = power
# n = number wins to smooth over
smooth_p <- 3
smooth_n <- 7

## Define number of permutations to perform
n_permutations <- 10000

## Get z-fst values for both comps

par_taig_zfst_auto <- z_fst_chr(par_taig_fst_auto)
pied_coll_zfst_auto <- z_fst_chr(pied_coll_fst_auto)


## Smooth observed values

par_taig_zfst_smooth <- smooth_windows(par_taig_zfst_auto, par_taig_fst_auto)
pied_coll_zfst_smooth <- smooth_windows(pied_coll_zfst_auto, pied_coll_fst_auto)

## Permute windows

par_taig_zfst_permute <- permute_zfst(par_taig_zfst_auto, par_taig_fst_auto)
pied_coll_zfst_permute <- permute_zfst(pied_coll_zfst_auto, pied_coll_fst_auto)

## Get pvalue

par_taig_zfst_pval <- pval_zfst(par_taig_zfst_permute, par_taig_zfst_smooth)
pied_coll_zfst_pval <- pval_zfst(pied_coll_zfst_permute, pied_coll_zfst_smooth)

## Write out significant windows
par_taig_sig_wins <- par_taig_fst_auto[par_taig_zfst_pval<=p_threshold,]
pied_coll_sig_wins <- pied_coll_fst_auto[pied_coll_zfst_pval<=p_threshold,]

write.table(par_taig_sig_wins, 'par_taig.sig_fst_peaks_auto.w200kb_s200kb.txt', col.names = T, row.names = F, quote = F, sep = '\t')
write.table(pied_coll_sig_wins , 'pied_coll.sig_fst_peaks_auto.w200kb_s200kb.txt', col.names = T, row.names = F, quote = F, sep = '\t')
```

After writing the set of Fst peaks to a new file, I merged adjacent windows using BEDTools with the following command, after converting the file to bed format.

```Bash
bedtools merge -i par_taig.sig_fst_peaks_auto.w200kb_s200kb.bed > par_taig.sig_fst_peaks_auto.w200kb_s200kb.adj_merged.bed
```

I then took the intersection of the Fst peaks for both species comparsions with BEDTools to identify Fst peaks that were shared between both comparisons.

```Bash
bedtools intersect -a par_taig.sig_fst_peaks_auto.w200kb_s200kb.adj_merged.bed -b pied_coll.sig_fst_peaks_auto.w200kb_s200kb.adj_merged.bed > pt_pc_sig_fst_peaks_auto_overlap.w200kb_s200kb.bed
```


## Plotting genomic landscape

Plots of each statistic in 200kb windows were created for collared flycatcher chromosomes in R.

```R
library(signal)
library(plyr)
library(dplyr)
library(ggplot2)
library(car)
library(circlize)


## The first couple hundred lines is just reading in all the separate datasets and combining into one dataset
## This means the windows that are kept for analysis in the end are only those for which there is data for all statistics
## This is a total of 4278 windows

#### Load pi data ####

pi_par <- read.table('par_pi.snps_only.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.txt.chrompos')
pi_taig <- read.table('taig_pi.snps_only.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.txt.chrompos')
pi_coll <- read.table('coll_pi.snps_only.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.txt.chrompos')
pi_pied <- read.table('pied_pi.snps_only.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.txt.chrompos')

colnames(pi_par) <- c('Chr', 'Start_par', 'Stop_par', 'SNPs_par', 'Pi_sum_par', 'Sites_par','Pi_par')
colnames(pi_taig) <- c('Chr', 'Start_taig', 'Stop_taig', 'SNPs_taig', 'Pi_sum_taig', 'Sites_taig', 'Pi_taig')
colnames(pi_coll) <- c('Chr', 'Start_coll', 'Stop_coll', 'SNPs_coll', 'Pi_sum_coll', 'Sites_coll', 'Pi_coll')
colnames(pi_pied) <- c('Chr', 'Start_pied', 'Stop_pied', 'SNPs_pied', 'Pi_sum_pied', 'Sites_pied', 'Pi_pied')

pi_par$Mid <- (pi_par$Start_par+pi_par$Stop_par)/2
pi_taig$Mid <- (pi_taig$Start_taig+pi_taig$Stop_taig)/2
pi_coll$Mid <- (pi_coll$Start_coll+pi_coll$Stop_coll)/2
pi_pied$Mid <- (pi_pied$Start_pied+pi_pied$Stop_pied)/2

pi_par <- pi_par[order(pi_par$Chr,pi_par$Mid),]
pi_taig <- pi_taig[order(pi_taig$Chr, pi_taig$Mid),]
pi_coll <- pi_coll[order(pi_coll$Chr, pi_coll$Mid),]
pi_pied <- pi_pied[order(pi_pied$Chr, pi_pied$Mid),]

pi_par <- pi_par[pi_par$Chr!='ChrZ' & pi_par$Chr!='ChrLGE22',]
pi_taig <- pi_taig[pi_taig$Chr!='ChrZ' & pi_taig$Chr!='ChrLGE22',]
pi_coll <- pi_coll[pi_coll$Chr!='ChrZ' & pi_coll$Chr!='ChrLGE22' & pi_coll$Chr!='ChrFal36',]
pi_pied <- pi_pied[pi_pied$Chr!='ChrZ',]

pi_all_comps_combined <- join(pi_pied, join(pi_coll, join(pi_par, pi_taig, by = c('Chr', 'Mid'), type = 'inner'), by = c('Chr', 'Mid'), type = 'inner'), by = c('Chr', 'Mid'), type = 'inner')


#### Load Dxy data ####

dxy_RT <- read.table('par_taig_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.txt.chrompos')
dxy_CP <- read.table('coll_pied_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.txt.chrompos')
dxy_CT <- read.table('coll_taig_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.txt.chrompos')
dxy_CR <- read.table('coll_par_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.txt.chrompos')
dxy_TP <- read.table('taig_pied_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.txt.chrompos')
dxy_RP <- read.table('par_pied_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.txt.chrompos')

colnames(dxy_RT) <- c("Chr", "Start_RT", "Stop_RT", "SNPs_RT", "dxy_sum_RT", "Sites_RT", "dxy_RT")
colnames(dxy_CP) <- c("Chr", "Start_CP", "Stop_CP", "SNPs_CP", "dxy_sum_CP", "Sites_CP", "dxy_CP")
colnames(dxy_CT) <- c("Chr", "Start_CT", "Stop_CT", "SNPs_CT", "dxy_sum_CT", "Sites_CT", "dxy_CT")
colnames(dxy_CR) <- c("Chr", "Start_CR", "Stop_CR", "SNPs_CR", "dxy_sum_CR", "Sites_CR", "dxy_CR")
colnames(dxy_TP) <- c("Chr", "Start_TP", "Stop_TP", "SNPs_TP", "dxy_sum_TP", "Sites_TP", "dxy_TP")
colnames(dxy_RP) <- c("Chr", "Start_RP", "Stop_RP", "SNPs_RP", "dxy_sum_RP", "Sites_RP", "dxy_RP")

dxy_RT$Mid <- (dxy_RT$Start_RT+dxy_RT$Stop_RT)/2
dxy_CP$Mid <- (dxy_CP$Start_CP+dxy_CP$Stop_CP)/2
dxy_CT$Mid <- (dxy_CT$Start_CT+dxy_CT$Stop_CT)/2
dxy_CR$Mid <- (dxy_CR$Start_CR+dxy_CR$Stop_CR)/2
dxy_TP$Mid <- (dxy_TP$Start_TP+dxy_TP$Stop_TP)/2
dxy_RP$Mid <- (dxy_RP$Start_RP+dxy_RP$Stop_RP)/2

dxy_RT <- dxy_RT[order(dxy_RT$Chr, dxy_RT$Mid),]
dxy_CP <- dxy_CP[order(dxy_CP$Chr, dxy_CP$Mid),]
dxy_CT <- dxy_CT[order(dxy_CT$Chr, dxy_CT$Mid),]
dxy_CR <- dxy_CR[order(dxy_CR$Chr, dxy_CR$Mid),]
dxy_TP <- dxy_TP[order(dxy_TP$Chr, dxy_TP$Mid),]
dxy_RP <- dxy_RP[order(dxy_RP$Chr, dxy_RP$Mid),]

dxy_RT <- dxy_RT[dxy_RT$Chr!='ChrZ' & dxy_RT$Chr!='ChrLGE22',]
dxy_CP <- dxy_CP[dxy_CP$Chr!='ChrZ',]
dxy_CT <- dxy_CT[dxy_CT$Chr!='ChrZ' & dxy_CT$Chr!='ChrLGE22',]
dxy_CR <- dxy_CR[dxy_CR$Chr!='ChrZ' & dxy_CR$Chr!='ChrLGE22',]
dxy_TP <- dxy_TP[dxy_TP$Chr!='ChrZ' & dxy_TP$Chr!='ChrLGE22',]
dxy_RP <- dxy_RP[dxy_RP$Chr!='ChrZ' & dxy_RP$Chr!='ChrLGE22',]

dxy_all_joined <- join(dxy_RP, join(dxy_TP, join(dxy_CR, join(dxy_CT, join(dxy_RT, dxy_CP, by = c('Chr', 'Mid'), type = 'inner'), by = c('Chr', 'Mid'), type = 'inner'), by = c('Chr', 'Mid'), type = 'inner'), by = c('Chr', 'Mid'), type = 'inner'), by = c('Chr', 'Mid'), type = 'inner')


#### Load Fst data ####

Fst_RT <- read.table('par_taig.minGQ.minDP5.maxDP200.biallelic.maxMiss10perc.femHetRmvd.RM.CM.fst_w200k_s200k_auto_scaffs_only.windowed.weir.fst.chrompos')
Fst_CP <- read.table('coll_pied.minGQ.minDP5.maxDP200.biallelic.maxMiss10perc.femHetRmvd.RM.CM.fst_w200k_s200k_auto_scaffs_only.windowed.weir.fst.chrompos')

colnames(Fst_RT) <- c('Chr', 'Start_fst_RT', 'Stop_fst_RT', 'Fst_sites_RT', 'Weight_fst_RT', 'Fst_RT')
colnames(Fst_CP) <- c('Chr', 'Start_fst_CP', 'Stop_fst_CP', 'Fst_sites_CP', 'Weight_fst_CP', 'Fst_CP')

Fst_RT$Mid <- (Fst_RT$Start_fst_RT+Fst_RT$Stop_fst_RT)/2
Fst_CP$Mid <- (Fst_CP$Start_fst_CP+Fst_CP$Stop_fst_CP)/2

Fst_RT <- Fst_RT[order(Fst_RT$Chr,Fst_RT$Mid),]
Fst_CP <- Fst_CP[order(Fst_CP$Chr, Fst_CP$Mid),]

Fst_RT <- Fst_RT[Fst_RT$Fst_sites_RT>=200,]
Fst_CP <- Fst_CP[Fst_CP$Fst_sites_CP>=200,]

Fst_RT <- Fst_RT[Fst_RT$Chr!='ChrLGE22',]
Fst_CP <- Fst_CP[Fst_CP$Chr!='ChrLGE22',]

fst_joined <- join(Fst_RT, Fst_CP, by = c('Chr', 'Mid'), type = 'inner')

#### Read in gene info ####

gene_info <- read.table('fAlb15_MTmask_mtfZan.RM.CM.exon.CNE.intron.intergenic.w200kb_s200kb.chrompos.with_peak_class.txt', header =T)
gene_info <- gene_info[gene_info$Chr!='Chr25' & gene_info$Chr!='ChrLGE22',]
gene_info$Mid <- (gene_info$Start+gene_info$Stop)/2

#### Read recombination rate in  ####

recom <- read.table('Chr.Rec.200kb.5kGap.txt.chrompos')
colnames(recom) <- c('Chr','Chr_recom', 'Start_recom', 'Stop_recom', 'cm.mb')
recom$Mid <- (recom$Start_recom + recom$Stop_recom)/2
recom <- recom[recom$Chr!='ChrZ' & recom$Chr!='ChrLGE22' & recom$Chr!='ChrFal34',]


#### Join all into a dataframe ####

all_stats_joined <- join(recom, join(gene_info, join(fst_joined, join(pi_all_comps_combined, dxy_all_joined, by = c('Chr', 'Mid'), type = 'inner'), by = c('Chr', 'Mid'), type = 'inner'), by = c('Chr', 'Mid'), type = 'inner'), by = c('Chr', 'Mid'), type = 'inner')
all_stats_joined$Sel_dens <- (all_stats_joined$Exon+all_stats_joined$CNE)/all_stats_joined$Total_sites


#### Making whole genome plot ####

# Set colours for plots

CP_dxy_col <- '#8B96F7'
RT_dxy_col <- '#FFA948'
pi_par_col <- "#FFCB48"
pi_taig_col <- "#FF6B48"
pi_pied_col <- "#B785F7"
pi_coll_col <- "#7EDBF5"

all_stats_joined <- all_stats_joined[all_stats_joined$Chr!='Chr26',]
all_stats_joined <- all_stats_joined[order(all_stats_joined$Chr, all_stats_joined$Mid),]
# List of chromosomes in order for plot
chr_list <- c('Chr1', 'Chr1A', 'Chr2', 'Chr3', 'Chr4', 'Chr4A', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10', 'Chr11', 'Chr12', 'Chr13', 'Chr14', 'Chr15', 'Chr17', 'Chr18', 'Chr19', 'Chr20', 'Chr21', 'Chr22', 'Chr23', 'Chr24', 'Chr27', 'Chr28')

# Get matrix of chromosome ranges
## Vector to hold the lengths of the chroms
chr_ends <- c()

for(chr in chr_list){
  chr_ends <- append(chr_ends, max(c(all_stats_joined$Stop_taig[all_stats_joined$Chr==chr],all_stats_joined$Start_taig[all_stats_joined$Chr==chr])))
}

## Matrix holding start and stop coordinates
chr_labs <- c('1', '1A', '2', '3', '4', '4A', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '17', '18', '19', '20', '21', '22', '23', '24', '27', '28')
chr_ranges <- matrix(c(rep(0,length(chr_list)), chr_ends), ncol=2)


# Set plotting parameters for circos

pdf('fig1_circos_plot.pdf')
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.par(gap.after=c(rep(1,26), 25))
circos.par(track.height=0.15)
circos.par(start.degree=90)

# Set up circos plot
circos.initialize(chr_list, xlim=chr_ranges)


# Add first track for pi plots
circos.track(ylim = c(0,0.0052))

# Plot pi values for each species along each chromosome
for(chr in chr_list){
  print(chr)
  circos.text(max(all_stats_joined$Mid[all_stats_joined$Chr==chr])/2, CELL_META$ylim[2] + mm_y(5), gsub("Chr", "", chr), chr)
  short_chrs <- c('Chr26', 'Chr27', 'Chr28')
  win_smooths <-ifelse(chr %in% short_chrs, 5, 7)
  pi_par_smooth <- sgolayfilt(all_stats_joined$Pi_par[all_stats_joined$Chr==chr], p = 3, n = win_smooths)
  pi_taig_smooth <- sgolayfilt(all_stats_joined$Pi_taig[all_stats_joined$Chr==chr], p = 3, n = win_smooths)
  pi_pied_smooth <- sgolayfilt(all_stats_joined$Pi_pied[all_stats_joined$Chr==chr], p = 3, n = win_smooths)
  pi_coll_smooth <- sgolayfilt(all_stats_joined$Pi_coll[all_stats_joined$Chr==chr], p = 3, n = win_smooths)

  circos.lines(all_stats_joined$Mid[all_stats_joined$Chr==chr], pi_coll_smooth, chr, col = pi_coll_col, lwd = 0.75)
  circos.lines(all_stats_joined$Mid[all_stats_joined$Chr==chr], pi_pied_smooth, chr, col = pi_pied_col, lwd = 0.75)
  circos.lines(all_stats_joined$Mid[all_stats_joined$Chr==chr], pi_par_smooth, chr, col = pi_par_col, lwd = 0.75)
  circos.lines(all_stats_joined$Mid[all_stats_joined$Chr==chr],  pi_taig_smooth, chr, col = pi_taig_col, lwd = 0.75)

}

# Add y axis
circos.yaxis('left', at = c(0,0.002,0.004), labels.cex = 0.75, sector.index = 'Chr1')


# Add second tract for dxy plots
circos.track(ylim=c(0,0.022))

for(chr in chr_list){
  short_chrs <- c('Chr26', 'Chr27', 'Chr28')
  win_smooths <-ifelse(chr %in% short_chrs, 5, 7)

  dxy_rt_smooth <- sgolayfilt(all_stats_joined$dxy_RT[all_stats_joined$Chr==chr], p = 3, n = win_smooths)
  dxy_cp_smooth <- sgolayfilt(all_stats_joined$dxy_CP[all_stats_joined$Chr==chr], p = 3, n = win_smooths)
  dxy_ct_smooth <- sgolayfilt(all_stats_joined$dxy_CT[all_stats_joined$Chr==chr], p = 3, n = win_smooths)

  circos.lines(all_stats_joined$Mid[all_stats_joined$Chr==chr], dxy_cp_smooth, chr, col = CP_dxy_col, lwd = 0.75)
  circos.lines(all_stats_joined$Mid[all_stats_joined$Chr==chr], dxy_rt_smooth, chr, col = RT_dxy_col, lwd = 0.75)
  circos.lines(all_stats_joined$Mid[all_stats_joined$Chr==chr],  dxy_ct_smooth, chr, col = 'magenta', lwd = 0.75)
}

circos.yaxis('left', at = c(0,0.01,0.02), labels.cex = 0.75, sector.index = 'Chr1')

# Add third tract for fst plots
circos.track(ylim=c(0,1))
for(chr in chr_list){
  short_chrs <- c('Chr26', 'Chr27', 'Chr28')
  win_smooths <-ifelse(chr %in% short_chrs, 5, 7)

  fst_cp_smooth <- sgolayfilt(all_stats_joined$Weight_fst_CP[all_stats_joined$Chr==chr], p = 3, n = win_smooths)
  fst_rt_smooth <- sgolayfilt(all_stats_joined$Weight_fst_RT[all_stats_joined$Chr==chr], p = 3, n = win_smooths)
  #circos.lines(all_stats_joined$Mid[all_stats_joined$Chr==chr], all_stats_joined$Weight_fst_CP[all_stats_joined$Chr==chr], chr, col = CP_dxy_col, lwd = 0.75)
  #circos.lines(all_stats_joined$Mid[all_stats_joined$Chr==chr], all_stats_joined$Weight_fst_RT[all_stats_joined$Chr==chr], chr, col = RT_dxy_col, lwd = 0.75)
  circos.lines(all_stats_joined$Mid[all_stats_joined$Chr==chr], fst_cp_smooth, chr, col = CP_dxy_col, lwd = 0.75)
  circos.lines(all_stats_joined$Mid[all_stats_joined$Chr==chr], fst_rt_smooth, chr, col = RT_dxy_col, lwd = 0.75)
}

circos.yaxis('left', at = c(0,0.5,1), labels.cex = 0.75, sector.index = 'Chr1')

# Add final track for recombination
circos.track(ylim=c(0,25))
for(chr in chr_list){
  circos.lines(all_stats_joined$Mid[all_stats_joined$Chr==chr], all_stats_joined$cm.mb[all_stats_joined$Chr==chr], chr, lwd = 0.75)
}

circos.yaxis('left', at = c(0,10,20), labels.cex = 0.75, sector.index = 'Chr1')

dev.off()

```

## Correlations among statistics

The correlations between diversity landscapes of different species, and differentiation and divergence landscapes between different species comparisons were performed in R.

```R
## Estimating correlations between pi for all species comparisons
cor.test(all_stats_joined$Pi_coll, all_stats_joined$Pi_pied)
cor.test(all_stats_joined$Pi_par, all_stats_joined$Pi_taig)
cor.test(all_stats_joined$Pi_coll, all_stats_joined$Pi_taig)
cor.test(all_stats_joined$Pi_coll, all_stats_joined$Pi_par)
cor.test(all_stats_joined$Pi_pied, all_stats_joined$Pi_par)
cor.test(all_stats_joined$Pi_pied, all_stats_joined$Pi_taig)

## Estimating correlations between fst and dxy for the two independent species comparisons
cor.test(all_stats_joined$Weight_fst_CP, all_stats_joined$Weight_fst_RT)
cor.test(all_stats_joined$dxy_CP, all_stats_joined$dxy_RT)

## Estimating correlation between mean pi and dxy for all species comparisons
CP_mean_pi <- rowMeans(cbind(all_stats_joined$Pi_pied,all_stats_joined$Pi_coll))
RT_mean_pi<- rowMeans(cbind(all_stats_joined$Pi_par,all_stats_joined$Pi_taig))
CR_mean_pi <- rowMeans(cbind(all_stats_joined$Pi_par,all_stats_joined$Pi_coll))
CT_mean_pi <- rowMeans(cbind(all_stats_joined$Pi_taig,all_stats_joined$Pi_coll))
PT_mean_pi <- rowMeans(cbind(all_stats_joined$Pi_pied,all_stats_joined$Pi_par))
PR_mean_pi <- rowMeans(cbind(all_stats_joined$Pi_pied,all_stats_joined$Pi_taig))

cor.test(CP_mean_pi, all_stats_joined$dxy_CP)
cor.test(RT_mean_pi, all_stats_joined$dxy_RT)
cor.test(CR_mean_pi, all_stats_joined$dxy_CR)
cor.test(CT_mean_pi, all_stats_joined$dxy_CT)
cor.test(PR_mean_pi, all_stats_joined$dxy_RP)
cor.test(PT_mean_pi, all_stats_joined$dxy_TP)
```


## Multiple linear regressions with recombination rate and functional density

Using the recombination map and gene and conserved noncoding element annotations from collared flycatcher, I performed multiple linear regressions to examine the relationship between these genomic features and diversity and differentiation measures. Analyses were performed in R in 50kb (regression with functional density only) and 200kb windows.

```R
#### Read in pi and Dxy with CNEs and exons masked ####

CR_neut_dxy <- read.table('coll_par_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.exon_cnes_masked.all_wins.txt.chrompos')
CP_neut_dxy <- read.table('coll_pied_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.exon_cnes_masked.all_wins.txt.chrompos')
CT_neut_dxy <- read.table('coll_taig_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.exon_cnes_masked.all_wins.txt.chrompos')
RP_neut_dxy <-  read.table('par_pied_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.exon_cnes_masked.all_wins.txt.chrompos')
RT_neut_dxy <-   read.table('par_taig_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.exon_cnes_masked.all_wins.txt.chrompos')
TP_neut_dxy <-  read.table('taig_pied_dxy.fixed_rmvd.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.exon_cnes_masked.all_wins.txt.chrompos')

coll_neut_pi<- read.table('coll_pi.snps_only.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.exon_cnes_masked.all_win.txt.chrompos')
par_neut_pi <- read.table('par_pi.snps_only.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.exon_cnes_masked.all_wins.txt.chrompos')
pied_neut_pi <- read.table('pied_pi.snps_only.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.exon_cnes_masked.all_wins.txt.chrompos')
taig_neut_pi <- read.table('taig_pi.snps_only.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.max_miss_10_perc.w200k_s200k.exon_cnes_masked.all_wins.txt.chrompos')

colnames(CR_neut_dxy) <- c('Chr', 'Start', 'Stop', 'CR_neut_dxy_seg_sites', 'CR_neut_dxy_sum', 'CR_neut_dxy_tot_sites', 'CR_neut_dxy')
colnames(CP_neut_dxy) <- c('Chr', 'Start', 'Stop', 'CP_neut_dxy_seg_sites', 'CP_neut_dxy_sum', 'CP_neut_dxy_tot_sites', 'CP_neut_dxy')
colnames(CT_neut_dxy) <- c('Chr', 'Start', 'Stop', 'CT_neut_dxy_seg_sites', 'CT_neut_dxy_sum', 'CT_neut_dxy_tot_sites', 'CT_neut_dxy')
colnames(RP_neut_dxy) <- c('Chr', 'Start', 'Stop', 'RP_neut_dxy_seg_sites', 'RP_neut_dxy_sum', 'RP_neut_dxy_tot_sites', 'RP_neut_dxy')
colnames(RT_neut_dxy) <- c('Chr', 'Start', 'Stop', 'RT_neut_dxy_seg_sites', 'RT_neut_dxy_sum', 'RT_neut_dxy_tot_sites', 'RT_neut_dxy')
colnames(TP_neut_dxy) <- c('Chr', 'Start', 'Stop', 'TP_neut_dxy_seg_sites', 'TP_neut_dxy_sum', 'TP_neut_dxy_tot_sites', 'TP_neut_dxy')

colnames(coll_neut_pi) <- c('Chr', 'Start', 'Stop', 'coll_neut_pi_seg_sites', 'coll_neut_pi_sum', 'coll_neut_pi_tot_sites', 'coll_neut_pi')
colnames(par_neut_pi) <- c('Chr', 'Start', 'Stop', 'par_neut_pi_seg_sites', 'par_neut_pi_sum', 'par_neut_pi_tot_sites', 'par_neut_pi')
colnames(pied_neut_pi) <- c('Chr', 'Start', 'Stop', 'pied_neut_pi_seg_sites', 'pied_neut_pi_sum', 'pied_neut_pi_tot_sites', 'pied_neut_pi')
colnames(taig_neut_pi) <- c('Chr', 'Start', 'Stop', 'taig_neut_pi_seg_sites', 'taig_neut_pi_sum', 'taig_neut_pi_tot_sites', 'taig_neut_pi')


neut_dxy_joined <- join(TP_neut_dxy, join(RT_neut_dxy, join(RP_neut_dxy, join(CT_neut_dxy, join(CR_neut_dxy, CP_neut_dxy, by = c('Chr', 'Start'), type = 'inner'), by = c('Chr', 'Start'), type = 'inner'), by = c('Chr', 'Start'), type = 'inner'), by = c('Chr', 'Start'), type = 'inner'), by = c('Chr', 'Start'), type = 'inner')
neut_dxy_joined <- neut_dxy_joined[,-8]

neut_pi_joined <- join(taig_neut_pi, join(pied_neut_pi, join(coll_neut_pi, par_neut_pi, by = c('Chr', 'Start'), type = 'inner'), by = c('Chr', 'Start'), type = 'inner'), by = c('Chr', 'Start'), type = 'inner')
neut_pi_joined <- neut_pi_joined[,-8]

neut_pi_dxy_joined <- join(neut_dxy_joined, neut_pi_joined, by = c('Chr', 'Start'), type = 'inner')
neut_pi_dxy_joined <- neut_pi_dxy_joined[,-28]
neut_pi_dxy_joined$Mid <- (neut_pi_dxy_joined$Start+neut_pi_dxy_joined$Stop)/2

all_stats_joined_neut_pi_dxy <- join(all_stats_joined, neut_pi_dxy_joined, by = c('Chr', 'Mid'), type = 'inner')


#### Read in neutral fst ####

CP_fst_neut <- read.table('coll_pied.minGQ.minDP5.maxDP200.biallelic.maxMiss10perc.femHetRmvd.RM.CM.fst_w200k_s200k_auto_scaffs_only.exon_cnes_masked.windowed.weir.fst.chrompos')
RT_fst_neut <- read.table('par_taig.minGQ.minDP5.maxDP200.biallelic.maxMiss10perc.femHetRmvd.RM.CM.fst_w200k_s200k_auto_scaffs_only.exon_cnes_masked.windowed.weir.fst.chrompos')

colnames(CP_fst_neut) <- c('Chr', 'Start', 'Stop', 'Fst_sites_cp_neut', 'Weight_fst_cp_neut', 'Fst_cp_neut')
colnames(RT_fst_neut) <- c('Chr', 'Start', 'Stop', 'Fst_sites_rt_neut', 'Weight_fst_rt_neut', 'Fst_rt_neut')

CP_fst_neut$Mid <- (CP_fst_neut$Start+CP_fst_neut$Stop)/2
RT_fst_neut$Mid <- (RT_fst_neut$Start+RT_fst_neut$Stop)/2

CP_fst_neut <- CP_fst_neut[order(CP_fst_neut$Chr, CP_fst_neut$Mid),]
RT_fst_neut <- RT_fst_neut[order(RT_fst_neut$Chr, RT_fst_neut$Mid),]

CP_fst_neut <- CP_fst_neut[CP_fst_neut$Fst_sites_cp_neut>=200,]
RT_fst_neut <- RT_fst_neut[RT_fst_neut$Fst_sites_rt_neut>=200,]

neut_fst_joined <- join(CP_fst_neut, RT_fst_neut, by = c("Chr", "Mid"), type = 'inner')
neut_fst_joined <- neut_fst_joined[,-c(8,9)]

all_stats_joined_neut_fst <- join(neut_fst_joined, all_stats_joined_neut_pi_dxy, by = c('Chr', 'Mid'), type = 'inner')

#### Relationship between neutral pi and genomic architecture ####

coll_neut_pi_gen_arch_lm_tfd <- lm(all_stats_joined_neut_fst$coll_neut_pi ~ log(all_stats_joined_neut_fst$cm.mb+1) + sqrt(all_stats_joined_neut_fst$Sel_dens))
par_neut_pi_gen_arch_lm_tfd <- lm(all_stats_joined_neut_fst$par_neut_pi ~ log(all_stats_joined_neut_fst$cm.mb+1) + sqrt(all_stats_joined_neut_fst$Sel_dens))
pied_neut_pi_gen_arch_lm_tfd <- lm(all_stats_joined_neut_fst$pied_neut_pi ~ log(all_stats_joined_neut_fst$cm.mb+1) + sqrt(all_stats_joined_neut_fst$Sel_dens))
taig_neut_pi_gen_arch_lm_tfd <- lm(all_stats_joined_neut_fst$taig_neut_pi ~ log(all_stats_joined_neut_fst$cm.mb+1) + sqrt(all_stats_joined_neut_fst$Sel_dens))

#### Relationship between neutral dxy and genomic architecture ####
CP_neut_dxy_gen_arch_lm_tfd <- lm(all_stats_joined_neut_fst$CP_neut_dxy ~ log(all_stats_joined_neut_fst$cm.mb+1) + sqrt(all_stats_joined_neut_fst$Sel_dens))
RT_neut_dxy_gen_arch_lm_tfd <- lm(all_stats_joined_neut_fst$RT_neut_dxy ~ log(all_stats_joined_neut_fst$cm.mb+1) + sqrt(all_stats_joined_neut_fst$Sel_dens))
CT_neut_dxy_gen_arch_lm_tfd <- lm(all_stats_joined_neut_fst$CT_neut_dxy ~ log(all_stats_joined_neut_fst$cm.mb+1) + sqrt(all_stats_joined_neut_fst$Sel_dens))
CR_neut_dxy_gen_arch_lm_tfd <- lm(all_stats_joined_neut_fst$CR_neut_dxy ~ log(all_stats_joined_neut_fst$cm.mb+1) + sqrt(all_stats_joined_neut_fst$Sel_dens))
PT_neut_dxy_gen_arch_lm_tfd <- lm(all_stats_joined_neut_fst$TP_neut_dxy ~ log(all_stats_joined_neut_fst$cm.mb+1) + sqrt(all_stats_joined_neut_fst$Sel_dens))
PR_neut_dxy_gen_arch_lm_tfd <- lm(all_stats_joined_neut_fst$RP_neut_dxy ~ log(all_stats_joined_neut_fst$cm.mb+1) + sqrt(all_stats_joined_neut_fst$Sel_dens))


#### Relationship between neutral fst and genomic architecture ####
CP_neut_fst_gen_arch_lm <- lm(all_stats_joined_neut_fst$Weight_fst_cp_neut ~ log(all_stats_joined_neut_fst$cm.mb+1) + all_stats_joined_neut_fst$Sel_dens)
RT_neut_fst_gen_arch_lm <- lm(all_stats_joined_neut_fst$Weight_fst_rt_neut ~ log(all_stats_joined_neut_fst$cm.mb+1) + all_stats_joined_neut_fst$Sel_dens)

```
