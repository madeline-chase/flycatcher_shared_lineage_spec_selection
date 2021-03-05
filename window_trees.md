# Window based trees

We analyzed phylogenetic discordance by estimating phylogenetic trees using [MVFtools](https://github.com/jbpease/mvftools).

## Converting VCF to MVF

To run MVFtools I first had to convert my VCF file to MVF format, which I did for each autosomal collared flycatcher scaffold separately with the following command, followed by estimating the phylogeny in 50kb and 200kb windows.

```Bash
for scaff in $scaffs
do

  vcftools --vcf par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.minDP5.minDP200.auto_scaffs.missing_samps_rmvd.recode.vcf --recode --recode-INFO-all --out par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.minDP5.minDP200.auto_scaffs.missing_samps_rmvd.$scaff --chr $scaff

  python mvftools.py ConvertVCF2MVF --vcf par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.minDP5.minDP200.auto_scaffs.missing_samps_rmvd.$scaff.recode.vcf --out par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.minDP5.minDP200.auto_scaffs.missing_samps_rmvd.$scaff.mvf

  python mvftools.py InferTree --mvf par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.minDP5.minDP200.auto_scaffs.missing_samps_rmvd.$scaff.mvf --windowsize 200000 --root-with F_hyperythra_X --raxml-path raxmlHPC-PTHREADS --raxml-opts "-T 10" --out mvf_win_tree.w200k.$scaff.pthreads10.missing_samps_rmvd

done
```


## Analyzing tree topologies

From the mvftools output, I then extracted only the topologies and used the script `analyze_window_trees.py` to root trees with *Ficedula hyperythra* and analyze the topologies.

```Bash
for scaff in $scaffs
do

  grep -v 'tree' mvf_win_tree.w200k.$scaff.pthreads10.missing_samps_rmvd | cut -f4 > mvf_win_tree.w200k.$scaff.pthreads10.missing_samps_rmvd.trees_only

  python analyze_window_trees_testing.py mvf_win_tree.w200k.$scaff.pthreads10.missing_samps_rmvd.trees_only $scaff > mvf_win_tree.w200k.$scaff.pthreads10.missing_samps_rmvd.topology_results.corrected.txt

done
```

Finally, the window tree outputs were combined with the other window summary statistic data (Fst, pi, dxy, recombination rate, functional density), so trees in windows with poor sequencing quality were removed. These were then used to count the number of windows showing the species tree topology, and where some splits were not resolved.

```R
#### Read in window tree data ####
win_trees_200k <- read.table('mvf_win_tree.w200k.scaffs_combined.missing_samps_rmvd.topology_results.corrected.txt.chrompos', sep = '\t')

colnames(win_trees_200k) <- c('Chr', 'Start', 'Stop', 'Scaff', 'Win_num' ,'Topo')
win_trees_200k$Mid <- (win_trees_200k$Start+win_trees_200k$Stop)/2

#### Combine window trees with all other statistics ####
all_stats_joined_neuts_trees <- join(all_stats_joined_neut_fst, win_trees_200k, by = c("Chr", 'Mid'), type = 'inner')
```
