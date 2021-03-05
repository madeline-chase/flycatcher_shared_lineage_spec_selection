# Identification of collapsed duplicated regions

I identified regions of the genome that showed signatures of being duplicates that were collapsed in the collared flycatcher assembly.

To identify these regions, I used data for collared and taiga flycatcher, after applying a minor allele frequency threshold of 10%, using VCFtools. This was performed separately for male and female samples, as some regions may contain collapsed duplicates associated with W sequences, which will only appear in female samples.

```Bash
vcftools --vcf $par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.all_chr.scaff_sorted.vcf --recode --recode-INFO-all --maf 0.1 --exclude-positions coll_missing_sites.max_miss_10_perc.txt --keep gotland.1993.2015.samples.txt --out coll_only.minDP5.maxDP200.minGQ.max_miss_10_perc.maf_10_perc
```

I then converted this from a vcf file to a .geno file using a script of Simon Martin's, parseVCF.py found [here](https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing).

```Bash
python2 parseVCF.py -i coll_only.minDP5.maxDP200.minGQ.max_miss_10_perc.maf_10_perc.recode.vcf > coll_only.minDP5.maxDP200.minGQ.max_miss_10_perc.maf_10_perc.geno
```

I then used the python script `GetHetProportions.py` to output the proportion of heterozygous genotypes for each site. I provide a list of samples to include, and separately provide the male and female samples.

```Bash
GetHetProportions.py -i coll_only.minDP5.maxDP200.minGQ.max_miss_10_perc.maf_10_perc.geno -g Coll_fem 93F24,93F26,93F30,93F32,93F34,93F35,93F42,93F44,93F45,93F47,93F54,93F56,93F59,93F74,93F75,93F77,93F82,93F88,93F89,93F90,93F92,93F93,93F94,15F129,15F130,15F131,15F135,15F142,15F143,15F145,15F149,15F151,15F17,15F18,15F21,15F22,15F23,15F24,15F25,15F29,15F447,15F448,15F450,15F453,15F457,15F459,15F460 > coll_fem.minDP5.maxDP200.minGQ.max_miss_10_perc.maf_10_perc.het_prop.txt
```

I then used the python script `collapsed_regions.py` to identify putatively collapsed regions. This command specifies that the excess heterozygosity threshold is 69%, that a collapsed region must contain at least 3 SNPs with heterozygosity above this threshold, that there can be no more than 5 SNPs below the excess heterozygosity threshold between and two of these SNPs, and that any two excessively heterozygous SNPs can be no more than 1000bp apart, and this adds 2000bp extension at the ends of the collapsed regions.

```Bash
collapsed_regions.py -i coll_fem.minDP5.maxDP200.minGQ.max_miss_10_perc.maf_10_perc.het_prop.txt -e 3 -n 5 -c 0.69 -d 1000 -x 2000 > coll.fem_collapsed.e69.e3.n5.d1000.x2000.bed
```

Because of the extension step, it's possibe that there are overlapping regions so I first merged this output using bedtools.

```Bash
bedtools merge -i coll.fem_collapsed.e69.e3.n5.d1000.x2000.bed > coll.fem_collapsed.e69.e3.n5.d1000.x2000.merged.bed
```

After running this script for both sexes and species, I then merged the regions that were identified as collapsed for everyone using BEDOPS, merging regions into one when they were less than 5000bp apart.

```Bash
bedops --merge coll.fem_collapsed.c69.e3.n5.d1000.x2000.merged.bed coll.male_collapsed.c69.e3.n5.d1000.x2000.merged.bed > coll.fem_male_collapsed_merge.c69.e3.n5.d1000.x2000.bed

bedops --merge coll.fem_male_collapsed_merge.c69.e3.n5.d1000.x2000.bed taig.fem_male_collapsed_merge.c69.e3.n5.d1000.x2000.bed > coll_taig.collapsed_merge.c69.e3.n5.d1000.x2000.bed
```

I then again used BEDTools to merge regions within this file that were less than 5000bp apart.

```Bash
bedtools merge -i coll_taig.collapsed_merge.c69.e3.n5.d1000.x2000.bed -d 5000 > coll_taig.collapsed_merge.c69.e3.n5.d1000.x2000.m5000.bed
```

This file was then used to filter SNPs located within collapsed regions.
