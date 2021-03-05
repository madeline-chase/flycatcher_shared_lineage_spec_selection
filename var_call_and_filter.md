# Variant calling and filtration

## Alignment

Red-breasted and taiga reads were aligned to the collared flycatcher assembly version xx using bwa-mem by NBIS staff.

```Bash
bwa mem $LOCAL_REFERENCE $R1path $R2path -t 5 -R '@RG\tID:170331_ST-E00216_0183_BHHTFCALXX.90_S11_L001\tPL:ILLUMINA\tLB:SX876_90.v1\tSM:Sample_90' | samtools view -h -T $LOCAL_REFERENCE -b -@ 5 - | samtools sort -O bam -T ${ALIGNMENT}.tmp -@ 5 - > $ALIGNMENT
```

The resulting BAM files were then merged and duplicates were marked in the same step, using Picard v2.0.1.

```Bash
java -jar $PICARD_HOME/picard.jar MarkDuplicates I=$LBAM1 I=$LBAM2 I=$LBAM3 I=$LBAM4  O=$ALLBAM METRICS_FILE=$METRICS REMOVE_DUPLICATES=false ASSUME_SORTED=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900
```


## Variant calling

Variant calling was performed in GATK v3.7 following best practices. First, BQSR was performed for red-breasted and taiga flycatchers by running HaplotypeCaller, followed by GenotypeGVCFs by NBIS staff.

```Bash
java -Xmx120g -Djava.io.tmpdir=$SNIC_TMP -jar $GATK -T HaplotypeCaller -R $LOCAL_REFERENCE -I $ALLBAM -o Sample_90.recal_0.g.vcf -ERC GVCF -nct 16
```

GenotypeGVCFs was performed with GATK v3.5.0
```Bash
java -Xmx256g -Djava.io.tmpdir=$SNIC_TMP -jar $GATK -T GenotypeGVCFs -R $LOCAL_REFERENCE -V SAMPLES_1-9.combined.recal_0.g.vcf -V SAMPLES_10-29.combined.recal_0.g.vcf -V SAMPLES_30-49.combined.recal_0.g.vcf -V SAMPLES_50-69.combined.recal_0.g.vcf -V SAMPLES_70-92.combined.recal_0.g.vcf -o $OUTFILE -nt 16
```

These initial variants were filtered according to GATK hard filtering thresholds (), and used to run BQSR.

```Bash
java -Xmx56g -Djava.io.tmpdir=$SNIC_TMP -jar $GATK -T BaseRecalibrator -R $LOCAL_REFERENCE -I $BAM -knownSites ALL_SAMPLES.merged.recal_0.filtered_indels.AC2.vcf -o Sample_90.recal_1.table -nct 8
```

PostBQSR, HaplotypeCaller was run a second time on recalibrated bam files for red-breasted and taiga flycatchers.

```Bash
java -Xmx120g -Djava.io.tmpdir=$SNIC_TMP -jar $GATK -T HaplotypeCaller -R $LOCAL_REFERENCE -I $ALLBAM -o Sample_90.recal_1.g.vcf -ERC GVCF -nct 16 -BQSR $RECAL_TABLE
```

Additionally, at this stage HaplotypeCaller was performed for pied and *Ficedula hyperythra*, using the retrieved postBQSR bam files.

```Bash
java -Xmx60g -jar $GATK -T HaplotypeCaller -R $LOCAL_REFERENCE -I $BAM -o SP_11_M.g.vcf -ERC gvcf -nct 10
```

After the second round of HaplotypeCaller, GenotypeGVCFs was performed with all samples, now including GVCF files for collared flycatcher that had previously undergone BQSR. For computational efficiency, variant calling was performed on multiple chunks of scaffolds at a time, with the vcf files then combined after.

```Bash
java -Xmx250g -Djava.io.tmpdir=$SNIC_TMP -jar $GATK -R $REFERENCE -T GenotypeGVCFs -V Taig.ANA.recal_1.HC.combined.gvcf -V Taig.YAK.recal_1.HC.combined.gvcf -V Taig.DOR.recal_1.HC.combined.gvcf -V Taig.IRK.recal_1.HC.combined.gvcf -V Taig.KAM.recal_1.HC.combined.gvcf -V Taig.MAG.recal_1.HC.combined.gvcf -V Taig.MAR.recal_1.HC.combined.gvcf -V Taig.small_pops.recal_1.HC.combined.gvcf -V Sample_93F24.g.vcf.gz -V Sample_93F26.g.vcf.gz -V Sample_93F30.g.vcf.gz -V Sample_93F32.g.vcf.gz -V Sample_93F34.g.vcf.gz -V Sample_93F35.g.vcf.gz -V Sample_93F42.g.vcf.gz -V Sample_93F44.g.vcf.gz -V Sample_93F45.g.vcf.gz -V Sample_93F47.g.vcf.gz -V Sample_93F54.g.vcf.gz -V Sample_93F56.g.vcf.gz -V Sample_93F59.g.vcf.gz -V Sample_93F74.g.vcf.gz -V Sample_93F75.g.vcf.gz -V Sample_93F77.g.vcf.gz -V Sample_93F82.g.vcf.gz -V Sample_93F88.g.vcf.gz -V Sample_93F89.g.vcf.gz -V Sample_93F90.g.vcf.gz -V Sample_93F92.g.vcf.gz -V Sample_93F93.g.vcf.gz -V Sample_93F94.g.vcf.gz -V Sample_93M25.g.vcf.gz -V Sample_93M27.g.vcf.gz -V Sample_93M28.g.vcf.gz -V Sample_93M29.g.vcf.gz -V Sample_93M36.g.vcf.gz -V Sample_93M38.g.vcf.gz -V Sample_93M39.g.vcf.gz -V Sample_93M40.g.vcf.gz -V Sample_93M41.g.vcf.gz -V Sample_93M46.g.vcf.gz -V Sample_93M53.g.vcf.gz -V Sample_93M55.g.vcf.gz -V Sample_93M58.g.vcf.gz -V Sample_93M71.g.vcf.gz -V Sample_93M72.g.vcf.gz -V Sample_93M73.g.vcf.gz -V Sample_93M78.g.vcf.gz -V Sample_93M79.g.vcf.gz -V Sample_93M80.g.vcf.gz -V Sample_93M81.g.vcf.gz -V Sample_93M83.g.vcf.gz -V Sample_93M84.g.vcf.gz -V Sample_93M86.g.vcf.gz -V Sample_93M91.g.vcf.gz -V Sample_15F129.g.vcf.gz -V Sample_15F130.g.vcf.gz -V Sample_15F131.g.vcf.gz -V Sample_15F135.g.vcf.gz -V Sample_15F142.g.vcf.gz -V Sample_15F143.g.vcf.gz -V Sample_15F145.g.vcf.gz -V Sample_15F149.g.vcf.gz -V Sample_15F151.g.vcf.gz -V Sample_15F17.g.vcf.gz -V Sample_15F18.g.vcf.gz -V Sample_15F21.g.vcf.gz -V Sample_15F22.g.vcf.gz -V Sample_15F23.g.vcf.gz -V Sample_15F24.g.vcf.gz -V Sample_15F25.g.vcf.gz -V Sample_15F29.g.vcf.gz -V Sample_15F447.g.vcf.gz -V Sample_15F448.g.vcf.gz -V Sample_15F450.g.vcf.gz -V Sample_15F453.g.vcf.gz -V Sample_15F457.g.vcf.gz -V Sample_15F459.g.vcf.gz -V Sample_15F460.g.vcf.gz -V Sample_15M153.g.vcf.gz -V Sample_15M155.g.vcf.gz -V Sample_15M158.g.vcf.gz -V Sample_15M160.g.vcf.gz -V Sample_15M161.g.vcf.gz -V Sample_15M162.g.vcf.gz -V Sample_15M163.g.vcf.gz -V Sample_15M201.g.vcf.gz -V Sample_15M202.g.vcf.gz -V Sample_15M203.g.vcf.gz -V Sample_15M204.g.vcf.gz -V Sample_15M207.g.vcf.gz -V Sample_15M468.g.vcf.gz -V Sample_15M469.g.vcf.gz -V Sample_15M475.g.vcf.gz -V Sample_15M477.g.vcf.gz -V Sample_15M49.g.vcf.gz -V Sample_15M537.g.vcf.gz -V Sample_15M568.g.vcf.gz -V Sample_15M571.g.vcf.gz -V Sample_15M573.g.vcf.gz -V Sample_15M589.g.vcf.gz -V Sample_15M684.g.vcf.gz -V Sample_15M724.g.vcf.gz -V Sample_31.recal_1.g.vcf -V Sample_32.recal_1.g.vcf -V Sample_4.recal_1.g.vcf -V Sample_48.recal_1.g.vcf -V Sample_49.recal_1.g.vcf -V Sample_51.recal_1.g.vcf -V Sample_52.recal_1.g.vcf -V Sample_53.recal_1.g.vcf -V Sample_54.recal_1.g.vcf -V Sample_55.recal_1.g.vcf -V Sample_56.recal_1.g.vcf -V Sample_77.recal_1.g.vcf -V Sample_82.recal_1.g.vcf -V Sample_83.recal_1.g.vcf -V Sample_92.recal_1.g.vcf -V SP_11_M.g.vcf -V SP_12_M.g.vcf -V SP_13_M.g.vcf -V SP_14_M.g.vcf -V SP_15_M.g.vcf -V SP_16_F.g.vcf -V SP_18-2_F.g.vcf -V SP_19_F.g.vcf -V SP_20-1_F.g.vcf -V SP_SP17-2_F.g.vcf -V SP_SV10_F.g.vcf -V SP_SV1_M.g.vcf -V SP_SV2_M.g.vcf -V SP_SV3_M.g.vcf -V SP_SV4_M.g.vcf -V SP_SV5_M.g.vcf -V SP_SV6_F.g.vcf -V SP_SV7_F.g.vcf -V SP_SV8_F.g.vcf -V SP_SV9_F.g.vcf -V F_hyperythra.g.vcf -L genotypeGVCF_4spec.region1.intervals -nt 1 -o $OUTPUT
```

## Variant filtration

After variant calling was finished, SNPs were extracted from the output and the same hard filtering thresholds were applied as in the initial variant calling for BQSR.

```Bash
vcftools --vcf par.taig.coll.pied.hyp.joint_genotype.combined_scaff.vcf --remove-indels --recode --recode-INFO-all --out par.taig.coll.pied.hyp.SNPs.unfilt
```

Applying hard filters.

```Bash
java -jar $GATK -T VariantFiltration -R $REFERENCE -V par.taig.coll.pied.hyp.SNPs.unfilt.recode.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "hardFilter" -o par.taig.coll.pied.hyp.SNPs.hard_filt.vcf
```

I then applied a number of more filters to the data using VCFtools. Minimum sequencing depth was set to 5, maximum depth was set to 200, and minimum genotype quality was set to 30. I removed SNPs that overlapped with annotated repeats and collapsed regions (see `collapsed_regions.md` for identification of collapsed regions).

```Bash
vcftools --vcf par.taig.coll.pied.hyp.SNPs.hard_filt.vcf --remove-filtered-all --recode --recode-INFO-all --out par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.all_chr.scaff_sorted --max-alleles 2 --minDP 5 --maxDP 200 --minGQ 30 --remove-filtered-geno-all --exclude-bed fAlb15.fa.out.bed
```

I then output the proportion of missing data per individual to identify if any individuals had high proportions of missing data.

```Bash
vcftools --vcf par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.all_chr.scaff_sorted.vcf --missing-indv --out par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.maxDP200
```

Based on these results, I removed samples with greater than 20% missing data, which included 9 pied flycatchers and 7 taiga flycatchers.

I then used VCFtools to output count data for each species, and identified sites with greater than 10% of samples with missing data, and wrote these to files to be excluded from analyses of those species.

```Bash
vcftools --vcf par.taig.coll.pied.hyp.SNPs.PASS.biallelic.minDP5.repeats_rmvd.minGQ.all_chr.scaff_sorted.vcf --counts --keep gotland.1993.2015.samples.txt --out coll_freq.minDP5.maxDP200.minGQ --exclude-positions missing_sites.CM.fem_het_Z.txt

awk '{if(NR>1) {if($4<171) print $1,$2}}' coll_freq.RM.CM.minDP5.maxDP200.minGQ.fem_het_rmvd.frq.count > missing_sites.CM.fem_het_Z.coll_max_miss_10_perc.txt
```

Although the final VCF includes all flycatcher scaffolds, all analyses only made use of scaffolds that map to autosomes according to the collared flycatcher genome assembly.
