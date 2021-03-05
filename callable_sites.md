# Getting callable sites

To accurately calculate pi and dxy, it is important to account for missing sites in the denominator. I estimated the number of callable sites for each species based on depth filters from alignment files, and by masking repeats and collapsed regions in the genome assemlby.


## Masking repeats and collapsed regions

Using BEDTools maskfasta, I set regions of the genome assembly annotated as repetitive elements or collapsed regions (see `collapsed_regions.md`) to 'Ns'.

```Bash
bedtools maskfasta -fi $REFERENCE -bed collapsed.repeats.merged.bed -fo fAblb15_rm_cm.fasta
```

## Identifying noncallable sites

To filter on depth, I ran the workflow in `Snakefile_coll_callable_sites` to get the read depth for each site with samtools, and then print only sites with a min of 5x and a max of 200x coverage to match SNP filters. I then used BEDTools to merge sites that were able to be called in a sample. I then used the `multiIntersectBed` command to intersect all callable sites across samples within a species, and extracted only the regions that were callable in 90% of the samples. I then used bedtools subtract, to subtract these regions from the reference assembly to get the regions that were not callable within a species.

Again using bedtools maskfasta, I used these noncallable regions to mask the assembly.

```Bash
bedtools maskfasta -fi fAblb15_rm_cm.fasta -bed coll.minDP5.maxDP200.mult_int.90_perc.noncallable.bed -fo fAblb15_rm_cm.noncall_mask_coll.auto_only.fasta
```

I then used the python script `count_non_N_sites_total.py` to count the total number of callable sites for each species.

```Bash
python3 count_non_N_sites_total.py fAblb15_rm_cm.noncall_mask_coll.auto_only.fasta
```

To get the callable sites in windows, I first used BEDTools getfasta to split the masked reference file into predefined genomic windows.

```Bash
bedtools getfasta -fi fAblb15_rm_cm.noncall_mask_coll.auto_only.fasta -bed fAlb15.w200kb_s200kb.bed -fo fAblb15_rm_cm.noncall_mask_coll.auto_only.w200k_s200k.fasta
```

I then used the python script `count_non_N_sites_windowed.py` to print the total number of callable sites in each window, and removed windows for which less than 50% of the window was callable.

```Bash
python3 count_non_N_sites_windowed.py fAblb15_rm_cm.noncall_mask_coll.auto_only.w200k_s200k.fasta > callable_sites.coll.w200kb_s200kb.bed

awk -v OFS='\t' '{if($4>($3-$2)/2) print}' callable_sites.coll.w200kb_s200kb.bed > callable_sites.coll.w200kb_s200kb.min_50_perc_window_called.bed
```


## Getting total sites with exons and CNEs masked

I also needed to get the total number of sites with exons and CNEs masked, which I did using BEDTools maskfasta.

```Bash
bedtools maskfasta -fi fAblb15_rm_cm.noncall_mask_coll.auto_only.fasta -bed exons.CNEs.merged.bed -fo fAblb15_rm_cm.noncall_mask_coll.auto_only.exons_cnes_masked.fasta
```

I then used the same python scripts to get the total number of callable sites, and the number of callable sites in genomic windows.
