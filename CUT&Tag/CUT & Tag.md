
This is the analysis pipeline of CUT & Tag applied in ME2 project.

In this file, only one sample is taken as an example. The program actually run in batch.


__1  Quality control of the fastq file__
```
fastqc H3K4me3_PBS_WT_rep2_R1.fastq.gz H3K4me3_PBS_WT_rep2_R2.fastq.gz
```

__2  Summarize the quality control reports__
```
multiqc .
```

__3  Build index of reference genome with Bowtie2__
```
bowtie2-build ../GRCm39_primary_assembly_genome.fa GRCm39
```

__4  Align the fastq file to the reference genome__

Remember to add the read group (RG) in the head of the BAM file, since __MarkDuplicates__ function in Picard asks the RG in the head of the BAM file
```
bowtie2 --local --very-sensitive --no-unal --no-mixed --no-discordant -p 40 -x /exec/wangyy/software/ChIP_seq/RefGenome/mouse/Bowtie2_index/GRCm39 -1 H3K4me3_PBS_WT_rep2_R1.fastq.gz -2 H3K4me3_PBS_WT_rep2_R2.fastq.gz --rg-id 1 --rg SM:H3K4me3_PBS_WT_rep2 --rg PL:ILLUMINA -S H3K4me3_PBS_WT_rep2.sam
```

__5  Convert SAM file to BAM file__
```
samtools view -bS --threads 40 -o H3K4me3_PBS_WT_rep2.bam H3K4me3_PBS_WT_rep2.sam
```

__6  Filter the multi-mapped reads__

BAM file don't have to be sorted or indexed
```
sambamba view -h -t 20 -f bam -F "[XS] == null" H3K4me3_PBS_WT_rep2.bam > H3K4me3_PBS_WT_rep2_filterhalf.bam
```

__7  Sort BAM file according to the name__

Since __fixmate__ funtion in samtools asks for the name-sorted BAM file
```
samtools sort --threads 20 -n H3K4me3_PBS_WT_rep2_filterhalf.bam -o H3K4me3_PBS_WT_rep2_filterhalf_namesorted.bam
```

__8  Mark the mates of the multi-mapped reads__

BAM file need to be name-sorted as documentation asked.

Actually, unsorted BAM files also make sense.

Samtools after 1.21 discarded this functionality.
```
samtools fixmate --threads 20 H3K4me3_PBS_WT_rep2_filterhalf_namesorted.bam H3K4me3_PBS_WT_rep2_filterhalf_namesorted_fixmate.bam
```

__9  Filter the mates of the multi-mapped reads__
```
samtools view -b -f 2 H3K4me3_PBS_WT_rep2_filterhalf_namesorted_fixmate.bam > H3K4me3_PBS_WT_rep2_filterpair_namesorted.bam
```

__10  Sort BAM files by coordinates__
```
samtools sort --threads 20 H3K4me3_PBS_WT_rep2_filterpair_namesorted.bam -o H3K4me3_PBS_WT_rep2_filterpair_sorted.bam
```

__11  Remove the duplicates__

BAM files must be sorted by coordinates but do not have to be indexed.
```
java -jar /exec/wangyy/software/ChIP_seq/Picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=H3K4me3_PBS_WT_rep2_filterpair_sorted.bam O=H3K4me3_PBS_WT_rep2_filtered_sorted.bam M=H3K4me3_PBS_WT_rep2_duplication_metrics.txt
```

__12  Generate the index__
```
samtools index H3K4me3_PBS_WT_rep2_filtered_sorted.bam
```

__13  Call peaks with MACS3__

--broad argument was added according to the peak shape
```
macs3 callpeak -t H3K4me3_PBS_WT_rep2_filtered_sorted.bam -f BAMPE -g mm --keep-dup all -n H3K4me3_PBS_WT_rep2 --outdir ./H3K4me3_PBS_WT_rep2_peaks
```

__14  Fix the blacklist from the ENCODE to merge overlapped regions__

The input (-i) file of merge function in bedtools must be sorted by chromatin.
```
bedtools sort -i ENCFF547MET.bed | bedtools merge -i - > ENCFF547MET_fixed.bed
```

__15  Convert BAM file to BigWig file, since the input of plotHeatmap must be BigWig file__

BAM files must be sorted and indexed
```
bamCoverage  -p 8 --blackListFileName /exec/wangyy/software/ChIP_seq/RefGenome/mouse/blacklist/ENCFF547MET_fixed.bed -b H3K4me3_PBS_WT_rep2_filtered_sorted.bam -o H3K4me3_PBS_WT_rep2_filtered_sorted.bw --binSize 10 --normalizeUsing RPKM --extendReads
```

__16  Compute matrix, gtf annotation file should be input here__

Reference-point mode was applied to peaks of H3K4me3.

Scale-regions mode was applied to peaks of H3K27me3, in which --regionBodyLength was set to 5000.

Only protein-coding genes were considered.
```
computeMatrix reference-point --referencePoint TSS -S H3K4me3_PBS_WT_rep2_filtered_sorted.bw -R /exec/wangyy/software/ChIP_seq/RefGenome/mouse/vM38_primary_assembly_protein_annotation.gtf --transcriptID gene --transcript_id_designator "gene_id" -p 40 -b 3000 -a 3000 -o H3K4me3_PBS_WT_rep2_gene_matrix.gz 
```


__17  Plot Heatmap and the profile plot__
```
plotHeatmap -m H3K4me3_PBS_WT_rep2_gene_matrix.gz -o H3K4me3_PBS_WT_rep2_heatmap.png --colorList white,blue --legendLocation none
```

