This is the analysis pipeline of RNA-seq applied in ME2 project.

In this tutorial, only one sample is taken as an example. The program actually ran in batches.

<br>

__1  Quality control of the fastq files__
```
fastqc OE_Me2_PBS_1_R1.fastq.gz OE_Me2_PBS_1_R2.fastq.gz
```
<br>

__2  Summarize the quality control reports__
```
multiqc .
```
<br>

__3  Build index of reference genome with Hisat2__
```
hisat2-build /path/to/genome/GRCm39_primary_assembly_genome.fa GRCm39
```
<br>

__4  Align the fastq file to the reference genome__
```
hisat2 --no-unal --no-mixed --no-discordant -p 40 -t -x /path/to/index/hisat2_index/GRCm39 -1 OE_Me2_PBS_1_R1.fastq.gz -2 OE_Me2_PBS_1_R2.fastq.gz -S OE_Me2_PBS_1.sam
```
<br>

__5  Convert SAM file to BAM file__
```
samtools view -bS --threads 40 -o OE_Me2_PBS_1.bam OE_Me2_PBS_1.sam
```
<br>

__6 Quantify the BAM file__
Not all annotation should be used. According to the points of view from 10X Genomics, we only preserve following annotation. 
```
gawk -F'\t' 'BEGIN{
	OFS="\t";
	split("IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_LV_gene IG_V_gene IG_V_pseudogene lncRNA protein_coding TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene TR_V_pseudogene",a," ");
	for(i in a) keep[a[i]]=1
}
/^#/ {print; next}
match($9,/gene_type "([^"]+)"/,m) && (m[1] in keep) {print}
' vM38_primary_assembly_annotation.gtf > vM38_primary_assembly_filtered_annotation.gtf
```

After filtering, remember to check the modification.
```
gawk -F'\t' '$3=="gene" {
    match($9, /gene_type "([^"]+)"/, arr)
    if (arr[1]!="") print arr[1]
}' vM38_primary_assembly_protein_annotation.gtf | sort | uniq
```

BAM files do not have to be sorted or indexed.
```
featureCounts -a /path/to/annotation/vM38_primary_assembly_filtered_annotation.gtf -o /output/dir -T 16 --extraAttributes gene_type,gene_name -p --countReadPairs /path/to/alignment/*bam
```
<br>

Following step were conducted in R.