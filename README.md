# evoclin
Codes and scripts used to analyse the evolution of P. aeruginosa in a CF patient

**Genomic DNA**

>**Trimming raw sequences, using TRIMMOMATIC**

```bash
trimmomatic-0.30.jar PE -threads N -phred33 sample_read1.fastq.gz \
sample_read2.fastq.gz output1_forward_paired.fq.gz output1_forward_unpaired.fq.gz \
output1_reverse_paired.fq.gz output1_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:20

```

>I intentionally left the trimming to be quite lenient, this is a personal preference.


>**Preliminary genome assembly of trimmed data using SPAdes**

```
spades.py --careful -k 21,71,91,101 -1 trimmed_forward_reads.fastq -2 trimmed_reverse_reads.fastq -o output_dir

```

>I trialled SPAdes using my own trimming and then using raw data to see how their internal QC worked, and it gave a comparable result. I would suggest to try both ways, different datasets might behave differently.

>**Validation of contigs using the A5 assembly pipeline**

```bash
a5_pipeline.pl --threads 20 trimmed_forward_reads.fastq trimmed_reverse_reads.fastq output_contigs

```

>Contigs output from both SPAdes and a5 were aligned to see if any contigs could be extended to give a cleaner alignment. NOTE* For this dataset a5 worked slightly better than SPAdes, however SPAdes was considerably quicker, and offers more customisability. I have a personal preference to use more than one program especially for assemblies.

>**Construction of 'core genome' of clinical isolates and reference strains using HARVEST**

```bash
parsnp -r draft_assembly_reference -d ~/dir/containing/all/assemblies/ -p 20 -o parsnp_output

```

>For parsnp and HARVEST, the draft assembly of S2239(16) was used as a reference for all steps. This is due to it being the oldest isolate obtained from the patient.

>Comparisons were visualised using gingr and a tree was exported using gingr that was visualised using MEGA 7.

**Variant analysis of patient isolates**

>**Mapping trimmed sequences to the S2239(16) draft assembly sequence**

```bash
bowtie2-build S2239_draft.fa outputname

samtools faidx reference.fa

bowtie2 -p N -t -x S2239_draft.fa -1 output1_forward_paired.fq.gz \
 -2 output1_reverse_paired.fq.gz -S output1.sam
``` 

>This is simply mapping the trimmed reads to the draft genome, using quite default parameters. There is quite a lot of changes to this that can be made to be more sensitive etc.

>**Processing mapped reads**

```bash
samtools view -b -S -o output1.bam output1.sam

samtools sort output1.bam > output1.sorted.bam

samtools index output1.sorted.bam

```
>**Calling variants using Freebayes**

```bash
freebayes.py -f S2239_draft.fa -p 1 output1.sorted.bam > output1.raw.vcf

```

>**Processing variants**

>***I personally found manually processing variants based on quality, coverage, strand bias etc. to be more effective than filtering the vcf files with tools like bcffilter or vcfutils.***

```
filtering based on
- QUAL >100
- DP (coverage) >20
- QR (phred score for variant) > 100
- AO (alternate allele count) < 10

```

**Alternate analysis of variants using breseq**

```bash
breseq -j 20 -r ~/path/to/S2239_draft.fasta -o sample1_out ~/path/to/reads/sample1_R1.fastq.gz \
~/path/to/reads/sample1_R2.fastq.gz

```

>Breseq (in this case) gave a comparable number of mutations, although there is no position information in respect to a reference genome, just scaffold number and base change. This makes associating gene specific changes more difficult. An alternative is to order the contigs in reference to PAO1 (using MAUVE), then perform an annotation on the ordered draft assembly using RAST or PROKKA. This will give approx position and gene information. To do these extra steps -

```bash
java -Xmx50000m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output ordered_S2239 -ref PAO1.gbk -draft S2239_draft.fasta

```

>This output fasta file can then be uploaded to RAST (http://rast.nmpdr.org/) or annotated using PROKKA. For PROKKA i created a P. aeruginosa specific database containing PA7, PA14, LESB58 and PAO1 for annotation. This lead to a more specific annotation of the draft S2239 assembly, and enabled more thorough analysis.

```bash
prokka --genus my_PSEUDOMONAS_db --outdir output_S2239_PROKKA_annotation --evalue 0.001

```

>After this, the annotated output can be put through breseq (as detailed above)

**RNAseq analysis**

>RNAseq analysis was done using bowtie2 and cuffdiff. There may be some contention amongst researchers as to the best tools to analyse RNAseq from bacteria. From experience however, there is not a whole lot of difference between methods, if something is 10x upregulated and significant with one method, it will still be significant with another method.

>RNAseq data was mapped using bowtie2 (similar as above) and then processed using SAMtools (as above) and then finally analysed for differential expression using cuffdiff.

```bash
cuffdiff -p 20 PAO1_reference.gtf sample1.1.sorted.bam,sample2.1.sorted.bam,sample3.1.sorted.bam \ sample1.2.sorted.bam,sample2.2.sorted.bam,sample3.2.sorted.bam -o outputRNA

```




