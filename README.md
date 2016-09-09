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

>I trialed SPAdes using my own trimming and then using raw data to see how their internal QC worked, and it gave a comparable result. I would suggest to try both ways, different datasets might behave differently.

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




