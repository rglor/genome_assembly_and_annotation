
Step 1: Preliminary QC & Quality Trimming
======
Preliminary QC of your sequences can be completed by applying [`fastqc`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to your fastq sequence files. [`fastqc`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a popular tool that "aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines." Carefully inspect the `html` output from [`fastqc`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) because this is the main way that you are going to make sure that nothing is seriously wrong with your data before you delve into the time consuming series of analyses discussed below.
```
fastqc -k 6 *
```

If your sequences look OK after preliminary QC, its time to get your sequences ready for downstream analyses by trimming adaptors and eliminating low quality sequences and basecalls. Here we do this operation on the reads generated from the short insert library. We are going to do this by using the function `cutadapt` to (1) trim Illumina adapter sequences, (2) discard reads <75 bp in length and (3) perform gentle trimming of low quality basecalls. This process should take around 12 hours to complete for a raw sequence file containing around 500 million 100-150 bp reads.
```
#PBS -N cutadapt_short.sh
#PBS -l nodes=1:ppn=1:avx,mem=16000m,walltime=48:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/Short_Insert
#PBS -j oe
#PBS -o cutadapterror_short

fastqc -k 6 /scratch/a499a400/anolis/deliv_Glor061715_raw/raw_fastq/Anolis_Genome_R1.fastq
fastqc -k 6 /scratch/a499a400/anolis/deliv_Glor061715_raw/raw_fastq/Anolis_Genome_R2.fastq
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 5 -m 25 -o Short_trimmed_R1.fastq.gz -p Short_trimmed_R2.fastq.gz /scratch/a499a400/anolis/deliv_Glor061715_raw/raw_fastq/Anolis_Genome_R1.fastq /scratch/a499a400/anolis/deliv_Glor061715_raw/raw_fastq/Anolis_Genome_R2.fastq > cutadapt.log
```
After trimming is complete, use `fastqc` on the resulting files to check that adapters have been trimmed and that the newly generated `fastq` files look good. 

Step 2: Assembly Free Estimates of Genome Coverage, Size and Heterozygosity
======
We will use an assembly-free k-mer counting method to estimate genome size, genome coverage and heterozygosity. Such k-mer counting approaches are widely used in genome studies. Although numerous applications for interpreting k-mer counts have been introduced (e.g., [GenomeScope](http://qb.cshl.edu/genomescope/), [estimate Genomesize.pl](http://josephryan.github.io/estimate_genome_size.pl/)), most rely on the rapid k-mer counting program [jellyfish](http://www.genome.umd.edu/jellyfish.html) to generate k-mer counts. Selecting a k-mer size for these analyses is non-trivial. Advice on selecting a k-mer size for various types of applications can be found in a few places, including the FAQ for the error-correction function [Quake](http://www.cbcb.umd.edu/software/quake/faq.html) and section 2 of the supplemental material for the program [KAT](https://doi.org/10.1093/bioinformatics/btw663). One approach that is widely-cited in recent genome papers follows the authors of the panda genome paper in relying on analyses of 17-mers ([Ruiqiang et al. 2010](http://dx.doi.org/10.1038/nature08696)), but other studies suggest that large k-mer sizes are appropriate for large complex genomes. However, large k-mer sizes can lead to longer runs and more intense memory demands. Running Jellyfish shouldn't take more than 12 hours when conducted on our short insert library including 183M paired reads.

```
#PBS -N jellyfish_short
#PBS -q bigm
#PBS -l nodes=1:ppn=16:avx,mem=20000m,walltime=12:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -j oe
#PBS -d /scratch/glor_lab/rich/distichus_genome/Jellyfish
#PBS -o jellyfish_short_error

work_dir=$(mktemp -d)
cat /scratch/glor_lab/rich/distichus_genome/Short_Insert/Short_trimmed_R1.fastq /scratch/glor_lab/rich/distichus_genome/Short_Insert/Short_trimmed_R2.fastq > $work_dir/Short_trimmed_merged.fastq
jellyfish count -m 17 -o fastq.counts -C $work_dir/Short_trimmed_merged.fastq -s 10000000000 -U 500 -t 16
rm $work_dir/Short_trimmed_merged.fastq
mv $work_dir/* /scratch/glor_lab/rich/distichus_genome/Jellyfish/
rm -rf $work_dir
```
We can then make a histogram from this data as follows.
```
jellyplot.pl fastq.counts > fastq.counts.histo
```
The resulting histogram file can then be uploaded and evaluated by [GenomeScope](http://qb.cshl.edu/genomescope/), an online tool capable of providing assembly free estimates of genome size, coverage, heterozygosity and other statistics.


Step 3: Correct Sequencing Errors
======
Although Illumina is generally regarded as a relatively error-free sequencing method, your sequences are still expected to include errors, most of which will be substitution errors. One way to fix these types of sequencing erros is to use a k-mer spectrum error correction framework such as EULER or Quake ([Kelley et al. 2010](http://genomebiology.com/2010/11/11/R116)). The basic idea of these approaches involves identification of particularly low frequency k-mers that are likely the result of sequencing errors. We have already seen these likely erroneous k-mers as the large peak near the Y-axis in the plots of k-mer coverage generated by GenomeScope at the end of Step 2. The [publication introducing Quake](http://genomebiology.com/2010/11/11/R116) has a nice introduction to the underlying algorithms. More recently Trowel has been introduced and has the advantage of incorporating information on quality scores.

Step 3a: Quake
-----
We use [Quake](http://www.cbcb.umd.edu/software/quake/index.html) here. Quake uses Jellyfish for k-mer counting. Running Quake is fairly straightforward, with simple instructions available via the [program's online manual](http://www.cbcb.umd.edu/software/quake/manual.html). If you are uncertain of what k-mer size to use, the [Quake FAQ](http://www.cbcb.umd.edu/software/quake/faq.html) suggests that you can determine an appropriate k-mer based on an estimate of genome size, where `optimal k-mer = ln(200 * Genome size in bases)/ln(4)`. We use a k-mer of 19 for Anolis distichus because we expect that anole genome is somewhere in the 1.5 Gb range based on our assembly-free analyses in GenomeScan. The `-q 33` tag denotes the manner in which your sequence quality scores are recorded and is correct for anole data. If you need to check for you data, you can use the simple solution offered during an [online discussion](https://www.biostars.org/p/63225/). Quake can't unzip zipped sequence files on the fly. For the second step of the Quake process 
```
#PBS -N quake_short
#PBS -q bigm -l nodes=1:ppn=24:avx,mem=512000,walltime=12:00:00,file=200gb
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/Quake
#PBS -j oe
#PBS -o quake_short_error

work_dir=$(mktemp -d)
cat /scratch/a499a400/anolis/contam_filtered/decontam_short_1.fastqc /scratch/a499a400/anolis/contam_filtered/decontam_short_2.fastqc | count-qmers -k 27 -q 33 > $work_dir/quake_27_counts
cov_model.py quake_27_counts
#cat /scratch/glor_lab/rich/distichus_genome/Short_Insert/Short_trimmed_R1.fastq /scratch/glor_lab/rich/distichus_genome/Short_Insert/Short_trimmed_R2.fastq | correct -k 27 -c 2 -m quake_19_counts -p 24 > $work_dir/
mv $work_dir/* /scratch/glor_lab/rich/distichus_genome/Quake
rm -rf $work_dir
```

Step 3b: Trowel
-----
For more on Trowel, see the publication responsible for this tool ([Lim et al. 2014](http://dx.doi.org/10.1093/bioinformatics/btu513)).
```
#PBS -N trowel_short
#PBS -q bigm -l nodes=1:ppn=24:avx,mem=512000m,walltime=148:00:00,file=300gb
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -j oe
#PBS -d /scratch/glor_lab/rich/distichus_genome/Trowel
#PBS -o trowel_short_error

work_dir=$(mktemp -d)
mkdir $work_dir
cp /scratch/glor_lab/rich/distichus_genome/Short_Insert/Short_trimmed* $work_dir
gunzip -c $work_dir/Short_trimmed_R*
trowel -k 17 -t 24 -f files_short.txt
gzip $work_dir/*fastq
rm $work_dir/Short_trimmed_R1.fastq.gz
rm $work_dir/Short_trimmed_R2.fastq.gz
mv $work_dir/* /scratch/glor_lab/rich/distichus_genome/Trowel
rm -rf $work_dir
```
Step 4:*De novo* Assembly of Small Insert Library with DISCOVAR
======
```
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=1450G
#SBATCH --partition=large-shared
#SBATCH -A kan110
#SBATCH --no-requeue
#SBATCH -t 2880 

export MALLOC_PER_THREAD=1
/home/aalexand/bin/discovardenovo-52488/bin/DiscovarDeNovo READS=/oasis/scratch/comet/aalexand/temp_project/combined_Anolis_Genome_R1.fastq,/oasis/scratch/comet/aalexand/temp_project/combined_Anolis_Genome_R2.fastq OUT_DIR=/oasis/scratch/comet/aalexand/temp_project/17Feb2016 NUM_THREADS=15 MAX_MEM_GB=1000 MEMORY_CHECK=True
```

Step 4: Mapping Reads to Genome
------
Mapping short reads back to the final genome is a basic way of assessing overall genome coverage. Here we map the short reads from the short insert library to the final Dovetail assembly.
```
#PBS -N mapping_dovetail_short
#PBS -q default -l nodes=1:ppn=24:avx,mem=50000m,walltime=168:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome
#PBS -j oe
#PBS -o mapping_dovetail_short_error

bowtie2-build /scratch/a499a400/anolis/dovetail/trunk_anole_19Jun2016_xkeD9.fasta Dovetail_bowtie_DB.fasta
bowtie2 --local --no-unal -x Dovetail_bowtie_DB.fasta -q -1 /scratch/glor_lab/rich/distichus_genome/Short_Insert/Short_trimmed_R1.fastq.gz -2 /scratch/glor_lab/rich/distichus_genome/Short_Insert/Short_trimmed_R2.fastq.gz | samtools view -Sb - | samtools sort -no - - > bowtie2.dovetail.nameSorted.bam
SAM_nameSorted_to_uniq_count_stats.pl bowtie2.dovetail.nameSorted.bam > bowtie_dovetail_mapping_summary

```
To get a summary of the mapping results, we can use a perl script from `Trinity`.
```
SAM_nameSorted_to_uniq_count_stats.pl bowtie2.nameSorted.bam

#read_type      count   pct
proper_pairs    30446172        94.39
improper_pairs  928085  2.88
left_only       699683  2.17
right_only      181046  0.56

Total aligned rnaseq fragments: 32254986
```

Step 5: Assess Genome Completeness by Identifying Coverage of Expected Proteins
======
One important benchmark for genome quality involves asking how many of the proteins that are expected to occur in our organism based on prior genomic studies are present in our assembled genome. 

Step 5a: Benchmarking Universal Single-Copy Orthologs (BUSCO)
------
"BUSCO v2 provides quantitative measures for the assessment of genome assembly, gene set, and transcriptome completeness, based on evolutionarily-informed expectations of gene content from near-universal single-copy orthologs selected from OrthoDB v9." For more about BUSCO visit the [project's website](http://busco.ezlab.org/) or read the paper reporting the method ([Simão et al. 2015](http://dx.doi.org/10.1093/bioinformatics/btv351)).
```
#PBS -N busco_dovetail.sh
#PBS -l nodes=1:ppn=24:avx,mem=64000m,walltime=72:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/
#PBS -j oe
#PBS -o busco_dovetail_error

BUSCO.py -o busco_dovetail -i /scratch/a499a400/anolis/dovetail/trunk_anole_19Jun2016_xkeD9.fasta -l /scratch/glor_lab/rich/distichus_genome_RNAseq/vertebrata_odb9 -m geno
```
After running this operation, you should find a simple summary of your data (`short_summary_busco_dewlap.txt`) in a folder called `run_busco_name` that will tell you how many complete BUSCOs were recovered in your dataset.
```
       XX
```

Step 5b: Assess Full-length Transcripts Relative to Reference Via BLAST+
------
We can also ask how many of the proteins in the Anolis carolinensis proteome are present in our genome.
```
#PBS -N blast_anole_dovetail.sh
#PBS -l nodes=1:ppn=24:avx,mem=200000m,walltime=140:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/Anole_BLAST
#PBS -j oe
#PBS -o blast_anole_dovetail_error

makeblastdb -in ASU_Acar_v2.1.prot.fa -dbtype prot
blastx -query /scratch/a499a400/anolis/dovetail/trunk_anole_19Jun2016_xkeD9.fasta -db ASU_Acar_v2.1.prot.fa -out blastx.outfmt6 -evalue 1e-20 -num_threads 24 -max_target_seqs 1 -outfmt 6
analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 /scratch/a499a400/anolis/dovetail/trunk_anole_19Jun2016_xkeD9.fasta ASU_Acar_v2.1.prot.fa > analyze_blastPlus_topHit_coverage.log
/tools/cluster/6.2/trinityrnaseq/2.2.0/util/misc/blast_outfmt6_group_segments.pl blastx.outfmt6 /scratch/a499a400/anolis/dovetail/trunk_anole_19Jun2016_xkeD9.fasta ASU_Acar_v2.1.prot.fa > blast.outfmt6.grouped
/tools/cluster/6.2/trinityrnaseq/2.2.0/util/misc/blast_outfmt6_group_segments.tophit_coverage.pl blast.outfmt6.grouped > analyze_groupsegments_topHit_coverage.log

```
Step 6: Assess Repetitive Content of Genome
======
We are going to use the function `RepeatMasker` to generate a basic assessment of repetitive content of our genome. We will use the basic function `RepeatMasker` to identify repetitive content from the RepeatMasker database. We will then use the function `RepeatModeler` to generate a library of repetitive elements for our species of interest. Similar processes are also integrated in Maker2.

Step 6a: RepeatMasker
------
```
#PBS -N repeatmasker_dovetail.sh
#PBS -q default -l nodes=1:ppn=24:avx,mem=64000m,walltime=148:00:00,file=200gb
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/RepeatMasker
#PBS -j oe
#PBS -o repeatmasker_dovetail_error

work_dir=$(mktemp -d)
mkdir $work_dir
cp /scratch/a499a400/anolis/dovetail/trunk_anole_19Jun2016_xkeD9.fasta $work_dir
RepeatMasker -e ncbi -species vertebrates $work_dir/trunk_anole_19Jun2016_xkeD9.fasta >> RepeatMasker_Dovetail.log
mv $work_dir/ /scratch/glor_lab/rich/distichus_genome/RepeatMasker
rm -rf $work_dir
```

Step 6b: RepeatModeler
------
An alternative to the homology-based assessment of repetitive content implemented above is to infer repetitive content from your genome using the program [RepeatModeler](http://www.repeatmasker.org/RepeatModeler.html), which uses two de novo repeat finding functions -- RECON and Repeatscout -- to "build, refine and classify consensus models of putative interspersed repeats." Once RepeatModeler has generated models for repeats in your genome, you can then use RepeatMasker to characterize the abundance of these repeats in your genome and to mask repetitive content. The first step is to run RepeatModeler on your assembled genome.
```
#PBS -N repeatmodeler_dovetail.sh
#PBS -l nodes=1:ppn=24:avx,mem=64000m,walltime=148:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/RepeatModeler
#PBS -j oe
#PBS -o repeatmodeler_dovetail_error

BuildDatabase -name distichus_dovetail -engine ncbi /scratch/a499a400/anolis/dovetail/scrubbed_genome.fasta
RepeatModeler -engine ncbi -pa 24 -database distichus_dovetail
```

Then run RepeatMasker with RepeatModeler output.
```
#PBS -N repeatmasker_modeler_dovetail.sh
#PBS -q long -l nodes=1:ppn=16:avx,mem=64000m,walltime=296:00:00,file=200gb
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/RepeatModeler
#PBS -j oe
#PBS -o repeatmasker_modelerdovetail_error

work_dir=$(mktemp -d)
cp /scratch/a499a400/anolis/dovetail/scrubbed_genome.fasta $work_dir
RepeatMasker -lib /scratch/glor_lab/rich/distichus_genome/RepeatModeler/RM_91293.WedDec210322572016/consensi.fa.classified $work_dir/scrubbed_genome.fasta >> $work_dir/RepeatMasker_Modeler_Dovetail.log
mv $work_dir/* /scratch/glor_lab/rich/distichus_genome/RepeatModeler/
rm -rf $work_dir
```
The output of this process will return both a repeat masked version of your genome as well as a table that summarizes repeat cotent.
```
==================================================
file name: scrubbed_genome.fasta
sequences:        878866
total length: 3184760166 bp  (3157651864 bp excl N/X-runs)
GC level:         41.38 %
bases masked: 1623293304 bp ( 50.97 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
SINEs:           654420    107661251 bp    3.38 %
      ALUs            0            0 bp    0.00 %
      MIRs       162403     21234576 bp    0.67 %

LINEs:           2291376    644222230 bp   20.23 %
      LINE1      201210     86424391 bp    2.71 %
      LINE2      857552    202300976 bp    6.35 %
      L3/CR1     657977    191925128 bp    6.03 %

LTR elements:    261422     80392966 bp    2.52 %
      ERVL        23709      3172739 bp    0.10 %
      ERVL-MaLRs      0            0 bp    0.00 %
      ERV_classI  45930     14408218 bp    0.45 %
      ERV_classII 17921      4979932 bp    0.16 %

DNA elements:    1956439    309567665 bp    9.72 %
     hAT-Charlie 306081     58862267 bp    1.85 %
     TcMar-Tigger342444     65377695 bp    2.05 %

Unclassified:    1823490    421033447 bp   13.22 %

Total interspersed repeats:1562877559 bp   49.07 %


Small RNA:         2350       456486 bp    0.01 %

Satellites:       18697      4300015 bp    0.14 %
Simple repeats:  932222     50330769 bp    1.58 %
Low complexity:   90897      9239884 bp    0.29 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element


The query species was assumed to be homo
RepeatMasker version open-4.0.5 , default mode

run with rmblastn version 2.2.27+
The query was compared to classified sequences in ".../consensi.fa.classified"
RepBase Update 20140131, RM database version 20140131
```

Step 7: *Ab initio* Gene Prediction With AUGUSTUS
======
We are going to run a preliminary ab initio assembly using [AUGUSTUS](http://dx.doi.org/10.1093/bioinformatics/btg1080).


Step 8: Genome Annotation in Braker1
======
Braker1 is a relatively new platform for *de novo* genome annotation. Braker requires as input two files: (1) a '.bam' file with transcriptome mapping information and (2) a repeat masked '.fasta' file for your assembled genome. We will generate the required '.bam' file using the spliced junction mapping tool tophat. This operation is considerably different from the mapping exercise discussed previously becuase it involves use of program that is able to map RNA reads across spliced portions of the genome. For the '.fasta' genome file, we will use the output from Repeatmasker.

First, we can use the following commands to generate the required '.bam' file.

```
#PBS -N tophat_dovetail_short
#PBS -q bigm -l nodes=1:ppn=24:avx,mem=500000m,walltime=296:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/Cufflinks
#PBS -j oe
#PBS -o tophat_dovetail_short_error

bowtie2-build /scratch/a499a400/anolis/dovetail/scrubbed_genome.fasta Dovetail_bowtie_DB.fasta
tophat -r 0 Dovetail_bowtie_DB.fasta /scratch/glor_lab/rich/distichus_genome_RNAseq/All_Tissues/decontam_RNA_1.fastq.gz /scratch/glor_lab/rich/distichus_genome_RNAseq/All_Tissues/decontam_RNA_2.fastq.gz
```

Before running Braker, you will need to sure that you have access to a number of important dependencies, including the general tools for dealing with BAM and SAM formatted files (BAMtools and SAMtools) and the *de novo* genome annotation tools (AUGUSTUS and Genemark). You will have an opportunity in the input file to tell Braker where these required dependencies are located.

```
#PBS -N braker_allRNA_dovetailScrubbed
#PBS -q bigm -l nodes=1:ppn=24:avx,mem=500000m,walltime=296:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/Braker1
#PBS -j oe
#PB. -o braker_allRNA_dovetailScrubbed_error

braker.pl --genome=/scratch/a499a400/anolis/dovetail/scrubbed_genome.fasta --bam=/scratch/glor_lab/rich/distichus_genome/Cufflinks/tophat_out/accepted_hits.bam --AUGUSTUS_CONFIG_PATH=/tools/cluster/6.2/augustus/3.2.2/config --GENEMARK_PATH=/tools/cluster/6.2/genemark-et/4.32 --BAMTOOLS_PATH=/tools/cluster/6.2/bamtools/2.3.0/bin --SAMTOOLS_PATH=/tools/cluster/6.2/samtools/1.2
```


Give Maker: assembled genome fasta file, transcripts, A. carolinensis proteome, repeat library from distichus, reference repeat library from vertebrates, softmasking
Use BUSCO output as initial parameters for Augustus/Maker. May need to use long flag with busco
Conflict over whether using bowtie output or raw transcripts does better in Maker
Check GMOD/Maker Wiki
Augustus does not require iterative training
SNAP may need to run multiple times
Run at first with parameters set to 1
Run again with re-trained SNAP file
AED (Annotation edit distance) smaller value = better evidence
