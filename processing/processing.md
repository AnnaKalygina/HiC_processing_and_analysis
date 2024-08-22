# From fasta to pairs
There are many ways to process the HiC data from the .fasta format to the final .bam pairs file both manually and using available pipelines. The ones that I would definitely recommend are the HiC-Pro and Mirny's lab pipeline. However, for some types of data, for example diploid HiC processing or Micro-C, many of these pipelines are either much harder to use or unavailable at all. Therefore, in this tutorial I will first go through step-by-step manual processing that is relevant for any kind of data. Make sure, the required packages are installed beforehand:

- samtools
- pairtools
- cooltools
- HiC-Pro

### Pre-processing
Before we start mapping the reads to the reference genome, we need to prepare several files that will be used alongside. 
1. The reference genome sequence in a .fasta format. Could be downloaded from the UCSC genome browser. The .fasta files for hg38, CHM13 and hg002 are available in this repositiry in the [data_example folder](HiC_processing_and_analysis/data_example/ref/).

2. Index file for the reference genome, only the main chromosomes should be used. The way to create an index file:
   ```
   samtools faidx <reference_genome.fasta>
   ```
   Seven files would be created and output in the target directory.
3. Genome file wiht chromnames and chromsizes. We can generate genome file using one of the index files:
   ```
   cut -f1,2 <reference_genome.fasta.fai> > <reference.genome>
   ```
4. bwa index file. As we are going to be using the Burrowss-Wheeler aligner, we need to create another index file compatible with bwa:

   ```
   bwa index <reference_genome.fasta>
   ```
### Read mapping 

<img width="697" alt="Screenshot 2024-08-21 at 22 18 54" src="https://github.com/user-attachments/assets/b27747e3-0633-4c8d-930b-53db0767b684">

After we prepared requred files beforehand we need to align the HiC reads library to a reference fasta file. If you download the data from someones experiments, sometimes they recommend using one algorithm over another. In our case, I would recommend using BWA-MEM alignment algorithm by default. In this case both reads are mapped together. For deep coverage data I would recommend requiesting a node for ~32 CPUs. Run the following command: 

```
bwa mem -5SP -T0 -t<n_cpus_requested> <ref.fasta> <HiC_R1.fastq> <HiC_R2.fastq> -o <aligned.sam>
```
Alternatively, if you have libraries from several papers, the multiple pairs could be aligned together:

```
bwa mem -5SP -T0 -t16 hg38.fasta <(zcat file1.R1.fastq.gz file2.R1.fastq.gz file3.R1.fastq.gz) <(zcat file1.R2.fastq.gz file2.R2.fastq.gz file3.R2.fastq.gz) -o aligned.sam
```

You can find more information on how to choose alignment parameters in [bwa documentation](https://bio-bwa.sourceforge.net/bwa.shtml). 

*Note:* if you are planning to be aligning diploid library on diploid genome, it is better to align reads separately. 

This step yields a single .sam format file with alignments.

### Recording ligation events
Now we need to filter aligned reads to .pairs file, which will consist of only valid biological interactions. For this purpose, I use a pairtools command `parse`. The great overview on how to parse pairs is in [this documentation](https://pairtools.readthedocs.io/en/latest/parsing.html).

A couple of terms to now before setting command parameters:

**Ligation junction**
Some kits, like Arima kit, introduce linkers to the ends if interacting pairs. Therefore, when sequenced, some reads will span ligation junctions introduced during experimental procedure. When these 'chimeric' single-end reads are mapped to the reference genome, and both 5' and 3' ends align to the sequence with a high mapping score, the 3' end portion that comes from the junction must be filtered out. Therefore, when a ligation event is identified in the alignment file the pairtools pipeline will record the outer-most (5’) aligned base pair.



**Walks**
It could happen, that during experimental procedure, more than 2 sequences get ligated together, yielding more than 2 hihg-quality alignments from only 2 reads, thwy are called walks. The most basic way to handle such walks is to disregard the middle portion with --walks-policy 5unique, however, if you want to save all high-quality mapping and consider them as valid combinations of pairs, the different walk policy could be used. 

<img width="736" alt="Screenshot 2024-08-21 at 22 25 06" src="https://github.com/user-attachments/assets/bdf60395-ffaa-4e1d-b4e4-f2077b6071d3">


**Alignment gaps**
As opposed to walks, some portions of reads could align only partially. If a part of a read doesn't map well to the reference genome, we call it a gap. Such gaps could be considered as an accidental insertion or a technical artifact, thus we assume that even with this missing part our reads were formed by one ligation event and therefore a pair must be reported. To set how big gaps the command must tolerate before reporting a pair we set --max-inter-align-gap. Traditionally this value is set to 30bp.

<img width="726" alt="Screenshot 2024-08-21 at 22 26 15" src="https://github.com/user-attachments/assets/87cc2acf-b6a3-4947-b177-fa179374ea44">


For optimal results run the following command:
```
pairtools parse --min-mapq 40 --walks-policy 5unique \
--max-inter-align-gap 30 --nproc-in <cores>\
--nproc-out <cores> --chroms-path <ref.genome> <aligned.sam> > <parsed.pairsam>

```
As a result, this step will record classification of each pair from the .sam entry to the .pairsam file. So, if both reads aligned uniquely to only one region in the genome, the type UU (Unique-Unique) will be assigned to the pair. Such high quality valid pairs are retained. 

### Removing PCR duplicates
Experimental protocols can create PCR duplicates that need to be filtered. First, sort the .pairsam file:

```
pairtools sort --nproc <cores> --tmpdir=<path/to/tmpdir> <parsed.pairsam> > <sorted.pairsam>
```
Next, remove duplicates:

```
pairtools dedup --nproc-in <cores> --nproc-out <cores> --mark-dups --output-stats <stats.txt> \
--output <dedup.pairsam> <sorted.pairsam>
```
pairtools dedup detects molecules that could be formed via PCR duplication and tags them as “DD” pair type. These pairs should be excluded from downstream analysis. Use the pairtools dedup command with the –output-stats option to save the dup stats into a text file.

*Note:* If you are alinging a diploid library, this step should be ignored, as every single read has a copy on both parental chromosomes. 

### Generating .pairs file
The pairtools split command is used to split the final .pairsam into two files: .sam (or .bam) and .pairs (.pairsam has two extra columns containing the alignments from which the Omni-C pair was extracted, these two columns are not included in .pairs files). 

```
pairtools split --nproc-in <cores> --nproc-out <cores> --output-pairs <mapped.pairs> \
--output-sam <unsorted.bam> <dedup.pairsam>
```
The .pairs file can be used for generating contact matrix.


# From pairs to matrices

### Convert .pairs (matrix) to .cool and .mcool
As we want to generate matrices for different resolutions I would recommend running several script instead of doing everything manually. The resulting pairs are processed using cooler cload, and then converted to npz using cool2npy.py script that was adopted from C.origami:

``` bash
#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --partition=normal,altemose,owners
#SBATCH --time=48:00:00
#SBATCH --chdir=/scratch/users/kalyanna/GM12878_microc/
#SBATCH --export=all
#SBATCH --requeue
#SBATCH --output=/scratch/users/kalyanna/GM12878_microc/sbatch_logs/download_%j.out

source activate hic

resolutions=(250 500 1000 2000 4000 8000 16000 32000 64000)
genome="/scratch/users/kalyanna/ref/hg38.genome"
pairs="/scratch/users/kalyanna/GM12878_micro_c/mapped_pairs_MicroC_800M.pairs.gz"
cools="/scratch/users/kalyanna/GM12878_micro_c/cools"
npz="/scratch/users/kalyanna/GM12878_micro_c/npz"

for res in "${resolutions[@]}"; do
    echo "Processing at $res resolution"
    cooler cload pairix -p 16 ${genome}:${res} $pairs $cools/GM12878_microc_${res}.cool
    

    echo "Balancing cool"
    cooler balance $cools/GM12878_microc_${res}.cool


    echo "Converting to npz"
    mkdir -p $npz/npz_${res}
    python /scratch/users/kalyanna/C.Origami/src/corigami/preprocessing/cool2npy.py \
    $cools/GM12878_microc_${res}.cool \
    $npz/npz_${res} \
    -r $res \
    --no-balance
done

```

If you want to conver .pairs into any other HiC formats, I recommend using [hicExplorer](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicConvertFormat.html) function hicConvert.

### Matrix correction
As you could notice in the script from matrix convertation, we use command cooler balance, which normalizes matrices. Matrix correction is necessary to remove biases like gc content or mappability.
In order to correct a matrix, it is assumed that if no biases were affecting the experiment, each bin should have equal “visibility” of contacts. This translates to an intuitive solution to matrix correction: transforming the matrix in such a way that the total number of contacts of every row and every column is the same. Such a procedure is called “matrix balancing”, and many algorithms have been described to achieve this for applications outside HiC data analysis. 
For HiC data, the most common ones are called Knight-Ruiz (KR), and Iterative Correction (ICE). Have a look at an amazing breakdown of [normalization methods](https://liorpachter.wordpress.com/2013/11/17/imakaev_explained/).

The cooler package has a funcrtion balance that performs iterative correction as it was developed in Imakaev 2012 [1]. Using this type of correction we filter bad based on MAD max (see explanation in the cooler balance -h). The "balancing weights" produced by cooler are the reciprocal of the "biases" as defined in Imakaev et al, 2012. By default, this function rescales the weights so that the corrected marginals sum to unity.

```
$ cooler balance /path/to/cool/file.cool
```

### Sum samples into one file
A common practice in HiC data is to sum biological replicate matrices in order to increase sequencing depth, and thus matrix resolution. This can be done after checking that the biological replicates are indeed similar. It is advised to also conduct downstream analyses separately on each replicate to assess differences at those levels.
```
hicSumMatrices -m replicate_1.cool replicate_2.cool -o merged_replicates.cool

```

# HiC-Pro Pipeline 
For a standard Hi-C procedure, HiC-Pro pipeline is probably the most helpful tool to use. It is a pain in the neck to set it up and sometimes it takes hours/days to debug, but when you get used to it, it becomes your best friend. The procedure is absolutely the same as described above, except the fact that it takes ~30 minutes to set it up and then it runs from .fasta to .matrix automatically on the dev node. Here I will outline a common procedure to set up the pipeline and major problems I have encountered. 

### Downloading HiC-Pro

### Setting up working directory

### Setting up configuration file

### Running script on interactive node

### HiC-Pro specificities




HiC-Pro is very robust in terms of filtering valid interaction pairs. For it to run it requires a list of possible fragments generated by the restriction enzymes in a mix. Having every possible restriction fragment, the pipeline assigns each aligned read to it. Only read that come from the same pair and span across different restriction fragments are considered to be valid interaction pair generated by the HiC protocol. Such procedure helps to filter out self circle pairs, singletons and multi-hits. Short range interactions within restriction fragment are also discarded. Next each pair is flagged according toits classification, only valid unique-unique mapping pairs make it to the next step. 

# Other pipelines
Depending on what experimental procedure you use, some different pipelines could be more straightforward or compatible. Here are the alternative tools you might want to consider:
- When preparing library with Arima HiC kit you can either follow their [recommended pipeline](https://github.com/ArimaGenomics/mapping_pipeline/blob/master/Arima_Mapping_UserGuide_A160156_v03.pdf) or use HiC-Pro
- When processing Micro-C / Omni-C library, the HiC-Pro would not work, so it is better to process it manually. Additional instruction could be found [here](https://micro-c.readthedocs.io/en/latest/index.html) and [here](https://omni-c.readthedocs.io/en/latest/).

# Read more
[1] Iterative correction of HiC matrices by Imakaev et al. : 10.1038/nmeth.2148

[2] Documentation on Micro-C processing that was adapted for this tutorial: https://micro-c.readthedocs.io/en/latest/index.html




