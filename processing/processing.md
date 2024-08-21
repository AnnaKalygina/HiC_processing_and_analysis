# From fasta to pairs
There are many ways to process the HiC data from the .fasta format to the final .bam pairs file both manually and using available pipelines. The ones that I would definitely recommend are the HiC-Pro and Mirny's lab pipeline. However, for some types of data, for example diploid HiC processing or Micro-C, many of these pipelines are either much harder to use or unavailable at all. Therefore, in this tutorial I will first go through step-by-step manual processing that is relevant for any kind of data. Make sure, the required packages are installed beforehand

- samtools
- 

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

### Recording ligation events
We use the parse module of the pairtools pipeline to find ligation junctions in Micro-C (and other proximity ligation) libraries. When a ligation event is identified in the alignment file the pairtools pipeline will record the outer-most (5’) aligned base pair and the strand of each one of the paired reads into .pairsam file (pairsam format captures SAM entries together with the Hi-C pair information). In addition, it will also asign a pair type for each event. e.g. if both reads aligned uniquely to only one region in the genome, the type UU (Unique-Unique) will be assigned to the pair. The following steps are necessary to identify the high quality valid pairs over low quality events (e.g. due to low mapping quality):


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
hicSumMatrices -m ZmEn_1_10k.cool ZmEn_2_10k.cool -o ZmEn_10k.cool
hicSumMatrices -m ZmMC_1_10k.cool ZmMC_2_10k.cool -o ZmMC_10k.cool
```

# HiC-Pro Pipeline 
For a standard Hi-C procedure, HiC-Pro pipeline is probably the most helpful tool to use. It is a pain in the neck to set it up and sometimes it takes hours/days to debug, but when you get used to it, it becomes your best friend. The procedure is absolutely the same as described above, except the fact that it takes ~30 minutes to set it up and then it runs from .fasta to .matrix automatically on the dev node. Here I will outline a common procedure to set up the pipeline and major problems I have encountered. 

### Downloading HiC-Pro

### Setting up working directory

### Setting up configuration file

### Running script on interactive node

# Other pipelines
Depending on what experimental procedure you use, some different pipelines could be more straightforward or compatible. Here are the alternative tools you might want to consider:
- When preparing library with Arima HiC kit you can either follow their [recommended pipeline](https://github.com/ArimaGenomics/mapping_pipeline/blob/master/Arima_Mapping_UserGuide_A160156_v03.pdf) or use HiC-Pro
- When processing Micro-C / Omni-C library, the HiC-Pro would not work, so it is better to process it manually. Additional instruction could be found [here](https://micro-c.readthedocs.io/en/latest/index.html) and [here](https://omni-c.readthedocs.io/en/latest/).

# Read more
[1] Iterative correction of HiC matrices by Imakaev et al. : 10.1038/nmeth.2148

[2] Documentation on Micro-C processing that was adapted for this tutorial: https://micro-c.readthedocs.io/en/latest/index.html




