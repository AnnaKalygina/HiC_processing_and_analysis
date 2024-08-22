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
The pairtools split command is used to split the final .pairsam into two files: .sam (or .bam) and .pairs (.pairsam has two extra columns containing the alignments from which the pair was extracted, these two columns are not included in .pairs files). 

```
pairtools split --nproc-in <cores> --nproc-out <cores> --output-pairs <mapped.pairs> \
--output-sam <unsorted.bam> <dedup.pairsam>
```
This way of generating the .pairs file will yield sorted but not indexed. The .pairs file can be used for generating contact matrix. Though it is not covered in this tutorial, the .bam file could be sorted and index to then be used for evaluating the library complexity


# From pairs to matrices

To generate a matrix .pairs file must be separated into bins, each containing ***n*** base pairs. The ***n*** base pairs in the bins, or the binsize, corresponds to the resolution of a matrix. For example, the map of 1kb resolution will have (len_genome // 1000) number of bins, each corresponding to a sequential 1000 base pair region in a genome. 

### Convert .pairs (matrix) to .cool and .mcool
We will use Cooler to generate .cool contact matrices.

***Note:*** Cooler supports two different ways of converting pairs to .cool files. The one is `cooler cload pairix` command, that requires pairix-indexed contact list file as input. To use this method, we first need to index and compress our pairs with the following command:

```
cooler cload pairix -p <cores> <ref.genome>:<bin_size_in_bp> <mapped.pairs.gz> <matrix.cool>
```
This will output compressed mapped.pairs.gz file that could be used with `cooler cload pairix` in the next step. 

However, in case you do not want to index and compress your pairs, there's another way to convert them to the .cool matrix by running:

```
cooler cload pairs --assembly <genome_assembly. e.g. hg38> -c1 <column with the 1st chromosome> -pos1 <column with position of the first read> -c2 <column with the 2nd chromosome> -p2 <coulmn with the second read position> <BINS> <PAIRS> <OUTDIR>
```
Here is the example of this command:
```
cooler cload pairs --assembly hg38 --chrom1 1 --pos1 2 --chrom2 3 --pos2 4 hg38.chrom.sizes:1024 mapped_pairs/concat_mapaq30_pairs.txt.gz K562_combined_1024.cool
```

When you are not sure what to choose, please consult [cooler documentation](https://readthedocs.org/projects/cooler/downloads/pdf/stable/).

**SLURM script**
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

If you want to convert .pairs into any other HiC formats, I recommend using [hicExplorer](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicConvertFormat.html) function hicConvert.

### Matrix correction
As you could notice in the script from matrix convertation, we use command cooler balance, which normalizes matrices. Matrix correction is necessary to remove biases like gc content or mappability.

In order to correct a matrix, it is assumed that if no biases were affecting the experiment, each bin should have equal “visibility” of contacts. This translates to an intuitive solution to matrix correction: transforming the matrix in such a way that the total number of contacts of every row and every column is the same. Such a procedure is called “matrix balancing”, and many algorithms have been described to achieve this for applications outside HiC data analysis. 
For HiC data, the most common ones are called Knight-Ruiz (KR), and Iterative Correction (ICE). Have a look at an amazing breakdown of [normalization methods](https://liorpachter.wordpress.com/2013/11/17/imakaev_explained/).

The cooler package has a function balance that performs iterative correction as it was developed in Imakaev 2012 [1]. Using this type of correction we filter bad based on MAD max (see explanation in the cooler balance -h). The "balancing weights" produced by cooler are the reciprocal of the "biases" as defined in Imakaev et al, 2012. By default, this function rescales the weights so that the corrected contact frequencies sum to unity.

```
$ cooler balance /path/to/cool/file.cool
```
The `hicExplorer` is another package that helps perform matrix correction in a more tailored and elaborated way. 

### Sum samples into one file
A common practice in HiC data is to sum biological replicate matrices in order to increase sequencing depth, and thus matrix resolution. This can be done after checking that the biological replicates are indeed similar. It is advised to also conduct downstream analyses separately on each replicate to assess differences at those levels.
```
hicSumMatrices -m replicate_1.cool replicate_2.cool -o merged_replicates.cool

```

# HiC-Pro Pipeline 
For a standard Hi-C procedure, HiC-Pro pipeline is probably the most helpful tool to use. It is a pain in the neck to set it up and sometimes it takes hours/days to debug, but when you get used to it, it becomes your best friend. The procedure is absolutely the same as described above, except the fact that it takes ~30 minutes to set it up and then it runs from .fasta to .matrix automatically on the dev node. Here I will outline a common procedure to set up the pipeline and major problems I have encountered. 

### Setting up environment
First, access HiC-Pro repository on [GitHub](https://github.com/nservant/HiC-Pro), additionaly they have a more detailed documentation in [pages](https://nservant.github.io/HiC-Pro/).
For working in this repository I recommend creating new environmnet with conda using .yml file attached:
```
conda env create -f hicpro_environment.yml
conda activate hicpro_environment
```
### Downloading HiC-Pro
Now, navigate to a directory where you want your HiC-Pro to be installed. To install HiC-Pro on a cluster in your personal folder, follow these steps:
```
# Clone repository from GitHub
git clone https://github.com/nservant/HiC-Pro.git
cd HiC-Pro

# Edit config-install.txt
nano config-install.txt
```
Change prefix that leads to your main directory, bowtie and samtools path to the packages that you downloaded with the environment:

``` bash
#########################################################################
## Paths and Settings  - Start editing here !
#########################################################################

PREFIX = /scratch/users/<your folder> #this could be a current directory or anywhere you want it to be installed
BOWTIE2_PATH = /home/groups/altemose/<your folder>/miniconda3/envs/hicpro_environment/bin/bowtie2
SAMTOOLS_PATH = /home/groups/altemose/<your folder>/miniconda3/envs/hicpro_environment/bin/samtools
R_PATH =
PYTHON_PATH =
CLUSTER_SYS = SLURM
```
Next, start installation:
```
make configure
make install
```
This command will create a new HiC-Pro folder at the directory that you put in PREFIX.
Now navigate to your $HOME directory and export path in .bashrc:

``` bash
nano .bashrc
```
Add the following line:
``` bash
export PATH=$PATH:"<your folder>"
### <your folder> must be the same directory that you put in PREFIX in config-install.txt
```
Now, try calling HiC-Pro -h. If it works, it works!

### Setting up working directory
Cool, now you have your tool ready, it's time to prepare your working directory. Everything need to be in order, nothing could be rearranged, otherwise, HiC-Pro wouldn't run.

```
# This is your working directory. I would call it by the name of a cell line
mkdir cell_line
cd cell_line


# Create directory for annotation files
mkdir annotation
cd annotation
```
**Annotation folder** must contain the following files:

- Reference genome that reads should be aligned to
- Chromsizes. This file could be downloaded from the UCSC Genome browser or generated in the text editor:

  ```
  chr1    249250621
  chr2    243199373
  chr3    198022430
  chr4    191154276
  chr5    180915260
  chr6    171115067
  chr7    159138663
  chr8    146364022
  chr9    141213431
  chr10   135534747
  (...)
  ```
- Restriction fragments file in a .bed format. Please consult your HiC protocol to find the names of the restriction enzymes that were used in the experiment. For example, [Arima kit](https://arimagenomics.com/faqs/#:~:text=The%20Arima%2DHiC%20chemistry%20uses,for%20your%20genome%20of%20interest.) uses a cocktail of enzymes that digest chromatin at ^GATC and G^ANTC, where N can be any of the 4 genomic bases. There is a specia script available in HiC-Pro package, called [digest_genome.py](https://github.com/nservant/HiC-Pro/tree/master/annotation) that generates the restriction fragments file for you:
  ```
  ## Digest the mm9 genome by HindIII
   HICPRO_PATH/bin/utils/digest_genome.py -r A^AGCTT -o mm9_hindiii.bed mm9.fasta

   ## The same ...
   HICPRO_PATH/bin/utils/digest_genome.py -r hindiii -o mm9_hindiii.bed mm9.fasta

   ## Double digestion, HindIII + DpnII
   HICPRO_PATH/bin/utils/digest_genome.py -r hindiii dpnii -o mm9_hindiii_dpnii.bed mm9.fasta
  ```
- Bowtie indexes for BWA alignment. To generate indexes run the folllowing:
  ```
  bowtie2-build --threads 8 path/to/reference_genome_fasta_file <ref_genome_prefix>
  ```
  This command only takes zipped fasta file as input. It will run for ~20 minutes and will ouput 6 index files with the [refix specified in the command. Make sure that the files do not contain .tmp suffic, otherwise, you will need to delete them and rerun the command with more memory allocation.

**Rawdata folder**
Now you have your annotation folder set up. Now create a `rawdata` folder in the main `cell_line` directory:
```
mkdir rawdata
```
Move all .fasta read files that you need to process into this directory. It is important to note, that when you are dealing with several samples, each sample must have its own folder in the `rawdata` directory. 

### Setting up configuration file
Before you run the main script it is very important to set the configuration file `config-hicpro.txt`, following [these instructions](https://nservant.github.io/HiC-Pro/MANUAL.html#setting-the-configuration-file). Copy `config-hicpro.txt` to your `cell_line` directory and edit it. Make sure, that in the section **ANNOTATION FILES** the REFERENCE_GENOME is set to the <ref_genome_prefix> that comes first in the names of your index files generated by bowtie2-build. 

### Running script on interactive node
To run script on interactive note, do the following:

``` bash
/scratch/users/kalyanna/HiC-Pro_3.1.0/bin/HiC-Pro \
-i /scratch/users/kalyanna/CHM13/by_hic_pro/rawdata/ \
-c /scratch/users/kalyanna/CHM13/by_hic_pro/config-hicpro.txt \
-o /scratch/users/kalyanna/CHM13/by_hic_pro/results/ \ # Set your own output directory
-s mapping -s proc_hic -s quality_checks -s merge_persample -s build_contact_maps -s ice_norm

```

As soon as the script strarts to run, the `results` folder in the `cell_line` directory is created. It will contain rawdata (the files are copied from your rawdat), bowtie_results, hic_results, logs and tmp folders.

### HiC-Pro specificities

Sometimes, you might want to run the commands in sequential mode specifying the steps yourself. Most often you will need to run several samples separately and then merge them. In this case, HiC-Pro will need to have .validPairs files from each samples to be moved to one sample results directory. 
It is going to look like this:

```
~/results_replica1/hic_results/data/merged_samples:

— replica1.validPairs

— replica2.validPairs
```

After the valid pairs are combined we run the following command, that will combine the pairs per sample and start building maps:

```python
/scratch/users/kalyanna/HiC-Pro_3.1.0/bin/HiC-Pro \
-c /scratch/users/kalyanna/CHM13/by_hic_pro/config-hicpro.txt \
-i /scratch/users/kalyanna/CHM13/by_hic_pro/results_replica1/hic_results/data/ \ # This is a directory where you moved all valid pairs from all samples
-o /scratch/users/kalyanna/CHM13/by_hic_pro/results_replica1/ \
-s merge_persample -s build_contact_maps -s ice_norm

```

This command will first output .allValidPairs in hic_results/data/merged_sample directory, and then it will use it to build and normalise matrices, output in hic_results/matrix.

**When not to use HiC-Pro**
HiC-Pro is very robust in terms of filtering valid interaction pairs. For it to run it requires a list of possible fragments generated by the restriction enzymes in a mix. Having every possible restriction fragment, the pipeline assigns each aligned read to it. Only read that come from the same pair and span across different restriction fragments are considered to be valid interaction pair generated by the HiC protocol. Such procedure helps to filter out self circle pairs, singletons and multi-hits. Short range interactions within restriction fragment are also discarded. Next each pair is flagged according toits classification, only valid unique-unique mapping pairs make it to the next step. 

# Other pipelines
Depending on what experimental procedure you use, some different pipelines could be more straightforward or compatible. Here are the alternative tools you might want to consider:
- When preparing library with Arima HiC kit you can either follow their [recommended pipeline](https://github.com/ArimaGenomics/mapping_pipeline/blob/master/Arima_Mapping_UserGuide_A160156_v03.pdf) or use HiC-Pro
- When processing Micro-C / Omni-C library, the HiC-Pro would not work, so it is better to process it manually. Additional instruction could be found [here](https://micro-c.readthedocs.io/en/latest/index.html) and [here](https://omni-c.readthedocs.io/en/latest/).

# Read more
[1] Iterative correction of HiC matrices by Imakaev et al. : 10.1038/nmeth.2148

[2] Documentation on Micro-C processing that was adapted for this tutorial: https://micro-c.readthedocs.io/en/latest/index.html




