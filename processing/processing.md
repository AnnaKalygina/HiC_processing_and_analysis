# From fasta to pairs

# From pairs to matrices

### Convert .pairs (matrix) to .cool and .mcool

### Convert .cool to npz

### Matrix correction
Matrix correction is necessary to remove biases like gc content or mappability.
In order to correct a matrix, it is assumed that if no biases were affecting the experiment, each bin should have equal “visibility” of contacts. This translates to an intuitive solution to matrix correction: transforming the matrix in such a way that the total number of contacts of every row and every column is the same. Such a procedure is called “matrix balancing”, and many algorithms have been described to achieve this for applications outside HiC data analysis. 
For HiC data, the most common ones are called Knight-Ruiz (KR), and Iterative Correction (ICE).https://liorpachter.wordpress.com/2013/11/17/imakaev_explained/

### Sum samples into one file
A common practice in HiC data is to sum biological replicate matrices in order to increase sequencing depth, and thus matrix resolution. This can be done after checking that the biological replicates are indeed similar. It is advised to also conduct downstream analyses separately on each replicate to assess differences at those levels.
```
hicSumMatrices -m ZmEn_1_10k.cool ZmEn_2_10k.cool -o ZmEn_10k.cool
hicSumMatrices -m ZmMC_1_10k.cool ZmMC_2_10k.cool -o ZmMC_10k.cool
```

# HiC-Pro Pipeline 

### Downloading HiC-Pro

### Setting up working directory

### Setting up configuration file

### Running script on interactive node
