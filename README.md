# HiC processing and analysis
In this GitHub repository, I will present a detailed, tutorial-style guide on how to process, visualize, and analyze Hi-C data with various tools. This includes everything from the initial data processing using the HiC-Pro pipeline to more advanced analyses such as calculating insulation scores to identify TAD boundaries, integrating epigenetic data, and exploring the functional implications of chromatin structure on genome regulation.

### How to use this repository?

To get started, please refer to the **Table of Contents** below, which outlines the sections available in this repository. Each section contains an `.md` file with step-by-step instructions for different stages of Hi-C data processing and analysis. If you'd like to follow along with the same data used in this tutorial, detailed instructions on how to download this data from the Sherlock cluster are available [here](data_example/data_example.md). Additionally, several data and annotation files are directly uploaded to this directory, with relevant sections linking to them individually.

# Table of contents
### 1. [Short introduction to HiC](background.md)
### 2. [HiC processing: from .fasta to .pairs](processing)
### 3. [HiC processing: from .pairs to maps](visualising)
### 4. [Map analysis: insulation and domains](metrics)
### 5. [Screening: cell-specific differences and perturbation](screening)
### 6. [Data available for processing and analysis](data_example)

# Miscellaneous tools
In addition to the main Hi-C workflow, this repository includes several tools that are useful for broader bioinformatics analysis:

### 7. [Motif analysis](motif_analysis)
Here, I upload documentation on how I analysed distribution of CTCF motifs across active &/alpha& satellites in CHM13 genome.

### 8. [Predicting HiC data using deep learning models](prediction)
Here, I upload python scripts that allow to [prepare data for training](prediction/tools_for_HiC_prediction) - aggregate problematic regions and filter them out of training dataset. And how to analyse [attention scores](prediction/transformer_attention.md) derived during model training and prediction to deduce the most important pieces of information that the model learns from. 
