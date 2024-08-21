# HiC processing and analysis
Tutorial-style documentation to how to process and analyse HiC data with different tools

### Backgorund
Chromatin, the complex of DNA and proteins that form chromosomes, is not just a static structure within the nucleus. It is highly dynamic, and its conformation plays a crucial role in regulating gene expression, replication, and other essential cellular processes. Understanding these conformational changes is key to deciphering the underlying mechanisms of gene regulation and chromosomal organization.

**Types of Chromatin Interactions**
Chromatin interactions can be observed at various scales, from the fine details of nucleosome positioning to large-scale inter-chromosomal contacts. Here are some key types of interactions:

- Promoter-Promoter and Enhancer-Promoter Interactions: These occur at the nucleosome level (~140 bp) and are vital for regulating gene expression. These focal points can be visualized as specific interaction spots in Hi-C maps, reflecting the direct physical interactions between regulatory elements and their target genes.

- CTCF Loops and Cohesin-Mediated Structures: At a slightly larger scale (~1kb - 40kb), chromatin forms loops mediated by the architectural proteins CTCF and cohesin. These loops are essential for organizing the chromatin into topologically associated domains (TADs) and ensuring proper gene regulation.

- Topologically Associated Domains (TADs): TADs represent functional regions of the genome where genes within the same domain interact more frequently with each other than with genes outside the domain. TADs are a fundamental unit of chromatin organization, contributing to the functional compartmentalization of the genome.

- Zygotic Fountains and Polycomb Domains: These interactions, observable at the scale of hundreds of kilobases, play a role in early developmental processes and maintaining gene silencing, respectively. Zygotic fountains are crucial during zygote activation, while Polycomb domains are associated with the regulation of developmental genes.

- A/B Chromosome Compartments: At the largest scale, chromatin is organized into active (A) and inactive (B) compartments, which correspond to large-scale regions of the genome that are either transcriptionally active or repressed. These compartments can be visualized through Hi-C as large, distinct domains on the chromatin interaction maps.

  ![image](https://github.com/user-attachments/assets/6f7411ce-2b56-4868-ae3c-0b7de4baf5b2)


**Capturing Chromatin Conformation: The Hi-C Technique**
To study these interactions, we use the Hi-C (High-throughput chromosome conformation capture) technique. Hi-C is a powerful method that allows us to capture and quantify the frequency of chromatin interactions across the entire genome. The protocol involves the following key steps:

- Crosslinking and Digestion: Cells are treated with a crosslinking agent, typically formaldehyde, which stabilizes interactions between DNA segments that are in close proximity. The chromatin is then digested with a restriction enzyme, cutting the DNA into smaller fragments.

- Ligation: DNA fragments that were in close proximity in the 3D space of the nucleus are ligated together, creating chimeric DNA molecules that reflect the original spatial organization of the chromatin.

- Sequencing and Mapping: The ligated DNA is then sequenced, and the resulting reads are mapped back to the reference genome to identify interacting regions.

- Data Analysis: The final step involves the analysis and visualization of the data, typically resulting in a Hi-C map. These maps provide a comprehensive view of chromatin interactions at various scales, from local interactions within TADs to global chromosome-wide interactions.

**Hi-C Data Processing and Analysis**
In this GitHub repository, I will present a detailed, step-by-step guide on how to process, visualize, and analyze Hi-C data. This includes everything from the initial data processing using the HiC-Pro pipeline to more advanced analyses such as calculating insulation scores to identify TAD boundaries, integrating epigenetic data, and exploring the functional implications of chromatin structure on genome regulation.

# [HiC processing: from .fasta to .pairs](HiC_processing_and_analysis/processing)
