### Calculating insulation score
Topological domains (TADs) are defined as genomic neighborhoods of highly interacting chromatin, with relatively more infrequent inter-domain interactions. 
Topological domains are demarcated by boundaries, i.e., genomic regions bound by insulators thus hampering DNA contacts across adjacent domains. 
For each genomic position, in a given resolution (typically 40 kb or less), we define a “boundary score” to quantify the insulation strength of this position. 
The higher the boundary score, the higher the insulation strength and the probability that this region actually acts as a boundary between adjacent domains. 

When we visualise the HiC matrices, we notice patterns of Topologically Associated Domains or (TADs) along the main diagonal. 
These TADs indicate increased frequency of interaction between regions encompassed by the boundaries around TADs. 
In other words, we would say that TADs are the regions of the highest insulation. 
The easy way to analyse the patterns dominating in certain chromatin regions is to calculate the insulation score that would reflect the how enclosed or insulated the regions is compared to its neighbours.
For this purpose, we would stride a diamond-shaped window along the main diagonal and calculate the score of compartmentalisation in each of these windows. 
There are many ways to calculate the score. In Crane et al. 2015, they simply adds up contacts per window, while 




In Lazaris et al., the insulation score is implemented as the ratio of maximum left and right region average intensity and the middle region intensity. 
We also added a pseudocount calculated from chromosome-wide average intensity to prevent division by zero in unmappable regions. 
Given that all the regions contain n interactions, the insulation score can be formulated as follows:

<img width="370" alt="Screenshot 2024-08-18 at 22 02 34" src="https://github.com/user-attachments/assets/e8604365-24e5-4a00-bde1-1c9f4b9ba0d0">



### Read more:
TAD calling: https://vaquerizaslab.github.io/fanc/fanc-executable/fanc-analyse-hic/domains.html
