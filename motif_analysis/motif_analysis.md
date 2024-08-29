# Ananlysing distribution of CTCT motifs in CHM13 and its impact on chromatin interactions

## Scanning for motifs in CHM13
Additionaly installed [HOMER software](http://homer.ucsd.edu/homer/index.html) is required for the following analysis.

The CTCF motifs were downloaded from JASPAR with the following ID: [MA0139.1](https://jaspar.elixir.no/matrix/MA0139.1). The file was converted into HOMER-compatible .matrix format. As it is instructed in the [HOMER documentation](http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html), the file must be tab-delimited and must contain a logg-odds threshold. The CTCF.motif file is avalable [in this repository](CTCF_homer_compatible.motif) or on cluster at `/oak/stanford/groups/altemose/kalyanna/ctcf_in_chm13/CTCF_homer_compatible.motif`.

Next, I pass the .motif file to scan the CHM13 genome for specified CTCF motifs using the `scanMotifGenomeWide.pl` from HOMER package:

```jsx
scanMotifGenomeWide.pl\
CTCF_homer_compatible.motif\
chm13v2.0.fa \
 -bed -p 10 > CTCF_in_CHM13_19bp_t3.bed
```
The threshold could be scpecified for more or less robust scanning.


## Analysing distribution of CTCF motifs in active $\alpha$ satellites 
The active $\alpha$ satellites are sampled from [Cen/Sat v2.1 annotation](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed) and are available at `/oak/stanford/groups/altemose/kalyanna/ctcf_in_chm13/	active_hor_in_CHM13.bed` , all CTCF motifs are available at `/oak/stanford/groups/altemose/kalyanna/ctcf_in_chm13/	CTCF_in_CHM13_19bp_t5.bed`. 

The CTCF motifs that occur in active HORs could be found by intersecting these two .bed files using `bedtools`:
```bash
bedtools intersect -a CTCF_in_CHM13_19bp_t5.bed -b active_hor_in_CHM13.bed > CTCF_hor_overlap.bed
```
The intersection file is available at `/oak/stanford/groups/altemose/kalyanna/ctcf_in_chm13/CTCF_hor_overlap.bed`.

Nextt, I analyse the distribution of CTCF motifs across centromeres in the CHM13. We can see that while the overall distribution of CTCF motifs across entire chromosomes appears nearly random and shows a positive correlation with chromosome size (graph 1), the distribution within $\alpha$ satellites is significantly uneven (graph 2).

**GRAPH 1:**
<img width="1061" alt="Screenshot 2024-08-28 at 11 56 17" src="https://github.com/user-attachments/assets/8fec7861-d2cf-4a42-9d35-b6dcaf87aeb5">

In particular, chromosomes 5, 11, 18, 19, and 20 shows a strong dominance of motifs in one orientation. I think it may suggest an expansion of a single motif, though further testing is needed. In contrast, chromosomes 1 and X show almost equal proportions of positively and negatively oriented motifs, which could indicate a more balanced expansion of motifs. (Even though the proportions are equal, the chi-square statistic, that I calculate for both chromosomes, indicates non-uniform distribution, p-value < 0.001 )

**GRAPH 2:**
<img width="1069" alt="Screenshot 2024-08-28 at 11 56 33" src="https://github.com/user-attachments/assets/86dc20e2-7f0c-4645-8036-928bfd577dea">


On chromosome 1, the CTCF motifs are densely packed and often occur in tandem (n_motifs = 1781 , graph 3). 

**GRAPH 3:**
<img width="1100" alt="Screenshot 2024-08-28 at 11 57 56" src="https://github.com/user-attachments/assets/f586569d-e7ce-40f1-aef4-42d8b607ed7f">

Conversely, on chromosome X, the motifs are much sparser (n_motifs = 78 motifs, graph 4). I would think, that the sparse distribution on chromosome X, combined with the presence of multiple convergent CTCF pairs, would suggest a higher potential for interesting 3D chromatin organisation, so it is a very promising candidate for predicting the Hi-C map of this region.

**GRAPH 4:***
<img width="1104" alt="Screenshot 2024-08-28 at 11 58 16" src="https://github.com/user-attachments/assets/205924d0-b53e-4cbc-b660-0934f7f67395">


Below (graph 5) are the satellites from Nick’s annotation (seems valid), below is the distribution of CTCF motifs (blue - “+”, orange - “-”). It seems like the distribution of + and - is somewhat uniform, but CTCF tends to avoid certain regions, like hsat2_1_6(A2, A1, B) and other hsat2/3, and several inactive hor. 

**GRAPH 5:**
<img width="1102" alt="Screenshot 2024-08-28 at 11 58 34" src="https://github.com/user-attachments/assets/aaa53b63-61bd-4b62-ae7f-96af1e3cd046">
