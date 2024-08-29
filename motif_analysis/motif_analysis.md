# Ananlysing distribution of CTCT motifs in CHM13 and its impact on chromatin interactions
### Scanning for motifs in CHM13
Additionaly installed [HOMER software](http://homer.ucsd.edu/homer/index.html) is required for the following analysis.

The CTCF motifs were downloaded from JASPAR with the following ID: [MA0139.1](https://jaspar.elixir.no/matrix/MA0139.1). The file was converted into HOMER-compatible .matrix format. As it is instructed in the [HOMER documentation](http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html), the file must be tab-delimited and must contain a logg-odds threshold. The CTCF.motif file is avalable in this repository or on cluster at `/oak/stanford/groups/altemose/kalyanna/ctcf_in_chm13/CTCF_homer_compatible.motif`.

Next, I pass the .motif file to scan the CHM13 genome for specified CTCF motifs using the `scanMotifGenomeWide.pl` from HOMER package:

```jsx
scanMotifGenomeWide.pl\
CTCF_homer_compatible.motif\
chm13v2.0.fa \
 -bed -p 10 > CTCF_in_CHM13_19bp_t3.bed
```
The threshold could be scpecified for more or less robust scanning.

### Analysing distribution of CTCF motifs in active $\alpha$ satellites 
<img width="1061" alt="Screenshot 2024-08-28 at 11 56 17" src="https://github.com/user-attachments/assets/8fec7861-d2cf-4a42-9d35-b6dcaf87aeb5">
<img width="1069" alt="Screenshot 2024-08-28 at 11 56 33" src="https://github.com/user-attachments/assets/86dc20e2-7f0c-4645-8036-928bfd577dea">

<img width="1100" alt="Screenshot 2024-08-28 at 11 57 56" src="https://github.com/user-attachments/assets/f586569d-e7ce-40f1-aef4-42d8b607ed7f">
<img width="1104" alt="Screenshot 2024-08-28 at 11 58 16" src="https://github.com/user-attachments/assets/205924d0-b53e-4cbc-b660-0934f7f67395">

<img width="1102" alt="Screenshot 2024-08-28 at 11 58 34" src="https://github.com/user-attachments/assets/aaa53b63-61bd-4b62-ae7f-96af1e3cd046">
