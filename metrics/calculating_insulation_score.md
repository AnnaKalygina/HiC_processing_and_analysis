### Calculating insulation score
Topological domains (TADs) are defined as genomic neighborhoods of highly interacting chromatin, with relatively more infrequent inter-domain interactions. 
Topological domains are demarcated by boundaries, i.e., genomic regions bound by insulators thus hampering DNA contacts across adjacent domains. 
For each genomic position, in a given resolution (typically 40 kb or less), we define a “boundary score” to quantify the insulation strength of this position. 
The higher the boundary score, the higher the insulation strength and the probability that this region actually acts as a boundary between adjacent domains. 

When we visualise the HiC matrices, we notice patterns of Topologically Associated Domains or (TADs) along the main diagonal. 
These TADs indicate increased frequency of interaction between regions encompassed by the boundaries around TADs. 
In other words, we would say that TADs are the regions of the highest insulation. 
The easy way to analyse the patterns dominating in certain chromatin regions is to calculate the insulation score that would reflect the how enclosed or insulated the regions is compared to its neighbours.


For this purpose, we would stride a diamond-shaped window surrounding the loci, along the main diagonal and calculate the insulation score in each of these windows. There are many ways to calculate the insulation score. 

- In Crane et al. 2015, they simply add up total contact count per window. Such method is implemented in the FAN-C package. [https://vaquerizaslab.github.io/fanc/fanc-executable/fanc-analyse-hic/domains.html]

- In Lazaris et al., they divide the maximum contact count between the right and the left regions from the loci by the contact count in the middle regions (see oicture below) . Given that all the regions contain n interactions, the insulation score can be formulated as follows:

$Insulation score = \frac{max(L, R)}{C}$

<img width="370" alt="Screenshot 2024-08-18 at 22 02 34" src="https://github.com/user-attachments/assets/e8604365-24e5-4a00-bde1-1c9f4b9ba0d0">

Here is the function adopted from the C.origami paper [https://doi.org/10.1038/s41587-022-01612-8] that we can use for insulation scores:

``` python
import numpy as np
import matplotlib.pyplot as plt

# Radius is the distance in base pairs from the locus to the corner of a window
# Pseudocount coefficient helps us avoid division by zero

def chr_score(matrix, res, radius, pseudocount_coeff = 10):
    pseudocount = matrix.mean() * pseudocount_coeff
    pixel_radius = int(radius / res)
    scores = []
    for loc_i, loc in enumerate(range(len(matrix))):
        scores.append(point_score(loc, pixel_radius, matrix, pseudocount))
    return scores

def point_score(locus, radius, matrix, pseudocount):
    l_edge = max(locus - radius, 0)
    r_edge = min(locus + radius, len(matrix))
    l_mask = matrix[l_edge : locus, l_edge : locus]
    r_mask = matrix[locus : r_edge, locus : r_edge]
    center_mask = matrix[l_edge : locus, locus : r_edge]
    score = (max(l_mask.mean(), r_mask.mean()) +  pseudocount) /\
            (center_mask.mean() + pseudocount)
    return score

```
### Plotting HiC insulation score
In the way I plot the insulation score graph I would like to rotate the HiC map 45 degrees and present it as a triangle lying on the main diagonal. It is quite easy to implement in R, but not that trivial in Python, but I will try. The plotting code is adopted from the cooltools documentation [https://cooltools.readthedocs.io/en/latest/notebooks/insulation_and_boundaries.html]:

``` python
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import EngFormatter

def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
    import itertools
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im


bp_formatter = EngFormatter('b')
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)
```
Now we can utilise the function for inuslation plotting:

``` python
# Get cooler from data
GM12878 = cooler.Cooler('/Users/tennisnyjmac/Downloads/merged_sample_GM12878_10000.cool')
region = ['chr1', 10000000, 12000000]
data = np.log1p(GM12878.matrix(balance=False).fetch(region))
resolution = 10000


# Calculate insulation score
GM12878_insulation_scores = chr_score(data, radius = 50000)

f, ax = plt.subplots(figsize=(20, 15))
im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, cmap='fall', vmax = 10)
ax.set_aspect(0.9)
ax.set_ylim(0, 200000)
format_ticks(ax, rotate=False)

# Hide the x-axis
ax.xaxis.set_visible(False)

# Add color bar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1, aspect=6)
plt.colorbar(im, cax=cax)

ins_ax = divider.append_axes("bottom", size="30%", pad=0., sharex=ax)
ins_ax.plot(np.arange(len(GM12878_insulation_scores)) * resolution + region[1], GM12878_insulation_scores, label='Insulation score')
ins_ax.set_ylabel('Insulation Score')
ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4)
format_ticks(ins_ax, y=False, rotate=False)
ax.set_xlim(region[1], region[2])
ax.set_title('Hi-C GM12878 (res = 10kb)'+ ' '+chrom+': '+str(start_bp/1000000)+'-'+ str(end_bp/1000000)+' Mb')

plt.tight_layout()
plt.show()
```
The figure that we get:
<img width="1186" alt="Screenshot 2024-08-20 at 13 09 43" src="https://github.com/user-attachments/assets/c364d9e1-b815-4261-ab90-23729bfa1e0a">

### Calling boundaries
From the picture above we see that some regions have higher insulation scores and some have lower scores. Using the insulation score we can identify the boundaries and calculate their strenght. Here is a function that would find local_maxima which correspond to the TAD boundary:

``` python
# I would advise putting threshold equal to the mean of all insulation scores, but it could be calibrated and replaced
def find_boundaries(insulation_scores, threshold=np.nanmean(insulation_scores)):
    insulation_scores = np.array(insulation_scores)
    local_maxima = signal.argrelextrema(insulation_scores, np.greater)[0]
    
    if threshold is not None:
        local_maxima = local_maxima[insulation_scores[local_maxima] > threshold]

    return local_maxima

```
Now let's plot the inuslation score graph with boundaries identified:

``` python
boundaries = find_boundaries(GM12878_insulation_scores, np.nanmean(GM12878_insulation_scores))
boundary_positions = boundaries * resolution + region[1]

# Experiment with no threshold
all_boundaries = find_boundaries(GM12878_insulation_scores)
all_boundary_positions = all_boundaries * resolution + region[1]

# Plot the insulation score and boundaries
region = ['chr1', 10000000, 12000000]
GM12878_matrix = GM12878.matrix(balance=False).fetch(region)
GM12878_insulation_scores = chr_score(data, radius = 50000)
resolution = 10000
start_bp = region[1]
end_bp = region[2]


plt.figure(figsize=(15, 5))
plt.plot(np.arange(len(GM12878_insulation_scores)) * resolution + region[1], GM12878_insulation_scores, label='Insulation Score')
plt.scatter(all_boundary_positions, np.array(GM12878_insulation_scores)[all_boundaries], color='green', label='All Boundaries')
plt.scatter(boundary_positions, np.array(GM12878_insulation_scores)[boundaries], color='red', label='Boundaries')


# Set ticks for genomic positions
ax = plt.gca() 
n_ticks = 10
ticks = ["{0:.2f}".format(x / 1000000) for x in np.arange(start_bp, end_bp, (end_bp - start_bp) / n_ticks)]
ax.set_xticklabels(ticks)

plt.xlabel('Genomic Position')
plt.ylabel('Insulation Score')
plt.title('Inuslation score with TAD boundaries \n GM12878 Hi-C (res = 10kb) chr1: 10Mb - 12Mb')
plt.legend()
plt.show()

```
<img width="1023" alt="Screenshot 2024-08-20 at 14 46 27" src="https://github.com/user-attachments/assets/93cedf03-a4af-4f70-84d7-cbae3f26b8b2">

### Other HiC features
At higher resolutions we can witness the focal points and fountains on the HiC maps.
Calling fountains: https://github.com/agalitsyna/fontanka


### Read more:
TAD calling: https://vaquerizaslab.github.io/fanc/fanc-executable/fanc-analyse-hic/domains.html

HiChew, method for TAD calling: https://github.com/encent/hichew

Cool tools, insulation score and TAD calling: https://cooltools.readthedocs.io/en/latest/notebooks/insulation_and_boundaries.html

