# Generating Contact Matrices
After processing our Hi-C data, we might want to visualize the contact matrix as a Hi-C map. While there are several sophisticated tools for visualization (e.g., Juicer Tools, HiGlass), a straightforward way to plot a matrix is by using the matplotlib package in Jupyter Notebook.

For this purpose, the only requirement is that the Hi-C data is in .cool or .mcool format. I am also providing a Jupyter Notebook for interactive follow-along data processing. Alternatively, you can use the provided Python scripts to run the analysis in a terminal.

### Extracting Matrix and Other Properties from Cooler
To load a .cool file into your notebook, we can use the cooler package. Loading data with cooler creates a new object of the Cooler class. This object contains attributes such as chromnames, chromsizes, binsize, and info. Several methods are available to retrieve the matrix, pixels, and bins tables for further manipulation:

``` python
region = ['chr8', 21000000, 23000000] # Alternatively ('chr8', 21000000, 23000000)

GM12878 = cooler.Cooler('path/to/cooler::resolution/10000')

GM12878_matrix = GM12878.matrix(balance=False, sparse=False) # Gives contact matrix of all chromosomes
GM12878_matrix = GM12878.matrix(balance=False, sparse=False).fetch(region) # Gives contact matrix of a specified region

GM12878_bins = GM12878.bins(balance=False, sparse=False)
GM12878_pixels = GM12878.pixels(balance=False, sparse=False)
```
Next, we will need matrix object for plotting HiC maps.

### Plotting single HiC map

``` python
# Set map attributes
region = ['chr8', 21000000, 23000000]
resolution = 10000
chrom = region[0]
start_bp = region[1]
end_bp = region[2]
start_matrix = start_bp // resolution
end_matrix = end_bp // resolution
n_ticks = 3

# How many pixlels are there in the window?
zoom = (end_bp - start_bp) / resolution
ticks = ["{0:.2f}".format(x / 1000000) for x in np.arange(start_bp, end_bp, (end_bp - start_bp) / n_ticks)]

# Set the color palette
color_map = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])

# Log transform the matrix 
GM12878_matrix = np.log1p(GM12878_matrix)

# Plot the map
fig, ax = plt.subplots(figsize=(5, 5))
cax = ax.imshow(GM12878_matrix, cmap = color_map)#, norm = norm)
cbar = fig.colorbar(cax, label='log(1 + contact frequency)')
ax.set_title('Hi-C GM12878 10kb\n'+ ' '+chrom+': '+str(start_bp/1000000)+'-'+ str(end_bp/1000000)+' Mb')

mb_formatter = EngFormatter(unit='b')
ax.xaxis.set_major_formatter(mb_formatter)
ax.yaxis.set_major_formatter(mb_formatter)
ax.set_xticks(np.arange(0, zoom, zoom / n_ticks))
ax.set_xticklabels(ticks)
ax.set_yticks(np.arange(0, zoom, zoom / n_ticks))
ax.set_yticklabels(ticks)
ax.set_ylabel('Genomic position (Mb)')

plt.show()

```
The following map is generated:


<img width="535" alt="Screenshot 2024-08-19 at 13 12 08" src="https://github.com/user-attachments/assets/faf13eb4-c985-4446-92fd-0dcc057e8fb4">



You can use the Python script plot_single_hic_map.py to generate a Hi-C plot in .png format for a specified region. If you want to generate maps for multiple windows, you can pass the following optional arguments:

--selected-windows with a .txt file containing a list of regions.
--stride-through with arguments [chr1, starting base pair, ending base pair, stride length] to generate continuous windows of a specified size.


### Plotting several maps side by side
Sometimes, it would be handy to compare the same regions across different cell lines. For this purpose you can use the following function:

``` python
def plot_2_maps(matrix1, matrix2, resolution, region):
    
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    norm = Normalize(vmin=0, vmax=np.max([matrix1.max(), matrix2.max()]))
    color_map = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])

    # Set ticks formatting and matrix sizing
    chrom = region[0]
    start_bp = region[1]
    end_bp = region[2]
    start_matrix = start_bp // resolution
    end_matrix = end_bp // resolution
    n_ticks = 5
    zoom = (end_bp - start_bp) / resolution
    ticks = ["{0:.2f}".format(x / 1000000) for x in np.arange(start_bp, end_bp, (end_bp - start_bp) / n_ticks)]

         
    # Set plotting options
    im1 = axs[0].imshow(matrix1, cmap=color_map, norm = norm)
    axs[0].set_title('GM12878\n'+ ' '+chrom+': '+str(start_bp/1000000)+'-'+ str(end_bp/1000000)+' Mb')
    axs[0].set_ylabel('Genomic position (Mb)')
    mb_formatter = EngFormatter(unit='b')
    axs[0].xaxis.set_major_formatter(mb_formatter)
    axs[0].yaxis.set_major_formatter(mb_formatter)
    axs[0].set_xticks(np.arange(0, zoom, zoom / n_ticks))
    axs[0].set_xticklabels(ticks)
    axs[0].set_yticks(np.arange(0, zoom, zoom / n_ticks))
    axs[0].set_yticklabels(ticks)
    
    im2 = axs[1].imshow(matrix2, cmap=color_map, norm = norm)
    axs[1].set_title('K562\n'+ ' '+chrom+': '+str(start_bp/1000000)+'-'+ str(end_bp/1000000)+' Mb')
    mb_formatter = EngFormatter(unit='b')
    axs[1].xaxis.set_major_formatter(mb_formatter)
    axs[1].yaxis.set_major_formatter(mb_formatter)
    axs[1].set_xticks(np.arange(0, zoom, zoom / n_ticks))
    axs[1].set_xticklabels(ticks)
    axs[1].set_yticks(np.arange(0, zoom, zoom / n_ticks))
    axs[1].set_yticklabels(ticks)

    for ax in axs:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(ax.images[0], cax=cax)

    plt.tight_layout()
    plt.show()
```
Now, prepare coolers for side by side comparison:

``` python
# Load matrices
GM12878 = cooler.Cooler('path/to/cooler::resolution/10000')
K562 = cooler.Cooler('path/to/cooler::resolution/10000')

#Get matrices and upscale the one with the lower coverage by the ratio between the mean coverages
region = ['chr1', 65000000, 67000000]
GM12878_matrix = GM12878.matrix().fetch(region)
K562_matrix = K562.matrix().fetch(region)

ratio = np.mean(np.diagonal(GM12878_matrix)) / np.mean(np.diagonal(K562_matrix))
K562_matrix = K562_matrix * ratio

# Log-transform the matrices 
GM12878_matrix = np.log1p(GM12878_matrix + 1)
K562_matrix = np.log1p(K562_matrix + 1)

plot_2_maps(GM12878_matrix, K562_matrix, 10000, region)
```


We get the following graph:

<img width="861" alt="Screenshot 2024-08-19 at 13 51 36" src="https://github.com/user-attachments/assets/4de7af86-d17a-47e8-b766-0ef5378b726b">

Alternatively, side-by-side maps can be generated using the Python script plot_2_hic_maps.py. To generate maps for multiple windows, you can pass optional arguments:

--selected-windows with a .txt file containing a list of regions.
--stride-through with arguments [chr1, starting base pair, ending base pair, stride length] to generate continuous windows of a specified size.

### Plotting expected vs. observed
The expected matrix is derived by calculating the average interaction frequency between loci at a given distance from the diagonal across the entire genome or within a specific chromosome. This matrix accounts for the inherent bias where contact frequency naturally decreases as the genomic distance between loci increases. In general terms, it models the general decay in contact frequency due to distance, without considering specific local interactions like loops or domains. The Observed vs. Expected (O/E) matrix is then generated by subtracting the expected contact frequencies from observed contact frequencies. Each cell in the expected matrix reflects the expected contact frequency based solely on distance, independent of structural features such as loops or TADs that might appear in the observed data. It also highlights regions with depletion in interaction frequencies.

To calculate the expected matrix, I average the interaction frequencies between bins at the same distance from each other. For instance, the first central diagonal represents the mean interaction frequency within the same bins, while the second diagonal reflects the interaction frequency between bins separated by one binsize.

``` python
def compute_expected(matrix):
    n = matrix.shape[0]
    expected = np.zeros_like(matrix)
    for d in range(n):
        diag_values = np.diag(matrix, k=d)
        mean_value = np.nanmean(diag_values)
        np.fill_diagonal(expected[d:], mean_value)
        np.fill_diagonal(expected[:, d:], mean_value)
    return expected
```

Next, I derive the Observed vs. Expected (OE) matrix via substraction:

``` python
def plot_OE_comparison(observed, expected, OE, resolution, region):
    
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))

    color_map = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])

    # Set ticks formatting and matrix sizing
    chrom = region[0]
    start_bp = region[1]
    end_bp = region[2]
    start_matrix = start_bp // resolution
    end_matrix = end_bp // resolution
    n_ticks = 5
    zoom = (end_bp - start_bp) / resolution
    ticks = ["{0:.2f}".format(x / 1000000) for x in np.arange(start_bp, end_bp, (end_bp - start_bp) / n_ticks)]

         
    # Set plotting options
    im1 = axs[0].imshow(observed, cmap=color_map)#, vmin = 0, vmax = 10)
    axs[0].set_title('Observed \n'+ ' '+chrom+': '+str(start_bp/1000000)+'-'+ str(end_bp/1000000)+' Mb')
    axs[0].set_ylabel('Genomic position (Mb)')
    mb_formatter = EngFormatter(unit='b')
    axs[0].xaxis.set_major_formatter(mb_formatter)
    axs[0].yaxis.set_major_formatter(mb_formatter)
    axs[0].set_xticks(np.arange(0, zoom, zoom / n_ticks))
    axs[0].set_xticklabels(ticks)
    axs[0].set_yticks(np.arange(0, zoom, zoom / n_ticks))
    axs[0].set_yticklabels(ticks)
    
    im2 = axs[1].imshow(expected, cmap=color_map)#, vmin = 0, vmax = 10)
    axs[1].set_title('Expected \n'+ ' '+chrom+': '+str(start_bp/1000000)+'-'+ str(end_bp/1000000)+' Mb')
    mb_formatter = EngFormatter(unit='b')
    axs[1].xaxis.set_major_formatter(mb_formatter)
    axs[1].yaxis.set_major_formatter(mb_formatter)
    axs[1].set_xticks(np.arange(0, zoom, zoom / n_ticks))
    axs[1].set_xticklabels(ticks)
    axs[1].set_yticks(np.arange(0, zoom, zoom / n_ticks))
    axs[1].set_yticklabels(ticks)

    im3 = axs[2].imshow(OE, cmap='bwr', vmin = -2, vmax = 2)
    axs[2].set_title('Observed vs. Expected\n'+ ' '+chrom+': '+str(start_bp/1000000)+'-'+ str(end_bp/1000000)+' Mb')
    mb_formatter = EngFormatter(unit='b')
    axs[2].xaxis.set_major_formatter(mb_formatter)
    axs[2].yaxis.set_major_formatter(mb_formatter)
    axs[2].set_xticks(np.arange(0, zoom, zoom / n_ticks))
    axs[2].set_xticklabels(ticks)
    axs[2].set_yticks(np.arange(0, zoom, zoom / n_ticks))
    axs[2].set_yticklabels(ticks)

    for ax in axs:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(ax.images[0], cax=cax)

    plt.tight_layout()
    plt.show()

# Set matrices

GM12878_matrix = GM12878.matrix(balance=False, sparse = False).fetch(region)
expected_matrix = compute_expected(GM12878_matrix)

GM12878_matrix = np.log1p(GM12878_matrix)
expected_matrix = np.log1p(expected_matrix)
oe_matrix = GM12878_matrix - expected_matrix

plot_OE_comparison(GM12878_matrix, expected_matrix , oe_matrix , 10000, region)
```
We get the following graph:

<img width="1138" alt="Screenshot 2024-08-19 at 16 41 40" src="https://github.com/user-attachments/assets/1af9054f-3d77-4427-8e9e-a02973e14719">


The way to interpret it is the following:
- Values > 0 indicate regions with more interactions than expected, suggesting structural features like chromatin loops.
- Values < 0 indicate regions with fewer interactions than expected, which might be indicative of boundaries or insulators.


# Read more:
[1] Cooler API dosumentation: https://cooler.readthedocs.io/en/latest/api.html#cooler.Cooler

[2] Source code for cooler.api: https://cooler.readthedocs.io/en/latest/_modules/cooler/api.html#Cooler.matrix

[3] Homer's guide on plotting: http://homer.ucsd.edu/homer/interactions/HiCmatrices.html



