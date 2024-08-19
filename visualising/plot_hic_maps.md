# Generating Contact Matrices
After we processed our data, we might want to visualise the contact matrix in a form of a HiC map. 
Even though there are several handy softwares used for visulalisation (Jucier Tools, HiGlass), the simplest way to plot a matrix is by using a matplotlib package in Jupyter Notebook.
For this purpose the only thing we need is the HiC file in a .cool or .mcool format. I am also attaching the jupyter notebook for interactive follow-along data processing. Alternatively, there are python scripts that could be run using Terminal.  

### Extracting matrix and other properties from cooler
The `.cool` file is loaded to the notebook using cooler package. When loading data with Cooler the new object of a cooler class is created. Cooler object contains attributes chromnames, chromsizes, binsize, and info. There are methods that could be used for retrieving matrix, pixels and bins tables for further manipulation:


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

```
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



You can use the python script `plot_single_hic_map.py` to generate a HiC plot .png of a specified region.

### Plotting several maps side by side

### Plotting expected vs. observed graph

# Read more:
[1] Cooler API dosumentation: https://cooler.readthedocs.io/en/latest/api.html#cooler.Cooler
[2] Source code for cooler.api: https://cooler.readthedocs.io/en/latest/_modules/cooler/api.html#Cooler.matrix


