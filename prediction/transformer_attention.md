# Recording attention of a model to specific patterns in the data

The model is based on transformer, the transformer has 8 heads. During training, the data fed to the transformer is attributed with specific weights, this attribution and tailoring of weight is called attention. 
Different models have different attention mechanisms, but mostly it is a black box for us to know what parts of the data the model will be attending to the most. 
However, we can record attention to different regions of DNA during training by extracting the weights from these 8 transformer heads. The validation dataset is based solely in chromosome ten (chrX), so the data always refers to this one chromosome.

### Extract TSS data from the SK1 yeast strain 

I had an assumption that the greatest attentnion would be associated with transcription start sites (TSS) and promoter regions, as they are the most biologically relevant. So, at first we need to extract data on promoter and TSS regions from SK1 yeast. 
The data is available in the file SK1_PacBio.all_feature_modified_2genes.gff uploaded in this folder. Load it to the Jupyter notebook:

``` python
gff_df = pd.read_csv('/Users/tennisnyjmac/Downloads/SK1_PacBio.all_feature_modified_2genes.gff', sep='\t', comment='#', header=None, 
                     names=['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

# Sort into different types of TSS (negative and positive strand)
tss_pos_array = gff_df.loc[(gff_df['type'] == 'mRNA') & (gff_df['chr'] == 'chrX') & (gff_df['strand'] == '+')]['start'].array
tss_neg_array = gff_df.loc[(gff_df['type'] == 'mRNA') & (gff_df['chr'] == 'chrX') & (gff_df['strand'] == '-')]['end'].array
all_tss = np.hstack((tss_pos_array, tss_neg_array))
end_pos_array = gff_df.loc[(gff_df['type'] == 'mRNA') & (gff_df['chr'] == 'chrX') & (gff_df['strand'] == '+')]['end'].array
end_neg_array = gff_df.loc[(gff_df['type'] == 'mRNA') & (gff_df['chr'] == 'chrX') & (gff_df['strand'] == '-')]['start'].array
all_end = np.hstack((end_pos_array, end_neg_array))
```
Some promoters could exert activation on both strand, they are called bidirectional promoters. The way to filter them:

```python
# Find the transcriprion start sites with bidirectional promoters
threshold_bp = 1000
bidirectional_promoters = []

for bp_pos in range(len(tss_pos_array)):
    for bp_neg in range(len(tss_neg_array)):
        if abs(tss_pos_array[bp_pos] - tss_neg_array[bp_neg]) <= 1000:
            bidirectional_promoters.append((tss_pos_array[bp_pos], tss_neg_array[bp_neg]))

pos_bidirectional_tss = [bidirectional_promoters[bp][0] for bp in range(len(bidirectional_promoters))]
neg_bidirectional_tss = [bidirectional_promoters[bp][1] for bp in range(len(bidirectional_promoters))]
```
### Extract weights for the windows the model was trained on
The weights extracted from transformer are continuous, however, we tested them on specific windows, so we need to extract weights for exactly those windows:

``` python
#Extract chromosome coordinates pairs from the data the model was tested on:
def extract_chrom_coor(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        pairs = []
        pattern = re.compile(r'features_(\d+)-(\d+)\.npz')
        for line in lines:
            match = pattern.search(line)
            if match:
                start = int(match.group(1))
                end = int(match.group(2))
                pairs.append((start, end))
                    
    return pairs
coord = extract_chrom_coor('/Users/tennisnyjmac/Downloads/datasheet_val.txt')


# Extract attention weights for these coordinates:
def attn_weights_into_dict(data, pairs_chrom_coor):
    base_pair_data = {}
    
    for idx_window in range(len(pairs_chrom_coor)):
        attn_weights = np.repeat(data[idx_window],8)
        start_bp = pairs_chrom_coor[idx_window][0]
        end_bp = pairs_chrom_coor[idx_window][1]
        for idx_bp in range(start_bp, end_bp):
            if idx_bp not in base_pair_data:
                base_pair_data[idx_bp] = []
            base_pair_data[idx_bp].append(attn_weights[idx_bp - start_bp])

    return base_pair_data

all_attn_scores = np.load('/Users/tennisnyjmac/Downloads/_attn_weights.npy', allow_pickle=True)
all_attn_scores_dict = attn_weights_into_dict(all_attn_scores, coord)


base_pairs = sorted(attn_weight_dict.keys())
attn_weights_means = [np.mean(attn_weight_dict[bp]) for bp in base_pairs]

plt.figure(figsize=(20,6))
plt.plot(base_pairs, attn_weights_means, label='attn_weight', color='blue')
for i in range(len(all_tss)):
    plt.axvline(all_tss[i], color='purple', alpha = 0.5)
for _ in range(len(all_end)):
    plt.axvline(all_end[_], color='orange', alpha = 0.5)

plt.xlabel('Base Pair')
plt.ylabel('Attention weight')
plt.xlim(0,53382)
plt.legend()
plt.show()
```
The following graph is produced. The vertical lines represent the specificities of the TSSs.

<img width="1126" alt="Screenshot 2024-08-21 at 13 47 38" src="https://github.com/user-attachments/assets/04538f8c-1111-49cc-951c-25d4c1c9f68c">


It would be cool to aggregate attention weights depending to their position relative to TSSs:
``` python
base_pairs = all_attn_scores_dict.keys()
attn_to_bp = [np.mean(all_attn_scores_dict[bp]) for bp in base_pairs]

def get_area_around_tss(tss_coord, left, right, attn_scores):
    return attn_scores[(tss_coord+left):(tss_coord+right+1)]

# Calculate average attention score for positive TSS
tss_pos_excl_array = [tss for tss in tss_pos_array if tss not in pos_bidirectional_tss]
pos_tss_attn_scores = []
for tss in tss_pos_excl_array:
    pos_tss_attn_scores.append(get_area_around_tss(tss, int(-5000), int(5000), attn_to_bp))
    
len(pos_tss_attn_scores) == len(tss_pos_excl_array)

mean_attn_scores_around_pos_tss = np.mean(pos_tss_attn_scores, axis = 0)
plt.figure(figsize=(20,6))
plt.plot(range(-5000,5001), mean_attn_scores_around_pos_tss,  label='attn_weights', color='blue')
plt.axvline(0, color='orange', alpha = 0.5)

plt.xlabel('Base Pair')
plt.title('Mean attention weights for positive strand tss (bidirecional promoters excluded)')
plt.ylabel('Mean Attention weights')
plt.xlim(-5000, 5000)
plt.legend()
plt.show()
```

<img width="1140" alt="Screenshot 2024-08-21 at 13 56 06" src="https://github.com/user-attachments/assets/1df1e2fe-a06b-4af2-aa1f-e5747194dedc">

Do the same procedure for exclusively negative promoters:

<img width="1121" alt="Screenshot 2024-08-21 at 13 56 48" src="https://github.com/user-attachments/assets/03bc14ef-67b6-4ad1-b852-9f53863d3222">

And now bidirectional promoters and combine all three in the same graph:

<img width="1136" alt="Screenshot 2024-08-21 at 13 57 19" src="https://github.com/user-attachments/assets/5e8deb83-864b-48f4-93ef-84464b5fa2e1">

Interesting results!


