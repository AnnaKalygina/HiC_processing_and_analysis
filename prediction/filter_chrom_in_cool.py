import pandas as pd
import numpy as np
import cooler
import cooltools.lib.plotting
from cooltools import insulation
import argparse
from pathlib import Path


def main(path, outdir):
    c = cooler.Cooler(path)
    
    primary_chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    
    bins = c.bins()[:]
    pixels = c.pixels()[:]
    
    primary_bins = bins[bins['chrom'].isin(primary_chromosomes)].reset_index(drop=True)
    
    old_to_new_bin_ids = {old: new for new, old in enumerate(primary_bins.index)}
    
    primary_pixels = pixels[
        (pixels['bin1_id'].isin(old_to_new_bin_ids)) & 
        (pixels['bin2_id'].isin(old_to_new_bin_ids))
    ].copy()
    
    primary_pixels['bin1_id'] = primary_pixels['bin1_id'].map(old_to_new_bin_ids)
    primary_pixels['bin2_id'] = primary_pixels['bin2_id'].map(old_to_new_bin_ids)
    
    
    cooler.create_cooler(
        outdir,
        bins=primary_bins,
        pixels=primary_pixels
    )
    
    print("Filtered cooler file created at:", outdir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help = 'Path to cool file to be filtered')
    parser.add_argument('outdir', help = 'Directory to save files to')
    args = parser.parse_args()
    outdir = Path(args.outdir)
    main(path=args.path, outdir=args.outdir)
    
