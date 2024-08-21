#!/usr/bin/env python

import pyranges as pr
import numpy as np
import pandas as pd
import argparse
from pathlib import Path
from tqdm import tqdm


def main(path_to_windows, path_to_ignore, save_path, total_threshold, middle_threshold):

    windows_of_interest = pr.read_bed(path_to_windows, as_df=True)
    windows_to_ignore = pr.read_bed(path_to_ignore, as_df=True)
    
    filtered_windows_of_interest = filter_windows(windows_of_interest, windows_to_ignore, total_threshold, middle_threshold)
    filtered_windows_of_interest = filtered_windows_of_interest[['Chromosome', 'Start', 'End']]
    
    filtered_windows_of_interest.to_csv(save_path, sep='\t', header=False, index=False)
    
    
def filter_windows(windows_df, ignore_dataset, total_threshold, middle_threshold):
    filtered_windows = []

    for _, window in tqdm(windows_df.iterrows()):
        chromosome = window['Chromosome']
        start = window['Start']
        end = window['End']
        window_length = end - start
        middle_start = start + window_length * 0.25
        middle_end = end - window_length * 0.25

        overlapping_regions = ignore_dataset[
            (ignore_dataset['Chromosome'] == chromosome) &
            (ignore_dataset['Start'] < end) &
            (ignore_dataset['End'] > start)]
        
        total_overlap_length = overlapping_regions.apply(lambda row: min(end, row['End']) - max(start, row['Start']), 
                                                         axis = 1).sum()

        middle_overlap_length = overlapping_regions.apply(lambda row: min(middle_end, row['End']) - max(middle_start, row['Start']),
                                                          axis=1).sum()

        if total_overlap_length <= total_threshold and middle_overlap_length <= middle_threshold:
            
            filtered_windows.append(window)

    return pd.DataFrame(filtered_windows)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter windows of interest against problematic regions')
    parser.add_argument('path_to_windows', help='Path to .bed file with windows to be processed')
    parser.add_argument('path_to_ignore', help='Path to .bed windows to be filtered against' )
    parser.add_argument('save_path', help='Directory to save files to. Will be created if need but not its parents')
    parser.add_argument('--total_threshold', type=int, default=50000,
                        help='Max length of total overlap between window of interest and a problematic regions [default: 50000]')
    parser.add_argument('--middle_threshold', type=int, default=5000,
                        help='Max length of overlap between middle 50% of window of interest and a problematic region [default: 5000]')
    argv = parser.parse_args()
    save_path = Path(argv.save_path)

    main(argv.path_to_windows, argv.path_to_ignore, save_path, total_threshold=argv.total_threshold, middle_threshold=argv.middle_threshold)
    