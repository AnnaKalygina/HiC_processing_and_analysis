#!/usr/bin/env python

import pyranges as pr
import numpy as np
import pandas as pd
import argparse
from pathlib import Path

'''
The dfs passed must contain three columns = ['Chromosome'. 'Start', 'End']

'''

def main(args):

    dfs = pd.concat([pr.read_bed(bed_file, as_df=True) for bed_file in args.bed_files])
    filtered_data = filter_dataframes(dfs, args.merge_threshold, args.drop_threshold)
    filtered_data = filtered_data[['Chromosome', 'Start', 'End']]
    
    filtered_data.to_csv(save_path, sep='\t', header=False, index=False)

### Pass the indecies of rows that are less than merge_threshold bases apart form each other. 
### This function gives the pairs of indecies that must be merged together. Two pointers function. 
def transform_to_intervals(array):
    intervals = []
    start = array[0]
    end = array[0]

    for i in range(1, len(array)):
        if array[i] == array[i-1] + 1:
            end = array[i]
        else:
            intervals.append((start, end + 1))
            start = array[i]
            end = array[i]

    intervals.append((start, end + 1))
    return intervals

### Assuming the rows are sorted. 
def merge_rows(df, first_index, second_index):
    start = min(df.loc[first_index, 'Start'], df.loc[second_index, 'Start'])
    end = max(df.loc[first_index, 'End'], df.loc[second_index, 'End'])
    final_row = [df.loc[first_index, 'Chromosome'], start, end]
    return final_row

def filter_dataframes(dfs, merge_threshold, drop_threshold):
    ignore_dataset = dfs.sort_values(by=['Chromosome', 'Start']).reset_index(drop=True)
    filtered_ignore_dataset = pd.DataFrame(columns=['Chromosome', 'Start', 'End'])

    for chromosome in ignore_dataset['Chromosome'].unique():
        df = ignore_dataset[ignore_dataset['Chromosome'] == chromosome].reset_index(drop=True)

        ### Do not process unconventional chromosomes
        if len(df) <= 1:
            filtered_ignore_dataset = pd.concat([filtered_ignore_dataset, df])
            continue
        ### Calculate the distance between intervals. Obtain indecies where the distance is lower than a threshold.
        starts = df["Start"].to_numpy()[1:]
        ends = df["End"].to_numpy()[:-1]
        difference = starts - ends
        difference_idx = np.where(difference < merge_threshold)[0].tolist()

        intervals_to_merge = transform_to_intervals(difference_idx)
        ### Obtain indecies of rows that were modified and drop them
        rows_to_drop = list(sum(intervals_to_merge, ()))

        rows_to_append = []
        for i in intervals_to_merge:
            merged_row = merge_rows(df, i[0], i[1])
            rows_to_append.append(merged_row)

        df_to_append = pd.DataFrame(rows_to_append, columns=['Chromosome', 'Start', 'End'])

        df = df.drop(rows_to_drop).reset_index(drop=True)
        df = pd.concat([df, df_to_append]).sort_values('Start').reset_index(drop=True)

        df['Length'] = df['End'] - df['Start']
        df = df.drop(df[df['Length'] < drop_threshold].index)

        filtered_ignore_dataset = pd.concat([filtered_ignore_dataset, df]).reset_index(drop=True)

    return filtered_ignore_dataset


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter windows of interest against problematic regions')
    parser.add_argument('--bed_files', nargs='+', help='Paths to .bed files to be processed')
    parser.add_argument('--save_path', help='Directory to save the filtered BED file')
    parser.add_argument('--merge_threshold', type=int, default=2500,
                        help='Max length between regions to be merged [default: 2500]')
    parser.add_argument('--drop_threshold', type=int, default=5000,
                        help='Min length of regions to be kept [default: 5000]')
    
    args = parser.parse_args()
    save_path = Path(args.save_path)

    main(args)
