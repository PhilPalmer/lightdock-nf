#!/usr/bin/env python

import argparse
import pandas as pd
from shutil import copyfile

def parse_args():
    parser = argparse.ArgumentParser(description='Get data for the top PDB docking structure',
                                     usage='get_top.py --scores_path rank_by_scoring.list --pdb_code 5mi0 --out_path top_ranked_by_scores.csv')
    parser.add_argument('-s', '--scores_path', action='store', type=str, default='rank_by_scoring.list',
                        help="Input file with the list of PDB IDs and their scores")
    parser.add_argument('-p', '--pdb_code', action='store', type=str,
                        help="PDB ID to get the top docking structure for")
    parser.add_argument('-o', '--out_path', action='store', type=str, default='top_ranked_by_scores.csv',
                        help="Output file with the top docking result including the score")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    # Load input arguments
    args = parse_args()
    scores_path = args.scores_path
    pdb_code = args.pdb_code[0:4]
    # Load input file
    scores_df = pd.read_csv(scores_path, sep='(?<!,)\s+')
    scores_df = scores_df.head(1)
    scores_df['pdb'] = pdb_code
    scores_df.to_csv(args.out_path, index=False, sep=',')
    # Rename PDB file
    top_pdb = scores_df.iloc[0]
    pdb_old_fname = f'swarm_{top_pdb.Swarm}/{top_pdb.PDB}'
    pdb_new_fname = f'{pdb_code}_top_lightdock.pdb'
    copyfile(pdb_old_fname, pdb_new_fname)
    print(f'{pdb_new_fname}')
