import argparse
import os
from sre_parse import parse

import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
from copy import deepcopy
import glob
import tqdm
from itertools import combinations
import pyranges

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="ChipTableInterset - Tool to add additional columns to chipsummary table by intersecting with bed-like files",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-input_table", help="Input table", type=str,
                        required=True)
    parser.add_argument("-intersect_table", help="Intersect table", type=str,
                        required=True)
    parser.add_argument("-output", help="Output table", type=str, required=True)

    parser.add_argument("--chr_col", help="Column with chromosome", type=str, default="Chromosome")
    parser.add_argument("--start_col", help="Column with start position", type=str, default="Start")
    parser.add_argument("--end_col", help="Column with end position", type=str, default="End")
    parser.add_argument("--match_col", help="Column with matching flag (default intersect table file name)", type=str, default="")
    parser.add_argument("--add_cols", help="List of columns that will be added in case of matching", type=str, default="")
    parser.add_argument("-verbose", help="Log level: 0 - error, 1 - info, 2 - debug  ", type=int,default=1)
    args = parser.parse_args()

    input_table = pd.read_csv(args.input_table, sep=",")
    intersect_table = pd.read_csv(args.intersect_table, sep=",")
    if args.match_col == "":
        args.match_col = os.path.basename(args.intersect_table)
    adds_cols = []
    if args.add_cols != "":
        adds_cols = args.add_cols.split(",")

    pr_input=pyranges.PyRanges(input_table)
    #rename in intersect table args.chr_col to Chromosome, args.start_col to Start, args.end_col to End
    intersect_table=intersect_table.rename(columns={args.chr_col: 'Chromosome', args.start_col: 'Start', args.end_col: 'End'})
    for col in adds_cols:
        intersect_table[f'{col}_added'] = intersect_table[col]

    pr_intersect=pyranges.PyRanges(intersect_table)
    dftmp = deepcopy(pr_input.df)

    for index,row in tqdm.tqdm(pr_input.df.iterrows(),total=len(pr_input.df),desc='Compute intersection for each input interval',disable=args.verbose!=1):
        interval = pyranges.from_dict({"Chromosome": [row['Chromosome']], "Start": [row['Start']], "End": [row['End']]})
        res = pr_intersect.intersect(interval)
        if len(res.df) > 0:
            for col in adds_cols:
                dftmp.at[index,f"{col}_added"] = res.df[f"{col}_added"].to_list()[0]
                list_of_str = []
                for i in res.df[f"{col}_added"].to_list():
                    list_of_str.append(str(i))
                dftmp.at[index, f"{col}_all"] = ",".join(list_of_str)
            dftmp.at[index, f"{args.match_col}_amount"] = len(res.df)
            dftmp.at[index,args.match_col] = 1
        else:
            dftmp.at[index,args.match_col] = 0
    dftmp.to_csv(args.output,index=False)
