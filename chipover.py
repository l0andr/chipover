import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import glob

import tqdm
from pybedtools import BedTool
from itertools import combinations
import pyranges

def get_combinations(input_list):
    all_combinations = []
    for r in range(1, len(input_list) + 1):
        for combo in combinations(input_list, r):
            all_combinations.append("_".join(combo))
    return all_combinations

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="ChipOver - Tool for compute common domain structure by overlap enchancers regions detected by LILY",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-indir", help="Directory with input bed files", type=str,
                        required=True)
    parser.add_argument("-outdir", help="Directory with output bed files", type=str,
                        required=True)
    parser.add_argument("-minoverlap", help="Minimum intersections", type=int,
                        default=0)
    parser.add_argument("-filter_regions", help="File that will be used for filtering results by there regions",
                        type=str,default="")
    parser.add_argument("-tissue_spec", help="if specified we will use only this tissue for analysis",
                        type=str, default="")

    parser.add_argument("-label", help="Additional output label",
                        type=str, default="")

    parser.add_argument("-metadata", help="sample-features csv table", type=str,default="")
    parser.add_argument("-verbose", help="Log level: 0 - error, 1 - info, 2 - debug  ", type=int,default=1)
    args = parser.parse_args()

    input_files = glob.glob(args.indir + '/*.bed')

    PyRangesObjects = {}
    domain_counts={}
    metadata = pd.read_csv(args.metadata)
    metadata['filename'] = metadata['Peaks'].apply(os.path.basename)
    if args.tissue_spec != "":
        metadata = metadata[metadata['Tissue'] == args.tissue_spec]

    for i in ['SampleID','Tissue','Condition']:
        if i not in metadata.columns:
            raise RuntimeError(f'column {i} expected in metadata file {args.metadata}')
    computed_domains = get_combinations(list(set(metadata['Condition'].to_list())))
    computed_labels = ['SE','enhancer','promoter']
    filter_bed = None
    if len(args.filter_regions) > 0:
        filter_bed = pyranges.read_bed(args.filter_regions)

    allowed_files = list(set(metadata['filename'].to_list()))
    if args.verbose > 0:
        print(f"Input directory:{args.indir}")
        print(f"Output directory:{args.outdir}")
        print(f"We use:")
        print(f"\tTissues and tissues combinations:{computed_domains}")
        print(f"\tRegion labels:{computed_labels}")
        print(f"\tBed files:{allowed_files}")
        print(f"\tMinimum regions per domain:{args.minoverlap+1}")
        if len(args.filter_regions) > 0:
            print(f"\tFilter_regions file:{args.filter_regions}")
        if len(args.label) > 0:
            print(f"\tLabel:{args.label}")

    for f in tqdm.tqdm(input_files, desc="Read files and construct PyRanges objects",disable= args.verbose!=1):
        if os.path.basename(f) not in allowed_files:
            continue
        rdf = pyranges.read_bed(f).df
        rdf['Tissue'] = metadata[metadata['filename']==os.path.basename(f)]['Tissue'].values[0]
        rdf['Condition'] = metadata[metadata['filename'] == os.path.basename(f)]['Condition'].values[0]
        rdf['SampleID'] = metadata[metadata['filename'] == os.path.basename(f)]['SampleID'].values[0]
        rdf['Samples'] = None
        rdf['Conditions'] = None
        rdf['Genes'] = None
        #rdf['Treatment'] = metadata[metadata['filename'] == os.path.basename(f)]['Treatment'].values[0]
        PyRangesObjects[f] = pyranges.PyRanges(rdf)
    stat = {}
    for f, bo in tqdm.tqdm(PyRangesObjects.items(), desc=f"Compute statistics for all input samples",
                           disable= args.verbose != 1 ):
        sample_id = metadata[metadata['filename'] == os.path.basename(f)]['SampleID'].values[0]
        stat_label = {}
        for label in computed_labels:
            if isinstance(label, str):
                bo_filtered = pyranges.PyRanges(bo.df[bo.df['Name'] == label])
            elif isinstance(label, list):
                bo_filtered = pyranges.PyRanges(bo.df[bo.df['Name'].isin(label)])
            else:
                raise RuntimeError(f"Incorrect type of label {type(label)}")
            if isinstance(label, str):
                label_str = label
            elif isinstance(label, list):
                label_str = "+".join(label)
            stat_label[label_str] = bo_filtered.__len__()
        stat[sample_id] = stat_label
        stat_df = pd.DataFrame(stat)
        stats_table_file_name = "stats.csv"
        pd.DataFrame(stat_df).to_csv(os.path.join(args.outdir, stats_table_file_name))

    PyRangeAll = None
    for f,obj in PyRangesObjects.items():
        if PyRangeAll is None:
            PyRangeAll = obj
        else:
            PyRangeAll = pyranges.PyRanges(pd.concat([PyRangeAll.df,obj.df]))
    PyRangeAllPatients = PyRangeAll #pyranges.PyRanges(PyRangeAll.df[PyRangeAll.df['Tissue'] == 'patient'])
    PyRangeSE = pyranges.PyRanges(PyRangeAll.df[PyRangeAll.df['Name'] == 'SE'])

    domain_counts = {}
    for label in computed_labels:
        common_domains = {}
        if isinstance(label, str):
            bo_filtered = pyranges.PyRanges(PyRangeAllPatients.df[PyRangeAllPatients.df['Name'] == label])
        elif isinstance(label, list):
            bo_filtered = pyranges.PyRanges(PyRangeAllPatients.df[PyRangeAllPatients.df['Name'].isin(label)])
        else:
            raise RuntimeError(f"Incorrect type of label {type(label)}")

        for cd in computed_domains:
            allowed_tissue = cd.split('_')
            common_domains[cd] = None
            bo_for_overlap = pyranges.PyRanges(bo_filtered.df[bo_filtered.df['Condition'].isin(allowed_tissue)])
            bo_for_exclude = pyranges.PyRanges(bo_filtered.df[~bo_filtered.df['Condition'].isin(allowed_tissue)])
            common_domains[cd] = bo_for_overlap.merge(count=True)
            if len(allowed_tissue) > 1:
                for one_of_allowed_tissue in allowed_tissue:
                    bo_for_intersect = pyranges.PyRanges(bo_filtered.df[bo_filtered.df['Condition'].isin([one_of_allowed_tissue])])
                    common_domains[cd] = common_domains[cd].intersect(bo_for_intersect, how='first', invert=False)
            common_domains[cd] = common_domains[cd].intersect(bo_for_exclude, how='first', invert=True)
            common_domains[cd] = pyranges.PyRanges(common_domains[cd].df[common_domains[cd].df['Count'] > args.minoverlap ])
            if filter_bed is not None:
                common_domains[cd] = common_domains[cd].intersect(filter_bed, how='first', invert=False)
            if isinstance(label, str):
                label_str = label
            elif isinstance(label, list):
                label_str = "+".join(label)
            outfilename = f"{label_str}-{cd}.bed"
            common_domains[cd].to_bed(path=os.path.join(args.outdir, outfilename),keep=True)
            if label_str not in domain_counts:
                domain_counts[label_str] = {}
            domain_counts[label_str][cd] = common_domains[cd].__len__()

    counts_table_file_name = "counts.csv"
    pd.DataFrame(domain_counts).to_csv(os.path.join(args.outdir, counts_table_file_name))

    #print(PyRangeAll)
