#!/usr/bin/env python
import argparse
import os

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

def get_combinations(input_list):
    all_combinations = []
    for r in range(1, len(input_list) + 1):
        for combo in combinations(input_list, r):
            all_combinations.append("_".join(combo))
    return all_combinations

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="ChipSummary - Tool for creation common summary table for ChIP-Seq results",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-indir", help="Directory with input bed files", type=str,
                        required=True)
    parser.add_argument("-file_mask", help="Input file mask. Default value '*.bed' ", type=str, default="*.bed")
    parser.add_argument("-out_perfix", help="Output file name perfix. Output filename format [out_perfix]_[type_of_regions]_report.csv", type=str, default="")

    parser.add_argument("-intersection_dir", help="Directory with bed files intersection with those will be shown in summary table", type=str,
                        required=False,default="")
    parser.add_argument("-names_filter",
                        help=f"Comma delimited list of names of regions in input bed files that "
                             f"should be taken into account. Default value: all regions from input files",
                        type=str, default="")
    parser.add_argument("-outdir", help="Directory with output files", type=str,
                        required=True)
    parser.add_argument("-filter_regions", help="File that will be used for filtering results by there regions",
                        type=str,default="")
    parser.add_argument("-tissue_spec", help="if specified we will use only this tissue for analysis",
                        type=str, default="")
    parser.add_argument("-genes_annotations", help="Path to gtf file with genes annotations",
                        type=str, default="Homo_sapiens.GRCh38.110.genes.gtf")
    parser.add_argument("-se_genes_range", help="Size of regione around enhancer where genes name will be writed out",
                        type=int, default=0)

    parser.add_argument("-genes_biotype", help="Genes biotype that will be taken into account (like protein_coding or so)"
                                               " By default all genes",
                        type=str, default="")
    parser.add_argument("-diffbind_regions",
                        help="Path to bed file with diffbind domains. ",
                        type=str, default="")


    parser.add_argument("-metadata", help="sample-features csv table in Diffbind metadata format", type=str,required=True)
    parser.add_argument("-cov_threshold", help="Minimal ratio of coverage at peak to differentiation. If set, "
                                               "peak in sample1 will be marked as diferential related to sample2"
                                               "even if is peak in sample 2, but mean coverage at this peak lower "
                                               "in 1/cov_threshold time than in sample 1. Default value 0 - and this feature will be disabled"
                                               " sample2_cov < sample1_cov*cov_threshold. Need column Coverage and coverage files",
                        type=float,default=0.0)

    parser.add_argument("-verbose", help="Log level: 0 - error, 1 - info, 2 - debug  ", type=int,default=1)
    args = parser.parse_args()

    input_files = glob.glob(args.indir + '/' + args.file_mask)

    PyRangesObjects = {}
    domain_counts={}
    metadata = pd.read_csv(args.metadata)
    metadata['filename'] = metadata['Peaks'].apply(os.path.basename)
    if args.tissue_spec != "":
        metadata = metadata[metadata['Tissue'] == args.tissue_spec]

    for i in ['SampleID','Tissue','Condition']:
        if i not in metadata.columns:
            raise RuntimeError(f'column {i} expected in metadata file {args.metadata}')
    if args.cov_threshold > 0.0 and 'Coverage' not in metadata.columns:
        raise RuntimeError(f'column Coverage expected in metadata file {args.metadata} because cov_threshold > 0.0')

    computed_domains = get_combinations(list(set(metadata['Condition'].to_list())))
    computed_labels = ['SE','enhancer','promoter']
    filter_bed = None
    if len(args.filter_regions) > 0:
        filter_bed = pyranges.read_bed(args.filter_regions)

    allowed_files = list(set(metadata['filename'].to_list()))

    coverage_files = []
    if args.cov_threshold > 0.0:
        coverage_files = metadata['Coverage'].to_list()
    if len(args.intersection_dir) > 0:
        intersection_files = glob.glob(args.intersection_dir + '/*.bed')
        if len(intersection_files) == 0:
            raise RuntimeError(f"No files found in {args.intersection_dir} with mask {args.file_mask}")
    if args.verbose > 0:
        print(f"Input directory:{args.indir}")
        print(f"Output directory:{args.outdir}")
        print(f"We use:")
        print(f"\tTissues and tissues combinations:{computed_domains}")
        if args.names_filter != "":
            print(f"\tRegion labels:{args.names_filter.split(',')}")
        else:
            print(f"\tRegion labels: all")
        print(f"\tBed files:{allowed_files}")
        if len(args.filter_regions) > 0:
            print(f"\tFilter_regions file:{args.filter_regions}")
        if len(coverage_files) > 0:
            print(f"\tCoverage files:{coverage_files}")
    diffbind_regions = None
    if args.diffbind_regions != "":
        diffbind_regions = pyranges.read_bed(args.diffbind_regions)
    if len(input_files) == 0:
        raise RuntimeError(f"No files found in {args.indir} with mask {args.file_mask}")
    if args.verbose > 1:
        for f in allowed_files:
            if f not in input_files:
                print(f"File {f} from metadata table are not found in {args.indir} with mask {args.file_mask}")

    for f in tqdm.tqdm(input_files, desc="Read files and construct PyRanges objects",disable= args.verbose!=1):
        if os.path.basename(f) not in allowed_files:
            if args.verbose > 1:
                print(f"Skip file {f} not in table with metadata")
            continue

        rdf = pyranges.read_bed(f).df
        if len(coverage_files) > 0:
            bwfn = metadata[metadata['filename']==os.path.basename(f)]['Coverage'].values[0]
            print(f"bwfn:{bwfn}")
            bw = pyranges.read_bigwig(bwfn)
            #iterate over all peaks in rdf and compute mean coverage in each peak and add it to rdf
            mean_coverage = []
            for index,row in rdf.iterrows():
                start_si = row['Start']
                end_si = row['End']
                super_interval = pyranges.from_dict({"Chromosome": [row['Chromosome']], "Start": [start_si], "End": [end_si]})
                cov_intersect = bw.intersect(super_interval, how='first', invert=False,nb_cpu=1)
                mean_coverage.append(cov_intersect.df['Value'].mean())
                print(f"cov_intersect:{cov_intersect.df['Value'].mean()}. {cov_intersect.df['Value']}")
            rdf['Coverage'] = mean_coverage
        rdf['Tissue'] = metadata[metadata['filename']==os.path.basename(f)]['Tissue'].values[0]
        rdf['Condition'] = metadata[metadata['filename'] == os.path.basename(f)]['Condition'].values[0]
        rdf['SampleID'] = metadata[metadata['filename'] == os.path.basename(f)]['SampleID'].values[0]
        rdf['Samples'] = None
        rdf['Conditions'] = None
        rdf['Genes'] = None
        rdf['nearest_gene'] = None
        if diffbind_regions is not None:
            rdf['db_score'] = None
            rdf['db_name'] = None
        #rdf['Treatment'] = metadata[metadata['filename'] == os.path.basename(f)]['Treatment'].values[0]
        PyRangesObjects[f] = pyranges.PyRanges(rdf)

    stat = {}
    PyRangeAll = None
    for f,obj in PyRangesObjects.items():
        if PyRangeAll is None:
            PyRangeAll = obj
        else:
            PyRangeAll = pyranges.PyRanges(pd.concat([PyRangeAll.df,obj.df]))
    list_of_allowed_names = args.names_filter.split(',')
    if list_of_allowed_names[0] == "":
        list_of_allowed_names = []
    if len(list_of_allowed_names) > 0:
        PyRangeSE = pyranges.PyRanges(PyRangeAll.df[PyRangeAll.df['Name'].isin(list_of_allowed_names)])
    else:
        PyRangeSE = PyRangeAll
    #show maximum columns in pandas
    #pd.set_option('display.max_columns', None)
    annot_gtf_file = args.genes_annotations
    gene_biotype = args.genes_biotype
    genes_range = args.se_genes_range
    annot = pyranges.read_gtf(annot_gtf_file)
    if len(gene_biotype) > 0:
        annot = pyranges.PyRanges(annot.df[annot.df['gene_biotype'] == gene_biotype])

    dftmp = deepcopy(PyRangeSE.df)
    i = 0
    MaxCPU = 1

    for index,row in tqdm.tqdm(PyRangeSE.df.iterrows(),total=len(PyRangeSE.df),desc='Compute samples and conditions for each interval',disable=args.verbose!=1):
        interval = pyranges.from_dict({"Chromosome": [row['Chromosome']], "Start": [row['Start']], "End": [row['End']]})
        if genes_range > 0:
            start_si = max([0,(row['Start'] + row['End']) // 2 - genes_range])
            end_si = (row['Start'] + row['End']) // 2 + genes_range
            super_interval = pyranges.from_dict({"Chromosome": [row['Chromosome']], "Start": [start_si],"End": [end_si]})
            annot_intersect = annot.intersect(super_interval, how='first', invert=False,nb_cpu=MaxCPU)
            #select nearest to (row['Start'] + row['End']) // 2 gene from annot_intersect
            if len(annot_intersect.df) > 0:
                list_of_genes = annot_intersect.df['gene_name'].dropna().unique().tolist()
                nearest_gene = annot_intersect.df.iloc[(annot_intersect.df['Start'] - (row['Start'] + row['End']) // 2).abs().argsort()[:1]]['gene_name'].tolist()[0]
            else:
                list_of_genes = []
                nearest_gene=""
        else:
            list_of_genes = []
            nearest_gene = ""
        intersect = PyRangeSE.intersect(interval, how='first', invert=False,nb_cpu=MaxCPU)

        list_of_samples = sorted(intersect.df['SampleID'].unique().tolist())
        if args.cov_threshold > 0.0:
            list_of_conditions=[]
            #If row['Coverage'] < intersected regions coverage * args.cov_threshold then we will add this condition to list_of_conditions
            for inter_index, inter_row in intersect.df.iterrows():
                if row['Coverage'] < inter_row['Coverage']*args.cov_threshold:
                    list_of_conditions.append(inter_row['Condition'])
            list_of_conditions = sorted(list(set(list_of_conditions)))
        else:
            list_of_condition = sorted(intersect.df['Condition'].unique().tolist())
        dftmp.loc[index,'Samples'] = ",".join(list_of_samples)
        dftmp.loc[index,'Conditions'] = ",".join(list_of_condition)
        dftmp.loc[index, 'Genes'] = ",".join(list_of_genes)
        dftmp.loc[index, 'nearest_gene'] = nearest_gene
        if diffbind_regions is not None:
            db_intersect = diffbind_regions.intersect(interval, how='first', invert=False,nb_cpu=MaxCPU)
            if len(db_intersect.df) > 0:
                dftmp.loc[index, 'db_score'] = db_intersect.df['Score'].max()
                dftmp.loc[index, 'db_name'] = db_intersect.df['Name'].max()
        i+=1
        #if i > 100:
        #    break
    #remove columns Strand and Tissue and Condition from dftmp

    clustered_pyrange = pyranges.PyRanges(dftmp).cluster(strand=False, count=True)
    dftmp = clustered_pyrange.df
    #sort dftmp by Descending by Cluster and score columns
    dftmp = dftmp.sort_values(by=['Cluster','Score'],ascending=[False,False])
    #group dftmp by Cluster column take minimum for Start column and maximum for End column,
    # take maximum for db_score and first values for all others
    if diffbind_regions is not None:
        dftmp = dftmp.groupby('Cluster').agg({'Chromosome':'first','Start':'min','End':'max','Score':'max',
                                              'SampleID': 'first','Samples':'first','Count':'first'
                                             ,'Conditions':'first','Genes':'first','nearest_gene':'first','db_score':'max','db_name':'first'}).reset_index()
    else:
        dftmp = dftmp.groupby('Cluster').agg(
            {'Chromosome': 'first', 'Start': 'min', 'End': 'max', 'Score': 'max', 'SampleID': 'first',
                  'Count':'first','Samples': 'first',
                 'Conditions': 'first', 'Genes': 'first','nearest_gene':'first'}).reset_index()
    dftmp = dftmp.drop(columns=['Cluster'])
    all_ids = set()
    dftmp['Samples'].fillna('', inplace=True)
    dftmp['Samples'].str.split(',').apply(all_ids.update)
    for id_ in all_ids:
        dftmp[id_] = dftmp['Samples'].str.contains(id_) #todo what if part of id is in another id?


    if len(args.intersection_dir) > 0:
        intersection_files = glob.glob(args.intersection_dir + '/*.bed')

        PyRangesIntersect = {}

        for f in tqdm.tqdm(intersection_files, desc="Read files and construct PyRanges objects",
                           disable=args.verbose != 1):
            PyRangesIntersect[f] = pyranges.read_bed(f)
        for f in intersection_files:
            dftmp[os.path.basename(f)] = 0
        for index, row in tqdm.tqdm(dftmp.iterrows(), total=len(dftmp),
                                    desc='Compute intersections of each SE with target samples',
                                    disable=args.verbose != 1):
            interval = pyranges.from_dict(
                {"Chromosome": [row['Chromosome']], "Start": [row['Start']], "End": [row['End']]})
            for f, pr in PyRangesIntersect.items():
                intersect = pr.intersect(interval, how='first', invert=False, nb_cpu=MaxCPU)
                dftmp.loc[index, os.path.basename(f)] = len(intersect.df)
    first_perfix = ""
    if args.out_perfix != "":
        first_perfix = args.out_perfix +"_"
    if args.names_filter != "":
        out_perfix = args.names_filter.replace(',','_').lower()
    else:
        out_perfix = 'all'
    dftmp.to_csv(os.path.join(args.outdir, first_perfix+out_perfix+'_report.csv'))




    #dftmp = dftmp.groupby('Cluster').first().reset_index()
