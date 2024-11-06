from copy import deepcopy

import numpy as np
import pandas as pd
import tqdm
from typing import List,Dict
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import argparse

def compute_specifity_statistics(df:pd.DataFrame)->Dict[str,Dict[int,int]]:
    """
    Compute specifity statistics.

    """
    #Get Unique values of 'Conditions' column
    conditions = df['Conditions'].unique()
    #iterate over conditions
    result = {}
    for condition in conditions:
        df_tmp = deepcopy(df[df['Conditions'] == condition])
        #for each row compute number of samples in samples list ['sample1','sample2','sample3'] in Samples column
        df_tmp.loc[:,'Number_of_samples'] = df_tmp['Samples'].apply(lambda x: len(x.split(',')))
        nos_unique = df_tmp['Number_of_samples'].unique()
        nos_list = {}
        for nos in nos_unique:
            nos_list[nos] = len(df_tmp[df_tmp['Number_of_samples'] == nos])
        result[condition] = nos_list

    return result

def compute_specifity_statistics_genes(df:pd.DataFrame)->Dict[str,Dict[int,int]]:
    """
    Compute specifity statistics.

    """
    #Get Unique values of 'Conditions' column
    conditions = df['Conditions'].unique()
    #iterate over conditions
    result = {}
    for condition in conditions:
        df_tmp = deepcopy(df[df['Conditions'] == condition])
        #for each row compute number of samples in samples list ['sample1','sample2','sample3'] in Samples column
        df_tmp['Number_of_samples'] = df_tmp['Samples'].apply(lambda x: len(x.split(',')))
        nos_unique = df_tmp['Number_of_samples'].unique()
        nos_list = {}
        for nos in nos_unique:
            df_tmp_tmp = df_tmp[df_tmp['Number_of_samples'] == nos]
            genes = []
            for index,row in df_tmp_tmp.iterrows():
                try:
                    genes+=row['Genes'].split(',')
                except AttributeError:
                    continue
            nos_list[nos] = len(set(genes))
        result[condition] = nos_list

    return result

def plot_specifity_statistics(specifity_statistics:Dict[str,Dict[int,int]],title=""):
    fig, ax = plt.subplots()
    #plot each conditions separately in different bar on the same bar plot
    i=0
    i_inc = 0.7/len(specifity_statistics)
    cond_list = list(specifity_statistics.keys())

    for cond in sorted(cond_list):
        values = specifity_statistics[cond]
        ax.bar([x-i_inc*1.3+i+i_inc for x in values.keys()], values.values(), width=i_inc*1.3,label=cond)
        i+=i_inc
    #set xticks with step 1
    ax.set_xticks(range(1,max([max(x.keys()) for x in specifity_statistics.values()])+1))
    ax.set_xlabel('Number of samples')
    ax.set_ylabel(f'Number of {title}')
    ax.grid()
    ax.legend()
    return fig,ax

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="ChipRnaCombiner - Tool for combine ChipSeq and RNASeq results to find genes that are controled by Super Enhancers\ Enhancers and\or have enrichment in promotor regions and are upregulated in RNASeq",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-indir", help="Directory with input files", type=str,
                        required=True)
    parser.add_argument("-chipseq_file", help="File with Super Enhancers", type=str,
                        required=True)
    parser.add_argument("-chipseq_promoter_file", help="File with Promoters", type=str,
                        default="")
    parser.add_argument("-chipseq_promoter_condition", help="Conditions of differencial enrichment", type=str,
                        default="tumor")
    parser.add_argument("-chipseq_promoter_condition_min_samples", help="Min number of samples that fit condition", type=int,
                        default=1)
    parser.add_argument("-chipseq_promoter_condition_max_samples", help="Max number of samples that fit condition", type=int,
                        default=5)

    parser.add_argument("-rna_enrich_file", help="File with RNASeq results", type=str,
                        default="")
    parser.add_argument("-rna_celllines_file", help="File with RNASeq results for cell lines", type=str,
                            default="")
    parser.add_argument("-lfc_thresh", help="Treshold odf absolute value of Log Fold Change (which genes are significant)", type=float,
                        default=1.0)
    parser.add_argument("-pvalue_thresh",
                        help="Treshold of pValue (which genes are significant)", type=float,
                        default=0.05)
    parser.add_argument("-gene_name_col",
                        help="Column name if RNA file with gene name", type=str,
                        required=True)

    parser.add_argument("-output_file", help="Output file", type=str,
                        required=True)
    parser.add_argument("-output_pdf", help="Output pdf file", type=str,
                        default="")
    args = parser.parse_args()
    dirname = args.indir
    chipseq_se_file = args.chipseq_file
    chipseq_promoter_file = args.chipseq_promoter_file
    promoter_condition = args.chipseq_promoter_condition
    rna_enrich_file = args.rna_enrich_file
    rna_celllines_file = args.rna_celllines_file.split(',')
    output_file = args.output_file
    output_pdf = args.output_pdf
    logfoldchangeCol = 'log2FoldChange'
    pvalueCol = 'padj'
    geneNameCol = args.gene_name_col
    if len(rna_enrich_file) == 0 and len(rna_celllines_file) == 0 and len(chipseq_promoter_file) == 0:
        raise RuntimeError("At least one of rna_enrich_file,rna_celllines_file,chipseq_promoter_file should be set")

    if len(output_pdf) > 0:
        pp = PdfPages(output_pdf)
    df_se = pd.read_csv(dirname + chipseq_se_file)
    spec_stat = compute_specifity_statistics(df_se)
    fig1,ax1=plot_specifity_statistics(spec_stat,title="Super Enhancers")
    if len(output_pdf) > 0:
        pp.savefig(fig1)
    spec_stat_genes= compute_specifity_statistics_genes(df_se)
    fig2,ax2=plot_specifity_statistics(spec_stat_genes,title="Super Enhancers controled Genes")
    if len(output_pdf) > 0:
        pp.savefig(fig2)
    if len(chipseq_promoter_file) > 0:
        df_promoter = pd.read_csv(dirname + chipseq_promoter_file)
        spec_stat = compute_specifity_statistics(df_promoter)
        fig3,ax3=plot_specifity_statistics(spec_stat,title="Promoters")
        if len(output_pdf) > 0:
            pp.savefig(fig3)
        enrichment_promotor = df_promoter['nearest_gene'].to_list()
        if len(promoter_condition)>0:
            enrichment_cond_promotor = {}
            counts_list = list(np.linspace(args.chipseq_promoter_condition_min_samples,args.chipseq_promoter_condition_max_samples,args.chipseq_promoter_condition_max_samples-args.chipseq_promoter_condition_min_samples+1,dtype=int))
            for counts in counts_list:
                enrichment_cond_promotor[counts] = \
                df_promoter[(df_promoter['Conditions'] == promoter_condition) & (df_promoter['Count'] >= counts)]['nearest_gene'].to_list()

    if len(rna_enrich_file) > 0:
        df_rna = pd.read_csv(dirname + rna_enrich_file)
        diff_expressed = df_rna[df_rna[pvalueCol] > args.pvalue_thresh]
        upreg = diff_expressed[np.abs(diff_expressed[logfoldchangeCol]) > args.lfc_thresh][geneNameCol].to_list()
    else:
        diff_expressed = None
    if len(rna_celllines_file) > 0:
        df_rna_celllines = {}
        cellines_upreg = {}
        for rna_cellline_file in rna_celllines_file:
            df_rna_celllines[rna_cellline_file] = pd.read_csv(dirname + rna_cellline_file)
            diff_expressed = df_rna_celllines[rna_cellline_file][df_rna_celllines[rna_cellline_file][pvalueCol] > args.pvalue_thresh]
            cellines_upreg[rna_cellline_file] = diff_expressed[np.abs(diff_expressed[logfoldchangeCol]) > args.lfc_thresh][geneNameCol].to_list()
    else:
        cellines_upreg = {}
    #print(df_rna)
    #add column genes_promotor to df_se
    df_se['Genes_promotor'] = None
    df_se['Genes_tumor_promotor'] = None
    df_se['Genes_upreg'] = None


    for index,row in tqdm.tqdm(df_se.iterrows(),desc="Iterating over ChipSeqs regions",total=df_se.shape[0],unit="SE(s)"):
        #print(row['Genes'])
        if row['Genes'] == [] or row['Genes'] is np.nan:
            continue
        list_of_genes = row['Genes'].split(',')
        # Find intersection of list_of_genes and enrichment_promotor
        genes_promotor = list(set(list_of_genes) & set(enrichment_promotor))
        df_se.at[index, f'Genes_promotor'] = ','.join(genes_promotor)
        df_se.at[index, f'Genes_upreg'] = ','.join(set(list_of_genes) & set(upreg))
        for rna_cellline_file in rna_celllines_file:
            df_se.at[index, f'Genes_{rna_cellline_file}'] = ','.join(set(list_of_genes) & set(cellines_upreg[rna_cellline_file]))
        for counts in counts_list:
            genes_cond_promotor = list(set(genes_promotor) & set(enrichment_cond_promotor[counts]))
            genes_cond_promotor_upreg = list(set(genes_cond_promotor) & set(upreg))
            df_se.at[index, f'Genes_{promoter_condition}_promotor_{counts}'] = ','.join(genes_cond_promotor)
            df_se.at[index, f'Genes_{promoter_condition}_promotor_upreg_{counts}'] = ','.join(genes_cond_promotor_upreg)
    gtp = {}
    gtpe = {}
    for counts in counts_list:
        gtp[counts] = len(df_se[f'Genes_{promoter_condition}_promotor_{counts}'].unique().tolist())
        gtpe[counts] = len(df_se[f'Genes_{promoter_condition}_promotor_upreg_{counts}'].unique().tolist())
    #bar plot of gtp and gtpe
    fig, ax = plt.subplots()
    ax.bar(gtp.keys(), gtp.values(), width=0.4, label=f'Genes with {promoter_condition} specified \n enrichment in promotor')
    ax.bar(gtpe.keys(), gtpe.values(), width=0.4, label=f'Genes with {promoter_condition} specified \n enrichment in promotor and upregulation')
    ax.set_xticks(counts_list)
    ax.set_xticklabels(counts_list)
    ax.set_xlabel('Number of samples')
    ax.set_ylabel('Number of genes')
    ax.grid()
    ax.legend()
    if len(output_pdf) > 0:
        pp.savefig(fig)
    df_se.to_csv(dirname + output_file, index=False)
    #print(df_se)
    #print(df_promoter)
    #print(df_rna)
    if len(output_pdf) > 0:
        pp.close()