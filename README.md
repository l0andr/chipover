# chipover (Version: 0.6.0)
A command-line toolkit for summarizing, analyzing, and interpreting ChIP-Seq and RNA-Seq experiments, designed for enhancer and super-enhancer region analysis in tumor-normal pairs, with applications in differential enrichment and multi-sample domain analysis. 

## Processing overview

<img src="img/analysis_graph.png">

## Features 

- **Combining results across samples**  
  Uses combinatorial logic to summarize overlapping peaks across a set of samples and conditions. Computes basic statistics on the number of peaks with different enrichment conditions across all samples, with the option to consider additional intersections (e.g., with cell lines). Outputs results in BED format, allowing users to create outputs such as: "BED file with regions differentially enriched in tumor samples and also present in a given cell line, marked as Super Enhancer," or for promoters, normal samples, HPV-related samples, and so on.  
  (**chipover tool**)

- **Creation of common summary tables**  
  Compiles ChIP-Seq results into a unified table, allowing easy comparison of peaks across multiple samples, integrating everything, including intersections with cell lines, DiffBind results, etc.  
  (**chipsummary tool**)

- **Integration with RNA-Seq data**  
  Provides joint analysis of ChIP-Seq and RNA-Seq data to correlate enrichment patterns with gene expression.  
  (**chiprnacombiner tool**)

- **Designed for super enhancer analysis in tumor-sample pairs**  
  Enables joint analysis of Super Enhancer regions, differentially enriched promoter regions, and gene expression data.  
  (**chiprnacombiner tool**)

- **Provide a lists of genes**  
  Generates a list of genes around differentially enriched regions, near differentially enriched promoters, and upregulated genes.  
  (**chipsummary tool**)

- **Compare results with other sets of regions and gene lists**  
  Allows comparison of derived differentially enriched domains with other sets of regions (e.g., differentially enriched regions in cell lines) and compares lists of genes with additional gene lists (e.g., genes overexpressed in cell lines).  
  (**chiprnacombiner tool**)
  
## Installation

Before installing, ensure you have the following:

- **Python 3.9 or later**: [Download Python](https://www.python.org/downloads/)
- **Git**: [Download Git](https://git-scm.com/downloads)
- **Bash-like shell**: A shell environment for running commands. This could be Git Bash (Windows), Terminal (macOS), or any Linux shell.

```bash
# Clone the repository
git clone https://github.com/l0andr/chipover.git
# Navigate to the directory
cd chipover
# Install dependencies
pip install -r requirements.txt
```

## Usage


<details>
<summary>Expected input file format and format of file with metadata </summary>

Expected format of input BED files:
| Chromosome | Start      | End        | Type | Score | Strand |
|------------|------------|------------|------|-------|--------|
| chr2       | 70023598   | 70143951   | SE   | 34323 | +      |
| chr10      | 61960248   | 62085351   | promotor| 34284 | +      |
| chr11      | 65403898   | 65455601   | enhancer| 28277 | +      |
| chr22      | 46034448   | 46088542   | other  | 27639 | +      |


Format of file with metadata:
| SampleID      | bamReads               | bamControl             | ControlID      | Peaks                  | Factor | Tissue   | Condition | Treatment | PeakCaller | Replicate |
|---------------|------------------------|------------------------|----------------|------------------------|--------|----------|-----------|-----------|------------|-----------|
| 30709t        | [path_to_file]         | [path_to_file]         | 30709t         | [path_to_file]         | h3k27ac| patient  | tumor     | noHPV     | bed        | 1         |
| 30866n        | [path_to_file]         | [path_to_file]         | 30866n         | [path_to_file]         | h3k27ac| patient  | normal    | HPV       | bed        | 1         |
| 30862n        | [path_to_file]         | [path_to_file]         | 30862n         | [path_to_file]         | h3k27ac| patient  | normal    | HPV       | bed        | 1         |
| hok16b        | [path_to_file]         | [path_to_file]         | hok16b         | [path_to_file]         | h3k27ac| cellline | normal    | HPV       | bed        | 1         |
| 30866t        | [path_to_file]         | [path_to_file]         | 30866t         | [path_to_file]         | h3k27ac| patient  | tumor     | HPV       | bed        | 1         |

Mandatory columns is SampleID,Peaks with paths to inpu BED files,Tissue,Condition. All other keep for compatability with DiffBind format
</details>

**Chipover**

ChipOver is a tool designed for analyzing and summarizing enhancer region overlaps across samples in ChIP-Seq data. It computes common domain structures by identifying overlaps in enhancer regions and generates output based on tissue-specific and sample-based metadata. It allows for flexible filtering of results and detailed logging options for various ChIP-Seq enrichment analyses.

Example:<br>
```
python chipover.py -indir $WORKDIR -metadata $WORKDIR/diffbind_enchancer.csv -outdir $WORKDIR/$OUTDIR -minoverlap 0
```

<details>
  <summary>chipover parameters description </summary>
  
- **`-indir`** (required):  
  Directory containing input BED files.

- **`-outdir`** (required):  
  Directory for saving output BED files.

- **`-minoverlap`** (optional, default = 0):  
  Minimum number of intersections for a domain to be considered.

- **`-filter_regions`** (optional):  
  BED file used for filtering results by specified regions.

- **`-tissue_spec`** (optional):  
  Specifies a tissue type for analysis; only data from this tissue will be used.

- **`-label`** (optional):  
  Additional label to append to output filenames.

- **`-metadata`** (optional):  
  CSV file containing sample metadata (required columns: `SampleID`, `Tissue`, `Condition`) - in general format of metadata file the same as used by DiffBind.

- **`-verbose`** (optional, default = 1):  
  Sets log level (0 = error, 1 = info, 2 = debug).
  
</details>

**Chipsummary**

ChipSummary is a tool for generating a comprehensive summary table from ChIP-Seq results. It identifies differentially enriched domains, extracts nearby genes, and genes on specified distance from domain,
also data about intersection with additional bed files and regions computed by DiffBind also can be added. 

Example:<br>
```
python chipsummary.py -indir $WORKDIR -intersection_dir $WORKDIR//celllines -metadata $WORKDIR/diffbind_enchancer.csv -outdir $WORKDIR -tissue_spec patient -genes_annotations $GENES_ANNOTATIONS -genes_biotype protein_coding -se_genes_range 1500000 -names_filter "SE"

or

python chipsummary.py -indir $WORKDIR -metadata $WORKDIR/diffbind_enchancer.csv -outdir $WORKDIR -tissue_spec patient -genes_annotations $GENES_ANNOTATIONS -genes_biotype protein_coding -se_genes_range 2000 -names_filter "promoter"

```

<details>
  <summary>chipsummary parameters description </summary>
  
- **`-indir`** (required):  
  Directory containing input BED files.

- **`-file_mask`** (optional, default = `'*.bed'`):  
  File pattern for input files.

- **`-out_perfix`** (optional):  
  Prefix for output filenames. The output format is `[out_perfix]_[type_of_regions]_report.csv`.

- **`-intersection_dir`** (optional):  
  Directory containing BED files for intersection, to be included in the summary.

- **`-names_filter`** (optional):  
  Comma-separated list of region names to include in the analysis. Default is all regions. 

- **`-outdir`** (required):  
  Directory for saving output files.

- **`-filter_regions`** (optional):  
  BED file for filtering results by specified regions.

- **`-tissue_spec`** (optional):  
  Restrict analysis to data from a specific tissue type.

- **`-genes_annotations`** (optional, default = `Homo_sapiens.GRCh38.110.genes.gtf`):  
  Path to a GTF file with gene annotations.

- **`-se_genes_range`** (optional, default = 0):  
  Size of the region around enhancers to extract nearby gene names.

- **`-genes_biotype`** (optional):  
  Biotype of genes to consider (e.g., `protein_coding`). Default includes all genes.

- **`-diffbind_regions`** (optional):  
  Path to a BED file containing DiffBind domains for analysis.

- **`-metadata`** (required):  
  Path to a CSV file in DiffBind metadata format with sample features.

- **`-cov_threshold`** (optional, default = 0.0):  
  Minimum coverage ratio at a peak for differential marking between samples. If specified, requires a "Coverage" column in the metadata.

- **`-verbose`** (optional, default = 1):  
  Log level: `0` for errors, `1` for info, `2` for debug.

  

</details>

**ChipRNAcombiner**

ChipRnaCombiner is a tool for next level generelesation or results and integration ChIP-Seq and RNA-Seq results to identify genes controlled by super-enhancers or enhancers, enriched in promoter regions, and showing upregulation in RNA-Seq data. 
It enables analysis of differential enrichment conditions across samples, providing a combined view of gene regulation. This tool take as input results of ChipSummary (up to two tables, second tables usually tables of promoters, but other usage scenarios also possible) and results of RNAseq analysis  

Example:<br>
```
python chiprnacombiner.py -indir $WORKDIR -chipseq_file se_report.csv -chipseq_promoter_file promoter_report.csv -rna_enrich_file prospective_DE_results005_1LFC.csv -rna_celllines_file "cell_line047_upreg.csv,cell_line090_upreg.csv,cell_line147_upreg.csv" -output_file "se_report_with_promoter_upreg.csv" -output_pdf "specifity_statistics.pdf" -gene_name_col "lbl"
```

<details>
  <summary>ChipRnaCombiner parameters description </summary>
  
- **`-indir`** (required):  
  Directory containing input files.

- **`-chipseq_file`** (required):  
  File with super-enhancer information.

- **`-chipseq_promoter_file`** (optional, default = ""):  
  File with promoter information.

- **`-chipseq_promoter_condition`** (optional, default = "tumor"):  
  Condition for differential enrichment analysis.

- **`-chipseq_promoter_condition_min_samples`** (optional, default = 1):  
  Minimum number of samples that must meet the specified condition.

- **`-chipseq_promoter_condition_max_samples`** (optional, default = 5):  
  Maximum number of samples that must meet the specified condition.

- **`-rna_enrich_file`** (optional, default = ""):  
  File containing RNA-Seq results.

- **`-rna_celllines_file`** (optional, default = ""):  
  File containing RNA-Seq results for cell lines.

- **`-lfc_thresh`** (optional, default = 1.0):  
  Threshold for the absolute value of log fold change to determine significant genes.

- **`-pvalue_thresh`** (optional, default = 0.05):  
  Threshold for p-value to determine significant genes.

- **`-gene_name_col`** (required):  
  Column name in the RNA-Seq file containing gene names.

- **`-output_file`** (required):  
  Output file path for the results.

- **`-output_pdf`** (optional, default = ""):  
  Output PDF file path for optional visualizations.

</details>

**ChipTableIntersect**

ChipTableIntersect is a tool designed to enhance the ChipSummary table by intersecting it with BED-like files. It enables users to add additional columns to the input table by computing intersections based on genomic regions and integrates relevant data from the intersecting table. This is useful for annotating or filtering ChIP-Seq results with data from other genomic studies.

Example:<br>
```bash
python chiptableintersect.py -input_table chip_summary.csv -intersect_table annotations.bed -output intersected_table.csv --chr_col "Chromosome" --start_col "Start" --end_col "End" --match_col "match_flag" --add_cols "GeneName,ExpressionLevel"
```
<details>
  <summary>ChipTableIntersect parameters description</summary>
  
- **`-input_table`** (required):  
  Path to the input table (e.g., ChipSummary table).

- **`-intersect_table`** (required):  
  Path to the table or file with genomic regions to intersect.

- **`-output`** (required):  
  Output file path for the resulting intersected table.

- **`--chr_col`** (optional, default = `"Chromosome"`):  
  Name of the column in both input and intersect tables containing chromosome information.

- **`--start_col`** (optional, default = `"Start"`):  
  Name of the column in both input and intersect tables containing the start position of regions.

- **`--end_col`** (optional, default = `"End"`):  
  Name of the column in both input and intersect tables containing the end position of regions.

- **`--match_col`** (optional, default = intersect table file name):  
  Name of the column to flag matches in the output table.

- **`--add_cols`** (optional, default = `""`):  
  List of columns (comma-separated) from the intersect table to add to the output when matches occur.

- **`-verbose`** (optional, default = `1`):  
  Log level:  
  - `0`: Error messages only  
  - `1`: Basic information messages (default)  
  - `2`: Debug-level messages  

</details>
