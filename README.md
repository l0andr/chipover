# chipover
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
