# chipover
A command-line toolkit for summarizing, analyzing, and interpreting ChIP-Seq and RNA-Seq experiments, designed for enhancer and super-enhancer region analysis in tumor-normal pairs, with applications in differential enrichment and multi-sample domain analysis. 

## Processing overview

<img src="img/analysis_graph.png">

## Features 

### 1. Creation of Common Summary Tables
- Compiles ChIP-Seq results into a unified summary table, allowing easy comparison of peaks across multiple samples.
(**chipsummary tool**)

### 2. Combining Results Across Samples
- Uses combinatorial logic to analyze and summarize overlapping and unique peaks across a set of samples and conditions 
(**chipsummary,chipover tools**)

### 3. Integration with RNA-Seq Data
- Provides joint analysis of ChIP-Seq and RNA-Seq data to correlate enrichment patterns with gene expression.
(**chiprnacombiner tool**)

### 4. Designed for Super Enhancer analysis in tumor-sampels pairs
- Provides joint analysis of Super Enhancer resions \ differently enriched promoter regions as well gene expression data
(**chiprnacombiner tool**)

### 5. Provide a list of genes  
- Provide list of genes around differentely enriched regions, near to differentely enriched promoters and upregulated
(**chipsummary tool**)

### 6. Compare results with other set of regions and list of genes   
- Allow compare of derived differentialy enriched domains with other set of regions (for example differentialy enriched regions in cell lines)
  and compare list of genes with additional list of genes (again, for example, list of genes pverexpressed in cell lines)
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
