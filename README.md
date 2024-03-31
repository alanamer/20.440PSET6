 # iPSC-derived vs. Primary Cardiomyocytes Transcriptome Pathway Analysis

This repo contains the code to reproduce all of the analyses in "Meta analysis of microbiome studies identifies shared and disease-specific patterns", Duvallet et al. 2017.

This repo contains the code to generate two graphs of differentially expressed genes in iPSC-derived vs. primary cardiomyocytes using transcriptome data from GSE146096 (Primary cardiomyocytes) and GSE226159 (iPSC-derived cardiomyocytes) from the NIH Gene Expression Omnibus (GEO) database. 

RNA count normalization and differential expression analysis was performed using the DESeq2 implementation in Python (PyDESeq2), and we identified differentially expressed genes based on the criteria of adjusted p-value (FDR) > 0.05 and |log2 fold change| > 2. These differentially expressed genes were then passed to the package PyKEGG (https://pypi.org/project/pykegg/) which allows visualization of KEGG information using a network approach. 

# Motivation

In vitro human organ models offer the potential for a more cost-effective and efficient drug development pipeline, yet there is still a critical lack of understanding of differences between cells used in these models derived from induced pluripotent stem cells (iPSCs) compared to primary cells. Elucidating the transcriptomic distinctions between cells of different developmental origins is imperative to identify the optimal method for mimicking primary cell responses, substantiate the equivalence of stem cell-derived cells, and inform future in vitro model development. There is a lack of in depth means of verifying stem cell derived lineages’ fidelity to primary tissues. Current validation methods include expensive functional assays and/or looking for expression of a few select markers which can be co-expressed across cell types of similar lineages.  We seek to provide an open source, standardized framework for transcriptomic analysis to allow researchers to validate their differentiation protocols more comprehensively against primary cell types. This repo offers one step in that process where differentially expressed genes can be viewed within a relevant pathway for the cell type.  

# Directory Structure and Reproducing analyses

The repo contains three folders within main: data_files, code, and figures.

data_files contains the tsv files of expression data used to generate the figures. code contains the jupyter notebook file that can be used to create the figures. figures contains the generated figures. 

## Installing

To re-make all of the analyses, you'll first need to install the required modules.

Please do this within a Python 3 (or latest) environment. You will need to import the following into your python file for the analysis.

```
import requests_cache
import numpy as np
from PIL import Image
import pykegg
import matplotlib as mpl
```

#### data

All data-related files are in `data/`:

#### code

All of the code is in the `code/` folder:

#### figures 

The generated figures are in the `figures` folder

