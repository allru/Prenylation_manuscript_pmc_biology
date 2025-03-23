# Uncovering Protein Prenylation in Th1 Cells


---

## Overview

This repository contains the code and data used in the analysis for the publication  
**"Uncovering Protein Prenylation in Th1 Cells: Novel Prenylation Sites and Insights into Statin and Farnesyltransferase Inhibition"**. The study maps known and novel prenylation events in human Th1 cells using proteomics and click-chemistry-based enrichment.

The computational analyses:
- Compare known and newly identified prenylated proteins
- Characterize canonical and non-canonical prenylation motifs
- Explore internal vs. C-terminal prenylation
- Examine the effects of farnesyltransferase inhibition (FTI)

---

## Repository Structure

- `0_functions`  
  Helper Python functions used across notebooks for motif extraction, formatting, and plotting.

- `1_Known_prenylated_proteins`  
  UniProt queries, literature-based curation, and motif annotation for proteins previously known to be prenylated.

- `2_Background_for_Exp1&Exp2`  
  Generates a reference background dataset of non-prenylated proteins. Used for motif discovery and enrichment analysis.

- `3_Union_Exp1&Exp2`  
  Combined results from Experiments 1 & 2; structural motif analysis and visualization.

- `4_Exp3_FTI_treatment`  
  Differential analysis of farnesylation after FTI treatment.

- `master_df.csv`  
  Master table summarizing proteins, annotations, motif classifications, and experiment-specific detection.

- `environment.yml`  
  Conda environment file listing dependencies (Python ≥3.10, pandas, biopython, etc.)

---

## Getting Started

Clone the repository:

```bash
git clone https://github.com/allru/Prenylation_manuscript_pmc_biology.git
cd Prenylation_manuscript_pmc_biology
```

Create the environment and activate:

```bash
conda env create -f environment.yml
conda activate prenylation
```

You can now run the notebooks inside each numbered subdirectory.

---

## Data Access

Mass spectrometry datasets on the MassIVE repository:
  - [Main dataset (MSV000095575)](ftp://massive.ucsd.edu/v08/MSV000095575/)
  - [FTI samples (MSV000096890)](ftp://massive.ucsd.edu/v08/MSV000096890/)

---

## Citation

If you use this code or data, please cite:

>  


---

## Contact

jana.koch@med.uni-augsburg.de  
alessandra.ruggia@library.ethz.ch

---


