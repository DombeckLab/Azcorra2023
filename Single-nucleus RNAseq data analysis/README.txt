"# Azcorra2023 - Single-nucleus RNAseq data analysis" 

Azcorra, Gaertner et al. Nat Neuro 2023 Custom R scripts.

In our manuscript (Azcorra, Gaertner et al. Nat Neuro 2023), we analyzed a large single-nucleus RNAseq dataset of midbrain dopamine neurons in order to establish subtypes, look at DEGs, and explore relationships between clusters. In the associated R script file, custom code necessary for the RNAseq analyses in the paper are provided.

Access to the datasets needed to generate the results is available at GEO (GSE222558). Libraries required to run the analysis are included at the top of the script file.

To get started, install R (multiple versions of R were used over the course of the analysis; versions used for each package/analysis step involved in the manuscript are listed in Key Resource table). In addition, we recommend that you please install RStudio to open the script file and run the analyses.

Featured in the scripts file are, in order:
1) A walkthrough of the standard Seurat pipeline, with the settings used in our analyses
2) Custom scripts for functions we developed, including: Bootstrapping the difference in correlations of gene expression between two pairs of clusters; making scatter plots for correlation of gene expression for a pair of clusters; Cluster stability calculations (several helper functions needed for this analysis are included).





