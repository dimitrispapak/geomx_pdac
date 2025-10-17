# Analysis of GeoMx dataset from PDAC patients in Gustave Roussy.
## main features
- Descriptive plots and unsupervized analysis like UMAP and PCA plots
- Enrichment of the annotations with spatial component from image analysis.
- Differential Expression based on various studied features in our dataset.
- Deconvolution of the macrophage RNA segments using a scRNA reference atlas. 
- Study of the correlation of proportions of macrophage subtypes and clinical features.

## Prerequisites
- dcc files located in the dccs folder
- GeoMx Human Whole Transcriptome Atlas "Hs_R_NGS_WTA_v1.0.pkc"
- Sample annotations in excel format.
- Folder with Images of slides.

## Steps:
1. add_proximity_score.py to calculate the proximity of each segment to each neighbors for each ROI based on the slide images.
2. geomx_qc.ipynb -> Quality Control, Filtering and Creating of the GeoMx dataset.
3. descriptive_plots.ipynb -> Overall descriptive plots of the dataset 
4. differential_expression.ipynb -> Do differential expression analysis based on two model Linear mixed models and a generalized log-linear model  
5. deconvolution_macrophages_segments.R -> Use BisquRNA library with scRNA atlas momac-verse for macrophages to calculate the proportions of macrophages subtypes in macrophages segments.
6. correlation_analysis.ipynb -> Aggregate the proximity scores and macrophage proportions and clinical data to search for significant correlatios among them. 
