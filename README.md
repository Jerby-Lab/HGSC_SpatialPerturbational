# Mapping ovarian cancer spatial organization uncovers immune evasion drivers at the genetic, cellular, and tissue level

This resource provides the code developed in the study of _Yeh, Aguirre, Laveroni _et al_. **_"Mapping ovarian cancer spatial organization uncovers immune evasion drivers at the genetic, cellular, and tissue level"_**. It provides the code to reproduce the main figures in this study and serves as a framework to perform integrated analyses with single cell spatial transcriptomics and perturb-seq data. 

## **Data**

The raw and processed data will be provided on the Gene Expression Omnibus and Zenodo **when the study is officially published**. The data will also be available for download and interactive exploration on the Single Cell Portal.

## **Quick start**
Follow the steps below to reproduce the figures and tables in Yeh, Aguirre, Laveroni _et al.: 

1. Clone this repository by executing `git clone https://github.com/Jerby-Lab/HGSC_SpatialPerturbational/`. 

2. Download _**HGSC_Data.zip**_, place and unzip it in the `HGSC_SpatialPerturbational` repository. This will result in two additional direcotries: (1) Data, and (2) Results.

3. In line 71 of `Code/HGSCgit_main.R` please replace `/path/to/your/local/clone/of/this/repo/` with the local path to this repository as it appears in your local machine. 

4. In RStudio or within a script, run: ```HGSC_main()```

You can find the resulting figure panels in _Figures_. The estimated run time is ~20 minutes on the platform aarch64-apple-darwin20 (64-bit) under macOS Sonoma 14.4.1 operating system with 10 cores. 

## Citation

Yeh, Aguirre, Laveroni _et al._ _**Mapping ovarian cancer spatial organization uncovers immune evasion drivers at the genetic, cellular, and tissue level**_.

## **Software Requirements**

* R (tested on version 4.3.2 -- "Eye Holes")
* R libraries (see [HGSC_lib.R](https://github.com/Jerby-Lab/HGSC_SpatialPerturbational/blob/main/Code/HGSC_lib.R)): survminer (v0.4.9), pheatmap (v1.0.12), stringr (v1.5.0), pROC (v1.18.4), argparse (v2.2.2), tidyr (v1.3.0), dplyr (v1.1.2), UpSetR (v1.4.0), RColorBrewer (v1.1-3), openxlsx (v4.2.5.2), e1071 (v1.7-13), heatmap3 (v1.1.9), devtools (v2.4.5), usethis (v2.2.2), lmerTest (v3.1-3), lme4 (v1.1-34), Rtsne (v0.16), Matrix (v1.5-4.1), plotrix (v3.8-2), reshape2 (v1.4.4), plyr (v1.8.8), EnhancedVolcano (v1.16.0), ggrepel (v0.9.3), ggpubr (v0.6.0), gplots (v3.1.3), tsne (v0.1-3.1), ROCR (v1.0-11), ppcor (v1.1), nnet (v7.3-19), ggplot2 (v3.4.2), MASS (v7.3-60), mixtools (v2.0.0), rms (v6.7-0), Hmisc (v5.1-0), survival (v3.5-5), EBImage (v4.40.1), SeuratObject (v4.1.3), Seurat (v4.3.0.1), cowplot (v1.1.1), beanplot (v1.3.1).
  

## License 

BSD 3-Clause License provided ([here](https://github.com/Jerby-Lab/HGSC_SpatialPerturbational/blob/main/LICENSE)).

## Contact 

For any inquiries on this repository please feel free to post these under "issues" or reach out to Christine Yiwen Yeh ([cyyeh@stanford.edu](cyyeh@stanford.edu)) and Livnat Jerby ([ljerby@stanford.edu](ljerby@stanford.edu)).
