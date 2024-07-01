# Mapping ovarian cancer spatial organization uncovers immune evasion drivers at the genetic, cellular, and tissue level

This resource provides the code developed in the study of Yeh, Aguirre, Laveroni _et al_. **"Mapping ovarian cancer spatial organization uncovers immune evasion drivers at the genetic, cellular, and tissue level"**. (2024). _**Nature Immunology (in press)**_. This code reproduces the main figures in this study and serves as a framework to perform integrated analyses with single cell spatial transcriptomics and perturb-seq data. 

## **Data**

Raw and processed data are provided from the Gene Expression Omnibus; processed and intermim .rds objects are available from [Zenodo](https://zenodo.org/records/11206564). The data is also available for download and interactive exploration on the [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP2640/hgsc-spatial-cohort-discovery-dataset#study-summary).

## **Quick start**
Follow the steps below to reproduce the figures and tables in Yeh, Aguirre, Laveroni _et al._: 

1. Clone this repository by executing `git clone https://github.com/Jerby-Lab/HGSC_SpatialPerturbational/`. 

2. Download _**HGSC_Data.zip**_ from [Zenodo](https://zenodo.org/records/11206564), place and unzip it in the `HGSC_SpatialPerturbational` repository. This will result in two additional direcotries: (1) Data, and (2) Results.

3. In line 71 of `Code/run/HGSCgit_main.R` please replace `/path/to/your/local/clone/of/this/repo/` with the local path to this repository as it appears in your local machine. 

4. In RStudio or within a script, run: ```HGSC_main()```

You can find the resulting figure panels in _Figures_. The estimated run time is ~20 minutes on the platform aarch64-apple-darwin20 (64-bit) under macOS 14.5 operating system. 

## Citation

Yeh, Aguirre, Laveroni et al. "Mapping ovarian cancer spatial organization uncovers immune evasion drivers at the genetic, cellular, and tissue level". (2024). Nature Immunology (in press).

## **Software Requirements**

* R (tested on version 4.3.2 -- "Eye Holes")
* R libraries (see [HGSC_lib.R](https://github.com/Jerby-Lab/HGSC_SpatialPerturbational/blob/main/Code/HGSC_lib.R)): patchwork (1.2.0), tibble (3.2.1), forestplot (3.1.3), abind (1.4-5), checkmate (2.3.1), tidyr (1.3.1), argparse (2.2.3), pROC (1.18.5), stringr (1.5.1), pheatmap (1.0.12), survminer (0.4.9), dplyr (1.1.4), UpSetR (1.4.0), RColorBrewer (1.1-3), openxlsx (4.2.5.2), e1071 (1.7-14), heatmap3 (1.1.9), devtools (2.4.5), usethis (2.2.3), lmerTest (3.1-3), lme4 (1.1-35.4), Rtsne (0.17), Matrix (1.6-5), plotrix (3.8-4), reshape2 (1.4.4), plyr (1.8.9), EnhancedVolcano (1.20.0), ggrepel (0.9.5), ggpubr (0.6.0), gplots (3.1.3.1), tsne (0.1-3.1), ROCR (1.0-11), ppcor (1.1), nnet (7.3-19), ggplot2 (3.5.1), MASS (7.3-60.0.1), mixtools (2.0.0), rms (6.8-1), Hmisc (5.1-3), survival (3.7-0), EBImage (4.44.0), Seurat (5.1.0), SeuratObject (5.0.2), sp (2.1-4), cowplot (1.1.3), beanplot (1.3.1)

## License 

BSD 3-Clause License provided ([here](https://github.com/Jerby-Lab/HGSC_SpatialPerturbational/blob/main/LICENSE)).

## Contact 

For any inquiries on this repository please feel free to post these under "issues" or reach out to Christine Yiwen Yeh ([cyyeh@stanford.edu](cyyeh@stanford.edu)) and Livnat Jerby ([ljerby@stanford.edu](ljerby@stanford.edu)).
