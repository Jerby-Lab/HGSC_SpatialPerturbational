#### Spatial mapping of tubo-ovarian cancer and Perturb-seq reveal immune evasion mechanisms ####
# This resource provides the code developed in the study of Yeh, Aguirre, Laveroni et al.
# "Spatial mapping of tubo-ovarian cancer and Perturb-seq reveal immune evasion mechanisms.".
# It provides the data and interim files (see below) to reproduces the main figures and tables in the study and serves as a framework to perform integrated analyses with single cell spatial transcriptomics and perturb-seq data.

### Data:
# The raw and processed data will be provided on the Gene Expression Omnibus and Zenodo when the study is officially published.

### Quick start:
# To reproduce the figures and tables in Yeh, Aguirre, Laveroni et al. download HGSC_Data.zip, follow the steps below:
# 1. Clone this repository by executing git clone https://github.com/Jerby-Lab/HGSC_SpatialPerturbational/.
# 2. Download HGSC_Data.zip, place and unzip it in the HGSC_SpatialPerturbational repository. This will result in two additional direcotries: (1) Data, and (2) Results.
# 3. In line 71 of Code/HGSCgit_main.R please replace /path/to/your/local/clone/of/this/repo/ with the local path to this repository as it appears in your local machine.
# 4. In RStudio or within a script, run: HGSC_main()

tictoc::tic()

HGSC_main()

tictoc::toc()


# You can find the resulting figure panels in Figures, and the tables in Tables.
# The estimated run time is ~20 minutes on the platform
# aarch64-apple-darwin20 (64-bit) under macOS Ventura 13.4.1 operating system with 10 cores.

### Software Requirements:
# R (tested on version 4.2.3 -- "Shortstop Beagle")
# R libraries: survminer (v0.4.9), pheatmap (v1.0.12), stringr (v1.5.0), pROC (v1.18.4), argparse (v2.2.2), tidyr (v1.3.0), dplyr (v1.1.2), UpSetR (v1.4.0), RColorBrewer (v1.1-3), openxlsx (v4.2.5.2), e1071 (v1.7-13), heatmap3 (v1.1.9), devtools (v2.4.5), usethis (v2.2.2), lmerTest (v3.1-3), lme4 (v1.1-34), Rtsne (v0.16), Matrix (v1.5-4.1), plotrix (v3.8-2), reshape2 (v1.4.4), plyr (v1.8.8), EnhancedVolcano (v1.16.0), ggrepel (v0.9.3), ggpubr (v0.6.0), gplots (v3.1.3), tsne (v0.1-3.1), ROCR (v1.0-11), ppcor (v1.1), nnet (v7.3-19), ggplot2 (v3.4.2), MASS (v7.3-60), mixtools (v2.0.0), rms (v6.7-0), Hmisc (v5.1-0), survival (v3.5-5), EBImage (v4.40.1), SeuratObject (v4.1.3), Seurat (v4.3.0.1), cowplot (v1.1.1), beanplot (v1.3.1).

### License: BSD 3-Clause License.

### Citation:
# Yeh, Aguirre, Laveroni et al. Spatial mapping of tubo-ovarian cancer and Perturb-seq reveal immune evasion mechanisms.




