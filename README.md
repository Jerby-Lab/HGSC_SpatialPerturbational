# Spatial mapping of tubo-ovarian cancer and Perturb-seq reveal immune evasion mechanisms

This resource provides the code developed in the study of Yeh, Aguirre, Laveroni _et al_. **_"Spatial mapping of tubo-ovarian cancer and Perturb-seq reveal immune evasion mechanisms."_**. It reproduces the main figures and tables in the study and serves as a framework to perform integreated analyses with single cell spatial transcriptomics and perturb-seq data. 

## **Software Requirements**

* R (tested on version 4.2.3 -- "Shortstop Beagle")
* R libraries: survminer (v0.4.9), pheatmap (v1.0.12), stringr (v1.5.0), pROC (v1.18.4), argparse (v2.2.2), tidyr (v1.3.0), dplyr (v1.1.2), UpSetR (v1.4.0), RColorBrewer (v1.1-3), openxlsx (v4.2.5.2), e1071 (v1.7-13), heatmap3 (v1.1.9), devtools (v2.4.5), usethis (v2.2.2), lmerTest (v3.1-3), lme4 (v1.1-34), Rtsne (v0.16), Matrix (v1.5-4.1), plotrix (v3.8-2), reshape2 (v1.4.4), plyr (v1.8.8), EnhancedVolcano (v1.16.0), ggrepel (v0.9.3), ggpubr (v0.6.0), gplots (v3.1.3), tsne (v0.1-3.1), ROCR (v1.0-11), ppcor (v1.1), nnet (v7.3-19), ggplot2 (v3.4.2), MASS (v7.3-60), mixtools (v2.0.0), rms (v6.7-0), Hmisc (v5.1-0), survival (v3.5-5), EBImage (v4.40.1), SeuratObject (v4.1.3), Seurat (v4.3.0.1), cowplot (v1.1.1), beanplot (v1.3.1)

## **Data**

The raw and processed data will be provided on the Gene Expression Omnibus and Zenodo when the study is officially published. 

## **Quick start**
To reproduce the figures and tables in Yeh, Aguirre, Laveroni _et al._ download _**HGSC_Data.zip**_, follow the steps below: 

1. Clone this repository by executing `git clone https://github.com/Jerby-Lab/HGSC_SpatialPerturbational/`. 

2. Download _**HGSC_Data.zip**_. Unzip the file and move the resulting files into the _Data_ directory within the `HGSC_SpatialPerturbational` repository.

3. Download _**HGSC_Results.zip**_. Unzip the file and move the resulting files into the _Results_ directory within the `HGSC_SpatialPerturbational` repository.

4. In Rstudio or within a script, run:
```
files <- list.files(paste0("/path/to/your/local/clone/of/this/repo/Code"),
                  include.dirs = F,
                  pattern = ".R",
                  full.names = T)
lapply(files, source)
HGSC_main()
```
5. You can find the resulting figure panels in _Figures_, and the tables in _Tables_. The estimated run time is ~20 minutes on a laptop. 

# General notes

For computational tractibility, intermediate results from this study have been provided in the ```Results``` directory. Therefore, the scope and goal of this repository is to provide higher-level analysis frameworks rather than code to reproduce the study from raw data. 

## Citation

Yeh, Aguirre, Laveroni _et al._ _**Spatial mapping of tubo-ovarian cancer and Perturb-seq reveal immune evasion mechanisms**_. (2023)

## Contact 

For any inquiries on this repository please feel free to reach out to Christine Yiwen Yeh ([cyyeh@stanford.edu](cyyeh@stanford.edu)) and/or Livnat Jerby ([ljerby@stanford.edu](ljerby@stanford.edu))

## License 


