#### Mapping ovarian cancer spatial organization uncovers immune evasion drivers at the genetic, cellular, and tissue level ####
# This resource provides the code developed in the study of Yeh, Aguirre, Laveroni et al.
# "Mapping ovarian cancer spatial organization uncovers immune evasion drivers at the genetic, cellular, and tissue level".
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
# R (tested on version 4.3.2 -- "Eye Holes")
### License: BSD 3-Clause License.

### Citation:
# Yeh, Aguirre, Laveroni et al. Mapping ovarian cancer spatial organization uncovers immune evasion drivers at the genetic, cellular, and tissue level.




