#### Mapping ovarian cancer spatial organization uncovers immune evasion drivers at the genetic, cellular, and tissue level ####
# This resource provides the code developed in the study of Yeh, Aguirre, Laveroni (2024) et al.
# "Mapping ovarian cancer spatial organization uncovers immune evasion at the genetic, cellular, and tissue level".
# It provides the data and interim files (see below) to reproduce the main figures and tables in the study and serves as a framework to perform integrated analyses with single cell spatial transcriptomics and perturb-seq data.

### Data:
# The raw and processed data will be provided on the Gene Expression Omnibus and Zenodo when the study is officially published.

### Quick start:
# Follow the steps below to reproduce the figures and tables in Yeh, Aguirre, Laveroni et al (2024).:
# 1. Clone this repository by executing git clone https://github.com/Jerby-Lab/HGSC_SpatialPerturbational/.
# 2. Download Yeh2024.zip from Zenodo, place and unzip it in the HGSC_SpatialPerturbational repository. This will result in two additional direcotries: (1) Data, and (2) Results.
# 3. In line 95 and 96 of Code/run/HGSCgit_main.R please replace /path/to/your/local/clone/of/this/repo/ with the local path to this repository as it appears in your local machine.
# 4. In RStudio or within a script, run: HGSC_main()
# You can find the resulting figure panels in Figures/. 

tictoc::tic()

HGSC_main()

tictoc::toc()

### Software Requirements:
# R (tested on version 4.3.2 -- "Eye Holes")
### License: BSD 3-Clause License.

### Citation:
# Yeh, Aguirre, Laveroni et al. Mapping ovarian cancer spatial organization uncovers immune evasion at the genetic, cellular, and tissue level (2024). 




