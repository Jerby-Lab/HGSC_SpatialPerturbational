# The code below provides a step-by-step guide
# to regenerate the figures of the manuscript
# in its final version as of June 2024. 

# Result Sections
#1 Data processing and annotations (Figure 1)
#2 Effector T Cells Infiltrate the Tumor (Figure 2)
#3 Malignant TIL (mTIL) program (Figure 3)
#4 mTIL program as predictor of clinical response (Figure 4)
#5 CNAs mapping to mTIL and TIL levels (Figure 5)
#6 Perturb-Seq meta-analyses (Figure 6)
#7 Perturb-Seq and genetic screen analyses (Figure 7)
#8 Validation of Genetic KO's (Figure 8)

HGSC_main<-function(){
  HGSC_main_download()
  HGSC_main_run()
  return()
}

HGSC_main_run<-function(){
  #Regenerate Figures from Yeh et al. 2024
  
  print("Figure 1")
  HGSC_Figure1_SpatiomolecularMapping(cell_2_rgb=cell_2_rgb)
  print("Figure 2")
  HGSC_Figure2_TIP_DF(r.smi = r.smi,rCD8.smi = rCD8.smi,r.xenium = r.xenium,
                      rTNK.xenium = rTNK.xenium,rCD8.xenium = rCD8.xenium,daf=daf)
  print("Figure 3")
  HGSC_Figure3_mTIL(r = r.smi,r1 = rMal.smi,Mtil.sig = Mtil.sig)
  print("Figure 4")
  HGSC_Figure4_ICB()
  print("Figure 5")
  HGSC_Figure5_CNAs(r1 = rCNA,rslts1 = rslts1_cna,rslts2 = rslts2_cna) 
  print("Figure 6")
  HGSC_Figure6_perturbMeta(rslts1 = rslts_prt1,rslts2 = rslts_prt2,rslts3 = rslts_prt3)
  print("Figure 7")
  HGSC_Figure7_perturbOC(r = r.prt,
                         rslts = rslts_prt4,
                         sig = Mtil.sig,
                         fitnessR = fitnessR)
  print("Figure 7")
  HGSC_Figure8_KOValidation()
  
  print("All Figures from Yeh et al., 2024 were reproduced and are avaliabile at the Figures directory.")
  return()
}


HGSC_main_download<-function(){
  overwrite.flag<<-T
  source("~/Desktop/R_code/6_Github/HGSC/HGSC_git_V2/Code/run/HGSCgit_main.R")
  if(!file.exists(get.file("Figures/Fig8c.pdf"))){HGSC_source()}
  
  #1 Download Presets 
  cell_2_rgb <<- readRDS(get.file("Data/Cell_Colors.rds"))
  
  #2 Download Discovery and Validation 1 datasets
  r.smi <<- readRDS(get.file("Data/ST_Discovery.rds"))
  r.xenium <<- readRDS(get.file("Data/ST_Validation1.rds"))
  r.smiT2 <<- readRDS(get.file("Data/ST_Test2_HGSC1.rds"))
  
  #3 Download T/NK and malignant spatial data.
  rCD8.smi <<- readRDS(get.file("Data/ST_Discovery_CD8.T.cells.rds"))
  rMal.smi <<- readRDS(get.file("Data/ST_Discovery_malignant.rds"))
  rTNK.xenium <<- readRDS(get.file("Data/ST_Validation1_TNK.cells.rds"))
  rCD8.xenium <<- readRDS(get.file("Data/ST_Validation1_CD8.T.cells.rds"))
  rMal.smiT2 <<-  readRDS(get.file("Data/ST_Test2_HGSC1_malignant.rds"))
  
  #4 Download Desmoplastic Fibroblast results
  daf <<- readRDS(get.file("Results/ST_DAF.rds"))
  
  #5 Download MTIL program.
  Mtil.sig <<- readRDS(get.file("/Results/ST_Mtil.sig.rds"))
  
  #6 Download CNA data and results. 
  rCNA<<-readRDS(get.file("Data/ST_CNAs.rds"))
  rslts1_cna<<-readRDS(get.file("Results/HGSC_CNAs.vs.TIL_ST.rds"))
  rslts2_cna<<-readRDS(get.file("Results/HGSC_CNAs.vs.TIL_TCGA.rds"))
  
  #6 Download Perturb-seq meta-analyses results.
  rslts_prt1<<-readRDS(get.file("Results/PerturbSeq_Mtil_K562_CRISPRa.rds"))
  rslts_prt2<<-readRDS(get.file("Results/PerturbSeq_Mtil_K562_CRISPRi.rds"))
  rslts_prt3<<-readRDS(get.file("Results/PerturbSeq_Mtil_RPE1_CRISPRi.rds"))
  
  #7 Download Perturb-Seq data and results.
  r.prt<<-readRDS(get.file("Data/PerturbSeq_TYKnu.rds"))
  fitnessR<<-readRDS(get.file("Results/SourceData_Fig7b.rds"))
  rslts_prt4<<-readRDS(get.file("Results/PerturbSeq_TYKnu_DEP.rds"))
  return()
  
}

get.file<-function(file1,code.flag = F){
  dir1 <- "/Volumes/ljerby/ljerby/HGSC/Zenodo/"
  # dir1 <- "/Volumes/ljerby/HGSC_Profiling/Data/Share/Zenodo_unzipped/HGSC_Data/"
  dir2 <- "~/Desktop/R_code/6_Github/HGSC/HGSC_git_V2/Code/"
  if(!code.flag){file1<-paste0(dir1,file1);return(file1)}
  file1<-paste0(dir2,file1)
  return(file1)
}

HGSC_source<-function(code.dir){
  print("Sourcing HGSC_SpatialPerturbational GitHub repository.")
  dir.create(get.file("Figures/"))
  files <- list.files(get.file("",code.flag = T),
                      include.dirs = F,
                      recursive = T,
                      pattern = ".R",
                      full.names = T)
  files<-files[!grepl("private",files)]
  lapply(files, source)
  return()
}
