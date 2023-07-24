
# The code below provides a step-by-step guide
# to regenerate all the figures of the manuscript
# in its final version as of 07/24/2023.

# Result Sections
#1 Data processing and annotations (Figure 1, Tables S1-S3)
#2 Malignant drift and TME reprogramming (Figure 2, Table S4)
#3 Immune Tumor Infiltration programs (Figure 3, Table S5)
#4 Malignant TIL (mTIL) program (Figure 4, Table S6)
#5 CNAs mapping to mTIL and TIL levels (Figure 5)
#6 Perturb-Seq meta-analyses (Figure 6, Table S6)
#7 Perturb-Seq and genetic screen analyses (Figure 7, Table S7)

HGSC_main<-function(path2git){
  files<-list.files(paste0(path2git,"/Code"),include.dirs = F,pattern = ".R",full.names = T)
  lapply(files, source)
  
  #1 Download Tumor Infiltration Programs (TIP)
  r.smi<-readRDS(get.file("Data/SMI_data.rds"))
  rTNK.smi<-HGSC_SMI.process.CD8()
  rmal.smi<-HGSC_SMI.process.mal()
  r.xenium<-readRDS(get.file("/Data/Xenium_data.rds"))
  rTNK.xenium<-HGSC_Xenium.process.CD8.NK(r = r.xenium)
  
  #2 Download Tumor Infiltration Programs (TIP)
  R<-TIP_find_all()
  
  #3 Download mTIL program
  rslts<-readRDS(get.file("Results/HGSC_mTIL_Malignant2env_TNK.cell.rds"))
  
  #4 Download Perturb-seq meta-analyses results.
  rslts1<-readRDS(get.file("Results/HGSC_mTIL_K562_CRISPRa.rds"))
  rslts2<-readRDS(get.file("Results/HGSC_mTIL_K562_CRISPRi.rds"))
  rslts3<-readRDS(get.file("Results/HGSC_mTIL_RPE1_CRISPRi.rds"))
  
  #5 Download Perturb-Seq data and results.
  r.prt<-readRDS(get.file("Data/PerturbSeq_TYKnuNK.rds"))
  rslts4<-readRDS(get.file("Results/PerturbSeq_TYKnuNK_DEGs.rds"))
  mTIL.sig<-readRDS(get.file("Results/HGSC_mTIL.rds"))
  fitnessR<-readRDS(get.file("Results/HGSC_CRISPR.rds"))
  
  #6 Regenerate all Figures.
  # HGSC_Figure1_map()
  # HGSC_Figure2_driftDAF()
  HGSC_Figure3_TIP(rTNK.smi,r.xenium,rTNK.xenium,R)
  HGSC_Figure4_mTIL(r = r.smi,r1 = rmal.smi,rslts = rslts)
  HGSC_Figure5_CNAs()
  HGSC_Figure6_perturbMeta(rslts = rslts,rslts1 = rslts1,
                           rslts2 = rslts2,rslts3 = rslts3)
  HGSC_Figure7_perturbOC(r = r.prt,rslts = rslts4,sig = mTIL.sig,fitnessR = fitnessR)
}

get.file<-function(file1){
  dir1<-"/Volumes/ljerby/HGSC_Profiling/Manuscript/GitHub/"
  return(paste0(dir1,file1))
}
