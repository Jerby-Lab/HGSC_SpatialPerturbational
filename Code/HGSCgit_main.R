# The code below provides a step-by-step guide
# to regenerate the figures of the manuscript
# in its final version as of September 2023.

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
  if(!file.exists(get.file("Figures/Fig10A.pdf"))){HGSC_source()}

  #1 Download Presets 
  cell_2_rgb <- readRDS(get.file("Data/Cell_Colors.rds"))

  #2 Download T/NK and malignant spatial data.
  rTNK.smi<-HGSC_SMI.process.CD8()
  rmal.smi<-HGSC_SMI.process.mal()
  rTNK.xenium<-HGSC_Xenium.process.CD8.NK(r = r.xenium)

  #3 Download Tumor Infiltration Programs (TIP).
  R<-TIP_find_all()

  #4 Download Desmoplastic Fibroblast results
  morph <- readRDS(get.file("Data/StromalMorphology.rds"))
  daf <- readRDS(get.file("Results/HGSC_DAF.rds"))

  #5 Download mTIL program.
  rslts<-readRDS(get.file("Results/HGSC_mTIL_Malignant2env_TNK.cell.rds"))

  #6 Download Perturb-seq meta-analyses results.
  rslts1<-readRDS(get.file("Results/HGSC_mTIL_K562_CRISPRa.rds"))
  rslts2<-readRDS(get.file("Results/HGSC_mTIL_K562_CRISPRi.rds"))
  rslts3<-readRDS(get.file("Results/HGSC_mTIL_RPE1_CRISPRi.rds"))

  #7 Download Perturb-Seq data and results.
  r.prt<-readRDS(get.file("Data/PerturbSeq_TYKnuNK.rds"))
  rslts4<-readRDS(get.file("Results/PerturbSeq_TYKnuNK_DEGs.rds"))
  mTIL.sig<-readRDS(get.file("Results/HGSC_mTIL.rds"))
  fitnessR<-readRDS(get.file("Results/HGSC_CRISPR.rds"))

  #8 Regenerate Main Figures
  # HGSC_Figure1_SpatiomolecularMapping(cell_2_rgb=cell_2_rgb)
  r.smi<-readRDS(get.file("Data/SMI_data.rds"))
  r.xenium<-readRDS(get.file("/Data/Xenium_data.rds"))
  HGSC_Figure2_TIP_DF(rTNK.smi,
                      r.xenium,
                      rTNK.xenium,R,
                      r = r.smi,
                      morph=morph, daf=daf,
                      cell_2_rgb=cell_2_rgb)
  HGSC_Figure3_mTIL(r = r.smi,
                    r1 = rmal.smi,
                    rslts = rslts)
  HGSC_Figure4_ICB()
  HGSC_Figure5_CNAs(r=r.smi)
  HGSC_Figure6_perturbMeta(rslts = rslts,rslts1 = rslts1,
                           rslts2 = rslts2,rslts3 = rslts3)
  HGSC_Figure7_perturbOC(r = r.prt,
                         rslts = rslts4,
                         sig = mTIL.sig,
                         fitnessR = fitnessR)
  HGSC_Figure8_KOValidation()
}

get.file<-function(file1){
  dir1 <- "~/Projects/HGSC_SpatialPerturbational/"
  return(paste0(dir1,file1))
}

HGSC_source<-function(){
  print("Sourcing HGSC_SpatialPerturbational GitHub repository.")
  files <- list.files(get.file("Code"),
                      include.dirs = F,
                      pattern = ".R",
                      full.names = T)
  lapply(files, source)
  return()
}
