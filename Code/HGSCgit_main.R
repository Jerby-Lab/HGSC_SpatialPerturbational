# The code below provides a step-by-step guide
# to regenerate the figures of the manuscript
# in its final version as of September 2023.

# Result Sections
#1 Data processing and annotations (Figure 1, Supplementary Tables 1-3)
#2 Effector T Cells Infiltrate the Tumor (Figure 2, Supplementary Tables 4-5)
#3 Malignant TIL (mTIL) program (Figure 3, Supplementary Table 6)
#4 CNAs mapping to mTIL and TIL levels (Figure 4, Supplementary Figure 5)
#5 Perturb-Seq meta-analyses (Figure 5, Supplementary Table 6)
#6 Perturb-Seq and genetic screen analyses (Figure 6, Supplementary Table 7)

HGSC_main<-function(){
  if(!file.exists(get.file("Figures/Fig3A.pdf"))){HGSC_source()}

  #1 Download Presets and Preprocessed Data for Spatiomolecular Profiling Analyses
  cell_2_rgb <- readRDS(get.file("Data/Cell_Colors.rds"))
  r.smi<-readRDS(get.file("Data/SMI_data.rds"))
  r.xenium<-readRDS(get.file("/Data/Xenium_data.rds"))
  r.merfish <- readRDS(get.file("/Data/MERFISH_data.rds"))

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

  #8 Regenerate main Figures and Tables.
  HGSC_Figure1_SpatiomolecularMapping(r = r.smi, q = r.xenium, s = r.merfish,
                                      cell_2_rgb=cell_2_rgb)
  # HGSC_Figure2_DriftDAF(r = r.smi,
  #                       sigs_cna = sigs_cna,sigs = sigs,
  #                       mal_umap = mal_umap, mal_drift = mal_drift,
  #                       mal_rf= mal_rf, lopo_rf = lopo_rf,
  #                       morph = morph, daf = daf,
  #                       cell_2_rgb = cell_2_rgb)
  HGSC_Figure3_TIP(rTNK.smi,
                   r.xenium,
                   rTNK.xenium,R)
  HGSC_Figure4_mTIL(r = r.smi,
                    r1 = rmal.smi,
                    rslts = rslts,
                    s = r.merfish)
  HGSC_Figure5_CNAs(r=r.smi)
  HGSC_Figure6_perturbMeta(rslts = rslts,rslts1 = rslts1,
                           rslts2 = rslts2,rslts3 = rslts3)
  HGSC_Figure7_perturbOC(r = r.prt,
                         rslts = rslts4,
                         sig = mTIL.sig,
                         fitnessR = fitnessR)
}

get.file<-function(file1){
  # dir1 <- "/path/to/your/local/clone/of/this/repo/"
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
