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

HGSC_main<-function(){
  #1 Download Presets and Preprocessed Data for Spatiomolecular Profiling Analyses
  cell_2_rgb <- readRDS(get.file("Data/Cell_Colors.rds"))
  r.smi<-readRDS(get.file("Data/SMI_data.rds"))
  r.xenium<-readRDS(get.file("/Data/Xenium_data.rds"))
  r.merfish <- readRDS(get.file("/Data/MERFISH_data.rds"))

  #2 Download Drift and DAF results
  sigs <- readRDS(get.file("Results/HGSC_NACTSites_pergene.rds"))
  sigs_cna <- readRDS(get.file("Results/HGSC_CNA_pergene.rds"))
  mal_umap <- readRDS(get.file("Results/HGSC_MalUMAP_30PCs.rds"))
  mal_rf <- readRDS(get.file("Results/HGSC_PatientsRF.rds"))
  lopo_rf <- readRDS(get.file("Results/HGSC_SitesRF.rds"))
  mal_drift <- readRDS(get.file("Results/HGSC_MalDrift.rds"))
  morph <- readRDS(get.file("Data/StromalMorphology.rds"))
  daf <- readRDS(get.file("Results/HGSC_DAF.rds"))

  #3 Download Tumor Infiltration Programs (TIP).
  rTNK.smi<-HGSC_SMI.process.CD8()
  rmal.smi<-HGSC_SMI.process.mal()
  rTNK.xenium<-HGSC_Xenium.process.CD8.NK(r = r.xenium)

  #4 Download Tumor Infiltration Programs (TIP).
  R<-TIP_find_all()

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

  #8 Regenerate all Figures and Tables.
  HGSC_Figure1_SpatiomolecularMapping(r = r.smi, q = r.xenium, s = r.merfish,
                                      cell_2_rgb=cell_2_rgb)
  HGSC_Figure2_DriftDAF(r = r.smi,
                        sigs_cna = sigs_cna,sigs = sigs,
                        mal_umap = mal_umap, mal_drift = mal_drift,
                        mal_rf= mal_rf, lopo_rf = lopo_rf,
                        morph = morph, daf = daf,
                        cell_2_rgb = cell_2_rgb)
  HGSC_Figure3_TIP(rTNK.smi,
                   r.xenium,
                   rTNK.xenium,R)
  HGSC_Figure4_mTIL(r = r.smi,
                    r1 = rmal.smi,
                    rslts = rslts,
                    s = r.xenium)
  HGSC_Figure5_CNAs()
  HGSC_Figure6_perturbMeta(rslts = rslts,rslts1 = rslts1,
                           rslts2 = rslts2,rslts3 = rslts3)
  HGSC_Figure7_perturbOC(r = r.prt,
                         rslts = rslts4,
                         sig = mTIL.sig,
                         fitnessR = fitnessR)
}

get.file<-function(file1){
  dir1 <- "~/Projects/HGSC_SpatialPerturbational/"
  return(paste0(dir1,file1))
}
