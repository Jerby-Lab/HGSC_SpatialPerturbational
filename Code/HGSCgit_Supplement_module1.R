# HGSC_Supplement_module1.R 
# author: cyyeh 
# -------------------------
# Code for Generating Supplement 
# r<-readRDS(get.file("Data/SMI_data.rds"))
# q<-readRDS(get.file("/Data/Xenium_data.rds"))
# s <- readRDS(get.file("/Data/MERFISH_data.rds"))
# t1 <- readRDS(get.file("Data/SMI6K_Data_wSubtypes.rds"))
# t2.1<- readRDS(get.file("Data/WT1_wSubtypes.rds"))
# t2.2 <- readRDS(get.file("Data/WT2_wSubtypes.rds"))
# t2.3 <- readRDS(get.file("Data/WT3_wSubtypes.rds"))
# t2.4 <- readRDS(get.file("Data/WT4_wSubtypes.rds"))

#2 Generate Table 1a, the specifications of ST and scRNA-seq platforms
print("Supplementary Table 1a")
HGSC_Table1a_write_specs(master)
#3 Generate Table 1b, genes of each dataset
print("Supplementary Table 1b")
HGSC_Table1b_write_genes(r=r,q=q,s=s)
#4 Generate Table 2a, key names of clinical data
print("Supplementary Table 2a")
HGSC_Table2a_write_colkeys()
#5 Generate Table 2b, colnames keys
print("Supplementary Table 2b")
HGSC_Table2b_write_fullmetadata(master)
#7 Generate Table 3a, cell type signatures from scRNA-seq
print("Supplementary Table 3a")
HGSC_Table3a_write_scRNAsigs()
#8 Generate Table 3b, cell type signatures from SMI
print("Supplementary Table 3b")

HGSC_Table1a_write_specs <- function(master){
  # per platform stats for ST
  per_platform <- master %>%
    group_by(platform) %>%
    summarize(n_patients = length(unique(patients)),
              n_sections = length(samples),
              n_cells = sum(n_cells)) %>%
    mutate(
      target_molecule = rep("RNA",3),
      n_genes = c(140, 960, 280),
      tissue_type = c("Whole Tissue", "FOV (0.6mm x 0.9um)",
                      "Tissue Cores (1.5mm x 1.5mm)"),
      resolution = rep("Subcellular",3),
      spatial = rep("Yes", 3))
  
  # aggregate stats for scRNA-seq
  files <- Sys.glob(get.file("Data/scRNA_*"))
  df <- do.call("rbind", mclapply(files, mc.cores = 10, function(x){
    r1 <- readRDS(x)
    out <- c(n_patients = length(unique(r1$patients)),
             n_sections = length(unique(r1$samples)),
             n_cells = length(r1$cells)
    )
    return(out)
  }))
  out <- apply(df, 2, sum)
  
  # combine
  final <- rbind(per_platform, c(platform = "scRNA",
                                 out,
                                 target_molecule = "RNA",
                                 n_genes = "Whole Transcriptome",
                                 tissue_type = "Dissociated Tissue",
                                 resolution = "Single Cell",
                                 spatial = "No"))
  row.names = final$platform
  final <- select(final, spatial, resolution, tissue_type, target_molecule, n_genes,
                  n_patients, n_sections, n_cells)
  final <- data.frame(t(final))
  colnames(final) <- row.names
  final <- cbind(field = row.names(final), final)
  
  # write to disk
  write.xlsx(final, file = get.file("Tables/Table1a.xlsx"),
             asTable = T)
}

HGSC_Table1b_write_genes <- function(r, q, s){
  out <- as.data.frame(list.2.mat(list(SMI = r$genes,
                                       ISS = q$genes,
                                       MERFISH = s$genes)))
  
  # write to disk
  write.xlsx(out, file = get.file("Tables/Table1b.xlsx"),
             asTable = T)
}

HGSC_Table2a_write_colkeys <- function(){
  colkeys <- readRDS(get.file("Data/ClinicalDataColumnKey.rds"))
  write.xlsx(colkeys, file = get.file("Tables/Table2a.xlsx"),
             asTable = T)
}

HGSC_Table2b_write_fullmetadata <- function(master){
  write.xlsx(master, file = get.file("Tables/Table2b.xlsx"),
             asTable = T)
}

HGSC_Table2b_mastertable <- function(r, q, s, t1, t2.1, t2.2, t2.3, t2.4){
  clin <- readRDS("Data/ClinicalUpdate.rds")
  muts <- readRDS(get.file("Data/BRCA_TP53_muts.rds"))
  
  # discovery dataset (SMI)
  smi <- data.frame(r[c("patients",
                        "samples",
                        "treatment",
                        "sites_binary")]) %>%
    group_by(samples, patients, sites_binary, treatment) %>%
    summarize(n_cells = length(treatment)) %>%
    mutate(platform = "SMI", 
           dataset = "Discovery")
  
  # validation dataset 2 (MERFISH)
  mer <- data.frame(s[c("patients", "samples", "sites")]) %>%
    group_by(samples, patients, sites) %>%
    mutate(sites = "Adnexa") %>%
    rename(sites_binary = sites) %>%
    summarize(n_cells = length(samples)) %>%
    mutate(platform = "MERFISH") %>% 
    mutate(dataset = "Validation 1")
  map <- data.frame(filter(clin, patients %in% mer$patients) %>% select(patients, specimen) %>%
                      mutate(treatment = c("Treated", "Untreated")[grepl("cyto", specimen) + 1]) %>%
                      select(patients, treatment))
  row.names(map) <- map$patients
  mer$treatment <- map[mer$patients,]$treatment
  mer <- mer %>% select(samples, patients, sites_binary, treatment, n_cells, platform, dataset)
  
  # validation dataset 1 (ISS)
  xen <- data.frame(q[c("patients", "samples", "treatment", "sites_binary")]) %>%
    group_by(samples, patients, sites_binary, treatment) %>%
    summarize(n_cells = length(treatment)) %>%
    mutate(platform = "ISS") %>% 
    mutate(dataset = "Validation 2")
  
  # test dataset 1 (SMI 6K)
  test1 <- data.frame(t1[c("patients", "samples", "treatment", "sites_binary")]) %>%
    group_by(samples, patients, sites_binary, treatment) %>%
    summarize(n_cells = length(sites_binary)) %>%
    mutate(platform = "SMI") %>% 
    mutate(dataset = "Test 1")
  
  # test datatset 2 (SMI WT)
  t2 <- bind_rows(lapply(list(t2.1, t2.2, t2.3, t2.4), function(x){
    df <- data.frame(x[c("patients", "treatment", "sites_binary")]) %>%
      mutate(samples = paste0("SMIWT_", gsub("HGSC", "", patients)), 
             treatment = c("Untreated", "Treated")[treatment]) %>%
      group_by(samples, patients, sites_binary, treatment) %>%
      summarize(n_cells = length(treatment)) %>%
      mutate(platform = "SMI") %>% 
      mutate(dataset = "Test 2") %>%
      select(samples, patients, sites_binary, treatment, n_cells, platform, dataset)
    
    return(df)
  }))
  
  # merge ST with clinical
  master <- bind_rows(list(smi, xen, mer, test1, t2))
  master <- merge(master, clin[!duplicated(clin$patients),], by.x = "patients", by.y = "patients", all.x = T)
  
  # merge dna with ST and clinical
  master <- merge(master, muts, by = "patients", all.x = T)
  
  # get TMB 
  molmaster <- read.csv("~/Projects/SpatialMapping_HGSC/Data/raw/xT/molecular_master_file.csv")
  tmb <- filter(molmaster, variant_type_code == "TMB",
                partner_patient_id %in% master$patients | partner_patient_id %in% master$patients) %>% 
    select(partner_patient_id, result) %>% 
    rename(patients = partner_patient_id, 
           tmb = result)
  master <- merge(master, tmb, by = "patients", all.x = T)
  
  # add transcripts quantification
  smi <- data.frame(r[c("comp.reads", "samples")]) %>%
    group_by(samples) %>%
    summarize(median_tpc = median(comp.reads, na.rm = T),
              mean_tpc = mean(comp.reads, na.rm = T))
  iss <- data.frame(q[c("comp.reads", "samples")]) %>%
    group_by(samples) %>%
    summarize(median_tpc = median(comp.reads, na.rm = T),
              mean_tpc = mean(comp.reads, na.rm = T))
  merfish <- data.frame(s[c("comp.reads", "samples")]) %>%
    group_by(samples) %>%
    summarize(median_tpc = median(comp.reads, na.rm = T),
              mean_tpc = mean(comp.reads, na.rm = T))
  test1 <- data.frame(t1[c("comp.reads", "samples")]) %>%
    group_by(samples) %>%
    summarize(median_tpc = median(comp.reads, na.rm = T),
              mean_tpc = mean(comp.reads, na.rm = T))
  t2 <- bind_rows(lapply(list(t2.1, t2.2, t2.3, t2.4), function(x){
    df <- data.frame(x[c("comp.reads", "patients")]) %>%
      mutate(samples = paste0("SMIWT_", gsub("HGSC", "", patients))) %>%
      group_by(samples) %>%
      summarize(median_tpc = median(comp.reads, na.rm = T),
                mean_tpc = mean(comp.reads, na.rm = T))
    return(df)
  }))
  
  transcripts <- bind_rows(list(smi, iss, merfish, test1, t2))
  master <- merge(master, transcripts, by = "samples", all.x = T)
  
  master$patients[is.na(master$patients)] <- "HGSC6K"
  master$sites_binary[is.na(master$sites_binary)] <- "Adnexa"
  
  master_out <- select(master, samples, dataset, platform, n_cells, median_tpc, mean_tpc, TMA,
                       patients, sites_binary, age, stage, 
                       treatment, "1L", parpi_recurrence, immunotherapy,
                       BRCA1_Somatic,
                       BRCA2_Somatic,
                       TP53_Somatic,
                       BRCA1_Germline,
                       BRCA2_Germline,
                       TP53_Germline, 
                       tmb, 
                       fu_time1, fu_time2, outcome, pfs4)
  
  master_out$TMA <- unlist(lapply(master_out$TMA, function(x) {
    split = strsplit(x, split = "\\/")[[1]]
    split[length(split)]
  }))
  
  master_out <- master_out %>% 
    mutate(dataset = factor(dataset, levels = c("Discovery", "Validation 1", 
                                                "Validation 2", 
                                                "Test 1", "Test 2"))) %>% 
    arrange(dataset, samples) %>% 
    rename(profile = samples) %>% 
    rename(parpi = parpi_recurrence) 
  #clean surv
  surv <- master_out
  surv$fu_time1 <- as.numeric(surv$fu_time1)
  surv$stage <- gsub("[^IV]", "", surv$stage)
  surv$stage <- gsub("IIII", "III", surv$stage)
  
  
  write.csv(surv, get.file("Tables/SupplementaryTable2B.csv"), quote = F, row.names = F)
  saveRDS(surv, get.file("Data/MasterProfileData.rds"))
}

HGSC_Table3a_write_scRNAsigs <- function(){
  cell.sig.scRNA <- readRDS(get.file("Results/cell.sigs.scRNA.combined.rds"))
  cell.sig.CellTypist <- readRDS(get.file("Results/cell.sigs.subtype.CellTypist.rds"))
  one <- cell.sig.scRNA[c("B.cell", "Endothelial", "Fibroblast",
                          "Myeloid", "Malignant", "TNK.cell", "Mast.cell")]
  names(one)[4] <- "Monocyte"
  two <- cell.sig.CellTypist[c("NK.cells", "Regulatory.T.cells",
                               "CD8.T.cell", "CD4.T.cell")]
  out <- c(one, two)
  
  write.xlsx(data.frame(list.2.mat(out)), file = get.file("Tables/Table3a.xlsx"),
             asTable = T)
}

HGSC_Table3b_write_SMIsigs <- function(r, so_smi, cell.sig.SMI.subtypes){
  so_smi <- readRDS(get.file("Results/HGSC_SMI_wUMAP.rds"))
  cell.sig.SMI.subtypes <- readRDS(get.file("Results/cell.sigs.SMI.subtypes.rds"))
  q <- subset_list(r, colnames(so_smi))
  q$cell.types <- sub("_LC", "", q$cell.types)
  rslts <- scRNA_denovo.cell.type.markers(q)
  out <- rslts$sig[-8]
  write.xlsx(data.frame(list.2.mat(c(out, cell.sig.SMI.subtypes))),
             file = get.file("Tables/Table3b.xlsx"),
             asTable = T)
}