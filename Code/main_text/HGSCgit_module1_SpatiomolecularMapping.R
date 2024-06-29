#### Results Section 1 ###
# Figure 1. Single cell spatial transcriptomics (ST) mapping of HGSC.
#
# Figure 1a. Spatial Cohort Design (generated in Figma)
# Figure 1b. CoMut Plot (output files in Results/ used for python comut code)
# Figure 1c. Cell Type UMAPs
# Figure 1d. Cell Types in situ 
# Figure 1e. Coembedding scross ST and scRNA-seq datasets
# Figure 1f. Cell Type Compositions
# Figure 1g. Hypergeometric Spatial Mapping
# Figure 1h. Cell Co-localization Quotient Analysis

#' Figure 1 Wrapper Function
#'
#' This function calls code to reproduce main text Figures 1b-h. 
#'
#' @param cell_2_rgb a named list of colors assigned to each of the cell types
#' @return this function returns nothing, but writes figures in .pdf format 
#' in the Figures/ folder. 
HGSC_Figure1_SpatiomolecularMapping <- function(cell_2_rgb){
  #1 Regenerate input files for CoMut plot in python
  master <- HGSC_Fig1b_make_comut()
  #2 Regenerate Figure 1c: Cell Type UMAPs
  HGSC_Fig1c_celltype_umaps(cell_2_rgb = cell_2_rgb)
  #3 Regenerate Figure 1d: Cell Types in situ
  HGSC_Fig1d_celltypes_insitu(cell_2_rgb=cell_2_rgb)
  #4 Regenerate Figure 1e: Coembedding scross ST and scRNA-seq datasets
  HGSC_Fig1e_coembedding(cell_2_rgb)
  #5 Regenerate Figure 1f: Cell Type Compositions
  HGSC_Fig1f_celltype_composition(cell_2_rgb = cell_2_rgb)
  #6 Regenerate Figure 1g
  HGSC_Fig1g_hg(cell_2_rgb = cell_2_rgb)
  #7 Regenerate Figure 1h: Cell Co-localization Quotient Analysis
  HGSC_Fig1h_clq(cell_2_rgb = cell_2_rgb)
  return()
}

#' Figure 1b. CoMut Files
#'
#' This function generates the .tsv files that are inputed into python "CoMut" 
#' code for generating a summary figure of multi-modal data for this HGSC cohort
#'
#' @return a data frame of simplified clinical meta data. 
HGSC_Fig1b_make_comut <- function(){
  master <- readRDS(get.file("Data/MasterProfileData.rds"))
  
  # fix immunotherapy text
  master$immunotherapy[grepl("embro", master$immunotherapy) | 
                         grepl( "varlilumab" , master$immunotherapy)] <- "Yes"
  master$immunotherapy[master$immunotherapy == "NA"] <- NA
  master$immunotherapy[grepl("placebo", master$immunotherapy)] <- NA
  master$immunotherapy[master$immunotherapy != "Yes"] <- "No"
  master$immunotherapy <- factor(master$immunotherapy, levels = c("No", "Yes"))
  
  # recode parpi inhibitor usage
  master$parpi <- master$parpi
  master$parpi[master$parpi == 3] <- NA
  master$parpi[master$parpi == 2 | master$parpi == 4] <- 1
  master$parpi <- factor(master$parpi, levels = c(0, 1))
  
  # fix text for whether patient receivved bevacizumab 
  master$`1L`[master$`1L` == 5] <- NA
  master$`1L`[master$`1L` == 2] <- NA
  master$`1L` <- c("No", "", "Yes")[master$`1L`]
  master$`1L` <- factor(master$`1L`, levels = c("No", "Yes"))
  
  # write clinical data to disk
  outdir <- paste0(get.file("Results/CoMut/"))
  if(!dir.exists(outdir)){dir.create(outdir, recursive = T)}
  df <- select(master,
               patients,
               age,
               stage,
               "1L", 
               parpi,
               treatment,
               immunotherapy,
               outcome,
               pfs4)
  
  # fix outcome 
  df$outcome <- unlist(lapply(df$outcome, function(x){
    if (is.na(x)) {
      return(NA)
    } else if (grepl("dead", x, ignore.case = T)) {
      return("Dead")
    } else if (x == "Alive") {
      return("Alive")
    } else {
      return(NA)
    }
  }))
  
  # format factors for data
  df <- df %>%
    mutate(patients = as.character(patients),
           age = c("<65", ">65")[(as.numeric(age) >= 65) + 1],
           stage = factor(stage, c("I", "II", "III", "IV")),
           treatment = c("Yes", "No")[(treatment == "Untreated") +1],
           beva = as.character(`1L`), 
           parpi = as.character(parpi),
           immunotherapy = as.character(immunotherapy),
           outcome = as.character(outcome),
           pfs4= c("<6m", ">6m")[(as.numeric(pfs4) >= 180) + 1]) %>%
    unique()
  
  # remove unwanted column. 
  df <- df[-4]
  colnames(df) <- c("Patient", "Age", "Stage", "PARPi", "NACT",
                    "Immunotherapy", "Outcome", "PFS","Bevacizumab")
  
  # write treatment history out 
  for (i in 2:9) {
    selected <- df[,c(1,i)]
    field <- colnames(df)[i]
    selected$categoory <- field
    out <- selected[,c(1, 3, 2)]
    colnames(out) <- c("sample", "category", "value")
    write.table(out, paste0(outdir, field,
                            "_CoMut_clin.tsv"), sep = "\t",
                row.names = F, quote = F)
  }
  
  # manage mutations (TP53 and BRCA1/2 mutations)
  ## tp53
  tp53<- master %>%
    select(patients, TP53_Somatic, TP53_Germline) %>% unique()
  tp53$somatic <- unlist(lapply(tp53$TP53_Somatic, function(x){
    if (is.na(x)) {
      return(NA)
    } else if (x == "P" | x == "LP") {
      return("Pathogenic Mutation (Somatic)")
    } else {
      return("No Pathogenic Mutation")
    }
  }))
  tp53$germline <- unlist(lapply(tp53$TP53_Germline, function(x){
    if (is.na(x)) {
      return(NA)
    } else if (x == "P" | x == "LP") {
      return("Pathogenic Mutation (Germline)")
    } else {
      return("No Pathogenic Mutation")
    }
  }))
  tp53 <- select(tp53,patients, somatic) %>%
    mutate(gene = "TP53") %>%
    select(patients, gene, somatic) %>%
    rename(mutation = somatic)
  ## brca1/2
  brca <- master %>%
    select(patients, BRCA1_Germline, BRCA2_Germline,
           BRCA1_Somatic, BRCA2_Somatic) %>%
    unique()
  
  brca$germline <- unlist(apply(brca, 1, function(row){
    if (is.na(row[2])) {
      return(NA)
    } else if (row[2] == "P" | row[3] == "P" ) {
      return("Pathogenic Mutation (Germline)")
    } else {
      return ("No Pathogenic Mutation")
    }
  }))
  brca$somatic <- unlist(apply(brca, 1, function(row){
    if (is.na(row[4])) {
      return(NA)
    } else if (row[4] == "P" | row[5] == "P" ) {
      return("Pathogenic Mutation (Somatic)")
    } else {
      return ("No Pathogenic Mutation")
    }
  }))
  blah <- filter(brca, grepl("Germline", germline),
                 grepl("Somatic", germline))
  brca <- rbind(filter(brca, grepl("Germline", germline)) %>%
                  select(patients, germline) %>%
                  rename(mutation = germline) %>%
                  mutate(gene = "BRCA1/2"),
                filter(brca, grepl("Somatic", somatic)) %>%
                  select(patients, somatic) %>%
                  rename(mutation = somatic) %>%
                  mutate(gene = "BRCA1/2"),
                filter(brca, grepl("No", germline), grepl("No", somatic)) %>%
                  select(patients, somatic) %>%
                  rename(mutation = somatic) %>%
                  mutate(gene = "BRCA1/2")) %>%
    select(patients, gene, mutation)
  
  ## combine mutations
  df_mut <- rbind(brca, tp53 %>% filter(!is.na(mutation)))
  write.table(df_mut, paste0(outdir, "CoMut_Muts.tsv"),
              sep = "\t", row.names = F, quote = F)
  
  # load tempus data to show which pts have tumor mutational burden data
  df_tmp <- data.frame(patients = unique(master$patients),
                       field = rep("Tempus xT",
                                   length(unique(master$patients))),
                       status = c("No", "Yes")[
                         (
                           unique(
                             master$patients) %in%
                             unique(filter(
                               read.csv(get.file("Data/molecular_master_file.csv")),
                               sample_type == "dna")$partner_patient_id
                             )
                         ) + 1]
  )
  write.table(df_tmp, paste0(outdir,"CoMut_Tempus.tsv"),
              sep = "\t", row.names = F, quote = F)
  
  # spatial transcriptomics data summarization
  df_st <- master %>%
    select(patients, dataset, sites_binary) %>%
    arrange(dataset) %>%
    unique()
  
  write.table(df_st, paste0(outdir,"CoMut_Spatial.tsv"),
              sep = "\t", row.names = F, quote = F)
  
  samples <- (df %>% arrange(Stage, Age))$Patient
  write.table(data.frame(samples = samples),
              file = paste0(outdir,"CoMut_patorder.tsv"),
              sep = "\t", row.names = F, quote = F)
  
  return(master)
}

#' Figure 1c. Cell Type UMAPs 
#'
#' This function loads Seurat Objects that store the UMAP embeddings for high-
#' confidence cells for each of the Discovery, Validation 1, and Validation 2
#' datasets. The function loads list objects with the UMAP embedding info for 
#' the Test 1 and Test 2 datasets. Finally this function plots the UMAP
#' embeddings for each dataset where the cells are colored by their assigned 
#' cell types. 
#'
#' @return this function returns nothing, but writes figures in .pdf format 
#' in the Figures/ folder. 
HGSC_Fig1c_celltype_umaps <- function(cell_2_rgb){
  # load data and results 
  r <- readRDS(get.file("Data/SMI_data.rds"))
  so_smi <- readRDS(get.file("Results/HGSC_SMI_wUMAP.rds"))
  so_iss <- readRDS(get.file("Results/HGSC_ISS_wUMAP.rds"))
  so_mer <- readRDS(get.file("Results/HGSC_MERFISH_wUMAP.rds"))
  so_t1 <- readRDS(get.file("Results/HGSC_SMI6K_UMAP_coord.rds"))
  so_t2 <- readRDS(get.file("Results/HGSC_SMIWT_wUMAP.rds"))
  
  # function for plotting cell type umaps. 
  plot_ct_umap <- function(so, 
                           cell_2_rgb, title, field = "cell.types.nonmal"){
    col_field <- names(
      cell_2_rgb
    )[names(cell_2_rgb) %in% unique(r$cell.types)]
    p <- DimPlot(so, group.by = field, pt.size = 0.1) +
      scale_color_manual(values = unlist(
        lapply(
          cell_2_rgb[col_field], rgb2hex
        )
      ), name = "Cell Types") +
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.title = element_text(size =12),
            legend.position = "none") +
      ggtitle(title)
    return(p)
  }
  
  # plot umap for Discovery, Validation 1, Validation 2
  a = plot_ct_umap(so_smi, cell_2_rgb, "Discovery Dataset")
  b = plot_ct_umap(so_iss, cell_2_rgb, "Validation Dataset  1")
  c = plot_ct_umap(so_mer, cell_2_rgb, "Validation Dataset 2")
  
  # plot umap for Test 1 and Test 2
  d = ggplot(so_t1, aes(x = umap_1, y = umap_2, col = ct)) + 
    geom_point(size = 0.1) + 
    theme_classic() +
    scale_color_manual(values = unlist(
      lapply(
        cell_2_rgb, rgb2hex
      )
    ), name = "Cell Types") +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.title = element_text(size =12),
          legend.position = "none") +
    coord_fixed()
  e = plot_ct_umap(so_t2, cell_2_rgb, "Test 2", field = "cell.types.conf")
  
  # format legend
  col_field <- names(
    cell_2_rgb
  )[names(cell_2_rgb) %in% unique(r$cell.types)]
  p <- DimPlot(so_smi, group.by = "cell.types.nonmal", pt.size = 0.1) +
    scale_color_manual(values = unlist(
      lapply(
        cell_2_rgb[col_field], rgb2hex
      )
    ), name = "Cell Types") +
    theme(legend.title = element_text(size =12))
  f = as_ggplot(get_legend(p))
  
  # write to disk 
  pdf(get.file("Figures/Fig1c.pdf"), width = 15, height = 10)
  p <- ggarrange(a,b,c,d,e,f)
  print(p)
  dev.off()
}

#' Figure 1d. Cell Types plotted in situ. 
#'
#' For smaller FOVs of tissue, we use segmentation maps to visualize the cell 
#' type calls in situ. For the Whole Tissue samples, we use 
#' 
#' @param cell_2_rgb a named list of colors assigned to each of the cell types
#' @return this function returns nothing, but writes figures in .pdf or .png 
#' format in the Figures/ folder. 
HGSC_Fig1d_celltypes_insitu<- function(cell_2_rgb){
  # plot 4 samples from the discovery dataset
  r<-readRDS(get.file("Data/SMI_data.rds"))
  samples <- c("SMI_T13_F001", "SMI_T10_F001", "SMI_T12_F016", "SMI_T10_F021")
  lapply(samples, function(x){
    segpath <- Sys.glob(paste0(get.file("Results/Segmentation/"),
                               x,
                               "*"))
    q <- subset_list(r, r$cells[r$samples == x])
    celltypes <- q$cell.types
    celltypes <- gsub("_LC", "", q$cell.types)
    spatial_sample_visualization(segpath,
                                 celltypes = celltypes,
                                 cell2rgb = cell_2_rgb,
                                 samplename = x,
                                 outfile = get.file(paste0("Figures/Fig1d_",
                                                           x,
                                                           ".png")))
  })
  remove(r) # for memory efficiency
  
  # plot the select sample from the validation 1 dataset (ISS)
  q<-readRDS(get.file("/Data/Xenium_data.rds"))
  x <- "XEN_T10_A2"
  spot = strsplit(x, split = "_")[[1]][3]
  segpath <- Sys.glob(paste0(get.file("Results/Segmentation/"),
                             spot,
                             "*"))
  q1 <- subset_list(q, q$cells[q$samples == x])
  celltypes <- q1$cell.types
  celltypes <- gsub("_LC", "", q1$cell.types)
  celltypes
  spatial_sample_visualization(segpath,
                               celltypes = celltypes,
                               cell2rgb = cell_2_rgb,
                               samplename = x,
                               outfile = get.file(paste0("Figures/Fig1d_",
                                                         "XEN_T10_", spot,
                                                         ".png")))
  remove(q) # for memory efficiency 
  
  # plot select sample from valifation 2 dataset (MERFISH)
  s <- readRDS(get.file("/Data/MERFISH_data.rds"))
  plot <- filter(data.frame(s[c("cell.types.nonmal", "samples")], s$coor),
                 samples == "MER_TB21361")
  p <- ggplot(plot, aes(x = x, y = y, col = cell.types.nonmal)) +
    geom_point(size = 0.1) +
    coord_fixed() +
    theme_classic() +
    scale_color_manual(values = lapply(cell_2_rgb, rgb2hex)) +
    theme(panel.background = element_rect(fill = "black"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "none")
  png(get.file("Figures/Fig1d_MER_TB21361.png"),
      height = 16, width = 16, units = "in", res = 500)
  print(p)
  dev.off()
  remove(s)
  
  # plot selected smaples from test 2
  t1 <- readRDS(get.file("Data/SMI6K_Data_wSubtypes.rds"))
  plot <- filter(data.frame(t1[c("cell.types.conf", 
                                 "samples", "patients", 
                                 "sites_binary")], t1$coor))
  samples = c("SMI6K_2_F00021", "SMI6K_2_F00027", 
              "SMI6K_2_F00019", "SMI6K_2_F00037")
  lapply(samples, function(x){
    print(x)
    fov = strsplit(x, split = "_")[[1]][3]
    segpath <- Sys.glob(paste0(get.file("Results/Segmentation/*"),
                               fov,
                               "*.csv"))
    q <- subset_list(t1, t1$cells[t1$samples == x])
    celltypes <- q$cell.types.conf
    celltypes <- gsub("_LC", "", q$cell.types.conf)
    cell_2_rgb <- list()
    cell_2_rgb[["B.cell"]] = c(255, 0, 0)
    cell_2_rgb[["Fibroblast"]] = c(34, 91, 224)
    cell_2_rgb[["Endothelial"]] = c(166, 128, 71)
    cell_2_rgb[["Monocyte"]] = c(255, 153, 0)
    cell_2_rgb[["Malignant"]] = c(7, 224, 0)
    cell_2_rgb[["TNK.cell"]] = c(255, 247, 0)
    spatial_sample_visualization(segpath,
                                 celltypes = celltypes,
                                 cell2rgb = cell_2_rgb,
                                 samplename = x,
                                 outfile = get.file(
                                   paste0("Figures/Fig1d_",
                                          x, "_", 
                                          unique(q$patients), 
                                          "_", 
                                          unique(q$sites_binary), 
                                          ".png")))
  })
  remove(t1) # for memory efficiency 
  
  # plot selected sample from test 2
  t2.3 <- readRDS(get.file("Data/WT3_wSubtypes.rds"))
  plot <- filter(data.frame(t2.3[c("cell.types.conf", "labels")], 
                            t2.3$coor.global))
  p <- ggplot(plot, aes(x = x_global_px, y = y_global_px, 
                        col = cell.types.conf)) +
    geom_point(size = 0.1) +
    coord_fixed() +
    theme_classic() +
    scale_color_manual(values = lapply(cell_2_rgb, rgb2hex)) +
    theme(panel.background = element_rect(fill = "black"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "none")
  png(get.file("Figures/Fig1d_WT.png"),
      height = 16, width = 16, units = "in", res = 500)
  print(p) 
  dev.off()
}

#' Figure 1e. Coembedding Validation with scRNA-seq
#'
#' For smaller FOVs of tissue, we use segmentation maps to visualize the cell 
#' type calls in situ. For the Whole Tissue samples, we use 
#' 
#' @param cell_2_rgb a named list of colors assigned to each of the cell types
#' @return this function returns nothing, but writes figures in .pdf or .png 
#' format in the Figures/ folder. 
HGSC_Fig1e_coembedding <- function(cell_2_rgb){
  c1<- readRDS(get.file("Results/HGSC_coembedding.rds"))
  c2 <- readRDS(get.file("Results/HGSC_coembeddingTest.rds"))
  coembedding <- rbind(c1, c2)
  
  # cell types figure
  a <- ggplot(coembedding[!is.na(coembedding$cell.types),],
              aes(x = UMAP_1, y = UMAP_2, col = cell.types) ) +
    geom_point(size = 0.05) +
    scale_color_manual(values = lapply(cell_2_rgb, opipes::rgb2hex)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),axis.ticks = element_blank(),
          legend.text = element_text(size = 18)) +
    guides(colour = guide_legend("" , override.aes = list(fill="white",
                                                          size = 2)))
  # dataset figure
  coembedding$dataset <- factor(coembedding$dataset, 
                                levels = c("SMI", "ISS", 
                                           "MERFISH",
                                           "Test 1", 
                                           "Test 2", 
                                           "Vazquez-García", "Geistlinger",
                                           "Olalekan", "Regner", "Qian",
                                           "Shih", "hornburg"))
  coembedding <- coembedding[sample(1:(dim(coembedding)[1])),]
  b <- ggplot(coembedding[!is.na(coembedding$cell.types),], 
              aes(x = UMAP_1, y = UMAP_2,
                  col = dataset) ) +
    geom_point(size = 0.05) +
    scale_color_manual(values = c(RColorBrewer::brewer.pal(9, "Set3")[1:3], 
                                  "#DDE78A", 
                                  "#C47AB3", 
                                  RColorBrewer::brewer.pal(9, "Set3")[4:9])) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),axis.ticks = element_blank(),
          legend.text = element_text(size = 18)) +
    guides(colour = guide_legend("" , override.aes = list(fill="white",
                                                          size = 2)))
  
  # print to disk
  pdf(get.file("Figures/Fig1e.pdf"),
      width = 6*2, height = 4.5*2)
  print(a);print(b)
  dev.off()
}

#' Figure 1f. Coembedding Validation with scRNA-seq
#'
#' For smaller FOVs of tissue, we use segmentation maps to visualize the cell 
#' type calls in situ. For the Whole Tissue samples, we use 
#' 
#' @param cell_2_rgb a named list of colors assigned to each of the cell types
#' @return this function returns nothing, but writes figures in .pdf or .png 
#' format in the Figures/ folder. 
HGSC_Fig1f_celltype_composition <- function(cell_2_rgb,
                                            ct = c("Malignant", "Fibroblast",
                                                   "Endothelial", "Monocyte",
                                                   "Mast.cell", "B.cell",
                                                   "TNK.cell", "Other")){
  # general composition function
  compute_composition <- function(r, out, ct) {
    out <- do.call("rbind", lapply(unique(r$samples), function(s){
      meta <- (filter(out, samples == s) %>%
                 select(-cells, -cell.types))[1:(length(ct)),] %>%
        cbind.data.frame(cell.types = ct) %>%
        cbind.data.frame(count =  unlist(
          unname(
            as.list(table(filter(out, samples == s)$cell.types))))) %>%
        mutate(proportion = count/(dim(filter(out, samples == s))[1]) )
      row.names(meta) <- NULL
      return(meta)
    }))
    return(out)
  }
  
  # compute composition for scRNA
  files = Sys.glob(get.file("Data/scRNA*"))
  if(!file.exists(get.file("Results/scRNA_Shih_r.rds"))){
    files = Sys.glob(get.file("Data/scRNA*"))
    f <- function(path){
      r <- readRDS(path)
      out <- data.frame(r[
        c("cells", "samples",
          "patients", "sites",
          "sites_binary", "cell.types",
          "platform", "dataset")
      ]) %>%
        mutate(
          cell.types = factor(cell.types,
                              levels = ct))
      out <- compute_composition(r, out, ct)
      split <- strsplit(path, split = "/")[[1]]
      l <- length(split)
      saveRDS(out,
              paste0(c(split[1:(l-2)], 
                       "/Results/", split[l]), 
                     collapse = "/"))
      return("")
    }
    tmp <-  lapply(files, f)
  }
  files = Sys.glob(get.file("Results/scRNA*"))
  df_sc <- do.call("rbind", lapply(files, readRDS))
  
  # composition of discovery dataset (SMI)
  r<-readRDS(get.file("Data/SMI_data.rds"))
  out <- data.frame(data.frame(r[c("cells", "samples", "patients",
                                   "sites", "sites_binary",
                                   "cell.types")]),
                    platform = "ST",
                    dataset = "SMI") %>%
    mutate(cell.types = gsub("_LC", "", cell.types)) %>%
    mutate(cell.types = factor(cell.types, levels = ct))
  df_smi <- compute_composition(r, out, ct)
  remove(r)
  
  # composition of validation datatset 1 (ISS)
  q<-readRDS(get.file("/Data/Xenium_data.rds"))
  out <- data.frame(data.frame(q[c("cells", "samples",
                                   "patients", "sites",
                                   "sites_binary", "cell.types")]),
                    platform = "ST",
                    dataset = "Xenium") %>%
    mutate(cell.types = gsub("_LC", "", cell.types)) %>%
    mutate(cell.types = factor(cell.types, levels = ct))
  df_xen <- compute_composition(q, out, ct)
  remove(q)
  
  # composition of validation datatset 2 (MERFISH)
  s <- readRDS(get.file("/Data/MERFISH_data.rds"))
  out <- data.frame(data.frame(s[c("cells", "samples", "patients",
                                   "sites", "cell.types.nonmal")])) %>%
    mutate(sites_binary = c("Adnexa", "Omentum")[(sites == "Omentum") + 1],
           cell.types = factor(cell.types.nonmal, levels = ct)) %>%
    select(-cell.types.nonmal) %>%
    mutate(platform = "ST",
           dataset = "MERFISH")
  df_mer <- compute_composition(s, out, ct)
  remove(s)
  
  # combine everying together and format data for plotting
  full <- rbind(df_sc, df_smi, df_xen, df_mer)
  full <- full %>% mutate(dataset = factor(dataset,
                                           levels = c("SMI", "Xenium",
                                                      "MERFISH",
                                                      "Vazquez-García",
                                                      "Geistlinger",
                                                      "Olalekan",
                                                      "Qian",
                                                      "Regner",
                                                      "Shih"))) %>%
    mutate(count = as.numeric(count),
           proportion = as.numeric(proportion))
  order <- (full %>% filter(
    cell.types == "Malignant") %>% arrange(dataset, proportion))$samples
  full <- full %>% mutate(samples = factor(samples, levels = order),
                          cell.types = factor(
                            cell.types,
                            levels = c("Malignant", "Fibroblast",
                                       "Endothelial","Monocyte",
                                       "Mast.cell", "TNK.cell",
                                       "B.cell", "Other")))
  # new color map
  cell_2_rgb <- list()
  cell_2_rgb[["B.cell"]] = c(255, 0, 0)
  cell_2_rgb[["Fibroblast"]] = c(34, 91, 224)
  cell_2_rgb[["Endothelial"]] = c(166, 128, 71)
  cell_2_rgb[["Monocyte"]] = c(255, 153, 0)
  cell_2_rgb[["Mast.cell"]] = c(157, 3, 252)
  cell_2_rgb[["Malignant"]] = c(7, 224, 0)
  cell_2_rgb[["TNK.cell"]] = c(255, 247, 0)
  cell_2_rgb[["Other"]] = c(153, 153, 153)
  
  # plot to disk
  pdf(get.file("Figures/Fig1f.pdf"),
      width = 18, height = 3)
  p <- ggplot(full, aes(x = samples, y = proportion, fill = cell.types)) +
    geom_bar(stat = "identity") +
    facet_grid(~dataset, scales = "free_x", space = "free", 
               labeller = as_labeller(labels)) +
    scale_fill_manual(values = lapply(cell_2_rgb, rgb2hex)) +
    theme(axis.text.x = element_blank(),
          strip.text = element_text(color = "lightgrey"))
  print(p)
  g <- ggplot(full, aes(x = samples, y = dataset, fill = dataset)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = RColorBrewer::brewer.pal(9, "Set3")) +
    theme(legend.position = "top")
  as_ggplot(get_legend(g))
  dev.off()
  
  # extract from test 1 
  t1 <- readRDS(get.file("Data/SMI6K_Data_wSubtypes.rds"))
  dfraw = data.frame(t1[c("samples", "cell.types.conf")])
  t1ncells = data.frame(t1[c("samples", "cell.types.conf")]) %>% 
    group_by(samples) %>% summarize(n = length(cell.types.conf))
  t1nct = bind_rows(lapply(unique(t1ncells$samples), function(s){
    dff = filter(dfraw, samples == s)
    dff = dff %>% mutate(cell.types = factor(cell.types.conf, 
                                             levels = names(cell_2_rgb)))
    c(samples = s, table(dff$cell.types))
  }))
  short = merge(t1nct, t1ncells, by = "samples", all.x = T) 
  short1 = t(apply(short, 1, function(row){
    as.numeric(row[2:9])/as.numeric(row[10])
  }))
  colnames(short1) = colnames(short)[2:9]
  dfshort = data.frame(samples = short$samples, short1)
  dflong = dfshort %>% gather(key = cell.types, 
                              value = proportion, -samples) %>% 
    mutate(dataset = "Test 1")
  remove(t1)
  
  # extract from test 2 
  t2.1<- readRDS(get.file("Data/WT1_wSubtypes.rds"))
  t2.2 <- readRDS(get.file("Data/WT2_wSubtypes.rds"))
  t2.3 <- readRDS(get.file("Data/WT3_wSubtypes.rds"))
  t2.4 <- readRDS(get.file("Data/WT4_wSubtypes.rds"))
  df1short = bind_rows(lapply(list(t2.1, t2.2, t2.3, t2.4), function(s){
    c(samples = unique(s$patients), 
      table(s$cell.types.conf)/length(s$cell.types.conf))
  }))
  df1long = df1short %>% 
    gather(key = cell.types, value = proportion, -samples) %>% 
    mutate(dataset = "Test 2")
  
  # combine test 1 and test 2 
  dfpatch = rbind(dflong, df1long)
  dfpatch <- dfpatch %>% mutate(dataset = factor(dataset,
                                                 levels = c("Test 1", 
                                                            "Test 2"))) %>%
    mutate(proportion = as.numeric(proportion))
  order <- (dfpatch %>% filter(
    cell.types == "Malignant") %>% arrange(proportion))$samples
  
  full = dfpatch %>% mutate(samples = factor(samples, levels = order),
                            cell.types = factor(
                              cell.types,
                              levels = c("Malignant", "Fibroblast",
                                         "Endothelial","Monocyte",
                                         "Mast.cell", "TNK.cell",
                                         "B.cell", "Other")))
  # make plots 
  pdf(get.file("Figures/Fig1f-p.pdf"),
      width = 4, height = 3)
  p <- ggplot(full, aes(x = samples, y = proportion, fill = cell.types)) +
    geom_bar(stat = "identity") +
    facet_grid(~dataset, scales = "free_x", space = "free", 
               labeller = as_labeller(labels)) +
    scale_fill_manual(values = lapply(cell_2_rgb, rgb2hex)) +
    theme(axis.text.x = element_blank(),
          strip.text = element_text(color = "lightgrey"))
  print(p)
  dev.off()
}

#' Figure 1g. Hypergeometric Tests for Spatial Data
#'
#' For the Discovery dataset, we perform hypergeometric tests to evaluate if 
#' cell types co-occur in spatial frames more often than expected by random. 
#' 
#' @param cell_2_rgb a named list of colors assigned to each of the cell types
#' @return this function returns nothing, but writes figures in .pdf format. 
HGSC_Fig1g_hg <- function(cell_2_rgb){
  # divide samples into grids
  r <-readRDS(get.file("Data/SMI_data.rds"))
  r1<-get.frames(r,n1=300) 
  row.names(r1$frames.metadata) <- r1$frames.metadata$Frame
  
  # clean up the neighborhood cell type counts
  r1$frames.metadata<-as.data.frame(table(r1$frames))
  colnames(r1$frames.metadata)<-c("Frame","No_cells")
  row.names(r1$frames.metadata) <- r1$frames.metadata$Frame
  r1$frames.metadata$samples <- unlist(
    lapply(r1$frames.metadata$Frame, 
           function(x){
             paste(strsplit(as.character(x), 
                            split = "_")[[1]][1:3], 
                   collapse = "_")}))
  r1$cell.types.all <- gsub("_LC", "", r1$cell.types)
  df <- data.frame(r1[c("frames", "cell.types.all")])
  df <- df %>% 
    dplyr::group_by(frames, cell.types.all) %>%
    dplyr::summarize(n = length(cell.types.all)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(key = cell.types.all, n)
  df[is.na(df)] <- 0
  tme <- t(apply(df[2:8], 1, function(x) x/sum(x)))
  row.names(tme) <- df$frames
  r1$frames.tme <- tme
  
  # do the hypergeometic test
  attr <- spatial.co.occur(r1, type = "Att")
  rep <- spatial.co.occur(r1, type = "Rep")
  attr[(attr < 0.05)] <- "High"
  attr[(attr != "High")] <- "NS"
  attr[(rep < 0.05)] <- "Low"
  
  # make the plotting data 
  plt <- data.frame(attr) 
  plt_long <- plt %>% mutate(pairs = row.names(plt)) %>% 
    tidyr::gather(value = hg, key = "samples", -pairs) %>% 
    dplyr::mutate(hg = factor(hg, levels = c("High","NS","Low"))) %>% 
    dplyr::mutate(pairs = gsub("_", ", ",pairs)) %>% 
    dplyr::mutate(pairs = gsub("\\.", " ",pairs)) %>% 
    dplyr::mutate(pairs = gsub("TNK", "T/NK", pairs))
  order <- (plt_long %>% 
              dplyr::group_by(pairs) %>% 
              dplyr::summarize(High = sum(hg == "High"),
                               NS = sum(hg == "NS"), 
                               Low = sum(hg == "Low")) %>% 
              dplyr::arrange(Low, desc(High)))$pairs
  plt_long <- plt_long %>% dplyr::mutate(pairs = factor(pairs, levels = order))
  
  # plot to disk 
  p <- ggplot(plt_long, aes(x = pairs, fill = hg)) + 
    geom_bar(stat = "count", position = "stack") +
    coord_flip() + 
    scale_fill_manual(values = c("darkred", "grey", "dodgerblue4"), 
                      name = "Co-localization\nFrequency") + 
    theme_classic() + 
    ylab("# Samples") + 
    xlab("Cell Type Pairs")
  pdf(get.file("Figures/Fig1g.pdf"),
      width = 8, height = 5)
  print(p)
  dev.off()
}

#' Figure 1h. Co-localization Quotient 
#'
#' For the Discovery dataset, we calculate co-localization quotient (CLQ, 
#' Methods) between fibroblasts and T/NK cells and malignant cells and T/NK 
#' cells to explore the relationship of infiltration in the two predominant 
#' spatial components of the tumor tissue. 
#' 
#' @param cell_2_rgb a named list of colors assigned to each of the cell types
#' @return this function returns nothing, but writes figures in .pdf format.
HGSC_Fig1h_clq <- function(cell_2_rgb){
  r <-readRDS(get.file("Data/SMI_data.rds"))
  
  # compute co-localization quotient
  clq <- function(a, b, rs, field) {
    C_b_a <- sum(rs$neighbors[rs[[field]] == a,b], na.rm = T)
    N_a <- sum(rs[[field]] == a, na.rm = T)
    N_b <- sum(rs[[field]] == b, na.rm = T)
    if(N_b == 0 | N_a == 0) {
      return(NA)
    }
    N <- length(rs$cells)
    clq <- (C_b_a/N_a)/(N_b/(N - 1))
    return(clq)
  }
  
  # fetch CLQ Values for hypothesis
  fetch_clq_for_fib_mal <- function(r){
    out <- do.call("rbind", mclapply(
      unique(r$samples),
      mc.cores = 10,
      function(s){
        print(s)
        rs <- subset_list(r, r$cells[r$samples == s])
        fib <-  clq("Fibroblast", "TNK.cell", rs, "cell.types")
        mal <- clq("Malignant", "TNK.cell", rs, "cell.types")
        out <- c(samples = s,
                 patients = unique(rs$patients),
                 sites_binary = unique(rs$sites_binary),
                 CLQ_Fib_TNK = fib,
                 CLQ_Mal_TNK = mal)
        return(out)
      }))
    return(out)
  }
  
  # plot clq wrapper fx
  plot_clq_boxplots <- function(clq_df, suffix = ""){
    # do ttest just fyi
    ttest <- t.test(clq_df$CLQ_Fib_TNK, clq_df$CLQ_Mal_TNK, paired = T)
    clq_df <- clq_df[complete.cases(clq_df),] %>%
      mutate(Fibroblast = log2(CLQ_Fib_TNK)) %>%
      mutate(Malignant = log2(CLQ_Mal_TNK))
    clq_df <- clq_df[is.finite(clq_df$Fibroblast) &
                       is.finite(clq_df$Malignant),]
    plt <-  dplyr::select(clq_df, -CLQ_Fib_TNK, -CLQ_Mal_TNK) %>%
      gather(pair, clq, -patients, -sites_binary, -samples)
    
    plt <- filter(plt, patients %!in% filter(plt, clq < 1.5)$patients)
    
    # plot with wilcoxon sum rank test
    pdf(get.file("Figures/Fig1h.pdf"), width = 4, height = 6)
    p <- ggpaired(plt, x = "pair", y = "clq",
                  fill = "pair", line.color = "gray", line.size = 0.4)+
      stat_compare_means(comparisons = list(c("Fibroblast", "Malignant")),
                         paired = TRUE, method = "t.test", label = "p.signif") +
      xlab("") +
      ylab("Log2 Colocalization with T/NK Cells") +
      scale_fill_manual(values = lapply(
        cell_2_rgb[c("Fibroblast", "Malignant"
        )], rgb2hex)
      ) +
      theme(legend.position = "none")
    print(p)
    dev.off()
  }
  
  # add neighbors
  r$neighbors <- readRDS(get.file("Results/HGSC_SMI_Neighbors.rds"))
  subset = r$cells[r$cells %in% row.names(r$neighbors)]
  r1 <- subset_list(r, subset)
  
  # calculate clq
  out <- fetch_clq_for_fib_mal(r1)
  
  # make clq_df
  clq_df <- data.frame(out) %>%
    mutate(CLQ_Fib_TNK = as.numeric(CLQ_Fib_TNK),
           CLQ_Mal_TNK = as.numeric(CLQ_Mal_TNK))
  
  # plot clq boxplots
  plot_clq_boxplots(clq_df)
}
