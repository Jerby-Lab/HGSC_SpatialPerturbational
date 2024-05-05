#### Results Section 1 ###
# Figure 1. Single cell spatial transcriptomics (ST) mapping of HGSC.
# Supplementary Table 1. Specifications of datasets collected and/or analyzed in this study. Associated with Figure 1.
# Supplementary Table 2. Metadata for each tissue section profiled by spatial transcriptomics. Associated with Figure 1.
# Supplementary Table 3. Cell Type Signature Genes. Includes both signature genes derived from scRNA-Seq and CellTypist Immune Encyclopedia (A) and from spatial transcriptomics (B). Associated with Figures 1, S1, S2.

# Figure 1a. Spatial Cohort Design (generated in Figma, with Lettie Mcguire)
# Figure 1b. CoMut Plot (output files in Results/ used for python comut code)
# Figure 1c. Cell Type UMAPs
# Figure 1d. Cell Types in situ (composite with H&E images made in Adobe Illustrator)
# Figure 1e. Coembedding scross ST and scRNA-seq datasets
# Figure 1f. Cell Type Compositions
# Figure 1g. Hypergeometric Spatial Mapping
# Figure 1h. Cell Co-localization Quotient Analysis

HGSC_Figure1_SpatiomolecularMapping <- function(d, v1, v2, t1, t2, cell_2_rgb){
  print("Section1")
  print("Fig 1a is generated outside of R")
  #1 Regenerate input files for CoMut plot in python
  print("Fig 1b")
  master <- HGSC_Fig1b_make_comut(d, v1, v2, t1, t2)
  #6 Regenerate Figure 1c: Cell Type UMAPs
  print("Fig 1c")
  HGSC_Fig1c_celltype_umaps(r=r,cell_2_rgb = cell_2_rgb)
  #9 Regenerate Figure 1d: Cell Types in situ
  print("Fig 1d")
  HGSC_Fig1d_celltypes_insitu(r=r, q=q, s=s, cell_2_rgb=cell_2_rgb)
  #10 Regenerate Figure 1e: Coembedding scross ST and scRNA-seq datasets
  print("Fig 1e")
  HGSC_Fig1e_coembedding(cell_2_rgb)
  #11 Regenerate Figure 1f: Cell Type Compositions
  print("Fig 1f")
  HGSC_Fig1f_celltype_composition(r=r, q=q, s=s, cell_2_rgb = cell_2_rgb)
  #12 Regenerate Figure 1g
  print("Fig 1g")
  # HGSC_Fig1G_hypergeometric(r=r)
  #13 Regenerate Figure 1h: Cell Co-localization Quotient Analysis
  print("Fig 1h")
  HGSC_Fig1h_clq(r=r, cell_2_rgb = cell_2_rgb)
  
  return()
}

HGSC_Fig1b_make_comut <- function(){
  master <- readRDS(get.file("Data/MasterProfileData.rds"))
  
  # fix immunotherapy
  master$immunotherapy[grepl("embro", master$immunotherapy) | grepl( "varlilumab" , master$immunotherapy)] <- "Yes"
  master$immunotherapy[master$immunotherapy == "NA"] <- NA
  master$immunotherapy[grepl("placebo", master$immunotherapy)] <- NA
  master$immunotherapy[master$immunotherapy != "Yes"] <- "No"
  master$immunotherapy <- factor(master$immunotherapy, levels = c("No", "Yes"))
  
  # parpi
  master$parpi <- master$parpi
  master$parpi[master$parpi == 3] <- NA
  master$parpi[master$parpi == 2 | master$parpi == 4] <- 1
  master$parpi <- factor(master$parpi, levels = c(0, 1))
  
  # fix bevacizumab 
  master$`1L`[master$`1L` == 5] <- NA
  master$`1L`[master$`1L` == 2] <- NA
  master$`1L` <- c("No", "", "Yes")[master$`1L`]
  master$`1L` <- factor(master$`1L`, levels = c("No", "Yes"))
  
  # CLINICAL
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
  
  df <- df[-4]
  colnames(df) <- c("Patient", "Age", "Stage", "PARPi", "NACT",
                    "Immunotherapy", "Outcome", "PFS","Bevacizumab")
  
  for (i in 2:9) {
    selected <- df[,c(1,i)]
    field <- colnames(df)[i]
    print(field)
    selected$categoory <- field
    out <- selected[,c(1, 3, 2)]
    colnames(out) <- c("sample", "category", "value")
    write.table(out, paste0(outdir, field,
                            "_CoMut_clin.tsv"), sep = "\t",
                row.names = F, quote = F)
  }
  
  # MUTATIONS
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
  
  # combine mutations
  df_mut <- rbind(brca, tp53 %>% filter(!is.na(mutation)))
  write.table(df_mut, paste0(outdir, "CoMut_Muts.tsv"),
              sep = "\t", row.names = F, quote = F)
  
  # TEMPUS
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
  
  # ST
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

HGSC_Fig1c_celltype_umaps <- function(r, cell_2_rgb){
  so_smi <- readRDS(get.file("Results/HGSC_SMI_wUMAP.rds"))
  so_iss <- readRDS(get.file("Results/HGSC_ISS_wUMAP.rds"))
  so_mer <- readRDS(get.file("Results/HGSC_MERFISH_wUMAP.rds"))
  
  plot_ct_umap <- function(r, so, cell_2_rgb, title){
    col_field <- names(
      cell_2_rgb
    )[names(cell_2_rgb) %in% unique(r$cell.types)]
    p <- DimPlot(so, group.by = "cell.types.nonmal", pt.size = 0.1) +
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
  
  a = plot_ct_umap(r, so_smi, cell_2_rgb, "Discovery Dataset")
  b = plot_ct_umap(r, so_iss, cell_2_rgb, "Validation Dataset 1")
  c = plot_ct_umap(r, so_mer, cell_2_rgb, "Validation Dataset 2")
  
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
  d = as_ggplot(get_legend(p))
  
  pdf(get.file("Figures/Fig1c.pdf"), width = 10, height = 10)
  p <- ggarrange(a,b,c,d)
  print(p)
  dev.off()
}

HGSC_Fig1d_celltypes_insitu<- function(r, q, s, cell_2_rgb){
  # plot 4 samples from the discovery dataset
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
  
  # plot the select sample from the validation 1 dataset (ISS)
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
  
  # plot select sample from valifation 2 dataset (MERFISH)
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
}

HGSC_Fig1e_coembedding <- function(cell_2_rgb){
  coembedding <- readRDS(get.file("Results/HGSC_coembedding.rds"))
  
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
  coembedding$dataset <- factor(coembedding$dataset, levels = c("SMI", "ISS", "MERFISH",
                                                                "Vazquez-García", "Geistlinger",
                                                                "Olalekan", "Regner", "Qian",
                                                                "Shih", "hornburg"))
  coembedding <- coembedding[sample(1:(dim(coembedding)[1])),]
  b <- ggplot(coembedding[!is.na(coembedding$cell.types),], aes(x = UMAP_1, y = UMAP_2,
                                                                col = dataset) ) +
    geom_point(size = 0.05) +
    scale_color_manual(values = RColorBrewer::brewer.pal(9, "Set3")) +
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
      width = 6, height = 4.5)
  print(a);print(b)
  dev.off()
}

HGSC_Fig1f_celltype_composition <- function(r, q, s, cell_2_rgb,
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
        cbind.data.frame(count =  unlist(unname(as.list(table(filter(out, samples == s)$cell.types))))) %>%
        mutate(proportion = count/(dim(filter(out, samples == s))[1]) )
      row.names(meta) <- NULL
      return(meta)
    }))
    return(out)
  }
  # compute composition for scRNA
  files = Sys.glob(get.file("Data/scRNA*"))
  df_sc <- do.call("rbind", mclapply(files, mc.cores = 10,
                                     function(path){
                                       print(path)
                                       start <- Sys.time()
                                       r <- readRDS(path)
                                       out <- data.frame(r[
                                         c("cells", "samples",
                                           "patients", "sites",
                                           "sites_binary", "cell.types",
                                           "platform", "dataset")
                                       ]) %>%
                                         mutate(
                                           cell.types = factor(cell.types,
                                                               levels = ct)
                                         )
                                       out <- compute_composition(r, out, ct)
                                       return(out)
                                       print(Sys.time() - start)
                                     }))
  
  # composition of discovery dataset (SMI)
  out <- data.frame(data.frame(r[c("cells", "samples", "patients",
                                   "sites", "sites_binary",
                                   "cell.types")]),
                    platform = "ST",
                    dataset = "SMI") %>%
    mutate(cell.types = gsub("_LC", "", cell.types)) %>%
    mutate(cell.types = factor(cell.types, levels = ct))
  df_smi <- compute_composition(r, out, ct)
  
  # composition of validation datatset 1 (ISS)
  out <- data.frame(data.frame(q[c("cells", "samples",
                                   "patients", "sites",
                                   "sites_binary", "cell.types")]),
                    platform = "ST",
                    dataset = "Xenium") %>%
    mutate(cell.types = gsub("_LC", "", cell.types)) %>%
    mutate(cell.types = factor(cell.types, levels = ct))
  df_xen <- compute_composition(q, out, ct)
  
  # composition of validation datatset 2 (MERFISH)
  out <- data.frame(data.frame(s[c("cells", "samples", "patients",
                                   "sites", "cell.types.nonmal")])) %>%
    mutate(sites_binary = c("Adnexa", "Omentum")[(sites == "Omentum") + 1],
           cell.types = factor(cell.types.nonmal, levels = ct)) %>%
    select(-cell.types.nonmal) %>%
    mutate(platform = "ST",
           dataset = "MERFISH")
  df_mer <- compute_composition(s, out, ct)
  
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
    facet_grid(~dataset, scales = "free_x", space = "free", labeller = as_labeller(labels)) +
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
}

HGSC_Fig1h_clq <- function(r, cell_2_rgb){
  # Compute Colocalization Quotient
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
    print(ttest)
    clq_df <- clq_df[complete.cases(clq_df),] %>%
      mutate(Fibroblast = log2(CLQ_Fib_TNK)) %>%
      mutate(Malignant = log2(CLQ_Mal_TNK))
    clq_df <- clq_df[is.finite(clq_df$Fibroblast) & is.finite(clq_df$Malignant),]
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
  subset <- row.names(r$neighbors[complete.cases(r$neighbors),])
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
