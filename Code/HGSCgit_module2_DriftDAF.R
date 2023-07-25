#### Results Section 2 ###
# Figure 2. Malignant drift and desmoplasia associated fibroblasts (DAFs).
# Table S4. Desmoplasia Associated Fibroblast (DAF) program (A) and its Gene
# Ontology enrichment analysis (B). Associated with Figures 2.

# Figure 2A. Mixed Effects variance association across cell types
# Figure 2B. Malignant-only UMAP
# Figure 2C. Random Forest Analyses
# Figure 2D. Malignant drift boxplot by patient
# Figure 2E. Patient-specific Malignant UMAPs
# Figure 2F. Fibroblast-only UMAP
# Figure 2G. Boxplots of Desmoplasia Status
# Figure 2H. Heatmap of DAF program
# Figure 2I. Gene Ontology Analysis of DAF
# Figure 2J. DAF program in situ (composite with H&E made in Adobe Illustrator)

HGSC_Figure2_DriftDAF <- function(r = r.smi,
                                  sigs_cna,sigs,
                                  mal_umap, mal_drift,
                                  mal_rf, lopo_rf,
                                  morph, daf,
                                  cell_2_rgb){
  print("Section2")
  #1 Regenerate Figure 2A. Mixed Effects variance association across cell types
  print("Fig2A")
  HGSC_Fig2A_mixed_effects(sigs_cna = sigs_cna, sigs = sigs)
  #2 Regenerate Figure 2B. Malignant-only UMAP
  print("Fig2B")
  HGSC_Fig2B_mal_umap(mal_umap = mal_umap)
  #3 Regenerate Figure 2C. Random Forest Analyses
  print("Fig2C")
  HGSC_Fig2C_rf(mal_rf = mal_rf, lopo_rf = lopo_rf, cell_2_rgb = cell_2_rgb)
  #4 Regenerate Figure 2D: Patient pair-wise expression drift
  print("Fig2D")
  HGSC_Fig2D_mal_drift(mal_drift = mal_drift)
  #5 Regenerate Figure 2E: Select patients' UMAPs showing drift
  print("Fig2E")
  HGSC_Fig2E_patient_mal_umap(r=r, mal_drift = mal_drift)
  #6 Regenerate Figure 2F: Fibroblast UMAP with Site and Morphology Annotations
  print("Fig2F")
  HGSC_Fig2F_fib_umap(r=r, morph = morph)
  #7 Write DAF genes to Table S4A
  print("TableS4A")
  HGSC_TableS4A_write_daf(daf = daf)
  #8 Regenerate Figure 2G: Boxplots associating DAF to desmoplasia
  print("Fig2G")
  HGSC_Fig2G_desmo_boxplots(morph = morph, daf= daf)
  #9 Regenerate Figure 2H: DAF heatmap
  print("Fig2H")
  HGSC_Fig2H_daf_heatmap(r=r, morph = morph, daf = daf)
  #10 Regenerate Figure 2I: Gene Ontology Analysis of DAF
  print("Fig2I")
  HGSC_Fig2I_daf_go(r=r, daf =daf)
  #11 Write DAF full gene set results to Table S4B
  print("TableS4B")
  HGSC_TableS4B_write_daf_gsea(r=r,daf=daf)
  #12 Regenerate PNG Images for Figure 2J: DAF in situ
  print("Fig2J")
  HGSC_Fig2J_daf_in_situ(r=r, daf = daf, cell_2_rgb = cell_2_rgb)
  return()
}

HGSC_Fig2A_mixed_effects <- function(sigs_cna, sigs){
  # determine order of cell.types
  order <- c("Malignant", "Monocyte", "TNK.cell", "B.cell", "Mast.cell", "Endothelial",
             "Fibroblast")

  # calculate proportions
  cna <- do.call("rbind.data.frame", mclapply(sigs_cna, function(df){
    out <- apply(df[,1:3], 2, function(x){p.adjust(x, method = "BH")}) < 0.05
    out <- (apply(out, 2, sum)/159)*100
    out <- c(cna = out[1], celltypes = unique(df[,4]))
  }))
  colnames(cna) <- c("cna", "cell.types")

  # wrangle significant genes from sites and treatment associations
  others <- do.call("rbind.data.frame", mclapply(sigs, function(df){
    df <- df[!grepl("Neg", row.names(df)), ]
    out <- apply(df[,1:2], 2, function(x){p.adjust(x, method = "BH")}) < 0.05
    out <- (apply(out, 2, sum)/960)*100
    out <- c(out, celltypes = unique(df[,3]))
  }))
  colnames(others) <- c("sites", "treatment", "cell.types")

  # make plotting data.frame
  df <- cbind.data.frame(cna = cna[,1], others) %>%
    mutate(cna = as.numeric(cna),
           sites = as.numeric(sites),
           treatment = as.numeric(treatment)) %>%
    dplyr::rename(CNA = cna, Sites = sites, NACT = treatment) %>%
    mutate(cell.types = factor(cell.types, levels = order)) %>%
    gather(factor, variable, -cell.types)

  # make plot
  p <- ggplot(df, aes(x = cell.types, y = variable, fill = factor)) +
    geom_bar(stat = "identity") +
    facet_wrap(~factor, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("#90a6d6", "lightgrey", "#f5af5b"),
                      name = "Effect",
                      labels = c("CNA", "NACT", "Site")) +
    theme_bw() +
    ylab("% of genes significantly associated") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))

  # plot to disk
  pdf(get.file("Figures/Fig2A.pdf"), width = 5, height = 4)
  print(p)
  dev.off()
}

HGSC_Fig2B_mal_umap <- function(mal_umap){
  # fetch embeddings from Seurat Object
  umap_malig <- mal_umap@reductions$umap@cell.embeddings

  # make plotting data frame
  plt <- cbind(umap_malig, mal_umap@meta.data) %>%
    mutate(patients = gsub("HGSC", "", patients))

  # make plot
  p <- ggplot() +
    geom_point(data = plt,
               aes(x = UMAP_1, y = UMAP_2, col = patients), size =0.05, alpha = 0.8) +
    theme_bw() +
    theme(axis.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=3)))

  # write to disk
  pdf(get.file("Figures/Fig2B.pdf"),
      width = 7,
      height = 5.5)
  print(p)
  dev.off()
}

HGSC_Fig2C_rf <- function(mal_rf, lopo_rf, cell_2_rgb){
  # read and compile data for each cell.type
  cell_types <- c("B.cell", "Endothelial", "Monocyte",
                  "TNK.cell", "Malignant", "Fibroblast")

  # extract roc statistics per model
  extract_statistics <- function(df){
    rocobj <- roc(df$sites, df$pred_score, smooth = F)
    roc <- data.frame(sensitivity =rocobj$sensitivities, specificity = rocobj$specificities)
    roc$celltype = rep(unique(df$celltype), dim(roc)[1])
    roc$auc = rep(rocobj$auc, dim(roc)[1])
    return(roc)
  }

  # extract statistics and compile into data.rame
  stats <- do.call("rbind", lapply(lopo_rf, extract_statistics))

  # make labels for leave-one-out cross validation
  labels <- select(stats, celltype, auc) %>% unique()
  row.names(labels) <- labels$celltype
  labels <- labels[names(cell_2_rgb),]
  labels$celltype <- names(cell_2_rgb)
  labels <- mutate(labels, label = paste0(celltype, " (AUC = ", signif(auc, 4), ")"))

  # plot the ROC for the leave one patient out analysis
  p1 <- ggplot(stats, aes(x = specificity, y = sensitivity, col = celltype)) +
    scale_x_reverse() +
    geom_line(lwd = 1)  +
    coord_fixed() +
    scale_color_manual(values = unlist(lapply(cell_2_rgb, rgb2hex)),
                       labels = labels$label, name = "") +
    theme_bw() +
    geom_abline(slope = 1, intercept = 1, col = "lightgrey") +
    theme(legend.position = c(0.75, 0.3)) +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          panel.border = element_rect(size = 1),
          legend.text = element_text(size = 10.5),
          legend.title = element_text(size = 1)) +
    ylab("Sensitivity") +
    xlab("Specificity")

  # generate colors for each patient in the patients rf analysis
  n <- 36
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal,
                             qual_col_pals$maxcolors,
                             rownames(qual_col_pals)))
  names(col_vector) <- names(mal_rf)

  # make plot for leave one out cross validation per patient
  p2 <- ggroc(mal_rf) +
    guides(color=guide_legend(title="Patients", ncol = 3)) +
    coord_fixed() +
    theme_bw() +
    geom_abline(slope = 1, intercept = 1, color = "lightgrey") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          panel.border = element_rect(size = 1),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12)) +
    scale_color_manual(values = col_vector) +
    theme(legend.position = c(0.6, 0.4))

  # plot to disk
  pdf(get.file("Figures/Fig2C.pdf"), width = 11, height = 5.5)
  print(p2 + p1)
  dev.off()
}

HGSC_Fig2D_mal_drift <- function(mal_drift){
  plt <- mal_drift
  plt$patients <- factor(plt$patients,
                         levels = (mal_drift %>%
                                     group_by(patients) %>%
                                     summarize(dist_med = median(pairwise_dist)) %>%
                                     arrange(desc(dist_med)))$patients)
  meds <- mal_drift %>%
             group_by(patients) %>%
             summarize(dist = median(pairwise_dist)) %>%
             arrange(desc(dist))

  plt2 <- filter(plt) %>%
    merge(meds,
          by = "patients", all.x = T)

  p <- ggplot(plt2, aes(x = patients, y = pairwise_dist, fill = dist)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size = 12)) +
    theme(axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14)) +
    xlab("Patients") +
    ylab("Drift Magnitude") +
    scale_fill_gradient2(low = '#009392',
                         mid = '#f6edbd',
                         high = '#A01C00', midpoint = median(plt2$dist),
                         name = "Drift\nMagnitude")

  pdf(get.file("Figures/Fig2D.pdf"),
      width = 9, height = 3)
  print(p)
  dev.off()
}

HGSC_Fig2E_patient_mal_umap <- function(r,
                                        mal_drift,
                                        pats = c("HGSC38",
                                                     "HGSC102",
                                                     "HGSC2",
                                                     "HGSC58")){
  # subset patients
  set.seed(1234)
  sub <- function(x){
    cells <- r$cells[r$samples == x & r$cell.types == "Malignant" &
                       as.numeric(r$cons$conf.score.cell.type) >= 0.95 &
                       !is.na(as.numeric(r$cons$conf.score.cell.type))]
    if (length(cells) < 500) {
      n <- length(cells)
    } else {
      n <- 500
    }
    cells <- sample(cells, n)
    return(cells)
  }
  cell_idx <- do.call("c", lapply(unique(r$samples), sub))
  q <- subset_list(r, cell_idx)

  # code to plot the umap per patient
  plot_pat <- function(pat){
    print(pat)
    s <- subset_list(q, q$cells[q$patients == pat])
    if (length(unique(s$samples)) > 2){
      s <- subset_list(s, s$cells[s$samples %in% samples])}
    so <- seuratify(s$tpm[s$genes[!grepl("Neg", s$genes)],], prj_name = "malig")
    so <- transfer_data_list_to_so(s, so, c("samples", "sites_binary",
                                            "patients", "cell.types.broad",
                                            "treatment"))
    so <- SCTransform(so, verbose = F)
    so <- FindVariableFeatures(so, nfeatures = 960, verbose = F)
    so <- RunPCA(so, npcs = 30, verbose = F)
    so <- RunUMAP(so, dims = 1:30, verbose = F)

    dist = median(filter(tmp, patients == pat)$pairwise_dist)
    sd = sd(filter(tmp, patients == pat)$pairwise_dist)

    cols = c("salmon", "dodgerblue")
    p = DimPlot(so, group.by = "sites_binary") +
      ggtitle(paste0(pat,
                     " (d = ", format(dist, digits = 3),
                     " +/- ", format(sd, digits = 2), ")")) +
      scale_color_manual(values = cols, na.value = "#DFDFDF") +
      guides(colour = guide_legend(override.aes = list(size=3)))
    print(p)
    return(pat)
  }

  # subsetting out only patients of interest
  tmp <- filter(mal_drift, patients %in% pats)
  patorder <- (tmp %>%
                 group_by(patients) %>%
                 summarize(dist = median(pairwise_dist)) %>%
                 arrange(desc(dist)))$patients

  #plot to disk
  pdf(get.file("Figures/Fig2E.pdf"),
      width = 7,
      height = 5)
  out <- lapply(patorder, plot_pat)
  dev.off()

}
HGSC_Fig2F_fib_umap<- function(r, morph){
  set.seed(1234)

  #subset out fibroblasts
  q <- subset_list(r, r$cells[r$cons$conf.score.cell.type > 0.95 &
                                r$cell.types == "Fibroblast"])
  meta <- morph[q$samples,]
  row.names(meta) <- q$cells

  # make umap
  so <- CreateSeuratObject(q$tpm[!grepl("NegPrb", r$genes), ], meta.data = meta)
  so <- ScaleData(so)
  so <- FindVariableFeatures(so, nfeatures = 920)
  so <- RunPCA(so, npcs = 30)
  so <- RunUMAP(so, dims = 1:30)

  # plot umap with sites annotations
  a <- DimPlot(so, group.by = "sites_binary",
               na.value = brewer.pal(n = 9, name = "Set3")[9]) +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank()) +
    scale_color_manual(values = c("dodgerblue", "salmon"),
                       na.value = brewer.pal(n = 9, name = "Set3")[9])
  b <- DimPlot(so, group.by = "type") +
    scale_color_manual(values = brewer.pal(n = 8,
                                           name = "Set3")[c(5,7:8)],
                       na.value = brewer.pal(n = 9, name = "Set3")[9]) +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank())

  pdf(get.file("Figures/Fig2F.pdf"), width = 10, height =5)
  print(a + b)
  dev.off()
}

HGSC_TableS4A_write_daf <- function(daf){
  ord1 <- order(union(daf$Fibroblast.paired.up, daf$Fibroblast.unpaired.up))
  up <- union(daf$Fibroblast.paired.up, daf$Fibroblast.unpaired.up)[ord1]
  ord2 <- order(union(daf$Fibroblast.paired.down, daf$Fibroblast.unpaired.down))
  down <- union(daf$Fibroblast.paired.down, daf$Fibroblast.unpaired.down)[ord2]

  daf_sig <- c()
  daf_sig$DAF_up <- up
  daf_sig$DAF_down <- down
  write.xlsx(x = as.data.frame(list.2.mat(daf_sig)),
                               file = get.file("Tables/TableS4A.xlsx"),
                               asTable = T)
}

HGSC_Fig2G_desmo_boxplots <- function(morph, daf){
  df <- morph
  deg <- daf$oe
  df$type = factor(df$type , levels = c("OvarianStroma+\nDesmoplasia-",
                                        "OvarianStroma+\nDesmoplasia+",
                                        "OvarianStroma-\nDesmoplasia+"))
  df <- deg %>% group_by(samples) %>% summarise_all(mean) %>%
    merge(df)

  df$desmoplasia <- factor(df$desmoplasia, levels = c(0, 1))
  a = ggboxplot(df[complete.cases(df),],
                x = "desmoplasia",
                y = "om",
                fill = "desmoplasia") +
    stat_compare_means(comparisons = list(c("0", "1")),
                       method = "wilcox.test",
                       method.args = list(alternative = "less"),
                       ref.group = "0",
                       label = "p.signif") +
    theme_classic() +
    theme(legend.position = "none") +
    stat_boxplot(geom = "errorbar") +
    geom_boxplot(data = df[complete.cases(df),], aes(fill = desmoplasia)) +
    xlab("Desmoplasia") +
    ggtitle("All Samples") +
    ylab("DAF Score") +
    theme(axis.text.x = element_text(vjust =0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "none")

  b <- ggboxplot(filter(df, sites_binary == "Adnexa"),
                 x = "type",
                 y = "om",
                 fill = "type") +
    stat_compare_means(comparisons = list(
      c("OvarianStroma-\nDesmoplasia+",
        "OvarianStroma+\nDesmoplasia+"),
      c("OvarianStroma-\nDesmoplasia+",
        "OvarianStroma+\nDesmoplasia-")),
      method = "wilcox.test",
      method.args = list(alternative = "greater"),
      ref.group = "OvarianStroma+\nDesmoplasia-",
      label = "p.signif") +
    stat_boxplot(geom = "errorbar") +
    geom_boxplot(aes(fill = type)) +
    scale_fill_manual(values = c("#FCCDE5","#B3DE69","#80B1D3")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust =0.5, size = 10)) +
    ylab("DAF Score") +
    xlab("") +
    ggtitle("Adnexa Samples Only") +
    theme(axis.text.x = element_text(angle = 90, vjust =0.5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "none")

  pdf(get.file("Figures/Fig2G.pdf"),
      width = 6, height = 5)
  print(a+b)
  dev.off()
}

HGSC_Fig2H_daf_heatmap <- function(r, morph, daf){
  subcells = r$cells[r$cell.types == "Fibroblast"]
  q <- subset_list(r, subcells)

  sig <- list("om.up" = daf$Fibroblast.paired.up,
              "om.down" = daf$Fibroblast.paired.down)
  sig <- list("om.up" = intersect(daf$Fibroblast.paired.up,
                                  daf$Fibroblast.unpaired.up),
              "om.down" = intersect(daf$Fibroblast.paired.down,
                                    daf$Fibroblast.unpaired.down))

  # calculate overall expression of the omentum signature in each cell
  q$oe <- get.OE(q, sig)
  q$oe <- cap_object(q$oe, 0.05)
  q$oe <- data.frame(q$oe)
  q$oe$samples <- q$samples

  # calculate mean overall expression
  oe_mean <- q$oe %>%
    group_by(samples) %>%
    summarize_all('mean') %>%
    arrange(desc(om))

  # scale, center, cap expression matrix
  q$zscores_fibroblast <- cap_object(q$tpm, 0.05)
  q$zscores_fibroblast <- scale_and_center(q$zscores_fibroblast, 1)

  mat_df <- data.frame(t(q$zscores_fibroblast), samples = q$samples) %>%
    group_by(samples) %>%
    summarize_all("mean")
  mat <- data.matrix(dplyr::select(mat_df, -samples))
  row.names(mat) <- mat_df$samples
  mat <- cap_object(mat, q = 0.05)

  # top correlated genes with omentum signature
  gene.cor<-spearman.cor(mat[oe_mean$samples,],oe_mean$om)
  sig1<-get.top.cor(gene.cor,q = 25,idx = "R")
  gene<-c(sig1$R.up,sig1$R.down)

  # format annotations
  df <- morph
  df[is.na(df)] <- "Unlabeled"
  ann <- unique(data.frame(samples = q$samples, sites = q$sites_binary,
                           morph = df[q$samples,]$type))
  ann <- merge(oe_mean, ann, by = "samples")
  ann[is.na(ann)] <- "Unlabeled"
  row.names(ann) <- ann$samples
  ann <- ann %>% dplyr::select(om, sites, morph)
  ann_colors = list(sites = c(Omentum = "salmon", Adnexa = "dodgerblue"),
                    direction = c(up = "#E10000", down = "#4782C9"),
                    morph = c(`OvarianStroma-\nDesmoplasia+` = "#80B1D3",
                              `OvarianStroma+\nDesmoplasia+`="#B3DE69",
                              `OvarianStroma+\nDesmoplasia-`="#FCCDE5",
                              `Unlabeled` = "#D9D9D9"))
  ann_col <- data.frame(direction = c(
    rep("up", length(sig1$R.up)), rep("down", length(sig1$R.down))
  ))
  row.names(ann_col) <- gene
  sample_order <- row.names(ann %>%
                              mutate(sites = factor(sites,
                                                    levels = c("Omentum",
                                                               "Adnexa"))) %>%
                              arrange(desc(om), sites))

  # plot to disk
  pdf(get.file("Figures/Fig2H.pdf"),
      width = 10, height = 5)
  plt <- pheatmap(mat[oe_mean$samples,gene],
                  cluster_rows = F,
                  cluster_cols = F,
                  annotation_row = ann,
                  show_rownames = F,
                  show_colnames = T,
                  annotation_colors = ann_colors,
                  border_color = NA,
                  annotation_col = ann_col)
  dev.off()
}

HGSC_Fig2I_daf_go <- function(r, daf){
  # calculate go enrichment for up direction.
  results1 <- (gprofiler2::gost(union(daf$Fibroblast.paired.up,
                                      daf$Fibroblast.unpaired.up),
                                custom_bg = r$genes[!grepl("NegPrb",
                                                           r$genes)],
                                user_threshold = 0.05,
                                correction_method = "fdr"))

  # filter results for plotting
  out <- filter(results1$result, grepl("GO:", source), term_size <300) %>%
    mutate(neglog10 = -log10(p_value)) %>%
    arrange(desc(neglog10)) %>%
    filter(term_id %!in% unlist(results1$result$parents)) %>%
    arrange(desc(intersection_size)) %>%
    filter(!grepl("development", term_name)) %>%
    filter(!grepl("sarcolemma", term_name))
  out <- out[1:8,] %>% arrange(neglog10)
  out <- out %>% mutate(out, term = paste0(term_name, "\n(", term_id, ")"))
  out$term <- factor(out$term, levels = c(out$term))
  out$intersection_size <- factor(unlist(lapply(out$intersection_size,
                                                function(x){
                                                  if(x < 5) {
                                                    return("<5")
                                                  } else if ( x < 10) {
                                                    return("6-10")
                                                  } else if (x < 20) {
                                                    return("11-20")
                                                  } else {
                                                    return(">21")
                                                  }
                                                })), levels = c("<5", "6-10", "11-20", ">21"))

  # make plot for up direction
  p1 <- ggplot(out, aes(x = neglog10, y = term, col = source)) +
    geom_point(aes(size = intersection_size)) +
    theme_classic() +
    guides(colour = guide_legend("GO Type" , override.aes = list(
      fill="white",
      size = 4,
      shape = 15)),
      size = guide_legend("No. Genes"))+
    xlab("-log10(p_value)") +
    ggtitle("Up in Omentum") +
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90),
          axis.title.x = element_text(size = 12))

  # calculate go enrichment for down direction.
  results2 <- (gprofiler2::gost(union(daf$Fibroblast.paired.down,
                                      daf$Fibroblast.unpaired.down),
                                custom_bg = r$genes[!grepl("NegPrb",
                                                           r$genes)],
                                user_threshold = 0.05,
                                correction_method = "fdr"))

  out <- filter(results2$result) %>%
    mutate(neglog10 = -log10(p_value)) %>%
    arrange(neglog10) %>%
    filter(term_id %!in% unlist(results2$result$parents),
           intersection_size >= 2,
           source != "TF")
  out <- out[c(1,3,5,6,10,11,12,22),]
  out <- out %>% mutate(term = paste0(term_name, "\n(", term_id, ")"))
  out <- out %>% arrange(neglog10)

  out$term <- factor(out$term, levels = c(out$term))
  out$intersection_size <- factor(unlist(lapply(out$intersection_size,
                                                function(x){
                                                  if(x < 5) {
                                                    return("<5")
                                                  } else if ( x < 10) {
                                                    return("6-10")
                                                  } else if (x < 20) {
                                                    return("11-20")
                                                  } else {
                                                    return(">21")
                                                  }
                                                })), levels = c("<5", "6-10", "11-20", ">21"))

  # make plot for down direction
  p2 <- ggplot(out, aes(x = neglog10, y = term, col = source)) +
    geom_point(aes(size = as.factor(intersection_size))) +
    theme_classic() +
    guides(
      size = guide_legend("No. Genes")) +
    scale_size_manual(values = c(1, 3, 5)) +
    xlab("-log10(p_value)") +
    ggtitle("Down in Omentum") +
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90),
          axis.title.x = element_text(size = 12))
  print(p2)

  # plot to disk
  pdf(get.file("Figures/Fig2I.pdf"), width =13, height = 4)
  print(p1 + p2)
  dev.off()
}

HGSC_TableS4B_write_daf_gsea <- function(r, daf){
  resultsout <- (gprofiler2::gost(union(daf$Fibroblast.paired.up,
                                        daf$Fibroblast.unpaired.up),
                                  organism = "gp__w8SG_tnEY_W4Q",
                                  custom_bg = r$genes[!grepl("NegPrb",
                                                             r$genes)],
                                  user_threshold = 0.05,
                                  correction_method = "fdr"))
  up <- dplyr::select(resultsout$result,
                term_id, query_size, precision, recall, p_value) %>%
    mutate(direction = "up")

  resultsout2 <- (gprofiler2::gost(union(daf$Fibroblast.paired.down,
                                         daf$Fibroblast.unpaired.down),
                                   organism = "gp__w8SG_tnEY_W4Q",
                                   custom_bg = r$genes[!grepl("NegPrb",
                                                              r$genes)],
                                   user_threshold = 0.05,
                                   correction_method = "fdr"))

  down <- dplyr::select(resultsout2$result,
                      term_id, query_size, precision, recall, p_value) %>%
    mutate(direction = "down")

  write.xlsx(x = rbind(up, down),
             file = get.file("Tables/TableS4B.xlsx"),
             asTable = T)
}

HGSC_Fig2J_daf_in_situ <- function(r, daf, cell_2_rgb, s = "SMI_T11_F019"){
    # subset cells for this sample's fibroblasts
    q <- subset_list(r, r$cells[r$samples == s &
                                  r$cell.types == "Fibroblast"])
    daf_scores <- daf$oe[q$cells, "om"]

    # calculate pvalue of daf in situ for this sample
    pval_df <- data.frame(daf_scores, morph = c("norm","desmo")[unlist(apply(q$coor, 1, function(row){
      row["y"] > (0.23*row["x"] + 1500)
    })) + 1], q$coor)
    out <- wilcox.test(filter(pval_df, morph == "norm")$daf_scores,
                       filter(pval_df, morph == "desmo")$daf_scores, alternative = "less")
    print(paste0("p-value for daf in situ: ", out$p.value))

    # make color bar legend
    names(daf_scores) <- q$cells
    daf_scores <- scale(daf_scores)
    daf_scores <- cap_object(daf_scores, 0.01)
    p <- ggplot(data.frame(daf_scores, q$coor), aes(x = x,
                                                    y =y,
                                                    col = daf_scores)) +
      geom_point() +
      scale_color_gradient2(low = '#009392', mid = '#f6edbd', high = '#A01C00',
                            midpoint = 0.8) +
      coord_fixed() +
      geom_abline(slope = 0.23, intercept =  1500) + theme_classic()
    daf_colors <- data.frame(ggplot_build(p)$data[[1]])[,1]
    names(daf_colors) <- row.names(daf_scores)
    leg <- ggpubr::as_ggplot(ggpubr::get_legend(p))

    # plot the
    seg_folder = get.file("Results/Segmentation/")
    seg_suffix = "_whole-cell_03.csv"
    cell2rgb <- list("Other" = c(200,200,200),
                     "Fibroblast" = c(34, 91, 224))

    # plot for each sample
    seg_path = paste0(seg_folder, s, seg_suffix)

    # subset out cells
    cells <- r$cells[r$samples == s]
    q <- subset_list(r, cells)
    celltypes <- q$cell.types
    names(celltypes) <- q$cells
    celltypes[celltypes != "Fibroblast"] <- "Other"

    # set up the colors of the continuous values
    contvals <- daf_colors[names(celltypes)]

    # run visualization for DAF
    spatial_sample_visualization(seg_path,
                                 celltypes,
                                 cell2rgb,
                                 s,
                                 background = "white",
                                 outpath = outpath,
                                 outfile = get.file("Figures/Fig2J-2.png"),
                                 cont_field = "Fibroblast",
                                 contvals = contvals)

    # run visualization for cell types for reference
    q <- subset_list(r, r$cells[r$samples == s])
    q$cell.types[grepl("_LC", q$cell.types)] <- "LC"
    cell_2_rgb$LC <- c(0,0,0)
    spatial_sample_visualization(seg_path,
                                 celltypes = q$cell.types,
                                 cell_2_rgb,
                                 s,
                                 background = "black",
                                 outpath = outpath,
                                 low_qc_color = 0,
                                 outfile = get.file("Figures/Fig2J-1.png"))

}

