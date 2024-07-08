#### Results Section 2 ###
# Figure 2. T/NK cell states reflect their tumor infiltration status.
#
# Figure 2a. Discovery CD8 T cell umaps
# Figure 2b. Discovery CD8 T cell TIP genes' association with tumor infiltration in different cell subtypes.
# Figure 2c. Validation 1 dataset spatial maps with CD8 TIP
# Figure 2d. Validation 1 CD8 T cell states
# Figure 2e. Discovery Heatmap of Desmoplastic Fibroblast Program
# Figure 2f. Discovery Desmoplastic Fibroblast Program in situ (composite with H&E made in Adobe Illustrator)
# Figure 2g. Spatial L-R Network (generated via Cytoscape)

#' HGSC Figure 2 T/NK cell states reflect their tumor infiltration status.
#' Generates various plots and analyses for reproducing Figures 2a-f
#' @param r.smi Discovery data.
#' @param rCD8.smi Discovery data, CD8 T cells.
#' @param r.xenium Validation 1 data.
#' @param rTNK.xenium Validation 1 data, T/NK cells.
#' @param rCD8.xenium Validation 1 data, CD8 T cells.
#' @param daf Desmoplastic Fibroblast results. 
#' @param cell_2_rgb Color mapping for cells.
#' @return None. This function prints plots to disk. 
HGSC_Figure2_TIP_DF<-function(r.smi ,rCD8.smi,r.xenium,rTNK.xenium,rCD8.xenium,daf){
  
  #1 Regenerate Figure 3a: UMAPs of CD8 T cells with TIP scores.
  TIP_Fig2a_umaps(r1 = rCD8.smi)
  
  #2 Regenerate Figure 2b: Dot-plot showing CD8 T cell TIP genes association with infiltration in different immune subsets.
  TIP_Fig2b_dotplot()
  
  #3 Regenerate Figure 2c: Spatial maps showing the CD8 TIP scores in Xenium data
  TIP_Fig2c_Xenium(r = r.xenium,r1 = rTNK.xenium)
  
  #4 Plot CD8 T cell states
  TIP_Fig2d_Xenium(r = rCD8.xenium)
  
  #5 Regenerate Figure 2e: DAF heatmap
  HGSC_Fig2e_daf_heatmap(r = r.smi, daf = daf)
  
  #6 Regenerate PNG Images for Figure 2f: DAF in situ
  HGSC_Fig2f_daf_in_situ(r = r.smi, daf = daf)
  return()
}

#' Figure 2a. CD8 T cell umaps
#' Generates UMAP plots for CD8 T cells with coloration corresponding to 
#' cell environment and TIP expression 
#' @param r1 Processed Discovery data specific to CD8 T cells (optional).
#' @return  None. This function prints plots to disk. 
#' @export
TIP_Fig2a_umaps<-function(r1){
  if(!overwrite.flag&file.exists(get.file("Figures/Fig2a.pdf"))){return()}
  if(missing(r1)){r1<-HGSC_SMI.CD8.process()}
  l1<-umap.ggplot(r1$umap,labels = r1$plot[,c("Malignant","Fibroblast","TIP")],main = "CD8 T cells (all genes)")
  l2<-umap.ggplot(r1$cts$umap,labels = cbind.data.frame(r1$plot[,c("Malignant","TIP")],Clusters = r1$cts$clusters),
                  main = "CD8 T cells (T cell genes)")

  pdf(get.file("Figures/Fig2a.pdf"))
  call.multiplot(c(l1,l2,"NULL"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)

  return()
}

# Figure 2b. Tumor infiltration genes in different immune cell subtypes.
#' Generates a dot-plot showing the association of CD8 T cell TIP genes with tumor infiltration in different immune subsets.
#' @param R Results from TIP analysis (optional).
#' @return None. Generates Fig2b.pdf
TIP_Fig2b_dotplot<-function(){
  if(!overwrite.flag&file.exists(get.file("Figures/Fig2b.pdf"))){return()}
  X0<-readRDS(get.file("/Data/SourceData_fig2b.rds"))
  p<-call.dotPlot(X0,cex = 8)
  pdf(get.file("Figures/Fig2b.pdf"))
  print(p)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

#' Figure 2c. Validation 1 dataset spatial maps with CD8 TIP
#' Generates spatial maps showing the CD8 TIP scores in Validation 1 dataset
#' @param r Processed Validation 1 dataset
#' @param r1 Processed Validation 1 dataset specific to CD8 T and NK cells.
#' @return None. Generates Fig2c.pdf
TIP_Fig2c_Xenium<-function(r,r1){
  if(!overwrite.flag&file.exists(get.file("Figures/Fig2c.pdf"))){return()}
  fun_color_range <- colorRampPalette(c("blue","grey","red"))
  f<-function(x){
    n1<-min(100,floor(length(x)/3))
    v<-call.discretize(x,n1)
    return(fun_color_range(n1)[v])}

  X<-apply(cap.mat(r1$scores[,1:2],cap = 0.01,MARGIN = 2),2,f)
  rownames(X)<-rownames(r1$scores)
  colnames(X)<-colnames(r1$scores[,1:2])

  # Expand to all cells, including non-malignant ones
  idx<-match(r$cells,r1$cells)
  X1<-X[idx,]
  X1[is.na(idx)&r$cell.types=="Malignant",]<-"grey" # T cells are marked in black
  # Non-malignant cells other than CD8 T cells and NK cells are not shown
  X1[is.na(idx)&r$cell.types!="Malignant",]<-"white"
  rownames(X1)<-r$cells
  r$cd8.b<-is.element(r$cells,r1$cells)
  r$b.plot<-r$cd8.b|r$cell.types=="Malignant"

  samples<-paste0("XEN_T10_",c("H4","C4","F2","F3","E5","G3","H1","I2"))

  pdf(get.file("Figures/Fig2c.pdf"))
  par(mfrow=c(2,2),oma = c(0, 0, 0, 1),xpd = T)
  for(x in samples){
    b1<-r$samples==x&r$b.plot;
    if(x=="XEN_T10_F3"){b1<-b1&r$coor[,1]>4000}
    if(x=="XEN_T10_G3"){b1<-b1&r$coor[,"y"]>3000}
    call.plot.plus(r$coor[b1,],labels = r$scores[b1,3],my.col = X1[b1,1],
                 b.top = r$cd8.b[b1],cex = ifelse(r$cd8.b[b1],0.4,0.2),
                 legend.flag = F,main = x,cex.axis = 1e-10)
  }
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

#' Figure 2c. Validation 1 dataset UMAPs
#' Generates Validation 1 dataset UMAPs showing CD8 cell states and TIP
#' @param r1 Processed Validation 1 dataset specific to CD8 T cells.
#' @return None. Generates Fig2d.pdf 
TIP_Fig2d_Xenium<-function(r1){
  if(!overwrite.flag&file.exists(get.file("Figures/Fig2d.pdf"))){return()}
  l1<-list(umap.ggplot(r1$umap,labels = r1$cell.types,main = "Validation 1 Dataset, CD8 T cells",labels.name = "Cell states"),
           umap.ggplot(r1$umap,labels = r1$tmeC[,"Malignant"],main = "Validation 1 Dataset, CD8 T cells",labels.name = "Malignant frq"))
  Export <- gridExtra::marrangeGrob(l1, nrow = 1, ncol = 1)
  ggsave(get.file("Figures/Fig2d.pdf"), width = 5, height = 4, units = "in",plot = Export)
  return()
}

#' Figure 2e. Heatmap of Desmoplastic Fibroblast Program
#' Generates a heatmap for the Desmoplastic Fibroblast (DAF) program.
#' @param r Processed Discovery dataset
#' @param daf Desmoplastic Fibroblast data.
#' @return None. This function is used for its side effect of creating a heatmap.
HGSC_Fig2e_daf_heatmap <- function(r, daf){
  if(!overwrite.flag&file.exists(get.file("Figures/Fig2e.pdf"))){return()}
  morph<-daf$morph
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
  # pdf(get.file("Figures/Fig2e.pdf"),
  #     width = 10, height = 5)
  plt <- pheatmap(mat[oe_mean$samples,gene],
                  cluster_rows = F,
                  cluster_cols = F,
                  annotation_row = ann,
                  show_rownames = F,
                  show_colnames = T,
                  annotation_colors = ann_colors,
                  border_color = NA,
                  annotation_col = ann_col, filename = get.file("Figures/Fig2e.pdf"))
  # dev.off()
}

#' Figure 2f. Desmoplastic Fibroblast Program in situ
#' Generates in situ visualizations for the Desmoplastic Fibroblast (DAF) program.
#' @param r Processed Discovery dataset.
#' @param daf Desmoplastic Fibroblast results.
#' @param cell_2_rgb Color mapping for cells.
#' @param s Sample name (default is "SMI_T11_F019").
#' @return None. This function is used for its side effect of creating in situ visualizations.
HGSC_Fig2f_daf_in_situ <- function(r, daf, s = "SMI_T11_F019"){
  if(!overwrite.flag&file.exists(get.file("Figures/Fig2f.pdf"))){return()}
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
                               outfile = get.file("Figures/Fig2f-2.png"),
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
                               outfile = get.file("Figures/Fig2f-1.png"))

}


