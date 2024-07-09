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
#' This function calls code to reproduce main text Figures 1b-h. 
#' @return this function returns nothing, but writes figures in .pdf format 
#' in the Figures/ folder. 
HGSC_Figure1_SpatiomolecularMapping <- function(cell_2_rgb, r, q){
  #1 Regenerate input files for CoMut plot in python
  master <- HGSC_Fig1b_make_comut()
  #2 Regenerate Figure 1c: Cell Type UMAPs
  HGSC_Fig1c_celltype_umaps(cell_2_rgb)
  #3 Regenerate Figure 1d: Cell Types in situ
  HGSC_Fig1d_celltypes_insitu(cell_2_rgb, r, q)
  #4 Regenerate Figure 1e: Coembedding across ST and scRNA-seq datasets
  HGSC_Fig1e_coembedding()
  #5 Regenerate Figure 1f: Cell Type Compositions
  HGSC_Fig1f_celltype_composition(r=r, q=q)
  #6 Regenerate Figure 1g
  HGSC_Fig1g_hg(r=r)
  #7 Regenerate Figure 1h: Cell Co-localization Quotient Analysis
  HGSC_Fig1h_clq(r=r)
  return()
}

#' Figure 1b. CoMut Files
#' This function generates the .tsv files that are inputed into python "CoMut" 
#' code for generating a summary figure of multi-modal data for this HGSC cohort
#' @return a data frame of simplified clinical meta data. 
#' Figure 1b. CoMut Files
#'
#' This function generates the .tsv files that are inputed into python "CoMut" 
#' code for generating a summary figure of multi-modal data for this HGSC cohort
#'
#' @return null, .tsv files can be found in Results/CoMut/
HGSC_Fig1b_make_comut <- function(){
  outdir <- paste0(get.file("Results/CoMut/"))
  if(!dir.exists(outdir)){dir.create(outdir, recursive = T)}
  sd <- readRDS(get.file("Data/SourceData_Fig1b.rds"))
  
  # write treatment data out. 
  df <- sd[[1]]
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
  
  ## write mutation data out 
  df_mut <- sd[[2]]
  write.table(df_mut, paste0(outdir, "CoMut_Muts.tsv"),
              sep = "\t", row.names = F, quote = F)
  
  # load tempus data to show which pts have tumor mutational burden data
  df_tmp <- sd[[3]]
  write.table(df_tmp, paste0(outdir,"CoMut_Tempus.tsv"),
              sep = "\t", row.names = F, quote = F)
  
  # spatial transcriptomics data summarization
  df_st <- sd[[4]]
  write.table(df_st, paste0(outdir,"CoMut_Spatial.tsv"),
              sep = "\t", row.names = F, quote = F)
  
  # write patient order! 
  samples <- (df %>% arrange(Stage, Age))$Patient
  write.table(data.frame(samples = samples),
              file = paste0(outdir,"CoMut_patorder.tsv"),
              sep = "\t", row.names = F, quote = F)
  
  return()
}

#' Figure 1c. Cell Type UMAPs 
#' This function loads Seurat Objects that store the UMAP embeddings for high-
#' confidence cells for each of the Discovery, Validation 1, and Validation 2
#' datasets. The function loads list objects with the UMAP embedding info for 
#' the Test 1 and Test 2 datasets. Finally this function plots the UMAP
#' embeddings for each dataset where the cells are colored by their assigned 
#' cell types. 
#' @return this function returns nothing, but writes figures in .pdf format 
#' in the Figures/ folder.
HGSC_Fig1c_celltype_umaps <- function(cell_2_rgb){
  df = readRDS(get.file("Results/SourceData_Fig1c.rds"))
  
  plts <- lapply(unique(df$dataset), function(dat){
    tmp <- filter(df, dataset == dat) 
    set.seed(12345)
    tmp <- tmp[sample(row.names(tmp), length(row.names(tmp))),]
    
    
    p = ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, col = cell.types)) + 
      geom_point(size = 0.1) + 
      theme_classic() +
      scale_color_manual(values = unlist(
        lapply(
          cell_2_rgb, rgb2hex
        )
      ), name = "Cell Types") +
      ggtitle(label = dat) +
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.title = element_text(size =12),
            legend.position = "none", 
            plot.title = element_text(size = 20, 
                                      face = "bold", 
                                      hjust = 0.5)) 
    return(p)
  })
  
  comb = plts[[1]] | (plts[[2]] + plts[[3]])/(plts[[4]] + plts[[5]])
  png(get.file("Figures/Fig1c.png"),
      height = 6, width = 16, units = "in", res = 500)
  print(comb) 
  dev.off()
}

#' Figure 1d. Cell Types plotted in situ. 
#' For smaller FOVs of tissue, we use segmentation maps to visualize the cell 
#' type calls in situ. For the Whole Tissue samples, we use 
#' @return this function returns nothing, but writes figures in .pdf or .png 
#' format in the Figures/ folder. 
HGSC_Fig1d_celltypes_insitu<- function(cell_2_rgb, r, q){
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
  s <- readRDS(get.file("/Data/ST_Validation2.rds"))
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
#' For smaller FOVs of tissue, we use segmentation maps to visualize the cell 
#' type calls in situ. For the Whole Tissue samples, we use 
#' @return this function returns nothing, but writes figures in .pdf or .png 
#' format in the Figures/ folder. 
HGSC_Fig1e_coembedding <- function(cell_2_rgb){
  coembedding <- readRDS(get.file("Results/SourceData_Fig1e.rds"))
  
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
#' For smaller FOVs of tissue, we use segmentation maps to visualize the cell 
#' type calls in situ. For the Whole Tissue samples, we use 
#' @return this function returns nothing, but writes figures in .pdf or .png 
#' format in the Figures/ folder. 
HGSC_Fig1f_celltype_composition <- function(ct = c("Malignant", "Fibroblast",
                                                   "Endothelial", "Monocyte",
                                                   "Mast.cell", "B.cell",
                                                   "TNK.cell", "Other"), r, q){
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
  
  # compute composition for all datasets
  files = Sys.glob(get.file("Data/scRNA*"))
  if(!file.exists(get.file("Results/SourceData_Fig1f.rds"))){
    
    # composition of scRNA-seq
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
    s <- readRDS(get.file("/Data/ST_Validation2.rds"))
    out <- data.frame(data.frame(s[c("cells", "samples", "patients",
                                     "sites", "cell.types.nonmal")])) %>%
      mutate(sites_binary = c("Adnexa", "Omentum")[(sites == "Omentum") + 1],
             cell.types = factor(cell.types.nonmal, levels = ct)) %>%
      select(-cell.types.nonmal) %>%
      mutate(platform = "ST",
             dataset = "MERFISH")
    df_mer <- compute_composition(s, out, ct)
    remove(s)
    
    # extract from test 1 
    t1 <- readRDS(get.file("Data/ST_Test1.rds"))
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
    t2.1<- readRDS(get.file("Data/ST_Test2_HGSC113.rds"))
    t2.2 <- readRDS(get.file("Data/ST_Test2_HGSC8.rds"))
    t2.3 <- readRDS(get.file("Data/ST_Test2_HGSC2.rds"))
    t2.4 <- readRDS(get.file("Data/ST_Test2_HGSC1.rds"))
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
  } else {
    sd = readRDS(get.file("Results/SourceData_Fig1f.rds"))
    df_sc <- sd[[1]];df_smi <- sd[[2]];df_xen <- sd[[3]];df_mer <- sd[[4]]
    full <- sd[[5]]
  }
  
  cell_2_rgb <- list()
  cell_2_rgb[["B.cell"]] = c(255, 0, 0)
  cell_2_rgb[["Fibroblast"]] = c(34, 91, 224)
  cell_2_rgb[["Endothelial"]] = c(166, 128, 71)
  cell_2_rgb[["Monocyte"]] = c(255, 153, 0)
  cell_2_rgb[["Mast.cell"]] = c(157, 3, 252)
  cell_2_rgb[["Malignant"]] = c(7, 224, 0)
  cell_2_rgb[["TNK.cell"]] = c(255, 247, 0)
  cell_2_rgb[["Other"]] = c(153, 153, 153)
  
  # combine scRNA-seq, Discovery and Validation datasets 
  full1 <- rbind(df_sc, df_smi, df_xen, df_mer)
  full1 <- full1 %>% mutate(dataset = factor(dataset,
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
  order <- (full1 %>% filter(
    cell.types == "Malignant") %>% arrange(dataset, proportion))$samples
  full1 <- full1 %>% mutate(samples = factor(samples, levels = order),
                          cell.types = factor(
                            cell.types,
                            levels = c("Malignant", "Fibroblast",
                                       "Endothelial","Monocyte",
                                       "Mast.cell", "TNK.cell",
                                       "B.cell", "Other")))
  
  # plot first plot to disk
  pdf(get.file("Figures/Fig1f.pdf"),
      width = 18, height = 3)
  p <- ggplot(full1, aes(x = samples, y = proportion, fill = cell.types)) +
    geom_bar(stat = "identity") +
    facet_grid(~dataset, scales = "free_x", space = "free", 
               labeller = as_labeller(labels)) +
    scale_fill_manual(values = lapply(cell_2_rgb, rgb2hex)) +
    theme(axis.text.x = element_blank(),
          strip.text = element_text(color = "lightgrey"))
  print(p)
  g <- ggplot(full1, aes(x = samples, y = dataset, fill = dataset)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = RColorBrewer::brewer.pal(9, "Set3")) +
    theme(legend.position = "top")
  as_ggplot(get_legend(g))
  dev.off()
  
  # make the patch plot for the test data 
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
#' For the Discovery dataset, we perform hypergeometric tests to evaluate if 
#' cell types co-occur in spatial frames more often than expected by random. 
#' @return this function returns nothing, but writes figures in .pdf format. 
HGSC_Fig1g_hg <- function(r){
  # divide samples into grids
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
#' For the Discovery dataset, we calculate co-localization quotient (CLQ, 
#' Methods) between fibroblasts and T/NK cells and malignant cells and T/NK 
#' cells to explore the relationship of infiltration in the two predominant 
#' spatial components of the tumor tissue. 
#' @return this function returns nothing, but writes figures in .pdf format.
HGSC_Fig1h_clq <- function(r=r){
  
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
  r$neighbors <- readRDS(get.file("Results/SourceData_Fig1h.rds"))
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
