#### Results Section 8 ###
# Figure 8. Knockout and Drug Treatment Validation 
# 8a. KO Validation in NK cells
# 8b. KO Validation in CD8 T cells
# 8c. PTPN1i Validation

#' Figure 8 Single Hit Validations 
#'
#' This function calls code to reproduce figures 8a-c
#'
#' @return this function returns nothing, but writes figures in .pdf format 
#' in the Figures/ folder. 
HGSC_Figure8_KOValidation <- function(){
  # Plot KO Validation in NK cells
  HGSC_Fig8a()
  # Plot KO Validation in CD8 T cells
  HGSC_Fig8b() 
  # Plot PTPN1i Validation
  HGSC_Fig8c()
}

#' Figure 8a. KO Validation in NK cells
#'
#' Wrapper around plot_caspase_final function to visualize caspase validation
#' experiment with NK cells. 
#' 
#' @return this function returns nothing, but prints a plot to device. 
HGSC_Fig8a <- function(){
  df1 <- read.csv(get.file("Data/NK_KO_Int.csv"))
  a = plot_caspase_final(df1)
  pdf(get.file("Figures/Fig8a.pdf"), width = 8, height = 6)
  print(a)
  dev.off()
}

#' Figure 8b. KO Validation in CD8 T cells
#'
#' Wrapper around plot_caspase_final function to visualize caspase validation
#' experiment with CD8 T cells
#' 
#' @return this function returns nothing, but prints a plot to device. 
HGSC_Fig8b <- function(){
  df2 <- read.csv(get.file("Data/T_KO_Int.csv"))
  b = plot_caspase_final(
    df2,title = "TYK-nu in CD8 T Co-culture vs. Monoculture", 
    subtitle = "Single Gene KOs", 
    colors = colors <- c("dodgerblue4", "#bcc6e8", 
                         "black", "#dbdbdb",
                         'darkred', "#f7cecb"),
    noise = T)
  pdf(get.file("Figures/Fig8b.pdf"), width = 8, height = 6)
  print(b)
  dev.off()
}

#' Figure 8c. KO Validation in CD8 T cells
#'
#' Uses data from cytotoxicity assasys with different cell lines to make 
#' unified plot of drug-specific cytotoxicity attributd to PTPN1i. 
#' 
#' @return this function returns nothing, but prints a plot to device. 
HGSC_Fig8c <- function(){
  # read in data from multiple experiments. 
  plate_ids <- as.matrix(read.csv(get.file("Data/0208_ptpni_tyknu.csv"), 
                                  header = F))
  plate_f <- as.matrix(read.csv(get.file("Data/0208_ptpni_tyknu_val.csv"), 
                                header = F))
  df <- data.frame(condition = c(plate_ids), f = c(plate_f))
  df$batch = "0208"
  plate_ids <- as.matrix(read.csv(get.file("Data/0313_ptpni_tyknu.csv"), 
                                  header = F))
  plate_f <- as.matrix(read.csv(get.file("Data/0313_ptpni_tyknu_val.csv"), 
                                header = F))
  df1 <- data.frame(condition = c(plate_ids), f = c(plate_f))
  df1$batch = "0313"
  
  # Correct for Media Baseline in Each Batch 
  bl<- mean(filter(df, condition == "Media")$f)
  df$fa = df$f - bl
  bl <- mean(filter(df1, condition == "Media")$f)
  df1$fa = df1$f - bl
  
  # Combine Both Datasets, remove Empty and Media
  df <- rbind(df, df1)
  df <- filter(df, condition != "Empty" & condition != "Media")
  
  # Parse Olivia's annotations
  split = data.frame(do.call("rbind", lapply(df$condition, function(cond){
    split = strsplit(cond, split = "Mono")[[1]]
    out = c(gsub("_", "", split[1]), "Monoculture", 
            gsub("uM", "", gsub("_", "", split[2])))
    if (length(split) == 1) {
      split = strsplit(cond, split = "Co-culture")[[1]]
      out = c(gsub("_", "", split[1]), "Co-culture", 
              gsub("uM", "", gsub("_", "", split[2])))
    }
    return(out)
  })))
  colnames(split) <- c("Treatment", "Condition", "Conc")
  plt <- cbind(df, split)
  
  # Filter out DMSO data we are NOT using 
  plt1 <- filter(plt, Treatment != "DMS0" & (Conc != 0 | Conc != 16))
  
  # Make 16uM DMSO controls
  dmso <- filter(plt1, Treatment == "DMSO", Conc == 16)
  controls = dmso %>% group_by(Condition, batch) %>% summarize(fa = mean(fa))
  
  ## use DMOS 16um as the 0 condition 
  plt3 <- filter(plt1, 
                 (Treatment == "PTPN1i" & Conc != "0") |
                   (Treatment == "DMSO" & Conc == "16"))
  b = grepl("DMSO", plt3$condition) & grepl("16uM", plt3$condition)
  plt3$Treatment[b] <- "PTPN1i"
  plt3$Conc[b] <- "0"
  plt4 <- merge(plt3, controls, by = c("Condition", "batch"), all.x = T)
  plt5 <- plt4 %>% 
    mutate(D = 1 - fa.x/fa.y) %>% 
    group_by(Condition, Conc) %>% 
    summarize(c = mean(D), 
              se = sd(D)/sqrt(length(D))) 
  plt5$Conc <- as.numeric(plt5$Conc) 
  
  # make TYK-nu plot
  cellline = "TYK-nu"
  a = ggplot(plt5, aes(x = Conc, y = c*100, col = Condition)) + 
    geom_point()+
    geom_line(aes(linetype = Condition)) +
    geom_errorbar(aes(ymin = c*100 - se*100, 
                      ymax = c*100+ se*100), width = 0.4) +
    theme_classic() +
    ggtitle(paste0(cellline, " in NK Co-culture vs. Monoculture"), 
            subtitle = "All Data") +
    xlab("[ABBV-CLS-484] (uM)") +
    ylab("ABBV-CLS-484 Specific Cytotoxicity (%)") +
    ylim(-10, 100) +
    scale_color_manual(values = c("darkred", "dodgerblue4")) +
    geom_hline(yintercept = 0, color = "grey")
  
  plate_ids <- as.matrix(
    read.csv(
      get.file("Data/0218_ptpni_ovcar.csv"), header = F))
  plate_f <- as.matrix(
    read.csv(
      get.file("Data/0218_ptpni_ovcar_val.csv"), header = F))
  df <- data.frame(condition = c(plate_ids), f = c(plate_f))
  df$batch = "0218"
  plate_ids <- as.matrix(
    read.csv(
      get.file("Data/0306_ptpni_ovcar.csv"), header = F))
  plate_f <- as.matrix(
    read.csv(
      get.file("Data/0306_ptpni_ovcar_val.csv"), header = F))
  df1 <- data.frame(condition = c(plate_ids), f = c(plate_f))
  df1$batch = "0306"
  
  # Correct for Media Baseline in Each Batch 
  bl<- mean(filter(df, condition == "Media Only")$f)
  df$fa = df$f - bl
  bl <- mean(filter(df1, condition == "Media Only")$f)
  df1$fa = df1$f - bl
  
  # Combine Both Datasets, remove Empty and Media
  df <- rbind(df, df1)
  df <- filter(df, condition != "Empty" & condition != "Media Only") %>% 
    mutate(condition = gsub("Co-Culture", "Co-culture", condition))
  
  # Parse Olivia's annotations
  split = data.frame(do.call("rbind", lapply(df$condition, function(cond){
    split = strsplit(cond, split = "Mono")[[1]]
    out = c(gsub("_", "", split[1]), 
            "Monoculture", 
            gsub("uM", "", gsub("_", "", split[2])))
    if (length(split) == 2 & grepl("NoDrug", split[2])) {
      out = c("", "Monoculture", "0")
    }
    if (length(split) == 1) {
      split = strsplit(cond, split = "Co-culture", )[[1]]
      out = c(gsub("_", "", split[1]), 
              "Co-culture", 
              gsub("uM", "",
                   gsub("_", "", split[2])))
    }
    return(out)
  })))
  colnames(split) <- c("Treatment", "Condition", "Conc")
  plt <- cbind(df, split)
  
  # Filter out DMSO data we are NOT using 
  plt1 <- filter(plt, Treatment != "DMS0" & (Conc != 0 | Conc != 16))
  
  # Make 16uM DMSO controls
  dmso <- filter(plt1, Treatment == "DMSO", Conc == 16)
  controls = dmso %>% group_by(Condition, batch) %>% summarize(fa = mean(fa))
  
  ## use DMOS 16um as the 0 condition 
  plt3 <- filter(plt1, (Treatment == "PTPN1i" & Conc != "0") | 
                   (Treatment == "DMSO" & Conc == "16"))
  b = grepl("DMSO", plt3$condition) & grepl("16uM", plt3$condition)
  plt3$Treatment[b] <- "PTPN1i"
  plt3$Conc[b] <- "0"
  plt4 <- merge(plt3, controls, by = c("Condition", "batch"), all.x = T)
  plt5 <- plt4 %>% 
    mutate(D = 1 - (fa.x/fa.y)) %>% 
    group_by(Condition, Conc) %>% 
    summarize(c = mean(D), 
              se = sd(D)/sqrt(length(D))) 
  plt5$Conc <- as.numeric(plt5$Conc) 
  
  # plot OVCAR3
  cellline = "OVCAR3"
  b = ggplot(plt5, aes(x = Conc, y = c*100, col = Condition)) + 
    geom_point() + 
    geom_line(aes(linetype = Condition)) +
    geom_errorbar(aes(ymin = c*100- se*100, 
                      ymax = c*100+ se*100), width = 0.4) +
    theme_classic() +
    ggtitle(paste0(cellline, " in NK Co-culture vs. Monoculture"), 
            subtitle = "All Data") +
    xlab("[ABBV-CLS-484] (uM)") +
    ylab("ABBV-CLS-484 Specific Cytotoxicity (%)") +
    ylim(-10, 100) +
    scale_color_manual(values = c("darkred", "dodgerblue4")) +
    geom_hline(yintercept = 0, color = "grey") 
  
  pdf(get.file("Figures/Fig8c.pdf"), width = 8, height = 4)
  print(a+b)
  dev.off()
}

#' Plotting function that is used for Figure 8a-b 
#'
#' Uses data from Incucyte readout to visualize cytotoxicity as function 
#' of single genetic perturbations vs. controls. 
#' 
#' @return this function returns nothing, but prints a plot to device. 
plot_caspase_final <- function(
    df, 
    title = "TYK-nu in NK Co-culture vs. Monoculture", 
    subtitle = "Single Gene KOs", 
    colors = colors <- c("dodgerblue4", "#bcc6e8", 
                         "#f5cc00","#fff9d9", 
                         "black", "#dbdbdb", 
                         'darkred', "#f7cecb"),
    noise = F){
  df <- df %>% gather(key = condition, value = i, -Elapsed)
  df$condition <- unlist(
    lapply(
      df$condition, function(x) strsplit(x, split = "\\.")[[1]][1]))
  df = df %>% separate(condition, into = c("Gene", "Condition"), sep = "_")
  
  # filter controls + add noise to intensities
  out <- filter(df, Condition != "Ctrl" & Condition != "CTRL")
  if (noise) {out = mutate(out, i = i/1000 + 1)} else {
    out = mutate(out, i = i/1000)}
  
  # get the monoculture controls
  plt_monos <- out %>% filter(Condition != "CC") %>%
    group_by(Elapsed, Gene, Condition) %>% 
    dplyr::summarize(c = mean(i), 
                     se = sd(i)/sqrt(length(i)))
  
  # do calculations
  plt = out %>%
    merge(plt_monos %>% select(-Condition), 
          by = c("Elapsed", "Gene"), all.x=T) %>% 
    mutate(cn = i - c) %>% 
    group_by(Elapsed, Gene, Condition) %>%
    dplyr::summarize(c = mean(cn),
                     se = sd(cn)/sqrt(length(cn)))  %>%
    mutate(cond = paste0(Gene, "_", Condition))
  plt$cond <- factor(plt$cond, unique(plt$cond))
  plt$data <- "Single Gene KOs"
  
  # make plot
  p1 = ggplot(plt, aes(x = Elapsed, y = c, col = cond, shape = Condition)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = c- se, ymax = c+ se), width = 0.3) +
    # facet_wrap(~name) +
    theme_pubclean() +
    scale_color_manual(values = colors) +
    ylab("Delta Caspase Intensity") +
    xlab("Time (Hours)") +
    ggtitle(title, subtitle = subtitle)+
    xlim(0, 16.5)
  
  return(p1)
}
