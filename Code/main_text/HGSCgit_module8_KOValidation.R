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
  source = readRDS(get.file("Data/SourceData_Fig8.rds"))
  # Plot KO Validation in NK cells
  HGSC_Fig8a(source[[1]])
  # Plot KO Validation in CD8 T cells
  HGSC_Fig8b(source[[2]]) 
  # Plot PTPN1i Validation
  HGSC_Fig8c(source[[3]])
}

#' Figure 8a. KO Validation in NK cells
#'
#' Wrapper around plot_caspase_final function to visualize caspase validation
#' experiment with NK cells. 
#' 
#' @return this function returns nothing, but prints a plot to device. 
HGSC_Fig8a <- function(df1){
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
HGSC_Fig8b <- function(df2){
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
HGSC_Fig8c <- function(df){
  plt4 <- filter(df, cellline == "Tyk-nu")
  plt5 <- plt4 %>% 
    mutate(D = 1 - f_adj/f_ctrl) %>% 
    group_by(Condition, Conc) %>% 
    dplyr::summarize(c = mean(D), 
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
  
  
  plt4 <- filter(df, cellline == "OVCAR3")
  plt5 <- plt4 %>% 
    mutate(D = 1 - f_adj/f_ctrl) %>% 
    group_by(Condition, Conc) %>% 
    dplyr::summarize(c = mean(D), 
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
  
  # do calculations   
  plt = df %>%
    mutate(cn = i - c) %>% 
    group_by(Elapsed, Gene, Condition) %>%
    dplyr::summarize(c = mean(cn),
                     se = sd(cn)/sqrt(length(cn)),
                     n = length(cn))  %>%
    mutate(cond = paste0(Gene, "_", Condition))
  plt$cond <- factor(plt$cond, unique(plt$cond))
  plt$data <- "Single Gene KOs"
  
  # make plot
  p1 = ggplot(plt, aes(x = Elapsed, y = c, col = cond, shape = Condition)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = c- se, ymax = c+ se), width = 0.3) +
    theme_pubclean() +
    scale_color_manual(values = colors) +
    ylab("Delta Caspase Intensity") +
    xlab("Time (Hours)") +
    ggtitle(title, subtitle = subtitle)+
    xlim(0, 16.5)
  
  return(p1)
}
