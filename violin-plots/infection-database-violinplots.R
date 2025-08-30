# cargamos librerias necesarias

# install.packages("ggplot2")
library(ggplot2)
# install.packages("dplyr")
library(dplyr)
# install.packages("survminer")
library(survminer)
# install.packages("pROC")
library(pROC)
# install.packages("foreign")
library(foreign)
# install.packages("haven")
library(haven)
# install.packages("rlang")
library(rlang)
# install.packages("patchwork")
library(patchwork)

# cargamos BDD
bddanciano <- 
  read_sav("C:/Users/tposc/master-in-bioinformatics/phd/BDD_SPSS_SYSTEMSBIOCOVID_08112023v18 - paper anciano.sav")

# View(bddanciano)
str(bddanciano)
dim(bddanciano)

# analizamos variable que divide a pacientes por edad
bddanciano$Edad_grupos_mínhasta70_71hastamáx
bddanciano$Edad_grupos_mínhasta70_71hastamáx <- as.factor(bddanciano$Edad_grupos_mínhasta70_71hastamáx)


# Hay variables que requieren ser factorizadas:
library(dplyr)
# Factorizamos las variables que están compuestas por niveles.
bddanciano <- bddanciano %>% 
  mutate(across (c(IL2_pgmL, GMCSF_pgmL, CCL2_pgmL, GRANA_pgmL, IFNg_pgmL, IL4_pgmL, 
                   TNFa_pgmL, B7H1_pgmL, EGF_pgmL, Endothelin1_pgmL, IFNa_pgmL, 
                   IL10_pgmL, IL7_pgmL, IL15_pgmL, CXCL10_pgmL, IL6_pgmL, DimerD_pgmL, 
                   MPO_pgmL, RANTESCCL5_pgmL, PENTRAXIN3TSG14_pgmL, Lipocalin2NGAL_pgmL, 
                   N1_copies_mL_plasm_1304, anti_N_título, 
                   anti_S_título), as.numeric))

# variables totales: "IL2_pgmL", "CCL2_pgmL", "GRANA_pgmL", "IFNg_pgmL", "IL4_pgmL", 
# "TNFa_pgmL", "B7H1_pgmL", "Endothelin1_pgmL", "IFNa_pgmL", "IL10_pgmL", "IL7_pgmL", 
# "IL15_pgmL", "CXCL10_pgmL", "IL6_pgmL", "DimerD_pgmL", "RANTESCCL5_pgmL", 
# "PENTRAXIN3TSG14_pgmL", "Lipocalin2NGAL_pgmL", "GMCSF_pgmL", "N1_copies_mL_plasm_1304", 
# "anti_N_título", "anti_S_título"



## violin plots:
## https://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization?title=ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
## https://ggplot2.tidyverse.org/reference/theme.html
## https://r-charts.com/ggplot2/themes/

resumen_vars <- function(df, vars) {
  lapply(df[vars], summary)
}

resumen_biomarc <- resumen_vars(
  bddanciano,
  c("IL2_pgmL", "GMCSF_pgmL", "CCL2_pgmL", "GRANA_pgmL", "IFNg_pgmL", "IL4_pgmL", 
    "TNFa_pgmL", "B7H1_pgmL", "EGF_pgmL", "Endothelin1_pgmL", "IFNa_pgmL", 
    "IL10_pgmL", "IL7_pgmL", "IL15_pgmL", "CXCL10_pgmL", "IL6_pgmL", "DimerD_pgmL", 
    "MPO_pgmL", "RANTESCCL5_pgmL", "PENTRAXIN3TSG14_pgmL", "Lipocalin2NGAL_pgmL", 
    "N1_copies_mL_plasm_1304", "anti_N_título", "anti_S_título")
)

resumen_biomarc


# ver si los NAs están bien:
# sapply(bddanciano[c("IL2_pgmL", "GMCSF_pgmL")], function(x) sum(is.na(x)))

# ver si las columnas se convirtieron a numeric
# sapply(bddanciano[c("IL2_pgmL", "GMCSF_pgmL")], class)


# bddanciano %>% ggplot(aes(Edad_grupos_mínhasta70_71hastamáx,
                         # TNFa, 
                         # color = Edad_grupos_mínhasta70_71hastamáx,
                         # xlab = "Grupo", 
                         # ylab = "TNFa (pg/mL)")) + 
                         # ggtitle("TNFa levels by age group") +
                  # geom_violin(trim=TRUE)

#################################################################################
### VAMOS A CAMBIAR EL NOMBRE A TODAS LAS VARIABLES Y AL GRUPO PARA REPRESENTACIÓN FINAL.
################################################################################

# View(bddanciano)

bddanciano$TNFa <- bddanciano$TNFa_pgmL
bddanciano$IL2 <- bddanciano$IL2_pgmL
bddanciano$GMCSF <- bddanciano$GMCSF_pgmL
bddanciano$CCL2 <- bddanciano$CCL2_pgmL
bddanciano$GranzymeA <- bddanciano$GRANA_pgmL
bddanciano$IFNg <- bddanciano$IFNg_pgmL
bddanciano$IFNa <- bddanciano$IFNa_pgmL
bddanciano$IL4 <- bddanciano$IL4_pgmL
bddanciano$PDL1 <- bddanciano$B7H1_pgmL
bddanciano$Endothelin1 <- bddanciano$Endothelin1_pgmL
bddanciano$IL6 <- bddanciano$IL6_pgmL
bddanciano$IL7 <- bddanciano$IL7_pgmL
bddanciano$IL15 <- bddanciano$IL15_pgmL
bddanciano$IL10 <- bddanciano$IL10_pgmL
bddanciano$CXCL10 <- bddanciano$CXCL10_pgmL
bddanciano$DDimer <- bddanciano$DimerD_pgmL
bddanciano$CCL5 <- bddanciano$RANTESCCL5_pgmL
bddanciano$Pentraxin3 <- bddanciano$PENTRAXIN3TSG14_pgmL
bddanciano$Lipocalin2 <- bddanciano$Lipocalin2NGAL_pgmL
bddanciano$N1copies <- bddanciano$N1_copies_mL_plasm_1304
bddanciano$IgGantiN <- bddanciano$anti_N_título
bddanciano$IgGantiS <- bddanciano$anti_S_título


bddanciano$Edad_grupos_mínhasta70_71hastamáx <- factor(bddanciano$Edad_grupos_mínhasta70_71hastamáx,
                           levels=c("0","1"),
                           labels=c("≤ 70 years old","> 70 years old"))




#################################################################################
#################################################################################
#################################################################################
################# >>>>>>>>>>> BUENO <<<<<<<<<<<<<<<<<<
################################################################################
################################################################################

# vamos a guardar el plot como .png:
png("violinp_il2_ccl2_ifng_grana.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
# svg("violin_plots1.svg", width=8, height=6)

plots_var <- function(df, var) {
  v <- df[[as_name(enquo(var))]]
  upper <- quantile(v, 0.999, na.rm = TRUE)   
  
  ggplot(df, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = {{var}})) +
    
    # violín → borde fuerte
    geom_violin(
      aes(colour = Edad_grupos_mínhasta70_71hastamáx),
      fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
    ) +
    
    geom_jitter(
      aes(fill = Edad_grupos_mínhasta70_71hastamáx, color = Edad_grupos_mínhasta70_71hastamáx),
      shape = 21, position = position_jitter(0.2), alpha = .2, size = 2
    ) +
    
    # boxplot → borde en tono distinto
    geom_boxplot(
      aes(fill = Edad_grupos_mínhasta70_71hastamáx),
      width = .05, outlier.shape = NA, color = "black", alpha = .7
    ) +
    
    # escalas distintas para borde (violín) y relleno (boxplot)
    
    scale_fill_manual(values = c("> 70 years old" = "firebrick1", 
                                 "≤ 70 years old" = "goldenrod1")) + # color relleno boxplot
    scale_colour_manual(values = c("> 70 years old" = "firebrick4", 
                                   "≤ 70 years old" = "goldenrod4")) + # color perfil violinplot
    
    
    coord_cartesian(ylim = c(0, upper)) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", colour = "grey50")
    ) +
    labs(
      # title = paste(as_name(enquo(var)), "levels by age group"),
      x = "Group by age",
      y = paste(as_name(enquo(var)), "(pg/mL)")
    )
}


vars <- c("IL2", "CCL2", "GranzymeA", "IFNg")

plots_list <- lapply(vars, function(v) {
  plots_var(bddanciano, !!sym(v))
})

wrap_plots(plots_list, ncol = 2)  # 2 columnas (ajusta a gusto)


## para meter más de 2x2 gráficos.
# library(gridExtra)
# do.call(grid.arrange, c(plots_list, ncol = 2))

# Cerramos la imagen .svg generada para el plot
dev.off()


#################################################################################
#################################################################################
################# >>>>>>>>>>> BUENO <<<<<<<<<<<<<<<<<< OTRAS VARIABLES

# vamos a guardar el plot como .png:
png("violinp_il4_tnfa_b7h1_endoth.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
# svg("violin_plots1.svg", width=8, height=6)

# variables totales: "IL2_pgmL", "CCL2_pgmL", "GRANA_pgmL", "IFNg_pgmL", "IL4_pgmL", 
# "TNFa_pgmL", "B7H1_pgmL", "Endothelin1_pgmL", "IFNa_pgmL", "IL10_pgmL", "IL7_pgmL", 
# "IL15_pgmL", "CXCL10_pgmL", "IL6_pgmL", "DimerD_pgmL", "RANTESCCL5_pgmL", 
# "PENTRAXIN3TSG14_pgmL", "Lipocalin2NGAL_pgmL", "GMCSF_pgmL", "N1_copies_mL_plasm_1304", 
# "anti_N_título", "anti_S_título"

vars2 <- c("IL4", "TNFa", "PDL1", "Endothelin1")

plots_list <- lapply(vars2, function(v) {
  plots_var(bddanciano, !!sym(v))
})

wrap_plots(plots_list, ncol = 2)  # 2 columnas (ajusta a gusto)


## para meter más de 2x2 gráficos.
# library(gridExtra)
# do.call(grid.arrange, c(plots_list, ncol = 2))

# Cerramos la imagen .svg generada para el plot
dev.off()


#################################################################################
#################################################################################
################# >>>>>>>>>>> BUENO <<<<<<<<<<<<<<<<<< OTRAS VARIABLES

# vamos a guardar el plot como .png:
png("violinp_cxcl10_il6_gmcsf_rantes.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
# svg("violin_plots1.svg", width=8, height=6)

# variables totales: "IL2_pgmL", "CCL2_pgmL", "GRANA_pgmL", "IFNg_pgmL", "IL4_pgmL", 
# "TNFa_pgmL", "B7H1_pgmL", "Endothelin1_pgmL", "IFNa_pgmL", "IL10_pgmL", "IL7_pgmL", 
# "IL15_pgmL", "CXCL10_pgmL", "IL6_pgmL", "DimerD_pgmL", "RANTESCCL5_pgmL", 
# "PENTRAXIN3TSG14_pgmL", "Lipocalin2NGAL_pgmL", "GMCSF_pgmL", "N1_copies_mL_plasm_1304", 
# "anti_N_título", "anti_S_título"

vars3 <- c("CXCL10", "IL6", "GMCSF", "CCL5")

plots_list <- lapply(vars3, function(v) {
  plots_var(bddanciano, !!sym(v))
})

wrap_plots(plots_list, ncol = 2)  # 2 columnas (ajusta a gusto)


## para meter más de 2x2 gráficos.
# library(gridExtra)
# do.call(grid.arrange, c(plots_list, ncol = 2))

# Cerramos la imagen .svg generada para el plot
dev.off()


#################################################################################
#################################################################################
################# >>>>>>>>>>> BUENO <<<<<<<<<<<<<<<<<< OTRAS VARIABLES

# vamos a guardar el plot como .png:
png("violinp_ptx3_lipocalin2_ddimer.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
# svg("violin_plots1.svg", width=8, height=6)

# variables totales: "IL2_pgmL", "CCL2_pgmL", "GRANA_pgmL", "IFNg_pgmL", "IL4_pgmL", 
# "TNFa_pgmL", "B7H1_pgmL", "Endothelin1_pgmL", "IFNa_pgmL", "IL10_pgmL", "IL7_pgmL", 
# "IL15_pgmL", "CXCL10_pgmL", "IL6_pgmL", "DimerD_pgmL", "RANTESCCL5_pgmL", 
# "PENTRAXIN3TSG14_pgmL", "Lipocalin2NGAL_pgmL", "GMCSF_pgmL", "N1_copies_mL_plasm_1304", 
# "anti_N_título", "anti_S_título"

vars4 <- c("Pentraxin3", "Lipocalin2", "DDimer")

plots_list <- lapply(vars4, function(v) {
  plots_var(bddanciano, !!sym(v))
})

wrap_plots(plots_list, ncol = 3)  # 3 columnas (ajusta a gusto)


## para meter más de 2x2 gráficos.
# library(gridExtra)
# do.call(grid.arrange, c(plots_list, ncol = 2))

# Cerramos la imagen .svg generada para el plot
dev.off()


#################################################################################
#################################################################################
################# >>>>>>>>>>> BUENO <<<<<<<<<<<<<<<<<< OTRAS VARIABLES

# vamos a guardar el plot como .png:
png("violinp_n1copies_antin_antis3.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
# svg("violin_plots1.svg", width=8, height=6)

plots_var <- function(df, var) {
  v <- df[[as_name(enquo(var))]]
  upper <- quantile(v, 0.99, na.rm = TRUE)   
  
  ggplot(df, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = {{var}})) +
    
    # violín → borde fuerte
    geom_violin(
      aes(colour = Edad_grupos_mínhasta70_71hastamáx),
      fill = "white", alpha = .5, trim = TRUE, linewidth = 0.5
    ) +
    
    geom_jitter(
      aes(fill = Edad_grupos_mínhasta70_71hastamáx, color = Edad_grupos_mínhasta70_71hastamáx),
      shape = 21, position = position_jitter(0.2), alpha = .2, size = 2
    ) +
    
    # boxplot → borde en tono distinto
    geom_boxplot(
      aes(fill = Edad_grupos_mínhasta70_71hastamáx),
      width = .05, outlier.shape = NA, color = "black", alpha = .7
    ) +
    
    # escalas distintas para borde (violín) y relleno (boxplot)
    
    scale_fill_manual(values = c("> 70 years old" = "firebrick1", 
                                 "≤ 70 years old" = "goldenrod1")) + # color relleno boxplot
    scale_colour_manual(values = c("> 70 years old" = "firebrick4", 
                                   "≤ 70 years old" = "goldenrod4")) + # color perfil violinplot
    
    
    coord_cartesian(ylim = c(0, upper)) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", colour = "grey50")
    ) +
    labs(
      # title = paste(as_name(enquo(var)), "levels by age group"),
      x = "Group by age",
      y = paste(as_name(enquo(var)), "units")
    )
}

# variables totales: "IL2_pgmL", "CCL2_pgmL", "GRANA_pgmL", "IFNg_pgmL", "IL4_pgmL", 
# "TNFa_pgmL", "B7H1_pgmL", "Endothelin1_pgmL", "IFNa_pgmL", "IL10_pgmL", "IL7_pgmL", 
# "IL15_pgmL", "CXCL10_pgmL", "IL6_pgmL", "DimerD_pgmL", "RANTESCCL5_pgmL", 
# "PENTRAXIN3TSG14_pgmL", "Lipocalin2NGAL_pgmL", "GMCSF_pgmL", "N1_copies_mL_plasm_1304", 
# "anti_N_título", "anti_S_título"

vars5 <- c("N1copies", "IgGantiN", "IgGantiS")

plots_list <- lapply(vars5, function(v) {
  plots_var(bddanciano, !!sym(v))
})

wrap_plots(plots_list, ncol = 3)  # 3 columnas (ajusta a gusto)


## para meter más de 2x2 gráficos.
# library(gridExtra)
# do.call(grid.arrange, c(plots_list, ncol = 2))

# Cerramos la imagen .svg generada para el plot
dev.off()



#################################################################################
#################################################################################
################# >>>>>>>>>>> BUENO <<<<<<<<<<<<<<<<<< MEJORA PARA VER VARIABLES

# vamos a guardar el plot como .png:
png("violinp_ddimer_n1copies_antis.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
# svg("violin_plots1.svg", width=8, height=6)

plots_var <- function(df, var) {
  v <- df[[as_name(enquo(var))]]
  upper <- quantile(v, 0.90, na.rm = TRUE)   
  
  ggplot(df, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = {{var}})) +
    
    # violín → borde fuerte
    # geom_violin(
    # aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    # fill = "white", alpha = .5, trim = TRUE, linewidth = 0.5
    # ) +
    
    geom_jitter(
      aes(fill = Edad_grupos_mínhasta70_71hastamáx, color = Edad_grupos_mínhasta70_71hastamáx),
      shape = 21, position = position_jitter(0.2), alpha = .2, size = 2
    ) +
    
    # boxplot → borde en tono distinto
    geom_boxplot(
      aes(fill = Edad_grupos_mínhasta70_71hastamáx),
      width = .05, outlier.shape = NA, color = "black", alpha = .7
    ) +
    
    # escalas distintas para borde (violín) y relleno (boxplot)
    
    scale_fill_manual(values = c("> 70 years old" = "firebrick1", 
                                 "≤ 70 years old" = "goldenrod1")) + # color relleno boxplot
    scale_colour_manual(values = c("> 70 years old" = "firebrick4", 
                                   "≤ 70 years old" = "goldenrod4")) + # color perfil violinplot
    
    
    coord_cartesian(ylim = c(0, upper)) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", colour = "grey50")
    ) +
    labs(
      # title = paste(as_name(enquo(var)), "levels by age group"),
      x = "Group by age",
      y = paste(as_name(enquo(var)), "(pg/mL)")
    )
}

# variables totales: "IL2_pgmL", "CCL2_pgmL", "GRANA_pgmL", "IFNg_pgmL", "IL4_pgmL", 
# "TNFa_pgmL", "B7H1_pgmL", "Endothelin1_pgmL", "IFNa_pgmL", "IL10_pgmL", "IL7_pgmL", 
# "IL15_pgmL", "CXCL10_pgmL", "IL6_pgmL", "DimerD_pgmL", "RANTESCCL5_pgmL", 
# "PENTRAXIN3TSG14_pgmL", "Lipocalin2NGAL_pgmL", "GMCSF_pgmL", "N1_copies_mL_plasm_1304", 
# "anti_N_título", "anti_S_título"

vars5 <- c("DDimer", "N1copies", "IgGantiS")

plots_list <- lapply(vars5, function(v) {
  plots_var(bddanciano, !!sym(v))
})

wrap_plots(plots_list, ncol = 3)  # 3 columnas (ajusta a gusto)


## para meter más de 2x2 gráficos.
# library(gridExtra)
# do.call(grid.arrange, c(plots_list, ncol = 2))

# Cerramos la imagen .svg generada para el plot
dev.off()