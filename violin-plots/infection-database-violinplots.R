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
bddanciano$Edad_grupos_mínhasta70_71hastamáx <- factor(bddanciano$Edad_grupos_mínhasta70_71hastamáx,
                                                       levels=c("0","1"),
                                                       labels=c("≤ 70 años","> 70 años"))


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

################################################################################
################################################################################
######### GRÁFICOS SEPARADOS PARA PONER EJE Y >>> BIEN <<< #####################

################################################################################
################################################################################
######################## ORDENADOS POR FUNCIÓN #################################
################################################################################
################################################################################


################################################################################
################################################################################
######### GRÁFICOS SEPARADOS PARA PONER EJE Y >>> BIEN <<< #####################
#### lipocalin2 - ptx3 - il6 - tnfa

# vamos a guardar el plot como .png:
png("violinp_lipocalin2_ptx3_il6_tnfa3.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
# svg("violin_plots1.svg", width=8, height=6)

# variables totales: "IL2_pgmL", "CCL2_pgmL", "GRANA_pgmL", "IFNg_pgmL", "IL4_pgmL", 
# "TNFa_pgmL", "B7H1_pgmL", "Endothelin1_pgmL", "IFNa_pgmL", "IL10_pgmL", "IL7_pgmL", 
# "IL15_pgmL", "CXCL10_pgmL", "IL6_pgmL", "DimerD_pgmL", "RANTESCCL5_pgmL", 
# "PENTRAXIN3TSG14_pgmL", "Lipocalin2NGAL_pgmL", "GMCSF_pgmL", "N1_copies_mL_plasm_1304", 
# "anti_N_título", "anti_S_título"


## --------- Plot 1: lipocalin2 ---------
upper_lipocalin2 <- max(bddanciano$Lipocalin2NGAL_pgmL, na.rm = TRUE) * 1.1 # na.rm = TRUE ignora los NAs

p1 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = Lipocalin2NGAL_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5 
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .2, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7 # width = .1 si: 4 en línea
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_lipocalin2)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "Lipocalina 2 (pg/mL)")


## --------- Plot 2: ptx3 ---------
upper_ptx3 <- max(bddanciano$PENTRAXIN3TSG14_pgmL, na.rm = TRUE) * 1.1

p2 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = PENTRAXIN3TSG14_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .3, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_ptx3)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "Pentraxina 3 (pg/mL)")


## --------- Plot 3: TNFa ---------
upper_TNFa <- max(bddanciano$TNFa_pgmL, na.rm = TRUE) * 1.1

p3 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = TNFa_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .25, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_TNFa)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "TNFa (pg/mL)")


## --------- Plot 4: IL6 ---------
upper_IL6 <- quantile(bddanciano$IL6_pgmL, 0.99, na.rm = TRUE)

p4 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = IL6_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .3, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_IL6)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "IL-6 (pg/mL)")


## --------- Juntar ---------
# (p1 | p2 | p3 | p4)
# wrap_plots(list(p1, p2, p3), ncol = 3)
wrap_plots(list(p1, p2, p3, p4), ncol = 2, nrow = 2)


## para meter más de 2x2 gráficos.
# library(gridExtra)
# do.call(grid.arrange, c(plots_list, ncol = 2))

# Cerramos la imagen .svg generada para el plot
dev.off()


################################################################################
################################################################################
######### GRÁFICOS SEPARADOS PARA PONER EJE Y >>> BIEN <<< #####################
#### IFNa - ifng - il7 - grana

# vamos a guardar el plot como .png:
png("violinp_ifna_ifng_il7_grana.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
# svg("violin_plots1.svg", width=8, height=6)

# variables totales: "IL2_pgmL", "CCL2_pgmL", "GRANA_pgmL", "IFNg_pgmL", "IL4_pgmL", 
# "TNFa_pgmL", "B7H1_pgmL", "Endothelin1_pgmL", "IFNa_pgmL", "IL10_pgmL", "IL7_pgmL", 
# "IL15_pgmL", "CXCL10_pgmL", "IL6_pgmL", "DimerD_pgmL", "RANTESCCL5_pgmL", 
# "PENTRAXIN3TSG14_pgmL", "Lipocalin2NGAL_pgmL", "GMCSF_pgmL", "N1_copies_mL_plasm_1304", 
# "anti_N_título", "anti_S_título"

## --------- Plot 1: IFNa ---------
upper_IFNa <- max(bddanciano$IFNa_pgmL, na.rm = TRUE) * 1.1 # na.rm = TRUE ignora los NAs

p1 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = IFNa_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5 
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .2, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7 # width = .1 si: 4 en línea
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_IFNa)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "IFNa (pg/mL)")

## --------- Plot 2: IFNg_pgmL ---------
upper_IFNg <- max(bddanciano$IFNg_pgmL, na.rm = TRUE) * 1.1

p2 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = IFNg_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .3, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_IFNg)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "IFNg (pg/mL)")


## --------- Plot 3: IL7 ---------
upper_IL7 <- max(bddanciano$IL7_pgmL, na.rm = TRUE) * 1.1

p3 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = IL7_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .3, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_IL7)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "IL-7 (pg/mL)")

## --------- Plot 4: GRANA ---------
upper_GRANA <- max(bddanciano$GRANA_pgmL, na.rm = TRUE) * 1.1

p4 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = GRANA_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .3, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_GRANA)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "Granzima A (pg/mL)")


## --------- Juntar ---------
# (p1 | p2 | p3 | p4)
# wrap_plots(list(p1, p2, p3), ncol = 3)
wrap_plots(list(p1, p2, p3, p4), ncol = 2, nrow = 2)


## para meter más de 2x2 gráficos.
# library(gridExtra)
# do.call(grid.arrange, c(plots_list, ncol = 2))

# Cerramos la imagen .svg generada para el plot
dev.off()

################################################################################
################################################################################
######### GRÁFICOS SEPARADOS PARA PONER EJE Y >>> BIEN <<< #####################
#### ddimer - endoth - il15 - ccl2

# vamos a guardar el plot como .png:
png("violinp_ddimerx2_endoth_il15.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
# svg("violin_plots1.svg", width=8, height=6)

# variables totales: "IL2_pgmL", "CCL2_pgmL", "GRANA_pgmL", "IFNg_pgmL", "IL4_pgmL", 
# "TNFa_pgmL", "B7H1_pgmL", "Endothelin1_pgmL", "IFNa_pgmL", "IL10_pgmL", "IL7_pgmL", 
# "IL15_pgmL", "CXCL10_pgmL", "IL6_pgmL", "DimerD_pgmL", "RANTESCCL5_pgmL", 
# "PENTRAXIN3TSG14_pgmL", "Lipocalin2NGAL_pgmL", "GMCSF_pgmL", "N1_copies_mL_plasm_1304", 
# "anti_N_título", "anti_S_título"


## --------- Plot 1: ddimer ---------
upper_ddimer <- quantile(bddanciano$DimerD_pgmL, 0.99, na.rm = TRUE)

p1 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = DimerD_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .3, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_ddimer)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "Dímero D (pg/mL)")


## --------- Plot 2: DDimer ---------
upper_DDimer <- quantile(bddanciano$DimerD_pgmL, 0.90, na.rm = TRUE) # na.rm = TRUE ignora los NAs

p2 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = DimerD_pgmL)) +
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .2, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .05, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_DDimer)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "Dímero D (pg/mL)")


## --------- Plot 3: Endothelin1 ---------
upper_Endothelin1 <- max(bddanciano$Endothelin1_pgmL, na.rm = TRUE) * 1.1

p3 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = Endothelin1_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .3, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_Endothelin1)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "Endothelina-1 (pg/mL)")


## --------- Plot 4: IL15 ---------
upper_IL15 <- max(bddanciano$IL15_pgmL, na.rm = TRUE) * 1.1

p4 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = IL15_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .3, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_IL15)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "IL-15 (pg/mL)")


## --------- Juntar ---------
# (p1 | p2 | p3 | p4)
# wrap_plots(list(p1, p2, p3), ncol = 3)
wrap_plots(list(p1, p2, p3, p4), ncol = 2, nrow = 2)


## para meter más de 2x2 gráficos.
# library(gridExtra)
# do.call(grid.arrange, c(plots_list, ncol = 2))

# Cerramos la imagen .svg generada para el plot
dev.off()


################################################################################
################################################################################
######### GRÁFICOS SEPARADOS PARA PONER EJE Y >>> BIEN <<< #####################
#### ccl2 - cxcl10 - il10 - b7h1

# vamos a guardar el plot como .png:
png("violinp_ccl2_cxcl10_il10_b7h12.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
# svg("violin_plots1.svg", width=8, height=6)

# variables totales: "IL2_pgmL", "CCL2_pgmL", "GRANA_pgmL", "IFNg_pgmL", "IL4_pgmL", 
# "TNFa_pgmL", "B7H1_pgmL", "Endothelin1_pgmL", "IFNa_pgmL", "IL10_pgmL", "IL7_pgmL", 
# "IL15_pgmL", "CXCL10_pgmL", "IL6_pgmL", "DimerD_pgmL", "RANTESCCL5_pgmL", 
# "PENTRAXIN3TSG14_pgmL", "Lipocalin2NGAL_pgmL", "GMCSF_pgmL", "N1_copies_mL_plasm_1304", 
# "anti_N_título", "anti_S_título"


## --------- Plot 1: CCL2 ---------
upper_CCL2 <- max(bddanciano$CCL2_pgmL, na.rm = TRUE)*1.1

p1 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = CCL2_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .25, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_CCL2)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "CCL2 (pg/mL)")

## --------- Plot 2: CXCL10 ---------
upper_CXCL10 <- max(bddanciano$CXCL10_pgmL, na.rm = TRUE)*1.1

p2 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = CXCL10_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .25, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_CXCL10)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "CXCL10 (pg/mL)")


## --------- Plot 3: IL10 ---------
upper_IL10 <- max(bddanciano$IL10_pgmL, na.rm = TRUE) * 1.1

p3 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = IL10_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .25, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_IL10)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "IL-10 (pg/mL)")


## --------- Plot 4: B7H1 ---------
upper_B7H1 <- max(bddanciano$B7H1_pgmL, na.rm = TRUE) * 1.1

p4 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = B7H1_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .3, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, 850)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "PD-L1 (pg/mL)")

## --------- Juntar ---------
# (p1 | p2 | p3 | p4)
# wrap_plots(list(p1, p2, p3), ncol = 3)
wrap_plots(list(p1, p2, p3, p4), ncol = 2, nrow = 2)


## para meter más de 2x2 gráficos.
# library(gridExtra)
# do.call(grid.arrange, c(plots_list, ncol = 2))

# Cerramos la imagen .svg generada para el plot
dev.off()

################################################################################
################################################################################
######### GRÁFICOS SEPARADOS PARA PONER EJE Y >>> BIEN <<< #####################
#### il2 - il4 - gmcsf - ccl5

# vamos a guardar el plot como .png:
png("violinp_il2_il4_gmcsf_ccl54.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
# svg("violin_plots1.svg", width=8, height=6)

# variables totales: "IL2_pgmL", "CCL2_pgmL", "GRANA_pgmL", "IFNg_pgmL", "IL4_pgmL", 
# "TNFa_pgmL", "B7H1_pgmL", "Endothelin1_pgmL", "IFNa_pgmL", "IL10_pgmL", "IL7_pgmL", 
# "IL15_pgmL", "CXCL10_pgmL", "IL6_pgmL", "DimerD_pgmL", "RANTESCCL5_pgmL", 
# "PENTRAXIN3TSG14_pgmL", "Lipocalin2NGAL_pgmL", "GMCSF_pgmL", "N1_copies_mL_plasm_1304", 
# "anti_N_título", "anti_S_título"

## --------- Plot 1: IL2 ---------
upper_IL2 <- quantile(bddanciano$IL2_pgmL, 0.99, na.rm = TRUE) # na.rm = TRUE ignora los NAs

p1 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = IL2_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5 
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .2, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7 # width = .1 si: 4 en línea
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_IL2)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "IL-2 (pg/mL)")

## --------- Plot 2: IL4 ---------
upper_IL4 <- quantile(bddanciano$IL4_pgmL, 1, na.rm = TRUE) # na.rm = TRUE ignora los NAs

p2 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = IL4_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5 
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .2, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7 # width = .1 si: 4 en línea
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_IL4)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "IL-4 (pg/mL)")

## --------- Plot 3: GMCSF ---------
upper_gmcsf <- quantile(bddanciano$GMCSF_pgmL, 0.95, na.rm = TRUE) # na.rm = TRUE ignora los NAs

p3 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = GMCSF_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5 
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .2, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7 # width = .1 si: 4 en línea
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_gmcsf)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "GM-CSF (pg/mL)")

## --------- Plot 4: ccl5 ---------
upper_ccl5 <- max(bddanciano$RANTESCCL5_pgmL, na.rm = TRUE) * 1.1

p4 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = RANTESCCL5_pgmL)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5 
  ) +
  
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .2, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7 # width = .1 si: 4 en línea
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_ccl5)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "CCL5 (pg/mL)")

## --------- Juntar ---------
# (p1 | p2 | p3 | p4)
# wrap_plots(list(p1, p2, p3), ncol = 3)
wrap_plots(list(p1, p2, p3, p4), ncol = 2, nrow = 2)


## para meter más de 2x2 gráficos.
# library(gridExtra)
# do.call(grid.arrange, c(plots_list, ncol = 2))

# Cerramos la imagen .svg generada para el plot
dev.off()


#################################################################################
#################################################################################
################# >>>>>>>>>>> BUENO <<<<<<<<<<<<<<<<<< variables por separado
########## n1copies - anti S - anti N


# vamos a guardar el plot como .png:
png("violinp_n1copies_antis_antin6.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
# svg("violin_plots1.svg", width=8, height=6)

## --------- Plot 1: N1copies ---------
upper_N1 <- max(bddanciano$N1_copies_mL_plasm_1304, na.rm = TRUE) * 1.1

p1 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = N1_copies_mL_plasm_1304)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5 
  ) +
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .25, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_N1)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "Viral RNA load in plasma")

## --------- Plot 2: N1copies ---------
upper_N1 <- quantile(bddanciano$N1_copies_mL_plasm_1304, 0.88, na.rm = TRUE)

p2 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = N1_copies_mL_plasm_1304)) +
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .25, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .05, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_N1)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "Viral RNA load in plasma")

## --------- Plot 3: IgGantiS ---------
upper_IgGs <- quantile(bddanciano$anti_S_título, 0.95, na.rm = TRUE)

p3 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = anti_S_título)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5 
  ) +
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .3, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_IgGs)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "Anti SARS-CoV-2 S1 IgG (UA/mL)")

## --------- Plot 4: IgGantiS ---------
upper_Iggn <- max(bddanciano$anti_N_título, na.rm = TRUE) * 1.3

p4 <- ggplot(bddanciano, aes(x = Edad_grupos_mínhasta70_71hastamáx, y = anti_N_título)) +
  geom_violin(
    aes(colour = Edad_grupos_mínhasta70_71hastamáx),
    fill = "white", alpha = .5, trim = FALSE, linewidth = 0.5 
  ) +
  geom_jitter(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx,
        color = Edad_grupos_mínhasta70_71hastamáx),
    shape = 21, position = position_jitter(0.2), alpha = .3, size = 2
  ) +
  geom_boxplot(
    aes(fill = Edad_grupos_mínhasta70_71hastamáx),
    width = .07, outlier.shape = NA, color = "black", alpha = .7
  ) +
  scale_fill_manual(values = c("> 70 años" = "firebrick1", 
                               "≤ 70 años" = "goldenrod1")) +
  scale_colour_manual(values = c("> 70 años" = "firebrick4", 
                                 "≤ 70 años" = "goldenrod4")) +
  coord_cartesian(ylim = c(0, upper_Iggn)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  labs(x = "Grupo por edad", y = "Anti SARS-CoV-2 N IgG (UA/mL)")


## --------- Juntar los tres ---------
# (p1 | p2 | p3)
# wrap_plots(list(p1, p2, p3), ncol = 3)
wrap_plots(list(p1, p2, p3, p4), ncol = 2, nrow = 2)

## para meter más de 2x2 gráficos.
# library(gridExtra)
# do.call(grid.arrange, c(plots_list, ncol = 2))

# Cerramos la imagen .svg generada para el plot
dev.off()
