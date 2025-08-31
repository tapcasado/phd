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
# install.packages("ggcorrplot")
library(ggcorrplot)
# install.packages("corrplot")
library(corrplot)
# install.packages("psych")
library(psych)
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


#################################################################################
### VAMOS A CAMBIAR EL NOMBRE A TODAS LAS VARIABLES Y AL GRUPO PARA REPRESENTACIÓN FINAL.
################################################################################

# View(bddanciano)

bddanciano$TNFa <- bddanciano$TNFa_pgmL
bddanciano$IL2 <- bddanciano$IL2_pgmL
bddanciano$GMCSF <- bddanciano$GMCSF_pgmL
bddanciano$CCL2 <- bddanciano$CCL2_pgmL
bddanciano$GranzimaA <- bddanciano$GRANA_pgmL
bddanciano$IFNg <- bddanciano$IFNg_pgmL
bddanciano$IFNa <- bddanciano$IFNa_pgmL
bddanciano$IL4 <- bddanciano$IL4_pgmL
bddanciano$PDL1 <- bddanciano$B7H1_pgmL
bddanciano$Endotelina1 <- bddanciano$Endothelin1_pgmL
bddanciano$IL6 <- bddanciano$IL6_pgmL
bddanciano$IL7 <- bddanciano$IL7_pgmL
bddanciano$IL15 <- bddanciano$IL15_pgmL
bddanciano$IL10 <- bddanciano$IL10_pgmL
bddanciano$CXCL10 <- bddanciano$CXCL10_pgmL
bddanciano$DimeroD <- bddanciano$DimerD_pgmL
bddanciano$CCL5 <- bddanciano$RANTESCCL5_pgmL
bddanciano$Pentraxina3 <- bddanciano$PENTRAXIN3TSG14_pgmL
bddanciano$Lipocalina2 <- bddanciano$Lipocalin2NGAL_pgmL
bddanciano$N1copias <- bddanciano$N1_copies_mL_plasm_1304
bddanciano$IgGantiN <- bddanciano$anti_N_título
bddanciano$IgGantiS <- bddanciano$anti_S_título

younger_group <- subset(bddanciano, 
                        bddanciano$Edad_grupos_mínhasta70_71hastamáx=="≤ 70 años")

older_group <- subset(bddanciano, 
                      bddanciano$Edad_grupos_mínhasta70_71hastamáx=="> 70 años")

################################################################################
############################### HEATMAPS #######################################
################################################################################
### INSTRUCCIONES:

variables_seleccionadas <- bddanciano[, c("Edad", "Pentraxina3", "Lipocalina2", "TNFa", "IL6", "IFNa", "IFNg", 
                                  "IL7", "GranzimaA", "DimeroD", "Endotelina1", "IL15", "CCL2", 
                                  "CXCL10", "IL10", "PDL1",  "IL2", "IL4", "GMCSF", "CCL5", "N1copias",
                                  "IgGantiS", "IgGantiN")]
                                           
matrix_age <- corr.test(variables_seleccionadas, method = "spearman", alpha = 0.05, 
                     use = "pairwise.complete.obs")   

matrix1 <- matrix_age$r

sig1 <- matrix_age$p

## best plots:

corrplot(matrix1, p.mat = sig1, method = 'color', type = 'lower', diag = FALSE,
         sig.level = c(0.01, 0.05), pch.cex = 0.6,
         insig = 'label_sig', tl.col = "black")

corrplot(matrix1, p.mat = sig1, method = 'color', type = 'lower', diag = FALSE,
         sig.level = c(0.01, 0.05), pch.cex = 0.6,
         insig = 'label_sig', tl.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(200))

corrplot(matrix1, p.mat = sig1, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', tl.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(200))


################################################################################
############ FUNCIÓN, ENTERA Y POR GRUPOS DE EDAD, corrplor upper
################################################################################

plot_hm_upper <- function(df) {
  variables_seleccionadas <- df[, c("Edad", "Pentraxina3", "Lipocalina2", "TNFa", "IL6", "IFNa", "IFNg", 
                                    "IL7", "GranzimaA", "DimeroD", "Endotelina1", "IL15", "CCL2", 
                                    "CXCL10", "IL10", "PDL1",  "IL2", "IL4", "GMCSF", "CCL5", "N1copias",
                                    "IgGantiS", "IgGantiN")]
  
  
  matrix <- corr.test(variables_seleccionadas, method = "spearman", alpha = 0.05, 
                          use = "pairwise.complete.obs")   
  matrix_r <- matrix$r
  matrix_p <- matrix$p
  corrplot(matrix_r, p.mat = matrix_p, method = 'color', diag = FALSE, type = 'upper',
           sig.level = c(0.01, 0.05), pch.cex = 0.9,
           insig = 'label_sig', tl.col = "black")
}

#png("plot-hm-whole-upp.png", width=1200, height=1000, res=150)
svg("plot-hm-whole-upp.svg", width=8, height=6)
plot_hm_upper(bddanciano)
dev.off()
# png("plot-hm-old-upp.png", width=1200, height=1000, res=150)
svg("plot-hm-old-upp.svg", width=8, height=6)
plot_hm_upper(older_group)
dev.off()
# png("plot-hm-young-upp.png", width=1200, height=1000, res=150)
svg("plot-hm-young-upp.svg", width=8, height=6)
plot_hm_upper(younger_group)
dev.off()



################################################################################
############ FUNCIÓN, ENTERA Y POR GRUPOS DE EDAD, corrplor lower
################################################################################

plot_hm_lower <- function(df) {
  variables_seleccionadas <- df[, c("Edad", "Pentraxina3", "Lipocalina2", "TNFa", "IL6", "IFNa", "IFNg", 
                                    "IL7", "GranzimaA", "DimeroD", "Endotelina1", "IL15", "CCL2", 
                                    "CXCL10", "IL10", "PDL1",  "IL2", "IL4", "GMCSF", "CCL5", "N1copias",
                                    "IgGantiS", "IgGantiN")]
  
  
  matrix <- corr.test(variables_seleccionadas, method = "spearman", alpha = 0.05, 
                      use = "pairwise.complete.obs")   
  matrix_r <- matrix$r
  matrix_p <- matrix$p
  corrplot(matrix_r, p.mat = matrix_p, method = 'color', diag = FALSE, type = 'lower',
           sig.level = c(0.01, 0.05), pch.cex = 0.9,
           insig = 'label_sig', tl.col = "black")
}

# png("plot-hm-whole-low.png", width=1200, height=1000, res=150)
svg("plot-hm-whole-low.svg", width=8, height=6)
plot_hm_lower(bddanciano)
dev.off()

# png("plot-hm-old-low.png", width=1200, height=1000, res=150)
svg("plot-hm-old-low.svg", width=8, height=6)
plot_hm_lower(older_group)
dev.off()

# png("plot-hm-young-low.png", width=1200, height=1000, res=150)
svg("plot-hm-young-low.svg", width=8, height=6)
plot_hm_lower(younger_group)
dev.off()


################################################################################
######## COMPROBACIONES: simetría de las matrices de los p valores #############
################################################################################

### tenemos valores faltantes en alguna de las variables y esto va a cambiar cómo
### va a calcularse el valor de las correlaciones y, por tanto, los p valores
### vamos a ver si als matrices son simétricas (igual número de valores, si no es 
## así hay que forzar simetría y quedarnos con los valores más restructivos para
## ser fieles a los resultados:

variables_seleccionadas <- bddanciano[, c("Edad", "Pentraxina3", "Lipocalina2", "TNFa", "IL6", "IFNa", "IFNg", 
                                          "IL7", "GranzimaA", "DimeroD", "Endotelina1", "IL15", "CCL2", 
                                          "CXCL10", "IL10", "PDL1",  "IL2", "IL4", "GMCSF", "CCL5", "N1copias",
                                          "IgGantiS", "IgGantiN")]
matrix <- corr.test(variables_seleccionadas, method = "spearman", alpha = 0.05, 
                    use = "pairwise.complete.obs")   
matrix_r <- matrix$r
matrix_p <- matrix$p

matrix_p["Edad", "Pentraxina3"]
matrix_p["Pentraxina3", "Edad"] # las matrices resultantes no son simétricas

# OPCIÓN 1:
matrix_p_sym <- (matrix_p + t(matrix_p)) / 2 # forzar simetría de las matrices

corrplot(matrix_r, p.mat = matrix_p_sym, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', tl.col = "black")

# OPCIÓN 2:
matrix_p_sym <- pmax(matrix_p, t(matrix_p)) # opción más conservadora

corrplot(matrix_r, p.mat = matrix_p_sym, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', tl.col = "black")
