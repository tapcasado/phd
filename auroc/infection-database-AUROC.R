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

# cargamos BDD
bddanciano <- 
  read.spss("C:/Users/tposc/master-in-bioinformatics/phd/BDD_SPSS_SYSTEMSBIOCOVID_08112023v18 - paper anciano.sav",
           rownames=FALSE, stringsAsFactors=TRUE, tolower=FALSE)
# View(bddanciano)
str(bddanciano)
dim(bddanciano)

# auroc curva
# https://xrobin.github.io/pROC/screenshots.html
# https://cran.r-project.org/web/packages/pROC/refman/pROC.html#pROC-package
# https://stackoverflow.com/questions/61105697/plotting-roc-problems-with-axes-limits-r-plots

# install.packages("pROC")
library(pROC)


################################################################################
############ APUNTES:

# 1. Calcular curva ROC
rocobj <- roc(
  var_categorica,  # Variable binaria de estado real (0/1). Ej.: mortalidad 0/1.
  var_continua,    # Variable continua que será el predictor (ej.: Edad).
  percent = TRUE,  # Devuelve sensibilidad/especificidad en porcentaje (0-100) en lugar de 0-1.
  plot = TRUE,     # Dibuja automáticamente la curva ROC al calcularla.
  asp = NA,        # Relación de aspecto del gráfico (NA = aspecto libre, no forzar 1:1).
  print.auc = TRUE # Muestra el valor del AUC (Área Bajo la Curva) en el gráfico.
)

################################################################################
############### variables definidas: 

var_categorica <- bddanciano$Mortalidad90días_mas_UCI_0_1
var_continua <- bddanciano$Edad

# paleta de colores

color_line <- "#1c61b6"
color_ic <- "#1c61b633"
color_bestpoint <- "red"


# Calculo curva ROC
rocobj <- roc(bddanciano$Mortalidad90días_mas_UCI_0_1, 
              bddanciano$Edad,
              percent=TRUE, ci=TRUE)

# vamos a guardar el plot como .png:
# png("curva_ROC_percent.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
svg("curva_ROC_percent.svg", width=8, height=6)

# Plot curva ROC principal
plot(rocobj,
     # main="Curva ROC con IC",
     col=color_line, lwd=3, lty=1,
     # legacy.axes: a logical indicating if the specificity axis (x axis) must be 
     # plotted as as decreasing “specificity” (FALSE, the default) or increasing 
     # “1 - specificity” (TRUE) as in most legacy software.
     legacy.axes=TRUE,
     xlim=c(100,0), ylim=c(0,100),
     asp = NA,
     xlab="100 - Especificidad (%)",   # etiqueta en castellano
     ylab="Sensibilidad (%)",         # etiqueta en castellano
     print.auc=TRUE, print.auc.x=50, print.auc.y=50)

# Añadir IC de la sensibilidad
ciobj <- ci.se(rocobj, specificities=seq(0, 100, 5))

plot(ciobj, type="shape", col=color_ic, add=TRUE, xlim=c(100,0), ylim=c(0,100))

# Dibujar la curva principal encima de nuevo
lines(rocobj, col=color_line, lwd=3)

# Calculo del mejor punto según Pitágoras (closest.topleft)
best_point <- coords(rocobj, "best",
                     best.method="closest.topleft",
                     ret=c("threshold","sensitivity","specificity"))

# Añadir el punto de corte en el gráfico
# Nota: en el gráfico, el eje X es 1 - especificidad
points(best_point["specificity"], 
       best_point["sensitivity"],
       col=color_bestpoint, pch=19, cex=1.5)

# añadir línea vertical y horizaontal del threshold de edad
abline(v=best_point["specificity"], col=color_bestpoint, lty=2, lwd=0.5)
abline(h=best_point["sensitivity"], col=color_bestpoint, lty=2, lwd=0.5)

# leyenda
# grid()
legend("bottomright",
       legend=c("Curva ROC", "IC sensibilidad", "Mejor corte (Pitágoras)"),
       col=c(color_line, color_ic, color_bestpoint),
       lwd=c(3, NA, 2),
       pch=c(NA, 15, 19),
       bty="n")

################################################################################
### TEXTO DEL MEJOR PUNTO: 

text(
  x = best_point["specificity"] + 17, # mover a la izq. del punto: eje x = 1 - specificity
  y = best_point["sensitivity"],  
  labels = paste0(
    "OOP = ", best_point["threshold"], " años\n",
    "Sensibilidad = ", round(best_point["sensitivity"],1), "%\n",
    "Especificidad = ", round(best_point["specificity"],1), "%"
  ),
  col = "#BA0900", cex = 0.9, pos = 3
)

# Cerramos la imagen .png generada para el plot
dev.off()
# Cerramos la imagen .svg generada para el plot
# dev.off()



###############################################################################
######### tanto por 1:

################################################################################
############### variables definidas: 

var_categorica <- bddanciano$Mortalidad90días_mas_UCI_0_1
var_continua <- bddanciano$Edad

# paleta de colores

color_line <- "#1c61b6"
color_ic <- "#1c61b633"
color_bestpoint <- "red"


# Calculo curva ROC
rocobj <- roc(bddanciano$Mortalidad90días_mas_UCI_0_1, 
              bddanciano$Edad,
              percent=FALSE, ci=TRUE)

# vamos a guardar el plot como .png:
# png("curva_ROC_perone.png", width=1200, height=1000, res=150)
# vamos a guardar el plot como SVG (vectorial)
svg("curva_ROC_perone.svg", width=8, height=6)

# Plot curva ROC principal
plot(rocobj,
     # main="Curva ROC con IC",
     col=color_line, lwd=3, lty=1,
     # legacy.axes: a logical indicating if the specificity axis (x axis) must be 
     # plotted as as decreasing “specificity” (FALSE, the default) or increasing 
     # “1 - specificity” (TRUE) as in most legacy software.
     legacy.axes=TRUE,
     xlim=c(1,0), ylim=c(0,1),
     asp = NA,
     xlab="1 - Especificidad (%)",   # etiqueta en castellano
     ylab="Sensibilidad (%)",         # etiqueta en castellano
     print.auc=TRUE, print.auc.x=0.50, print.auc.y=0.50)

# Añadir IC de la sensibilidad
ciobj <- ci.se(rocobj, specificities=seq(0, 1, 0.05))

plot(ciobj, type="shape", col=color_ic, add=TRUE, xlim=c(1,0), ylim=c(0,1))

# Dibujar la curva principal encima de nuevo
lines(rocobj, col=color_line, lwd=3)

# Calculo del mejor punto según Pitágoras (closest.topleft)
best_point <- coords(rocobj, "best",
                     best.method="closest.topleft",
                     ret=c("threshold","sensitivity","specificity"))

# Añadir el punto de corte en el gráfico
# Nota: en el gráfico, el eje X es 1 - especificidad
points(best_point["specificity"], 
       best_point["sensitivity"],
       col=color_bestpoint, pch=19, cex=1.5)

# añadir línea vertical y horizaontal del threshold de edad
abline(v=best_point["specificity"], col=color_bestpoint, lty=2, lwd=0.5)
abline(h=best_point["sensitivity"], col=color_bestpoint, lty=2, lwd=0.5)

# leyenda
# grid()
legend("bottomright",
       legend=c("Curva ROC", "IC sensibilidad", "Mejor corte (Pitágoras)"),
       col=c(color_line, color_ic, color_bestpoint),
       lwd=c(3, NA, 2),
       pch=c(NA, 15, 19),
       bty="n")

################################################################################
### TEXTO DEL MEJOR PUNTO: 

text(
  x = best_point["specificity"] + 0.17, # mover a la izq. del punto: eje x = 1 - specificity
  y = best_point["sensitivity"],  
  labels = paste0(
    "OOP = ", best_point["threshold"], " años\n",
    "Sensibilidad = ", 100 * round(best_point["sensitivity"],3), "%\n",
    "Especificidad = ", 100 * round(best_point["specificity"],3), "%"
  ),
  col = "#BA0900", cex = 0.9, pos = 3
)

# Cerramos la imagen .png generada para el plot
dev.off()
# Cerramos la imagen .svg generada para el plot
# dev.off()
