# install.packages("DiagrammeR")
# install.packages("DiagrammeRsvg")
# install.packages("rsvg")

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

g2 <- grViz("
digraph xgboost_shap {
  graph [rankdir = TB, size = 8]
  node [fontname = Helvetica]

  # --- XGBoost ---
  A [label = 'XGBoost\\nFunción objetivo = Pérdida + Regularización',
     shape = box, style = filled, color = lightblue4, fillcolor = lightblue]

  B [label = 'Árboles de decisión\\n(entrenados con G y H)',
     shape = ellipse, color = lightgrey]

  C [label = 'Predicción del modelo f(x)\\n(en log-odds o probabilidad)',
     shape = box, style = filled, color = lemonchiffon3, fillcolor = lemonchiffon2]

  # --- SHAP ---
  D [label = 'Valor base E[f(x)]',
     shape = ellipse, color = lightgrey]

  E [label = 'phi_j(x): contribución de cada atributo',
     shape = ellipse, color = lightgreen]

  F [label = 'Descomposición aditiva:\\n f(x) = E[f(x)] + Σ phi_j(x)',
     shape = box, style = filled, color = lightpink3, fillcolor = lightpink]

  # Conexiones
  A -> B
  B -> C
  C -> F
  D -> F
  E -> F
}
")

# --- Exportar ---
# Guardar como SVG
export_svg(g2) %>% charToRaw() %>% writeBin("xgboost_shap_esquema.svg")

# Guardar como PNG en alta resolución
rsvg_png(charToRaw(export_svg(g2)),
         file = "xgboost_shap_esquema.png",
         width = 1200, height = 1600)  # ajusta tamaño según lo necesites
