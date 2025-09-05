
# install.packages("DiagrammeR")
# install.packages("DiagrammeRsvg")
# install.packages("rsvg")

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)


# 1) Esquema SHAP

g_shap <- grViz("
digraph shap_esquema {
  graph [rankdir = TB, size = 8]
  node [fontname = Helvetica]

  A [label = 'Predicción del modelo f(x)',
     shape = box, style = filled, color = lightblue4, fillcolor = lightblue]

  B [label = 'Valor base E[f(x)]\\n(predicción promedio)',
     shape = ellipse, style = filled, color = grey49, fillcolor = lightgrey]

  C [label = 'phi_1(x): contribución atributo 1',
     shape = ellipse, color = grey49]

  D [label = 'phi_2(x): contribución atributo 2',
     shape = ellipse, color = grey49]

  E [label = 'phi_j(x): contribución atributo j',
     shape = ellipse, color = grey49]

  F [label = 'f(x) = E[f(x)] + Σ phi_j(x)',
     shape = box, style = filled, color = lemonchiffon3, fillcolor = lemonchiffon2]

  # Conexiones
  B -> F
  C -> F
  D -> F
  E -> F
  F -> A
}
")

# Guardar Esquema SHAP
export_svg(g_shap) %>% charToRaw() %>% writeBin("shap_esquema.svg")
rsvg_png(charToRaw(export_svg(g_shap)),
         file = "shap_esquema.png",
         width = 1200, height = 1600)

# 2) Cálculo SHAP

g_calc3 <- grViz("
digraph shap_calculo_esquema {
  graph [rankdir = TB, size = 8]
  node [fontname = Helvetica]

  A [label = 'Instancia x (valores de atributos)',
     shape = box, style = filled, color = lightblue4, fillcolor = lightblue]

  B [label = 'Conjunto de atributos F = {1,...,p}',
     shape = box, color = lightgrey]

  C [label = 'Todos los subconjuntos S de F\\n(sin el atributo j)',
     shape = box, style = filled, color = lemonchiffon3, fillcolor = lemonchiffon2]

  D [label = 'Contribución marginal:\\n f(S + j) - f(S)',
     shape = box, color = lightgrey]

  E [label = 'Promediar sobre todos los subconjuntos\\n(ponderado por nº de combinaciones)',
      shape = box, style = filled, color = hotpink3, fillcolor = lightpink]

  F [label = 'phi_j(x): valor SHAP\\ncontribución del atributo j',
     shape = box, style = filled, color = grey49, fillcolor = lightgrey]

  # Conexiones
  A -> B
  B -> C
  C -> D
  D -> E
  E -> F
}
")

# Guardar Cálculo SHAP
export_svg(g_calc3) %>% charToRaw() %>% writeBin("shap_calculo_esquema.svg")
rsvg_png(charToRaw(export_svg(g_calc3)),
         file = "shap_calculo_esquema.png",
         width = 1200, height = 1600)
