## Diagramas para explicar algoritmo de machine learning.
# install.packages("DiagrammeR")

library(DiagrammeR)

grViz("
digraph xgboost_esquema {
  
  graph [rankdir = TB, size = 8]

  # Estilos de nodos
  node [fontname = Helvetica]

  A [label = 'Función objetivo = Pérdida + Regularización',
     shape = box,
     style = filled,
     color = lightblue4,
     fillcolor = lightblue]

  B [label = 'Pérdida (log-loss)\nMide el error del modelo',
     shape = ellipse]

  C [label = 'Regularización (λ, α)\nControla complejidad y sobreajuste',
     shape = ellipse]

  D [label = 'Expansión de 2º orden (Taylor)',
     shape = box,
     style = filled,
     color = grey49,
     fillcolor = lightgrey]

  E [label = 'Gradiente (G):\nmedida del error',
     shape = ellipse]

  F [label = 'Hessiana (H):\nconfianza de la predicción',
     shape = ellipse]

  G [label = 'Derivar función objetivo respecto a w\n→ calcular peso óptimo w*',
     shape = box,
     style = filled,
     color = lemonchiffon3,
     fillcolor = lemonchiffon2]

  H [label = 'Función objetivo por hoja\n(pérdida + regularización (λ, α))',
     shape = ellipse]

  I [label = 'Calcular reducción de pérdida\n(con w* incluido)',
     shape = box,
     style = filled,
     color = lemonchiffon3,
     fillcolor = lemonchiffon2]

  J [label = 'Comparar:\n(Hijas) - (Madre)',
     shape = box,
     style = filled,
     color = hotpink3,
     fillcolor = lightpink]

  K [label = '¿Ganancia > 0\n(reducción > γ)?',
     shape = diamond,
     style = filled,
     color = darkorange2,
     fillcolor = tan1]

  L1 [label = 'Sí → Hacer Split',
      shape = box,
      style = filled,
      color = forestgreen,
      fillcolor = palegreen]

  L2 [label = 'No → No dividir',
      shape = box,
      style = filled,
      color = firebrick3,
      fillcolor = salmon]

  # Conexiones
  A -> B
  A -> C
  B -> D
  D -> E
  D -> F
  E -> G
  F -> G
  G -> H
  H -> I
  I -> J
  J -> K
  K -> L1 [label = 'Sí']
  K -> L2 [label = 'No']
}
")



