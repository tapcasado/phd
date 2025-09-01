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
library(dplyr)
library(tidyr)
library(ggplot2)
# install.packages("ggnewscale")
library(ggnewscale)
library(dplyr)
library(tidyr)   # para pivot_longer
library(dplyr)
library(tidyr)
library(ggplot2)

# Datos directamente en R
df <- data.frame(
  Virus = c("SARS-CoV-2", "MERS-CoV", "SARS-CoV-1", "Influenza A (H1N1) 2009",
            "Avian Influenza A (H7N9)", "Other Influenza type A and B", 
            "Adenovirus", "Rhinovirus", "RSV", "Multiple viruses"),
  Mortality_NF = c(4,1,0,0,2,0,0,0,0,0),
  Critical_NF = c(4,0,0,2,0,0,0,0,0,0),
  Respiratory_NF = c(3,0,1,3,0,0,0,0,0,0),
  Extrapulmonary_NF = c(0,0,0,0,0,0,0,0,0,0),
  Hospital_NF = c(0,0,0,0,0,0,0,0,0,0),
  SeverityScores_NF = c(0,0,0,0,1,1,0,0,0,0),
  SeverityScores_F = c(7,1,1,1,0,0,0,0,0,0),
  Hospital_F = c(4,0,0,0,0,0,0,0,0,0),
  Extrapulmonary_F = c(8,0,1,0,0,0,0,0,0,0),
  Respiratory_F = c(18,0,1,5,0,0,0,0,1,0),
  Critical_F = c(65,2,2,3,1,0,1,0,1,1),
  Mortality_F = c(51,3,1,4,1,0,2,1,1,1)
)

# Convertir a formato largo para "No encontrada"
df_nf <- df %>%
  select(Virus, ends_with("_NF")) %>%
  pivot_longer(-Virus, names_to = "Severity", values_to = "Count") %>%
  mutate(Association = "No encontrada",
         Severity = gsub("_NF", "", Severity))

# Convertir a formato largo para "Encontrada"
df_f <- df %>%
  select(Virus, ends_with("_F")) %>%
  pivot_longer(-Virus, names_to = "Severity", values_to = "Count") %>%
  mutate(Association = "Encontrada",
         Severity = gsub("_F", "", Severity))

# Combinar ambos
df_long <- bind_rows(df_nf, df_f)

# Traducir niveles de gravedad al castellano
df_long <- df_long %>%
  mutate(
    Severity = recode(Severity,
                      "Mortality" = "Mortalidad",
                      "Critical" = "Estado crítico",
                      "Respiratory" = "Afectación respiratoria",
                      "Extrapulmonary" = "Fallo extrapulmonar",
                      "Hospital" = "Ingreso hospitalario",
                      "SeverityScores" = "Scores de gravedad")
  )

# Ordenar virus para que aparezcan de arriba a abajo según tu lista
virus_order <- c(
  "SARS-CoV-2",
  "MERS-CoV",
  "SARS-CoV-1",
  "Influenza A (H1N1) 2009",
  "Avian Influenza A (H7N9)",
  "Other Influenza type A and B",
  "Adenovirus",
  "Rhinovirus",
  "RSV",
  "Multiple viruses"
)

df_long <- df_long %>%
  mutate(
    Virus = factor(Virus, levels = rev(virus_order)),  # invertir para ggplot
    Severity = factor(Severity, levels = c("Mortalidad", "Estado crítico", "Afectación respiratoria",
                                           "Fallo extrapulmonar", "Ingreso hospitalario", "Scores de gravedad")),
    Association = factor(Association, levels = c("No encontrada", "Encontrada"))
  )

################################################################################
### eliminamos los puntos de papers = 0
################################################################################
df_long_plot <- df_long %>%
 filter(Count > 0)
################################################################################

# Crear dot plot estilo scRNA-seq
# Crear dot plot con tamaños crecientes
ggplot(df_long_plot, aes(x = Severity, y = Virus)) +
  geom_point(
    data = df_long %>% filter(Association == "Encontrada" & Count == 0),
    aes(size = 1), shape = 3, color = "gray80" # "gray95" para que no se vea el fondo, antes: "gray80"
  ) +
  # Dibujar una cruz hueco donde no hay estudios
  geom_point(aes(size = Count, color = Association)) +
  scale_color_manual(values = c("No encontrada" = "blue", "Encontrada" = "red")) +
  scale_size_continuous(
    range = c(2, 10),       # tamaño mínimo y máximo de los puntos
    breaks = c(1,2,3,4,5,10,20,40,60)  # puntos de referencia para la leyenda
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(size = "Número de estudios", color = "Asociación",
       x = " ", y = " ")


################################################################################
##### pruebas de colores:

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Crear gradientes para cada tipo de asociación
grad_found <- scales::col_numeric(
  palette = c("mistyrose", "red"), 
  domain = c(1, max(df_long$Count, na.rm = TRUE))
)

grad_notfound <- scales::col_numeric(
  palette = c("lightblue", "blue"), 
  domain = c(1, max(df_long$Count, na.rm = TRUE))
)

# Añadir columna de colores según Association y Count
df_long <- df_long %>%
  mutate(
    Color = case_when(
      Count > 0 & Association == "Encontrada"   ~ grad_found(Count),
      Count > 0 & Association == "No encontrada" ~ grad_notfound(Count),
      TRUE ~ "gray95" # para Count == 0
    )
  )

# Gráfico
ggplot(df_long, aes(x = Severity, y = Virus)) +
  geom_point(
    data = df_long %>% filter(Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  geom_point(
    data = df_long %>% filter(Count > 0),
    aes(size = Count, color = Color)
  ) +
  scale_size_continuous(
    range = c(2, 8),
    breaks = c(1,2,3,4,5,10,20,40,60)
  ) +
  scale_color_identity() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    size = "Número de estudios",
    x = " ", y = " ",
    title = " "
  )

################################################################################
##### pruebas de colores:

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)


# Gráfico con intensidad por número de estudios
ggplot(df_long, aes(x = Severity, y = Virus)) +
  # Cruces grises para Count = 0
  geom_point(
    data = df_long %>% filter(Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos coloreados por Association, intensidad con alpha según Count
  geom_point(
    data = df_long %>% filter(Count > 0),
    aes(size = Count, color = Association, alpha = Count)
  ) +
  # Escala de colores fija para Association
  scale_color_manual(
    values = c("No encontrada" = "blue", "Encontrada" = "red"),
    name = "Asociación"
  ) +
  # Escala de tamaños
  scale_size_continuous(
    range = c(2, 8),
    breaks = c(1,2,3,5,10,20,40,60),
    name = "Número de estudios"
  ) +
  # Escala de transparencia (intensidad)
  scale_alpha_continuous(
    range = c(0.3, 1),   # 0.3 = más claro, 1 = más intenso
    guide = "none"       # ocultamos la leyenda de alpha porque ya está en size
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Nivel de gravedad",
    y = "Virus",
    title = "Asociaciones encontradas y no encontradas con intensidad según número de estudios"
  )


################################################################################
################################################################################
##### pruebas de colores (la mejor de momento):

ggplot(df_long, aes(x = Severity, y = Virus)) +
  # Cruces grises para Count = 0
  geom_point(
    data = df_long %>% filter(Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos coloreados con alpha según Count
  geom_point(
    data = df_long %>% filter(Count > 0),
    aes(size = Count, color = Association, alpha = Count)
  ) +
  # Colores base para asociación
  scale_color_manual(
    values = c("No encontrada" = "blue", "Encontrada" = "red"),
    name = "Asociación"
  ) +
  # Escala de tamaños
  scale_size_continuous(
    range = c(2, 8),
    breaks = c(1,2,3,5,10,20,40,60),
    name = "Tamaño (nº estudios)"
  ) +
  # Escala de alpha
  scale_alpha_continuous(
    range = c(0.3, 1),
    breaks = c(1,2,3,5,10,20,40,60),
    name = "Intensidad (nº estudios)"
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 0.8, size = 5)), # leyenda colores en sólido
    alpha = guide_legend(order = 2),  # leyenda separada para intensidad
    size  = guide_legend(order = 3)   # leyenda separada para tamaño
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = " ",
    y = " ",
    title = " "
  )

################################################################################
################################################################################
##### pruebas de colores (no es la mejor):



# Crear columna combinada de color según asociación y número de estudios
df_long <- df_long %>%
  mutate(ColorValue = case_when(
    Association == "Encontrada" ~ Count,
    Association == "No encontrada" ~ -Count  # negativos para separarlos visualmente
  ))

# Gráfico
ggplot(df_long %>% filter(Count > 0), aes(x = Severity, y = Virus)) +
  geom_point(
    data = df_long %>% filter(Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  geom_point(
    aes(size = Count, color = ColorValue)
  ) +
  scale_size_continuous(
    range = c(2, 8),
    breaks = c(1,2,3,5,10,20,40,60),
    name = "Número de estudios"
  ) +
  # Escala de color: azules para negativos, rojos para positivos
  scale_color_gradient2(
    low = "blue", mid = "gray90", high = "red",
    midpoint = 0,
    name = "Asociación\n(intensidad = nº estudios)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = " ",
    y = " ",
    title = " "
  )

################################################################################
################################################################################
##### pruebas de colores (buena graduación en leyenda):

# Crear variable combinada para color
df_long <- df_long %>%
  mutate(ColorValue = case_when(
    Association == "Encontrada" ~ Count,
    Association == "No encontrada" ~ -Count,
    TRUE ~ 0
  ))

# Gráfico
ggplot(df_long, aes(x = Severity, y = Virus)) +
  # Cruces grises para Count = 0
  geom_point(
    data = df_long %>% filter(Count == 0),
    shape = 3, color = "gray85", size = 2
  ) +
  # Puntos coloreados según ColorValue
  geom_point(
    data = df_long %>% filter(Count > 0),
    aes(size = Count, color = ColorValue)
  ) +
  # Escala de color continua combinada
  scale_color_gradient2(
    low = "blue",    # No encontrada
    mid = "gray90",  # 0
    high = "red",    # Encontrada
    midpoint = 0,
    limits = c(-39, 60),
    breaks = c(-39, -27, -15, 0, 15, 30, 45, 60),
    labels = c("4 (No encontrada)", "3", "2", "0", "15", "30", "45", "60 (Encontrada)"),
    name = "Número de estudios\n(azul = No encontrada,\n rojo = Encontrada)"
  ) +
  scale_size_continuous(
    range = c(2,8),
    name = "Tamaño (nº estudios)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Nivel de gravedad",
    y = "Virus",
    title = "Asociaciones encontradas y no encontradas con gradación continua de color"
  )

################################################################################
################################################################################
##### pruebas de colores (buena graduación en leyenda) version 2:

# Crear variable combinada para color
df_long <- df_long %>%
  mutate(ColorValue = case_when(
    Association == "Encontrada" ~ Count,
    Association == "No encontrada" ~ -Count,
    TRUE ~ 0
  ))

# Gráfico
ggplot(df_long, aes(x = Severity, y = Virus)) +
  # Cruces grises para Count = 0
  geom_point(
    data = df_long %>% filter(Count == 0),
    shape = 3, color = "gray85", size = 2
  ) +
  # Puntos coloreados según ColorValue
  geom_point(
    data = df_long %>% filter(Count > 0),
    aes(size = Count, color = ColorValue)
  ) +
  # Escala de color continua combinada
  scale_color_gradient2(
    low = "blue",    # No encontrada
    mid = "gray90",  # 0
    high = "red",    # Encontrada
    midpoint = 0,
    limits = c(-39, 60),
    breaks = c(-39, -27, -15, 0, 15, 30, 45, 60),
    labels = c("4 (No encontrada)", "3", "2", "0", "15", "30", "45", "60 (Encontrada)"),
    name = "Número de estudios\n(azul = No encontrada,\n rojo = Encontrada)"
  ) +
  scale_size_continuous(
    range = c(2,8),
    breaks = c(1,2,3,5,10,20,40,60),
    name = "Tamaño (nº estudios)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = " ",
    y = " ",
    title = " "
  )

################################################################################
################################################################################
################################################################################
##### pruebas de colores:

ggplot(df_long, aes(x = Severity, y = Virus)) +
  # Cruces grises para Count = 0
  geom_point(
    data = df_long %>% filter(Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos coloreados con alpha según Count
  geom_point(
    data = df_long %>% filter(Count > 0),
    aes(size = Count, color = Association, alpha = Count)
  ) +
  # Colores base para asociación
  scale_color_manual(
    values = c("No encontrada" = "blue", "Encontrada" = "red"),
    name = "Asociación"
  ) +
  # Escala de tamaños
  scale_size_continuous(
    range = c(2, 8),
    breaks = c(1,2,3,5,10,20,40,60),
    name = "Tamaño (nº estudios)"
  ) +
  # Escala de alpha
  scale_alpha_continuous(
    range = c(0.3, 1),
    breaks = c(1,2,3,5,10,20,40,60),
    name = "Intensidad (nº estudios)"
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 0.8, size = 5)), # leyenda colores en sólido
    alpha = guide_legend(order = 2),  # leyenda separada para intensidad
    size  = guide_legend(order = 3)   # leyenda separada para tamaño
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = " ",
    y = " ",
    title = " "
  )


################################################################################
###### POR GRUPO:
################################################################################

################################################################################
################################################################################
## SOLUCIONES PARA VER LAS LÍNEAS DE VIRUS DONDE NO HAY ESTUDIOS:
##
## 1 geom_point de los puntos donde no hay casos:
## geom_point(
## data = df_long %>% filter(Association == "Encontrada" & Count == 0),
## aes(size = 1), shape = 3, color = "gray95" # "gray95" para que no se vea el fondo
## ) +
## luego geom_point del resto de casos:
## geom_point(aes(size = Count, color = Association)) +
################################################################################

# Ordenar virus para que aparezcan de arriba a abajo según tu lista
virus_order <- c(
  "SARS-CoV-2",
  "MERS-CoV",
  "SARS-CoV-1",
  "Influenza A (H1N1) 2009",
  "Avian Influenza A (H7N9)",
  "Other Influenza type A and B",
  "Adenovirus",
  "Rhinovirus",
  "RSV",
  "Multiple viruses"
)


### POR GRUPOS:

# Filtrar solo "Encontrada"
df_encontrada <- df_long_plot %>%
  filter(Association == "Encontrada")

# Gráfico "Encontrada"
ggplot(df_encontrada, aes(x = Severity, y = Virus)) +
  # Dibujar una cruz donde no hay estudios
  geom_point(
    data = df_long %>% filter(Association == "Encontrada" & Count == 0),
    aes(size = 1), shape = 3, color = "gray95" # "gray95" para que no se vea el fondo
  ) +
  geom_point(aes(size = Count, color = Association)) +
  scale_color_manual(values = c("Encontrada" = "red")) +
  scale_size_continuous(range = c(2,10), breaks = c(1,2,3,4,5,10,20,40,60)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(size = "Número de estudios", color = "Asociación",
       x = " ", y = "",
       title = "Asociación Encontrada")

# Filtrar solo "No encontrada"
df_noencontrada <- df_long_plot %>%
  filter(Association == "No encontrada")

# Gráfico "No encontrada"
ggplot(df_noencontrada, aes(x = Severity, y = Virus)) +
  # Dibujar una cruz donde no hay estudios
  geom_point(
    data = df_long %>% filter(Association == "Encontrada" & Count == 0),
    aes(size = 1), shape = 3, color = "gray95" # "gray95" para que no se vea el fondo
  ) +
  geom_point(aes(size = Count, color = Association)) +
  scale_color_manual(values = c("No encontrada" = "blue")) +
  scale_size_continuous(range = c(2,5), breaks = c(1,2,3,4,5)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(size = "Número de estudios", color = "Asociación",
       x = " ", y = " ",
       title = "Asociación No Encontrada")

################################################################################
## POR GRUPOS EN UNA MISMA IMAGEN:


# Gráfico con gradiente de color por Association y Count
ggplot(df_long, aes(x = Severity, y = Virus)) +
  # Puntos cuando Count == 0 (cruz gris, no entra en escala)
  geom_point(
    data = df_long %>% filter(Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos con estudios > 0 (con gradiente por asociación y Count)
  geom_point(
    data = df_long %>% filter(Count > 0),
    aes(size = Count, color = Count)
  ) +
  scale_size_continuous(
    range = c(2, 8),
    breaks = c(1,2,3,4,5,10,20,40,60)
  ) +
  scale_color_gradientn(
    colours = c("gray80", "red"),
    values = scales::rescale(c(1, max(df_long$Count, na.rm = TRUE))),
    guide = "colorbar",
    name = "Número de estudios"
  ) +
  facet_wrap(~Association, ncol = 1) +  # separar Encontrada / No encontrada
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    size = "Número de estudios",
    x = " ", y = " ",
    title = " "
  )

#################################################################################
################################################################################


#################################################################################
## POR GRUPOS: grupo encontrada.

df_encontrada <- df_long_plot %>%
  filter(Association == "Encontrada",)

# Gráfico con gradiente de color por Association y Count
ggplot(df_encontrada, aes(x = Severity, y = Virus)) +
  # Puntos cuando Count == 0 (cruz gris, no entra en escala)
  geom_point(
    data = df_long %>% filter(Association == "Encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos con estudios > 0 (con gradiente por asociación y Count)
  geom_point(
    data = df_encontrada %>% filter(Count > 0),
    aes(size = Count, color = Count)
  ) +
  scale_size_continuous(
    range = c(2, 10),
    breaks = c(1,2,3,4,5,10,20,40,60)
  ) +
  scale_color_gradientn(
    colours = c("red", "red4"),
    values = scales::rescale(c(1, max(df_long$Count, na.rm = TRUE))),
    guide = "colorbar",
    name = "Número de estudios"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    size = "Número de estudios",
    x = " ", y = " ",
    title = " "
  )

#################################################################################
## POR GRUPOS: grupo no encontrada.

df_noencontrada <- df_long_plot %>%
  filter(Association == "No encontrada")

# Gráfico con gradiente de color por Association y Count
ggplot(df_noencontrada, aes(x = Severity, y = Virus)) +
  # Puntos cuando Count == 0 (cruz gris, no entra en escala)
  geom_point(
    data = df_long %>% filter(Association == "No encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos con estudios > 0 (con gradiente por asociación y Count)
  geom_point(
    data = df_noencontrada %>% filter(Count > 0),
    aes(size = Count, color = Count)
  ) +
  scale_size_continuous(
    range = c(2, 5),
    breaks = c(1,2,3,4)
  ) +
  scale_color_gradientn(
    colours = c("blue", "blue4"),
    values = scales::rescale(c(1, max(df_long$Count, na.rm = TRUE))),
    guide = "colorbar",
    name = "Número de estudios"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    size = "Número de estudios",
    x = " ", y = " ",
    title = " "
  )


################################################################################
####

library(patchwork)

# Gráfico encontrada
p_encontrada <- ggplot(df_encontrada, aes(x = Severity, y = Virus)) +
  geom_point(
    data = df_long %>% filter(Association == "Encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  geom_point(
    data = df_encontrada %>% filter(Count > 0),
    aes(size = Count, color = Count)
  ) +
  scale_size_continuous(range = c(2, 10), breaks = c(1,2,3,4,5,10,20,40,60)) +
  scale_color_gradientn(colours = c("red", "red4"),
                        values = scales::rescale(c(1, max(df_long$Count, na.rm = TRUE))),
                        guide = "colorbar",
                        name = "Número de estudios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Encontrada", x = " ", y = " ")

# Gráfico no encontrada
p_noencontrada <- ggplot(df_noencontrada, aes(x = Severity, y = Virus)) +
  geom_point(
    data = df_long %>% filter(Association == "No encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  geom_point(
    data = df_noencontrada %>% filter(Count > 0),
    aes(size = Count, color = Count)
  ) +
  scale_size_continuous(range = c(2, 10), breaks = c(1,2,3,4,5,10,20,40,60)) +
  scale_color_gradientn(colours = c("blue", "blue4"),
                        values = scales::rescale(c(1, max(df_long$Count, na.rm = TRUE))),
                        guide = "colorbar",
                        name = "Número de estudios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "No encontrada", x = " ", y = " ")

# Combinar con patchwork
library(patchwork)
p_encontrada / p_noencontrada  # uno arriba (Encontrada) y otro abajo (No encontrada)

##################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# Crear variable continua combinada para color
df_long_plot <- df_long %>%
  mutate(
    ColorValue = case_when(
      Association == "Encontrada" & Count > 0 ~ Count,
      Association == "No encontrada" & Count > 0 ~ -Count,
      TRUE ~ 0
    )
  )

# Gráfico encontrada
p_encontrada <- ggplot(df_long_plot, aes(x = Severity, y = Virus)) +
  # Cruces grises para Count = 0
  geom_point(
    data = df_long_plot %>% filter(Association == "Encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos con gradiente rojo
  geom_point(
    data = df_long_plot %>% filter(Association == "Encontrada" & Count > 0),
    aes(size = Count, color = ColorValue)
  ) +
  scale_size_continuous(range = c(2,10), breaks = c(1,2,3,4,5,10,20,40,60), name = "Número de estudios") +
  scale_color_gradient2(
    low = "blue",   # No encontrada
    mid = "gray90", # 0
    high = "red",   # Encontrada
    midpoint = 0,
    limits = c(-max(df_long$Count), max(df_long$Count)),
    name = "Asociación y número de estudios"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Encontrada", x = " ", y = " ")

# Gráfico no encontrada
p_noencontrada <- ggplot(df_long_plot, aes(x = Severity, y = Virus)) +
  # Cruces grises para Count = 0
  geom_point(
    data = df_long_plot %>% filter(Association == "No encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos con gradiente azul
  geom_point(
    data = df_long_plot %>% filter(Association == "No encontrada" & Count > 0),
    aes(size = Count, color = ColorValue)
  ) +
  scale_size_continuous(range = c(2,10), breaks = c(1,2,3,4,5,10,20,40,60), name = "Número de estudios") +
  scale_color_gradient2(
    low = "blue",
    mid = "gray90",
    high = "red",
    midpoint = 0,
    limits = c(-max(df_long$Count), max(df_long$Count)),
    name = "Asociación y número de estudios"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "No encontrada", x = " ", y = " ")

# Combinar gráficos con patchwork y una sola leyenda de tamaño
p_encontrada / p_noencontrada + plot_layout(guides = "collect")

#################################################################################


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# Crear variable continua combinada para color
df_long_plot <- df_long %>%
  mutate(
    ColorValue = case_when(
      Association == "Encontrada" & Count > 0 ~ Count,
      Association == "No encontrada" & Count > 0 ~ -Count,
      TRUE ~ 0
    )
  )

# Gráfico encontrada
p_encontrada <- ggplot(df_long_plot, aes(x = Severity, y = Virus)) +
  # Cruces grises para Count = 0
  geom_point(
    data = df_long_plot %>% filter(Association == "Encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos con gradiente rojo/azul según ColorValue
  geom_point(
    data = df_long_plot %>% filter(Association == "Encontrada" & Count > 0),
    aes(size = Count, color = ColorValue)
  ) +
  scale_size_continuous(
    range = c(2,8),
    breaks = c(1,2,3,5,10,20,40,60),
    name = "Tamaño (nº estudios)"
  ) +
  # Escala de color continua combinada
  scale_color_gradient2(
    low = "blue3",    # No encontrada
    mid = "gray90",  # 0
    high = "red3",    # Encontrada
    midpoint = 0,
    limits = c(-39, 60),
    breaks = c(-39, -27, -15, 0, 15, 30, 45, 60),
    labels = c("4 (No encontrada)", "3", "2", "0", "15", "30", "45", "60 (Encontrada)"),
    name = "Número de estudios\n(azul = No encontrada,\n rojo = Encontrada)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Encontrada", x = " ", y = " ",)

# Gráfico no encontrada
p_noencontrada <- ggplot(df_long_plot, aes(x = Severity, y = Virus)) +
  # Cruces grises para Count = 0
  geom_point(
    data = df_long_plot %>% filter(Association == "No encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos con gradiente rojo/azul según ColorValue
  geom_point(
    data = df_long_plot %>% filter(Association == "No encontrada" & Count > 0),
    aes(size = Count, color = ColorValue)
  ) +
  scale_size_continuous(
    range = c(2,8),
    breaks = c(1,2,3,5,10,20,40,60),
    name = "Tamaño (nº estudios)"
  ) +
  # Escala de color continua combinada
  scale_color_gradient2(
    low = "blue3",    # No encontrada
    mid = "gray90",  # 0
    high = "red3",    # Encontrada
    midpoint = 0,
    limits = c(-39, 60),
    breaks = c(-39, -27, -15, 0, 15, 30, 45, 60),
    labels = c("4 (No encontrada)", "3", "2", "0", "15", "30", "45", "60 (Encontrada)"),
    name = "Número de estudios\n(azul = No encontrada,\n rojo = Encontrada)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "No encontrada", x = " ", y = " ",)

# Combinar gráficos con patchwork y una sola leyenda de tamaño
p_encontrada / p_noencontrada + plot_layout(guides = "collect")

################################################################################
#### MÁS PRUEBAS:

# Crear variable continua combinada para color
df_long_plot <- df_long %>%
  mutate(
    ColorValue = case_when(
      Association == "Encontrada" & Count > 0 ~ Count,
      Association == "No encontrada" & Count > 0 ~ -Count,
      TRUE ~ 0
    )
  )

df_encontrada <- df_long_plot %>%
  filter(Association == "Encontrada",)

# Gráfico con gradiente de color por Association y Count
p_encontrada <- ggplot(df_encontrada, aes(x = Severity, y = Virus)) +
  # Puntos cuando Count == 0 (cruz gris, no entra en escala)
  geom_point(
    data = df_long_plot %>% filter(Association == "Encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos con estudios > 0 (con gradiente por asociación y Count)
  geom_point(
    data = df_long_plot %>% filter(Association == "Encontrada" & Count > 0),
    aes(size = Count, color = ColorValue)
  ) +
  scale_size_continuous(
    range = c(2, 8),
    breaks = c(1,2,3,4,5,10,20,40,60)
  ) +
  # Escala de color continua combinada
  scale_color_gradient2(
    low = "blue3",    # No encontrada
    mid = "gray90",  # 0
    high = "red3",    # Encontrada
    midpoint = 0,
    limits = c(-4, 60),
    breaks = c(-4, -3, -2, -1, 0, 15, 30, 45, 60),
    labels = c("4 (No encontrada)", "3", "2", "1", "0", "15", "30", "45", "60 (Encontrada)"),
    name = "Número de estudios\n(azul = No encontrada,\n rojo = Encontrada)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    size = "Número de estudios",
    x = " ", y = " ",
    title = " "
  )

#################################################################################
## POR GRUPOS: grupo no encontrada.

df_noencontrada <- df_long_plot %>%
  filter(Association == "No encontrada")

# Gráfico con gradiente de color por Association y Count
p_noencontrada <- ggplot(df_noencontrada, aes(x = Severity, y = Virus)) +
  # Puntos cuando Count == 0 (cruz gris, no entra en escala)
  geom_point(
    data = df_long_plot %>% filter(Association == "No encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos con estudios > 0 (con gradiente por asociación y Count)
  geom_point(
    data = df_long_plot %>% filter(Association == "No encontrada" & Count > 0),
    aes(size = Count, color = ColorValue)
  ) +
  scale_size_continuous(
    range = c(2, 4),
    breaks = c(1,2,3,4,5,10,20,40,60)
  ) +
  # Escala de color continua combinada
  scale_color_gradient2(
    low = "blue3",    # No encontrada
    mid = "gray90",  # 0
    high = "red3",    # Encontrada
    midpoint = 0,
    limits = c(-4, 60),
    breaks = c(-4, -3, -2, -1, 0, 15, 30, 45, 60),
    labels = c("4 (No encontrada)", "3", "2", "1", "0", "15", "30", "45", "60 (Encontrada)"),
    name = "Número de estudios\n(azul = No encontrada,\n rojo = Encontrada)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    size = "Número de estudios",
    x = " ", y = " ",
    title = " "
  )

# Combinar gráficos con patchwork y una sola leyenda de tamaño
p_encontrada / p_noencontrada + plot_layout(guides = "collect")

###################################################################################
###################### CAMBIO DE BREAKS EN LEGEND

# Crear variable continua combinada para color
df_long_plot <- df_long %>%
  mutate(
    ColorValue = case_when(
      Association == "Encontrada" & Count > 0 ~ Count,
      Association == "No encontrada" & Count > 0 ~ -Count,
      TRUE ~ 0
    )
  )

df_encontrada <- df_long_plot %>%
  filter(Association == "Encontrada",)

# Gráfico con gradiente de color por Association y Count
p_encontrada <- ggplot(df_encontrada, aes(x = Severity, y = Virus)) +
  # Puntos cuando Count == 0 (cruz gris, no entra en escala)
  geom_point(
    data = df_long_plot %>% filter(Association == "Encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos con estudios > 0 (con gradiente por asociación y Count)
  geom_point(
    data = df_long_plot %>% filter(Association == "Encontrada" & Count > 0),
    aes(size = Count, color = ColorValue)
  ) +
  scale_size_continuous(
    range = c(2, 8),
    breaks = c(1,2,3,4,5,10,20,40,60)
  ) +
  # Escala de color continua combinada
  scale_color_gradient2(
    low = "blue3",    # No encontrada
    mid = "gray90",  # 0; antes = mid = "gray90"
    high = "red3",    # Encontrada
    midpoint = 0,
    limits = c(-4, 60),
    breaks = c(-4, 0, 15, 30, 45, 60),
    labels = c("4 (No encontrada)", "0", "15", "30", "45", "60 (Encontrada)"),
    name = "Número de estudios\n(azul = No encontrada,\n rojo = Encontrada)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    size = "Número de estudios",
    x = " ", y = " ",
    title = " "
  )

#################################################################################
## POR GRUPOS: grupo no encontrada.

df_noencontrada <- df_long_plot %>%
  filter(Association == "No encontrada")

# Gráfico con gradiente de color por Association y Count
p_noencontrada <- ggplot(df_noencontrada, aes(x = Severity, y = Virus)) +
  # Puntos cuando Count == 0 (cruz gris, no entra en escala)
  geom_point(
    data = df_long_plot %>% filter(Association == "No encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  # Puntos con estudios > 0 (con gradiente por asociación y Count)
  geom_point(
    data = df_long_plot %>% filter(Association == "No encontrada" & Count > 0),
    aes(size = Count, color = ColorValue)
  ) +
  scale_size_continuous(
    range = c(2, 4),
    breaks = c(1,2,3,4,5,10,20,40,60)
  ) +
  # Escala de color continua combinada
  scale_color_gradient2(
    low = "blue3",    # No encontrada
    mid = "gray90" ,  # 0
    high = "red3",    # Encontrada
    midpoint = 0,
    limits = c(-4, 60),
    breaks = c(-4, 0, 15, 30, 45, 60),
    labels = c("4 (No encontrada)", "0", "15", "30", "45", "60 (Encontrada)"),
    name = "Número de estudios\n(azul = No encontrada,\n rojo = Encontrada)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    size = "Número de estudios",
    x = " ", y = " ",
    title = " "
  )

# Combinar gráficos con patchwork y una sola leyenda de tamaño
p_encontrada / p_noencontrada + plot_layout(guides = "collect")


################################################################################

library(patchwork)

# Gráfico encontrada
p_encontrada <- ggplot(df_encontrada, aes(x = Severity, y = Virus)) +
  geom_point(
    data = df_long %>% filter(Association == "Encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  geom_point(
    data = df_encontrada %>% filter(Count > 0),
    aes(size = Count, color = Count)
  ) +
  scale_size_continuous(range = c(2, 10), breaks = c(1,2,3,4,5,10,20,40,60)) +
  scale_color_gradientn(colours = c("red", "red4"),
                        values = scales::rescale(c(1, max(df_long$Count, na.rm = TRUE))),
                        guide = "colorbar",
                        name = "Número de estudios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Encontrada", x = " ", y = " ")

# Gráfico no encontrada
p_noencontrada <- ggplot(df_noencontrada, aes(x = Severity, y = Virus)) +
  geom_point(
    data = df_long %>% filter(Association == "No encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  geom_point(
    data = df_noencontrada %>% filter(Count > 0),
    aes(size = Count, color = Count)
  ) +
  scale_size_continuous(range = c(2, 10), breaks = c(1,2,3,4,5,10,20,40,60)) +
  scale_color_gradientn(colours = c("blue", "blue4"),
                        values = scales::rescale(c(1, max(df_long$Count, na.rm = TRUE))),
                        guide = "colorbar",
                        name = "Número de estudios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "No encontrada", x = " ", y = " ")

# Combinar con patchwork
library(patchwork)
p_encontrada / p_noencontrada  # uno arriba (Encontrada) y otro abajo (No encontrada)



#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
### ÉL GRÁFICO QUE MÁS ME CONVENCE 1:

library(patchwork)

png("bubble-plot-viremia2.png", width=1350, height=1500, res=150)
# svg("plot-hm-whole-upp.svg", width=8, height=6)

# Gráfico encontrada
p_encontrada <- ggplot(df_encontrada, aes(x = Severity, y = Virus)) +
  geom_point(
    data = df_long %>% filter(Association == "Encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  geom_point(
    data = df_encontrada %>% filter(Count > 0),
    aes(size = Count, color = Count)
  ) +
  scale_size_continuous(range = c(2, 8), breaks = c(1,2,3,4,5,10,20,40,60)) +
  scale_color_gradientn(
    colours = c("#FFCCCC", "#FF6363", "#FF0000", "#990000"),  # varios tonos de rojo
    values = scales::rescale(c(1, 5, 20, max(df_long$Count, na.rm = TRUE))), 
    guide = "colorbar",
    name = "Número de estudios que\nencuentran asociación"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = " ", x = " ", y = " ")

# Gráfico no encontrada
p_noencontrada <- ggplot(df_noencontrada, aes(x = Severity, y = Virus)) +
  geom_point(
    data = df_long %>% filter(Association == "No encontrada" & Count == 0),
    shape = 3, color = "gray95", size = 2
  ) +
  geom_point(
    data = df_long_plot %>% filter(Association == "No encontrada" & Count > 0),
    aes(size = Count, color = Count)
  ) +
  scale_size_continuous(guide = "none") + # oculta solo la leyenda de TAMAÑO
  # scale_color_continuous(name = "Número de estudios")  
  scale_size_continuous(
    range = c(2, 4),
    breaks = c(1,2,3,4,5,10,20,40,60),
    guide = "none"
  ) +
  scale_color_gradientn(colours = c("cornflowerblue", "blue4"),
                        values = scales::rescale(c(1, max(df_long$Count, na.rm = TRUE))),
                        guide = "colorbar",
                        name = "Número de estudios que\nno encuentran asociación") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = " ", x = " ", y = " ")

# Combinar con patchwork
library(patchwork)
p_encontrada / p_noencontrada + plot_layout(guides = "collect")  # uno arriba (Encontrada) y otro abajo (No encontrada)

dev.off()


#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
### ÉL GRÁFICO QUE MÁS ME CONVENCE 2:

library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)


# png("bubble-plot-viremia7.png", width=1350, height=1500, res=150)
svg("bubble-plot-viremia9.svg", width=8, height=10)

# Gráfico encontrada
p_encontrada <- ggplot(df_encontrada, aes(x = Severity, y = Virus)) +
  geom_point(
    data = df_long %>% filter(Association == "Encontrada" & Count == 0),
    aes(x = Severity, y = Virus),
    shape = 3, color = "gray95", size = 2
  ) +
  geom_point(
    data = df_encontrada %>% filter(Count > 0),
    aes(x = Severity, y = Virus, size = Count, color = Count)
  ) +
  scale_size_continuous(
    range = c(2, 12),
    breaks = c(1,2,3,4,5,10,20,40,60),
    name = "Número de estudios total",
    guide = guide_legend(order = 2)  # leyenda de tamaño en segundo lugar
  ) +
  scale_color_gradientn(
    colours = c("#FFCCCC", "#FF6363", "#FF0000", "#990000"),
    values = scales::rescale(c(1,5,20,max(df_long$Count, na.rm = TRUE))),
    guide = guide_colorbar(order = 1), # leyenda de color rojo primero
    name = "Número de estudios que\nencuentran asociación"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = " ", x = " ", y = " ")

# Gráfico no encontrada
p_noencontrada <- ggplot(df_noencontrada, aes(x = Severity, y = Virus)) +
  geom_point(
    data = df_long %>% filter(Association == "No encontrada" & Count == 0),
    aes(x = Severity, y = Virus),
    shape = 3, color = "gray95", size = 2
  ) +
  geom_point(
    data = df_long_plot %>% filter(Association == "No encontrada" & Count > 0),
    aes(x = Severity, y = Virus, size = Count, color = Count)
  ) +
  scale_size_continuous(
    range = c(2, 4),
    breaks = c(1,2,3,4,5,10,20,40,60),
    guide = "none"  # ocultamos leyenda de tamaño duplicada
  ) +
  scale_color_gradientn(
    colours = c("cornflowerblue", "blue4"),
    values = scales::rescale(c(1, max(df_long$Count, na.rm = TRUE))),
    guide = guide_colorbar(order = 3), # leyenda de color azul al final
    name = "Número de estudios que\nno encuentran asociación"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = " ", x = " ", y = " ")

# Combinar gráficos con patchwork, coleccionando las leyendas
(p_encontrada / p_noencontrada) + plot_layout(guides = "collect") & 
  theme(legend.position = "right")

dev.off()
