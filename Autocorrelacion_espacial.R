setwd("F:/8.Magister_Recursos_Hidricos/1-Tesis/deteccion_Velocidad_desplazamiento/Selecciones/Raster/")
getwd()

# Librerías
library(terra)
library(spdep)
library(ggplot2)
library(patchwork)

# ---- Parámetros & archivo ----
in_path <- "G_Mosco.tif"   # <- ajusta ruta
nb_type  <- "queen"          # "queen" o "rook"
alpha    <- 0.05              # umbral LISA (dos colas)

# ---- 1) Leer raster y preparar vectores ----
r <- rast(in_path)
stopifnot(nlyr(r) == 1)
vals_all <- values(r)              # incluye NA
keep <- !is.na(vals_all)
y <- vals_all[keep]
stopifnot(length(y) > 10)

# ---- 2) Vecindad y pesos (solo celdas válidas) ----
nb_all <- cell2nb(nrow = nrow(r), ncol = ncol(r), type = nb_type)
nb <- subset.nb(nb_all, subset = keep)               # vecinos alineados con y
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# ---- 3) Moran global ----
m_global <- moran.test(y, listw = lw, zero.policy = TRUE)
cat("Moran global I =", round(m_global$estimate[1], 4),
    " p =", signif(m_global$p.value, 4), "\n\n")

# ---- 4) Moran local (LISA) ----
lisa_m <- localmoran(y, listw = lw, zero.policy = TRUE)    # cols: Ii, E.Ii, Var.Ii, Z, Pr(z>0)
Ii <- lisa_m[,1]
Z  <- lisa_m[,4]
p2 <- 2 * pnorm(-abs(Z))        # p-value (dos colas)

# centrados para clasificación
y_c  <- as.numeric(scale(y, center = TRUE, scale = FALSE))
Wy_c <- as.numeric(scale(lag.listw(lw, y, zero.policy = TRUE), center = TRUE, scale = FALSE))

# clasificar HH/LL/HL/LH/NS
cluster <- rep("NS", length(y))
sig <- p2 < alpha
cluster[sig & y_c > 0 & Wy_c > 0] <- "HH"
cluster[sig & y_c < 0 & Wy_c < 0] <- "LL"
cluster[sig & y_c > 0 & Wy_c < 0] <- "HL"
cluster[sig & y_c < 0 & Wy_c > 0] <- "LH"

# ---- 5) Preparar dataframes para graficar ----
df_var <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(df_var)[3] <- "var"
df_cl  <- df_var
df_cl$cluster <- NA_character_
df_cl$cluster[!is.na(df_var$var)] <- cluster
df_cl$cluster <- factor(df_cl$cluster, levels = c("HH","LL","HL","LH","NS"))

df_scat <- data.frame(y_c = y_c, Wy_c = Wy_c)

# ---- 6) Plots: mapa var, mapa clúster, scatterplot ----
p_map_var <- ggplot(df_var, aes(x, y, fill = var)) +
  geom_raster() + coord_equal() +
  scale_fill_viridis_c(option = "C", na.value = "grey95") +
  labs(title = "Variable (raster)", fill = "Valor") + theme_minimal()

p_map_var


p_map_cl <- ggplot(df_cl, aes(x, y, fill = cluster)) +
  geom_raster() + coord_equal() +
  scale_fill_manual(values = c(HH="#b2182b", LL="#2166ac", HL="#ef8a62", LH="#67a9cf", NS="grey90"),
                    na.value = "grey95", name = "Cluster") +
  labs(title = paste0("LISA (α=", alpha, ", dos colas)")) + theme_minimal()

p_map_cl

p_scat <- ggplot(df_scat, aes(x = y_c, y = Wy_c)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_point(alpha = 0.45) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = paste0("Moran scatterplot — I = ", round(m_global$estimate[1], 3),
                      " (p=", signif(m_global$p.value,3), ")"),
       x = "y (centrado)", y = "W·y (centrado)") + theme_minimal()

p_scat


# combinar (dos mapas arriba, scatter abajo)
p_final <- (p_map_var | p_map_cl) / p_scat + plot_annotation(title = "Moran global y LISA")

# mostrar y guardar
print(p_final)
ggsave("moran_summary.png", p_final, width = 12, height = 8, dpi = 300)

# ---- 7) Resultados resumidos (texto y tabla) ----
cat("\n--- Moran global (resumen) ---\n"); print(m_global)
cat("\n--- LISA (primeras 6 filas: Ii, Z, p2) ---\n")
print(head(data.frame(Ii = Ii, Z = Z, p2 = p2), 6))
cat("\n--- Conteo clústeres LISA ---\n"); print(table(cluster))
