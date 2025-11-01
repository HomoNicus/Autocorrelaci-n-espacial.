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



# Simulaciones ------------------------------------------------------------



# === SIM_MIN_MOSCO: simulaciones breves usando tus datos reales (y, keep, r) ===
suppressPackageStartupMessages({ library(spdep); library(ggplot2); library(dplyr) })
set.seed(2025)

# Utilidad: Moran I con permutaciones -> I, p, IC (compatibilidad)
moran_perm_ci <- function(y, lw, nsim = 499) {
  mt <- spdep::moran.test(y, listw = lw, zero.policy = TRUE)
  mc <- spdep::moran.mc(y,  listw = lw, nsim = nsim, zero.policy = TRUE)
  I  <- unname(mt$estimate[["Moran I statistic"]])
  ci <- quantile(mc$res, c(0.025, 0.975), na.rm = TRUE)
  tibble(I = I, p = mc$p.value, ci_lo = ci[1], ci_hi = ci[2])
}

# --- Vecindarios sobre TU raster (solo celdas válidas) ---
nr <- nrow(r); nc <- ncol(r)
# queen/rook en la grilla completa y luego SUBSET por celdas válidas
nb_q_all  <- spdep::cell2nb(nr, nc, type = "queen")
nb_rk_all <- spdep::cell2nb(nr, nc, type = "rook")
nb_q  <- spdep::subset.nb(nb_q_all,  subset = keep)  # keep es lógico (TRUE = válida)
nb_rk <- spdep::subset.nb(nb_rk_all, subset = keep)

# kNN en coordenadas de celdas válidas
df_xy <- as.data.frame(r, xy = TRUE, na.rm = FALSE)[keep, c("x","y")]
xy_v  <- as.matrix(df_xy)
nb_k6 <- spdep::knn2nb(spdep::knearneigh(xy_v, k = 6))

lw_q  <- spdep::nb2listw(nb_q,  style = "W", zero.policy = TRUE)
lw_rk <- spdep::nb2listw(nb_rk, style = "W", zero.policy = TRUE)
lw_k6 <- spdep::nb2listw(nb_k6, style = "W", zero.policy = TRUE)

# --- (A) Sensibilidad a W en tus datos: queen vs rook vs kNN-6 ---
sensW_mosco <- bind_rows(
  moran_perm_ci(y, lw_q)  |> mutate(W = "queen"),
  moran_perm_ci(y, lw_rk) |> mutate(W = "rook"),
  moran_perm_ci(y, lw_k6) |> mutate(W = "kNN-6")
) |> mutate(p_FDR = p.adjust(p, method = "BH"))

readr::write_csv(sensW_mosco, "mosco_sensW.csv")
print(sensW_mosco)

# --- (B) Efecto del tamaño muestral en tus datos (submuestreo de celdas válidas) ---
n_valid <- length(y)
sizes0  <- c(10, 20, n_valid)    # se debe cambiar con los valores de n a evaluar
sizes   <- sort(unique(pmin(sizes0, n_valid)))
#sizes   <- sizes[sizes >= 30]  # asegura un mínimo razonable de n

make_mask <- function(idx, n_total) { m <- rep(FALSE, n_total); m[idx] <- TRUE; m }

S2_mosco <- lapply(sizes, function(n){
  idx  <- sample.int(n_valid, n)                 # índices SOBRE celdas válidas
  mask <- make_mask(idx, n_valid)                # lógico para subset.nb
  nb_sub <- spdep::subset.nb(nb_q, subset = mask)  # nb_q ya está en válidas
  lw_sub <- spdep::nb2listw(nb_sub, style = "W", zero.policy = TRUE)
  moran_perm_ci(y[idx], lw_sub) |> mutate(n = n)
}) |> bind_rows()

readr::write_csv(S2_mosco, "mosco_tamano_n.csv")
print(S2_mosco)

# --- Figura única: p (perm) vs n en tus datos ---
p_n_mosco <- ggplot(S2_mosco, aes(n, p)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = sizes) +
  labs(title = "Caso MOSCO: efecto del tamaño muestral",
       x = "n (celdas válidas)", y = "p-valor (permutaciones)") +
  theme_minimal(base_size = 11)
ggsave("mosco_p_vs_n.png", p_n_mosco, width = 4.6, height = 3.4, dpi = 300)

cat("\n[SIM_MIN_MOSCO] OK -> Figura: mosco_p_vs_n.png | Tablas: mosco_sensW.csv, mosco_tamano_n.csv\n")
# === FIN SIM_MIN_MOSCO ===



