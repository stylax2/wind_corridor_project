# ==============================================================================
# Statistical Validation of the Least-cost Corridor Model
# Description: 
#   1. Spatial feature extraction (LCC, Elevation, Slope, TPI) using exactextractr
#   2. Panel Mixed-effects Model (LMM) construction to control for spatiotemporal autocorrelation
#   3. Permutation test (n=100) and publication-ready visualization
# ==============================================================================

# 1. Environment Setup & Package Loading
# ------------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(terra, sf, dplyr, lme4, exactextractr, ggplot2)

# Set seed for reproducibility (crucial for permutation tests)
set.seed(1234)

# Define global parameters
n_perm <- 100 # Number of permutations for the validation test

# Define file paths (Adjust to local environment)
path_lcc     <- "02Output/Final_WindCorridor_Heatmap.tif"
path_domains <- "02Output/anisotropic_domains_final.gpkg"
path_dem     <- "01data/inside_150/dem.tif"
path_slope   <- "01data/inside_150/slope.tif"
path_tpi     <- "01data/inside_150/tpi.tif"
path_wind    <- "01data/wind.rds"

# 2. Data Loading and Coordinate Reference System (CRS) Alignment
# ------------------------------------------------------------------------------
message("Step 1: Loading raster and vector datasets...")
lcc     <- rast(path_lcc)
r_dem   <- rast(path_dem)
r_slope <- rast(path_slope)
r_tpi   <- rast(path_tpi)

# Load spatial domains and align CRS with the raster data
domains <- st_read(path_domains, quiet = TRUE) %>% st_transform(crs(lcc))

# Load daily meteorological time-series data (409,180 obs)
df <- readRDS(path_wind) %>% 
  mutate(site = as.character(site), Date = as.Date(Date))

# 3. Spatial Feature Extraction (Zonal Statistics)
# ------------------------------------------------------------------------------
message("Step 2: Extracting spatial statistics for each domain...")

# (1) Extract 90th percentile of LCC density (Core wind corridor index)
domains$lcc_p90 <- exact_extract(lcc, domains, function(x, f, ...) {
  val <- x[is.finite(x)]
  if(length(val) == 0) return(NA)
  quantile(val, 0.9, na.rm = TRUE, names = FALSE)
})

# (2) Extract topographic covariates (Elevation, Slope, TPI)
domains$topo_elev_mean <- exact_extract(r_dem, domains, 'mean')
domains$topo_slope_sd  <- exact_extract(r_slope, domains, 'stdev')
domains$topo_tpi_range <- exact_extract(r_tpi, domains, function(x, f, ...) {
  val <- x[is.finite(x)]
  if(length(val) == 0) return(NA)
  max(val, na.rm = TRUE) - min(val, na.rm = TRUE)
})

# 4. Panel Data Construction & Standardization
# ------------------------------------------------------------------------------
# Drop geometry and merge spatial covariates with meteorological time-series
site_cov <- domains %>% 
  st_drop_geometry() %>% 
  select(site, lcc_p90, topo_elev_mean, topo_slope_sd, topo_tpi_range)

df_panel <- df %>%
  left_join(site_cov, by = "site") %>%
  filter(is.finite(ws), !is.na(lcc_p90)) %>%
  # Z-score standardization for model convergence and relative effect comparison
  mutate(across(c(ws, lcc_p90, topo_elev_mean, topo_slope_sd, topo_tpi_range), 
                ~ as.numeric(scale(.x))))

# 5. Panel Mixed-effects Model (LMM) Analysis
# ------------------------------------------------------------------------------
message("Step 3: Fitting Panel Mixed-effects Models...")

# Base Model: Topographic covariates only
form_base <- ws ~ topo_elev_mean + topo_slope_sd + topo_tpi_range + (1|site) + (1|Date)
model_base <- lmer(form_base, data = df_panel, REML = FALSE, 
                   control = lmerControl(optimizer = "bobyqa"))

# Full Model: Topography + LCC density index
form_full <- ws ~ lcc_p90 + topo_elev_mean + topo_slope_sd + topo_tpi_range + (1|site) + (1|Date)
model_full <- lmer(form_full, data = df_panel, REML = FALSE, 
                   control = lmerControl(optimizer = "bobyqa"))

# Output statistical results
cat("\n--- Likelihood Ratio Test (Model Comparison) ---\n")
print(anova(model_base, model_full))

cat("\n--- Fixed Effects Coefficients (Full Model) ---\n")
print(summary(model_full)$coefficients)

# 6. Spatial Robustness Validation (Permutation Test)
# ------------------------------------------------------------------------------
message(sprintf("Step 4: Running permutation test (n=%d)...", n_perm))
real_beta <- fixef(model_full)["lcc_p90"]
perm_betas <- numeric(n_perm)

for(i in 1:n_perm) {
  # Spatially shuffle the LCC index across sites
  site_perm <- site_cov %>% mutate(lcc_p90 = sample(lcc_p90))
  
  df_temp <- df %>% left_join(site_perm, by = "site") %>% 
    mutate(across(where(is.numeric), ~as.numeric(scale(.x))))
  
  m_perm <- lmer(form_full, data = df_temp, REML = FALSE, 
                 control = lmerControl(optimizer = "bobyqa"))
  perm_betas[i] <- fixef(m_perm)["lcc_p90"]
}

# Calculate permutation p-value
perm_p <- mean(abs(perm_betas) >= abs(real_beta))
cat(sprintf("\n✅ Permutation P-value: %.4f\n", perm_p))

# 7. Publication-Ready Visualization (Figure Generation)
# ------------------------------------------------------------------------------
message("Step 5: Generating and saving permutation histogram...")

df_perm <- data.frame(beta = perm_betas)

p_perm <- ggplot(df_perm, aes(x = beta)) +
  geom_histogram(fill = "gray70", color = "black", bins = 20, alpha = 0.8) +
  geom_vline(xintercept = real_beta, color = "red", linetype = "dashed", linewidth = 1.2) +
  annotate("text", x = real_beta - 0.01, y = max(table(cut(df_perm$beta, 20))) * 0.8, 
           label = paste0("Observed Effect\n(\u03b2 = ", round(real_beta, 3), ")"), 
           color = "red", fontface = "bold", hjust = 1, size = 5) +
  labs(
    x = expression(paste("Permuted Fixed Effect Coefficients (", beta, ") for LCC Index")),
    y = "Frequency"
  ) +
  theme_classic(base_size = 14) +
  theme(
    # Titles are typically removed for journals (handled via Figure Captions in manuscript)
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black", size = 12),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

# Display and save the high-resolution figure
print(p_perm)
ggsave("Figure_Permutation_Test.png", plot = p_perm, width = 8, height = 5, dpi = 300)
message("Validation complete. Figure saved successfully.")
