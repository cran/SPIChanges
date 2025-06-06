## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----monte carlo 1, message=FALSE, warning=FALSE, eval=FALSE------------------
# library(gamlss)
# library(gamlss.dist)
# library(MuMIn)
# library(spsUtil)
# # Designing the selection function based on the AICc
# select.model <- function(rain.week) {
#   best <- matrix(NA, 1, 2)
#   colnames(best) <- c("NonLinear", "Linear")
#   t.gam <- quiet(gamlss(
#     rain.week ~ 1,
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns10 <- quiet(gamlss(
#     rain.week ~ time,
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns01 <- quiet(gamlss(
#     rain.week ~ 1,
#     sigma.formula = ~time,
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns11 <- quiet(gamlss(
#     rain.week ~ time,
#     sigma.formula = ~time,
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns20 <- quiet(gamlss(
#     rain.week ~ time + I(time^2),
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns02 <- quiet(gamlss(
#     rain.week ~ 1,
#     family = GA,
#     mu.link = "log",
#     sigma.formula = ~ time +
#       I(time^2),
#     sigma.link = "log"
#   ))
#   t.gam.ns21 <- quiet(gamlss(
#     rain.week ~ time + I(time^2),
#     sigma.formula = ~time,
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns12 <- quiet(gamlss(
#     rain.week ~ time,
#     sigma.formula = ~ time + I(time^2),
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns22 <- quiet(gamlss(
#     rain.week ~ time + I(time^2),
#     sigma.formula = ~ time + I(time^2),
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns30 <- quiet(gamlss(
#     rain.week ~ time +
#       I(time^2) +
#       I(time^3),
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns03 <- quiet(gamlss(
#     rain.week ~ 1,
#     family = GA,
#     mu.link = "log",
#     sigma.formula = ~ time +
#       I(time^2) +
#       I(time^3),
#     sigma.link = "log"
#   ))
#   t.gam.ns31 <- quiet(gamlss(
#     rain.week ~ time +
#       I(time^2) +
#       I(time^3),
#     sigma.formula = ~time,
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns13 <- quiet(gamlss(
#     rain.week ~ time,
#     sigma.formula = ~ time +
#       I(time^2) +
#       I(time^3),
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns32 <- quiet(gamlss(
#     rain.week ~ time +
#       I(time^2) +
#       I(time^3),
#     sigma.formula = ~ time +
#       I(time^2),
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns23 <- quiet(gamlss(
#     rain.week ~ time +
#       I(time^2),
#     sigma.formula = ~ time +
#       I(time^2) +
#       I(time^3),
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   t.gam.ns33 <- quiet(gamlss(
#     rain.week ~ time +
#       I(time^2) +
#       I(time^3),
#     sigma.formula = ~ time +
#       I(time^2) +
#       I(time^3),
#     family = GA,
#     mu.link = "log",
#     sigma.link = "log"
#   ))
#   # Selection among all 16 models
#   best[1, 1] <- which.min(list(
#     MuMIn::AICc(t.gam, k = k),
#     MuMIn::AICc(t.gam.ns10, k = k),
#     MuMIn::AICc(t.gam.ns01, k = k),
#     MuMIn::AICc(t.gam.ns11, k = k),
#     MuMIn::AICc(t.gam.ns20, k = k),
#     MuMIn::AICc(t.gam.ns02, k = k),
#     MuMIn::AICc(t.gam.ns21, k = k),
#     MuMIn::AICc(t.gam.ns12, k = k),
#     MuMIn::AICc(t.gam.ns22, k = k),
#     MuMIn::AICc(t.gam.ns30, k = k),
#     MuMIn::AICc(t.gam.ns03, k = k),
#     MuMIn::AICc(t.gam.ns31, k = k),
#     MuMIn::AICc(t.gam.ns13, k = k),
#     MuMIn::AICc(t.gam.ns32, k = k),
#     MuMIn::AICc(t.gam.ns23, k = k),
#     MuMIn::AICc(t.gam.ns33, k = k)
#   ))
#   # Selection among the stationary and linear models (models 1 to 4)
#   best[1, 2] <- which.min(list(
#     MuMIn::AICc(t.gam, k = k),
#     MuMIn::AICc(t.gam.ns10, k = k),
#     MuMIn::AICc(t.gam.ns01, k = k),
#     MuMIn::AICc(t.gam.ns11, k = k)
#   ))
# 
#   return(best)
# }
# library(gamlss)
# library(gamlss.dist)
# library(MuMIn)
# library(spsUtil)
# # Monte Carlo Simulation with three Sample Sizes
# sample_sizes <- c(30, 60, 90) # Define the sample sizes
# ns <- 10000 # Number of simulations
# k.seq <- seq(from = 2, to = 7, by = 0.2)
# penalty <- length(k.seq)
# 
# results_list <- list()
# selected_nonlinear_list <- list()
# selected_linear_list <- list()
# 
# for (n in sample_sizes) {
#   time <- 1L:n
#   rain.week <- matrix(NA, nrow = n, ncol = ns)
# 
#   sigma0 <- 0.22
#   mu0 <- 809
#   slope.mu <- 0 # trend-free
#   slope.sigma <- 0 # trend-free
#   sigma <- sigma0 + slope.sigma * time
#   mu <- mu0 + slope.mu * time
# 
#   # Simulating precipitation data
#   rain.week <- replicate(ns, sapply(1:n, function(t) rGA(1, mu = mu[t], sigma = sigma[t])))
# 
#   result <- matrix(NA, nrow = ns, ncol = 2 * penalty)
#   selected.nonlinear <- matrix(NA, 16, penalty)
#   selected.linear <- matrix(NA, 4, penalty)
# 
#   # Calculating loop for each k
#   for (j in seq_along(k.seq)) {
#     k <- k.seq[j]
#     result1 <- apply(rain.week, 2, select.model)
#     result[, (2 * j - 1):(2 * j)] <- t(result1)
# 
#     for (model in 1:16) {
#       selected.nonlinear[model, j] <- 100 * (sum(result[, 2 * j - 1] == model) / ns)
#       if (model <= 4) {
#         selected.linear[model, j] <- 100 * (sum(result[, 2 * j] == model) / ns)
#       }
#     }
#   }
#   rownames(selected.nonlinear) <- paste0("Model", 1:16)
#   rownames(selected.linear) <- paste0("Model", 1:4)
#   colnames(selected.nonlinear) <- as.character(k.seq)
#   colnames(selected.linear) <- as.character(k.seq)
# 
#   results_list[[as.character(n)]] <- result
#   selected_nonlinear_list[[as.character(n)]] <- selected.nonlinear
#   selected_linear_list[[as.character(n)]] <- selected.linear
# }
# 
# for (n in sample_sizes) {
#   cat("\nSample Size:", n, "\n")
#   cat("\nSelected Nonlinear Models:\n")
#   print(selected_nonlinear_list[[as.character(n)]])
#   cat("\nSelected Linear Models:\n")
#   print(selected_linear_list[[as.character(n)]])
# }

## ----monte carlo 2, message=FALSE, warning=FALSE, eval=FALSE------------------
# ## Monte Carlo Simulation
# n <- 90
# ns <- 10000
# time <- 1L:n
# 
# rain.week <- matrix(NA, nrow = n, ncol = ns)
# sigma0 <- 0.22
# mu0 <- 809
# slope.mu <- (-0.004 * mu0) # super-imposed trend
# slope.sigma <- 0 # trend-free
# sigma <- sigma0 + slope.sigma * time
# mu <- mu0 + slope.mu * time
# k.seq <- seq(from = 2, to = 7, by = 0.2)
# penalty <- length(k.seq)
# 
# result <- matrix(NA, nrow = ns, ncol = 2 * penalty)
# selected.nonlinear <- matrix(NA, 16, penalty)
# selected.linear <- matrix(NA, 4, penalty)
# # simulating again
# rain.week <- replicate(ns, sapply(1:n, function(t) rGA(1, mu = mu[t], sigma = sigma[t])))
# 
# for (j in seq_along(k.seq)) {
#   k <- k.seq[j]
#   result1 <- apply(rain.week, 2, select.model)
#   result[, (2 * j - 1):(2 * j)] <- t(result1)
# 
#   for (model in 1:16) {
#     selected.nonlinear[model, j] <- 100 * (sum(result[, 2 * j - 1] == model) / ns)
#     if (model <= 4) {
#       selected.linear[model, j] <- 100 * (sum(result[, 2 * j] == model) / ns)
#     }
#   }
# }
# 
# rownames(selected.nonlinear) <- paste0("Model", 1:16)
# rownames(selected.linear) <- paste0("Model", 1:4)
# colnames(selected.nonlinear) <- as.character(k.seq)
# colnames(selected.linear) <- as.character(k.seq)
# 
# print(selected.nonlinear)
# print(selected.linear)

## ----Case Study 1-------------------------------------------------------------
library(SPIChanges)
data(OxfordRain)
rainTS4 <- TSaggreg(daily.rain = OxfordRain, start.date = "1827-01-01", TS = 4)

### Fit All Models: Stationary, Linear and Nonlinear

Oxford.Changes <- SPIChanges(rain.at.TS = rainTS4, only.linear = "No")

### Fit Only the Stationary and Linear Models

Oxford.Changes.linear <- SPIChanges(rain.at.TS = rainTS4, only.linear = "Yes")

## ----selected models, fig.width=8, fig.height=5-------------------------------
eixo.x.Month <- Oxford.Changes$model.selection[, 1]
eixo.y.quasiWeek <- Oxford.Changes$model.selection[, 2]
valores <- Oxford.Changes$model.selection[, 3]
dados <- cbind(eixo.x.Month, eixo.y.quasiWeek, valores)
df <- as.data.frame(dados)

library(ggplot2)
df$valores <- as.factor(df$valores)
fig1 <- ggplot(df) +
  aes(x = eixo.x.Month, y = eixo.y.quasiWeek, fill = valores) +
  geom_tile() +
  scale_fill_manual(values = c(
    "1" = "#D3D3D3",
    "2" = "#ADD8E6",
    "3" = "#87CEEB",
    "4" = "#4682B4",
    "6" = "#FA8072",
    "11" = "#FF8C00"
  )) +
  theme_bw() +
  labs(x = "Months", y = "quasi-week", fill = NULL, title = "only.linear = 'No'", caption = ("Figure 1. Gamma-based models selected by the SPIChanges package. The plot corresponds to setting \n the `only.linear`  argument of the SPIChanges() function to `No`. Monthly precipitation series recorded \n at the Radcliffe Observatory site  in Oxford-UK (January 1827 to January 2020).")) +
  scale_x_continuous(breaks = 1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  geom_text(aes(label = valores), color = "black", size = 3, fontface = "bold") +
  theme(
    strip.text = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", color = "black"),
    axis.title.y = element_text(size = 10, face = "bold", color = "black"),
    plot.title = element_text(face = "bold", size = 14, color = "black"),
    plot.caption = element_text(hjust = 0, size = 10, color = "black", lineheight = 1.0),
    plot.margin = margin(10, 10, 20, 10)
  )

fig1

## ----year-to-year changes, fig.width=10, fig.height=5-------------------------
library(dplyr)
library(ggplot2)
library(patchwork)
library(purrr)

# Define constants
eixo.x.years <- 1827L:2019L

# Define conditions
conditions <- list(
  jan1 = list(Month = 1, quasiWeek = 1, col = 4),
  jan2 = list(Month = 1, quasiWeek = 2, col = 4),
  jan3 = list(Month = 1, quasiWeek = 3, col = 4),
  jan4 = list(Month = 1, quasiWeek = 4, col = 4),
  jul3 = list(Month = 7, quasiWeek = 3, col = 4),
  jul4 = list(Month = 7, quasiWeek = 4, col = 4),
  dec4mu = list(Month = 12, quasiWeek = 4, col = 4),
  jun4 = list(Month = 6, quasiWeek = 4, col = 5),
  jul1 = list(Month = 7, quasiWeek = 1, col = 5),
  jul2 = list(Month = 7, quasiWeek = 2, col = 5),
  sep3 = list(Month = 9, quasiWeek = 3, col = 5),
  dec2 = list(Month = 12, quasiWeek = 2, col = 5),
  dec4 = list(Month = 12, quasiWeek = 4, col = 5),
  mar2 = list(Month = 3, quasiWeek = 2, col = 5),
  mar3 = list(Month = 3, quasiWeek = 3, col = 5),
  sep1 = list(Month = 9, quasiWeek = 1, col = 5)
)

# Extract data using purrr::map
results <- map(names(conditions), ~ {
  cond <- conditions[[.x]]
  Oxford.Changes$Statistics %>%
    filter(Month == cond$Month, quasiWeek == cond$quasiWeek) %>%
    pull(cond$col)
}) %>% set_names(names(conditions))

# Assign results to variables
list2env(results, envir = .GlobalEnv)

# Adjust jan4
jan4 <- head(jan4, -1)

# Combine data into a single dataframe
df2 <- data.frame(
  ano = rep(eixo.x.years, length.out = length(unlist(results))),
  valores = unlist(results),
  line = rep(names(results), lengths(results)),
  letra = rep(c("a", "b", "c"), times = c(
    sum(lengths(results[1:7])),
    sum(lengths(results[8:13])),
    sum(lengths(results[14:16]))
  ))
)

# Recode `line` for better readability
df2 <- df2 %>%
  mutate(line = recode(line,
    "jan1" = "Jan (1st quasi-week)",
    "jan2" = "Jan (2nd quasi-week)",
    "jan3" = "Jan (3rd quasi-week)",
    "jan4" = "Jan (4th quasi-week)",
    "jul1" = "Jul (1st quasi-week)",
    "jul2" = "Jul (2nd quasi-week)",
    "jul3" = "Jul (3rd quasi-week)",
    "jul4" = "Jul (4th quasi-week)",
    "jun4" = "Jun (4th quasi-week)",
    "sep1" = "Sep (1st quasi-week)",
    "sep3" = "Sep (3rd quasi-week)",
    "dec2" = "Dec (2nd quasi-week)",
    "dec4" = "Dec (4th quasi-week)",
    "mar2" = "Mar (2nd quasi-week)",
    "mar3" = "Mar (3rd quasi-week)"
  ))

# Define colors
cores <- c(
  "Jan (1st quasi-week)" = "blue",
  "Jan (2nd quasi-week)" = "gold",
  "Jan (3rd quasi-week)" = "firebrick",
  "Jan (4th quasi-week)" = "black",
  "Jul (3rd quasi-week)" = "cyan",
  "Jul (4th quasi-week)" = "red",
  "Dec (4th quasi-week)" = "#1f77b4",
  "Jun (4th quasi-week)" = "#7f7f7f",
  "Jul (1st quasi-week)" = "#008080",
  "Jul (2nd quasi-week)" = "#8B4513",
  "Sep (3rd quasi-week)" = "#4B0082",
  "Dec (2nd quasi-week)" = "#ffbb78",
  "Mar (2nd quasi-week)" = "#9467bd",
  "Mar (3rd quasi-week)" = "#bcbd22",
  "Sep (1st quasi-week)" = "#c49c94"
)

# Function to create plots (updated with legend.position.inside)
create_plot <- function(data, y_limits, y_breaks, y_lab, title) {
  ggplot(data) +
    aes(x = ano, y = valores, colour = line) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = cores) +
    theme_bw() +
    theme(
      legend.position.inside = c(0.9, 0.3),  # Updated here
      legend.title = element_blank(),
      legend.key.size = unit(0.1, "cm"),
      legend.box.spacing = unit(0.1, "cm"),
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.background = element_blank(),
      legend.text = element_text(size = 9),
      legend.key.width = unit(1, "cm"),
      axis.title.y = element_text(size = 9, face = "bold"),
      strip.text = element_text(face = "bold", size = 12),
      plot.margin = margin(3, 3, 3, 3),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    scale_x_continuous(breaks = seq(1827, 2019, by = 25)) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +
    ylab(y_lab) +
    ggtitle(title) +
    guides(colour = guide_legend(ncol = 2))
}

# Create plots
plot_a <- create_plot(df2 %>% filter(letra == "a"), c(0, 100), seq(0, 100, by = 20), "mean parameter \n (μ) GAMLSS", "A")
plot_b <- create_plot(df2 %>% filter(letra == "b"), c(0, 1), seq(0, 1, by = 0.2), "dispersion parameter \n (δ) GAMLSS", "B")
plot_c <- create_plot(df2 %>% filter(letra == "c"), c(0, 2), seq(0, 2, by = 0.4), "dispersion parameter \n (δ) GAMLSS", "C") +
  labs(caption  = "Figure 2. Year-to-year changes in the dispersion parameter of gamma-based models estimated from
       monthly precipitation records of the Radcliffe Observatory site in Oxford-UK.
       (A) linear changes in the μ parameter
       (B) linear changes in the δ parameter
       (C) nonlinear changes in the δ parameter") +
  theme(axis.title.x = element_text(size = 10, face = "bold"), plot.caption = element_text(hjust = 0, size = 10, color = "black", lineheight = 1.0))

# Combine plots
fig2 <- plot_a + plot_b + plot_c + plot_layout(ncol = 1, heights = c(1.3, 1.3, 1))
fig2

## ----changes in drought frequency, fig.width=8, fig.height=7------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
Oxford.Changes$Changes.Freq.Drought <- as.data.frame(Oxford.Changes$Changes.Freq.Drought)
# Define quasi-weeks and labels
quasi_weeks <- list(
  "1st - Jan" = c(1,1), "2nd - Jan" = c(1,2), "3rd - Jan" = c(1,3), "4th - Jan" = c(1,4),
  "2nd - Mar" = c(3,2), "3rd - Mar" = c(3,3), "4th - Jun" = c(6,4),
  "1st - Jul" = c(7,1), "2nd - Jul" = c(7,2), "3rd - Jul" = c(7,3), "4th - Jul" = c(7,4),
  "1st - Sep" = c(9,1), "3rd - Sep" = c(9,3), "2nd - Dec" = c(12,2), "3rd - Dec" = c(12,4)
)

# Function to extract data
drought_data <- function(index) {
  sapply(quasi_weeks, function(qw) {
    Oxford.Changes$Changes.Freq.Drought %>% 
      filter(Month == qw[1], quasiWeek == qw[2]) %>%
      pull(index)
  })
}

# Create data frame
df3 <- data.frame(
  x = names(quasi_weeks),
  moderate = drought_data(7),
  severe = drought_data(8),
  extreme = drought_data(9),
  stat = drought_data(5),
  nonstat = drought_data(6)
) %>% pivot_longer(-x, names_to = "category", values_to = "value")

# Define colors
color_map <- c(
  "moderate" = "#00243F", "severe" = "#517C9D", "extreme" = "#36A7C1", 
  "stat" = "#1c340a", "nonstat" = "#3ba551"
)

# Generate plot function
generate_plot <- function(data, y_label, title, limits, colors, legend_labels) {
  ggplot(data, aes(x = x, y = value, fill = category)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = colors, labels = legend_labels) +
    scale_y_continuous(limits = limits) +
    labs(title = title, x = NULL, y = y_label) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 10, face = "bold", color = "black"),
      axis.title.y = element_text(size = 9, face = "bold", color = "black"),
      legend.position = "top", legend.title = element_blank(),
      legend.text = element_text(size = 10, color = "black"),
      legend.key.width = unit(1, "cm"), legend.key.size = unit(0.4, "cm"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Create plots
plot_3a <- generate_plot(
  df3 %>% filter(category %in% c("moderate", "severe", "extreme")),
  "Percentual change (%)", "A", c(-20, 20), 
  color_map[c("moderate", "severe", "extreme")],
  c("Change Moderate", "Change Severe", "Change Extreme")
)

plot_3b <- generate_plot(
  df3 %>% filter(category %in% c("stat", "nonstat")),
  "mm", "B", c(0, 100), 
  color_map[c("stat", "nonstat")],
  c("Stationary Normal Rain", "NonStationary Normal Rain")
) +
  labs(caption = "Figure 3. Year-to-year changes in drought frequency estimated from 
       quasi-weekly precipitation records of the Radcliffe Observatory site in Oxford-UK. 
       (A) Percentual changes in the frequency of moderate, severe, and extreme droughts.
       (B) Changes in the normal precipitation amount for stationary and nonstationary normality assumptions.") +
  theme(axis.title.x = element_text(size = 10, face = "bold"), plot.caption = element_text(hjust = 0, size = 10, color = "black", lineheight = 1.0))

# Combine plots
Fig3 <- plot_3a + plot_3b + plot_layout(ncol = 1, heights = c(1.5, 1.2))
Fig3

## ----Monitoring drought in Oxford-UK, 2019------------------------------------
Table1 <- Oxford.Changes$data.week[which(Oxford.Changes$data.week[, 1] == 2019), ]
# Table 1. Output data field /$data.week obtained by applying the SPIChanges package at the 4-quasi week time scale to the precipitation series recorded at the Radcliffe Observatory site in Oxford-UK (January 1827 to January 2020): rain.at.TS is the precipitation amount accumulated at the time scale specified by the users, SPI is the Standardized Precipitation Index, Exp.Acum.Prob and Actual.Acum.Prob are the expected frequency of the SPI estimates calculated under the stationary and under non-stationary approaches, respectively and, ChangeFreq is the changes in the frequency of dry events (SPI < 0). It is the difference between Actual.Acum.Prob and Exp.Acum.Prob.
Table1

## ----selected linear models, fig.width=8, fig.height=5------------------------
eixo.x.Month <- Oxford.Changes.linear$model.selection[, 1]
eixo.y.quasiWeek <- Oxford.Changes.linear$model.selection[, 2]
valores <- Oxford.Changes.linear$model.selection[, 3]
dados <- cbind(eixo.x.Month, eixo.y.quasiWeek, valores)
df <- as.data.frame(dados)

library(ggplot2)
df$valores <- as.factor(df$valores)
fig4 <- ggplot(df) +
  aes(x = eixo.x.Month, y = eixo.y.quasiWeek, fill = valores) +
  geom_tile() +
  scale_fill_manual(values = c(
    "1" = "#D3D3D3",
    "2" = "#ADD8E6",
    "3" = "#87CEEB",
    "4" = "#4682B4"
  )) +
  theme_bw() +
  labs(x = "Months", y = "quasi-week", fill = NULL, title = "only.linear = 'Yes'", caption = "Figure 4. Gamma-based models selected by the SPIChanges package. The plot corresponds to \n setting the `only.linear`  argument of the SPIChanges() function to `Yes`. Monthly precipitation series \n recorded at the  Radcliffe Observatory site  in Oxford-UK (January 1827 to January 2020).") +
  scale_x_continuous(breaks = 1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  geom_text(aes(label = valores), color = "black", size = 3, fontface = "bold") +
  theme(
    strip.text = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", color = "black"),
    axis.title.y = element_text(size = 10, face = "bold", color = "black"),
    plot.title = element_text(face = "bold", size = 14, color = "black"),
    plot.caption = element_text(hjust = 0, size = 10, color = "black", lineheight = 1.0),
    plot.margin = margin(10, 10, 20, 10)
  )

fig4

## ----Monitoring drought in Oxford-UK, 2019, only linear models----------------
Table2 <- Oxford.Changes.linear$data.week[which(Oxford.Changes.linear$data.week[, 1] == 2019), ]
# Table 2. Output data field /$data.week obtained by applying the SPIChanges package at the 4-quasi week time scale to the precipitation series recorded at the Radcliffe Observatory site in Oxford-UK (January 1827 to January 2020): rain.at.TS is the precipitation amount accumulated at the time scale specified by the users, SPI is the Standardized Precipitation Index, Exp.Acum.Prob and Actual.Acum.Prob are the expected frequency of the SPI estimates calculated under the stationary and under non-stationary approaches, respectively and, ChangeFreq is the changes in the frequency of dry events (SPI < 0). It is the difference between Actual.Acum.Prob and Exp.Acum.Prob.
Table2

## ----get-CPC-rain, eval=FALSE-------------------------------------------------
# # this chunk may take a long time to get daily precipitation data for entire Brazil (1980-2024)
# # Load necessary libraries
# library(ncdf4)
# library(curl)
# library(SPIChanges)
# 
# n <- nrow(lonlat)
# # Configure curl with a longer timeout
# h <- new_handle()
# handle_setopt(h, timeout = 10800L)
# 
# download_and_process <- function(year, lonlat) {
#   # Define file paths and URLs
#   rain_file <- file.path(tempdir(), paste0("cpcdata", year, ".nc"))
#   download_url <- paste0("https://psl.noaa.gov/thredds/fileServer/Datasets/cpc_global_precip/precip.", year, ".nc")
#   # Download the files
#   curl::curl_download(
#     url = download_url,
#     destfile = rain_file,
#     mode = "wb",
#     quiet = FALSE,
#     handle = h
#   )
# 
#   # Open and process the netCDF file
#   cep.prec <- nc_open(rain_file)
#   lons <- ncvar_get(cep.prec, "lon")
#   lats <- ncvar_get(cep.prec, "lat")
#   prate <- ncvar_get(cep.prec, "precip") # mm/day
#   NN <- dim(prate)[3]
#   # Flip latitude and precipitation order
#   prate <- prate[, rev(seq_len(dim(prate)[2])), 1:NN]
#   lats <- rev(lats)
#   # Initialize matrix for the year
#   pdata <- matrix(NA, NN, n)
#   for (i in seq_len(n)) {
#     longitude <- lonlat[i, 1] + 360 # range 0, 360
#     latitude <- lonlat[i, 2] # range -90, 90
#     pdata[, i] <- prate[
#       which.min((lons - longitude)^2),
#       which.min((lats - latitude)^2),
#     ]
#   }
# 
#   nc_close(cep.prec)
#   return(pdata)
# }
# 
# # Process data for multiple years
# years <- 1980:2024
# pdata_list <- lapply(years, download_and_process, lonlat = lonlat)
# 
# # Combine all years' data
# pdata1 <- do.call(rbind, pdata_list)

## ----aggregate-CPC-rain, eval=FALSE-------------------------------------------
# library(SPIChanges)
# # This chunk takes a long time to aggregate precipitation data for entire Brazil (1980-2024)
# pdata1[is.na(pdata1)] <- 0 # replacing NA by zero.
# TS12 <- TSaggreg(daily.rain = pdata1[, 1], start.date = "1980-01-01", TS = 12)
# N <- nrow(TS12)
# rain.at.TS <- matrix(NA, n, 5)
# rainTS12.bank <- matrix(NA, N, (n + 3))
# rainTS12.bank[, 1:4] <- as.matrix(TS12)
# for (i in 2:n) {
#   TS12 <- TSaggreg(daily.rain = pdata1[, i], start.date = "1980-01-01", TS = 12)
#   rainTS12.bank[, (i + 3)] <- TS12[, 4]
# }

## ----SPIChanges for entire Brazil, eval=FALSE---------------------------------
# # trying to speed up the things
# library(foreach)
# library(doParallel)
# library(SPIChanges)
# numCores <- detectCores() - 1 # Use one less core than the total available
# cl <- makeCluster(numCores)
# registerDoParallel(cl)
# rain <- matrix(NA, N, 4)
# rain[, 1:3] <- as.matrix(rainTS12.bank[, 1:3])
# final <- ncol(rainTS12.bank)
# n <- nrow(lonlat)
# Map <- matrix(NA, n, 18)
# model <- matrix(NA, n, 3)
# Map[, 1:2] <- as.matrix(lonlat)
# model[, 1:2] <- as.matrix(lonlat)
# changes <- matrix(NA, 1, 16)
# # Perform parallel processing
# results <- foreach(i = 4:final, .combine = rbind, .packages = "SPIChanges") %dopar% {
#   rain[, 4] <- rainTS12.bank[, i]
#   Local.Changes <- SPIChanges(rain.at.TS = rain, only.linear = "no")
#   changes[1, 1:3] <- Local.Changes$Changes.Freq.Drought[which(Local.Changes$Changes.Freq.Drought[, 1] == 2 &
#     Local.Changes$Changes.Freq.Drought[, 2] == 4), 7:9]
#   changes[1, 4] <- Local.Changes$model.selection[which(Local.Changes$model.selection[, 1] == 2 &
#     Local.Changes$model.selection[, 2] == 4), 3]
#   changes[1, 5:7] <- Local.Changes$Changes.Freq.Drought[which(Local.Changes$Changes.Freq.Drought[, 1] == 5 &
#     Local.Changes$Changes.Freq.Drought[, 2] == 4), 7:9]
#   changes[1, 8] <- Local.Changes$model.selection[which(Local.Changes$model.selection[, 1] == 5 &
#     Local.Changes$model.selection[, 2] == 4), 3]
#   changes[1, 9:11] <- Local.Changes$Changes.Freq.Drought[which(Local.Changes$Changes.Freq.Drought[, 1] == 8 &
#     Local.Changes$Changes.Freq.Drought[, 2] == 4), 7:9]
#   changes[1, 12] <- Local.Changes$model.selection[which(Local.Changes$model.selection[, 1] == 8 &
#     Local.Changes$model.selection[, 2] == 4), 3]
#   changes[1, 13:15] <- Local.Changes$Changes.Freq.Drought[which(Local.Changes$Changes.Freq.Drought[, 1] == 11 &
#     Local.Changes$Changes.Freq.Drought[, 2] == 4), 7:9]
#   changes[1, 16] <- Local.Changes$model.selection[which(Local.Changes$model.selection[, 1] == 11 &
#     Local.Changes$model.selection[, 2] == 4), 3]
#   return(changes)
# }
# # Combine results into the Map matrix
# a <- 1
# for (i in 1:nrow(results)) {
#   Map[a, 3:18] <- results[i, ]
#   a <- a + 1
# }
# 
# # Assign column names and save the results
# colnames(Map) <- c(
#   "lon", "lat", "SummerModerate", "SummerSevere", "SummerExtreme", "SummerModel",
#   "AutomnModerate", "AutomnSevere", "AutomnExtreme", "AutomnModel",
#   "WinterModerate", "WinterSevere", "WinterExtreme", "WinterModel",
#   "SpringModerate", "SpringSevere", "SpringExtreme", "SpringModel"
# )
# # Stop the cluster
# stopCluster(cl)
# head(Map)

## ----map, fig.width=8, fig.height=5-------------------------------------------
library(ggplot2)
library(sf)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(archive)
library(SPIChanges)
# Loading data
rar_url <- "https://github.com/gabrielblain/Grid_Brazil/raw/main/regioes_2010.rar"
temp_file <- tempfile(fileext = ".rar")
temp_dir <- tempdir()
download.file(rar_url, temp_file, mode = "wb")
archive_extract(temp_file, dir = temp_dir)
shp_path <- file.path(temp_dir, "regioes_2010.shp") # Adjust if files are extracted into a subdirectory
limitere <- st_read(shp_path)
# Tidying the data
df <- Map %>%
  select(lat, lon, SpringSevere, AutomnSevere, WinterSevere, SummerSevere) %>%
  pivot_longer(
    cols = c(SpringSevere, AutomnSevere, WinterSevere, SummerSevere),
    names_to = "Station", values_to = "Change"
  )
df_ok <- df %>%
  mutate(
    Station = recode(Station,
      "SpringSevere" = "Spring - Severe",
      "SummerSevere" = "Summer - Severe",
      "AutomnSevere" = "Autumm - Severe",
      "WinterSevere" = "Winter - Severe"
    )
  )
# Turning into sf
limitere_sf <- st_as_sf(limitere)
dados.sf <- st_as_sf(df_ok, coords = c("lon", "lat"), crs = 4326)
# Fixing possible errors
limitere_sf <- st_make_valid(limitere_sf)
# Finding centroids
limitere_centroids <- st_centroid(limitere_sf)
# Adding siglas (annacron)
limitere_centroids$sigla <- limitere$sigla
# Finding coordinates
coords <- st_coordinates(limitere_centroids)
limitere_centroids$lon <- coords[, 1]
limitere_centroids$lat <- coords[, 2]
dados.sf$Station <- factor(
  dados.sf$Station,
  levels = c(
    "Spring - Severe",
    "Summer - Severe",
    "Autumm - Severe",
    "Winter - Severe"
  )
)
cores_roxo_branco_vermelho <- colorRampPalette(c("purple", "white", "red"))(n = 100)
cores5 <- colorRampPalette(c("purple", "white", "red"))(7)
# The Map
map <- ggplot() +
  geom_sf(data = dados.sf, aes(color = Change), shape = 15, size = 3) +
  geom_sf(data = limitere_sf, fill = NA, color = "black", size = 1) +
  scale_colour_gradient2(
    low = "#6B238E", mid = "white", high = "red",
    midpoint = 0, limits = c(-30, 80)
  ) +
  theme_light() +
  labs(color = "Change (%)") +
  facet_wrap(vars(Station)) +
  theme(
    title = element_text(size = 12, face = "bold", color = "black"),
    text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.height = unit(.6, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 11),
    axis.title.x = element_blank(), # Remove o rotulo do eixo x
    axis.title.y = element_blank(), # Remove o rotulo do eixo y
    strip.text = element_text(size = 11, face = "bold", color = "black")
  )
map <- map +
  geom_text(
    data = limitere_centroids, aes(x = lon, y = lat, label = sigla),
    color = "black", size = 3,
    fontface = "bold",
    position = position_nudge(y = 0.01), show.legend = FALSE
  )

print(map)

