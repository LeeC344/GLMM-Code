
#### Packages ####

library(dplyr)
library(readr)
library(lubridate)
library(janitor)
library(glmmTMB)
library(gt)
library(performance)
library(MuMIn)
library(tibble)
library(DHARMa)
library(tidyverse)
library(lme4)

rm(list=ls())

#Load in Data
combined_detections<- readRDS("combined_actel_detections.rds")

env_model<- read.csv("env_model_daily.csv") %>% 
  mutate(date = as.Date(date))


#Filter Detections to Sal/Temp Logger deployments
GLMMDets <- combined_detections %>%
  mutate(
    date = as.Date(Timestamp)
  ) %>%
  filter(
    date >= as.Date("2025-06-12"),
    date <= as.Date("2026-04-16")
  )

#attaching associated RKMs 
rec_km <- tibble::tribble(
  ~Array, ~rkm,
  "LA10", 0,
  "LA9", 14.4,
  "LA8", 31.5,
  "LA7", 55.5,
  "LA6", 64.2,
  "LA5", 70.9,
  "LA4", 94.8,
  "LA3", 105.6,
  "LA2", 112.5,
  "LA1", 116.5,
  "CL8", 120.1,
  "CL7", 130.1,
  "CL6", 144.5,
  "CL5", 159,
  "CL4", 169.1,
  "CL3", 188,
  "CL2", 199.5,
  "CL1", 212.3,
  "MBL1", -9.6,
  "MBL2", -23.7,
  "MBL3", -39.9,
  "MBL4", -58.6,
  "MBL5", -63,
  "TNSW1", -16.7,
  "MDL1", -27,
  "TNSW4", -67.9,
  "BMB1", -77.8,
  "TNSW5", -69.7,
  "BSN1", -43.2,
  "RFT1", -70.8,
  "BBC1", -50.5,
  "TNSW2", -31.3,
  "TNSW3", -32.9,
  "TB1", 24.3,
  "TB2", 4.3,
  "MFLN1", -34.5
)


GLMMDets <- GLMMDets %>%
  mutate(
    Array=as.character(Array)
  ) %>%
  left_join(rec_km, by = "Array")



# Create fish-day detection dataset

fish_day <- GLMMDets %>%
  group_by(fish_id, date) %>%
  summarise(
    delta_use = as.integer(any(Section == "Delta")),
    n_detections = n(),
    mean_rkm = mean(rkm, na.rm = TRUE),
    min_rkm = min(rkm, na.rm = TRUE),
    max_rkm = max(rkm, na.rm = TRUE),
    .groups = "drop"
  )


fish_day_model <- fish_day %>%
  inner_join(env_model, by = "date")


fish_day_model %>%
  summarise(
    n_fish = n_distinct(fish_id),
    n_fish_days = n(),
    n_delta_days = sum(delta_use == 1),
    n_non_delta_days = sum(delta_use == 0),
    first_date = min(date),
    last_date = max(date)
  )


#-----------------------------------------------------
#Exploring this fish_day_model data frame
#-----------------------------------------------------
length(unique(fish_day_model$fish_id)) #130 unique fish
length(unique(fish_day_model$date)) #284 unique days

tmp = fish_day_model %>% 
  group_by(fish_id, date) %>% 
  summarise(count=n()) %>% 
  ungroup()

#Ok, so it appears as though there is only one observation per fish per day

#--------------------------------------------------------------------------
#Jonathon's glmm
#--------------------------------------------------------------------------

mGlobal = glmer(delta_use ~ sal_logger1_z + sal_logger3_z +
                  temp_mean_z + stage_ft_z + (1|fish_id),
                family=binomial, data=fish_day_model)



# 1. Run simulations on your fitted glmer model
simRes <- simulateResiduals(fittedModel = mGlobal)

# 2. Plot the main diagnostic panel (QQ plot and residuals vs. predicted)
plot(simRes)


# Test for temporal autocorrelation
testTemporalAutocorrelation(simulationOutput = simRes, time = date)




############### GLMM ######################


m_global <- glmmTMB(
  delta_use ~ sal_logger1_z + sal_logger3_z +
    temp_mean_z+
    stage_ft_z +
    (1 | fish_id),
  family = binomial,
  data = fish_day_model,
  na.action = na.fail
)

summary(m_global)

#Multi-Collinearity
check_collinearity(m_global)


#Autocorrelation
sim_res <- simulateResiduals(m_global, n = 1000)

res_df <- fish_day_model %>%
  arrange(fish_id, date) %>%
  mutate(dharma_resid = residuals(sim_res))

acf_by_fish <- res_df %>%
  group_by(fish_id) %>%
  arrange(date, .by_group = TRUE) %>%
  filter(n() >= 10) %>%
  summarise(
    n_days = n(),
    lag1_acf = acf(dharma_resid, plot = FALSE)$acf[2],
    .groups = "drop"
  )

acf_by_fish %>%
  summarise(
    n_fish_tested = n(),
    mean_lag1_acf = mean(lag1_acf, na.rm = TRUE),
    median_lag1_acf = median(lag1_acf, na.rm = TRUE),
    max_abs_lag1_acf = max(abs(lag1_acf), na.rm = TRUE)
  )

#### DREDGE ####

dredge_results<- dredge(m_global)

dredge_results

dredge_table<- as.data.frame(dredge_results) %>%
  tibble::rownames_to_column("Model1")

dredge_table

dredge_table %>%
  gt() %>%
  tab_header(
    title = "Model-selection results for Delta use"
  ) %>%
  fmt_number(
    columns = c(logLik,AICc,delta,weight),
    decimals =3
    
  )

Valid_AIC<-subset(dredge_results, delta <2)


# Convert dredge results to clean data frame


clean_aic_table <- as.data.frame(Valid_AIC) %>%
  rownames_to_column("Model") %>%
  rename(
    `Upstream Salinity Logger` = `cond(sal_logger1_z)`,
    `Downstream Salinity Logger` = `cond(sal_logger3_z)`,
    `Gage Height` = `cond(stage_ft_z)`,
    `Mean Temperature` = `cond(temp_mean_z)`,
    K = df,
    `Delta AICc` = delta,
    `AICc Weight` = weight
  ) %>%
  mutate(
    `Upstream Salinity Logger` = ifelse(is.na(`Upstream Salinity Logger`), "", "Yes"),
    `Downstream Salinity Logger` = ifelse(is.na(`Downstream Salinity Logger`), "", "Yes"),
    `Gage Height` = ifelse(is.na(`Gage Height`), "", "Yes"),
    `Mean Temperature` = ifelse(is.na(`Mean Temperature`), "", "Yes")
  ) %>%
  select(
    Model,
    `Upstream Salinity Logger`,
    `Downstream Salinity Logger`,
    `Gage Height`,
    `Mean Temperature`,
    K,
    logLik,
    AICc,
    `Delta AICc`,
    `AICc Weight`
  )


clean_aic_table %>%
  gt() %>%
  tab_header(
    title = "Model-selection results for Delta use"
  ) %>%
  cols_label(
    Model = "Model",
    `Upstream Salinity Logger` = "Upstream salinity",
    `Downstream Salinity Logger` = "Downstream salinity",
    `Gage Height` = "Gage height",
    `Mean Temperature` = "Mean temperature",
    K = "K",
    logLik = "logLik",
    AICc = "AICc",
    `Delta AICc` = "\u0394AICc",
    `AICc Weight` = "Weight"
  ) %>%
  fmt_number(
    columns = c(logLik, AICc, `Delta AICc`, `AICc Weight`),
    decimals = 3
  ) %>%
  tab_options(
    table.font.names = "Times New Roman"
  )

clean_aic_table

coef_table <- summary(m_global)$coefficients$cond %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Predictor") %>%
  rename(
    Beta = Estimate,
    SE = `Std. Error`,
    z = `z value`,
    p = `Pr(>|z|)`
  ) %>%
  mutate(
    Odds_Ratio = exp(Beta)
  )

coef_table

