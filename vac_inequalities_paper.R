################################
### Vaccination inequalities ###
############ Paper #############
################################

# Libraries
library(viridisLite)
library(broom.mixed)
library(data.table)
library(corrplot)
library(ggplot2)
library(broom)
library(dplyr)
library(MASS)
library(lme4)
library(sf)

# Set seed for reproducibility
set.seed(69420)  


## Get data into an analysis ready format ##

# Load ansd tidy influenza data
influenza_dat <- fread("./VaccinationUptake_Study2 wtJ20J21 2025-06-18.csv")
influenza_dat$agegroup[influenza_dat$agegroup == "85\\+"] <- "85+" # Change age group label
influenza_dat$pc_vax <- (influenza_dat$vacc_seasonMSOAAgegroup / influenza_dat$totalpop_seasonMSOAAgegroup) * 100 # Calculate percentage vaccinated
#influenza_dat$hospILI_seasonMSOAAgegroup[is.na(influenza_dat$hospILI_seasonMSOAAgegroup)] <- 0 # NAs should be 0
influenza_dat$hospSARI2_seasonMSOAAgegroup[is.na(influenza_dat$hospSARI2_seasonMSOAAgegroup)] <- sample(0:5, sum(is.na(influenza_dat$hospSARI2_seasonMSOAAgegroup)), replace = TRUE)
names(influenza_dat)[names(influenza_dat) == "hospSARI2_seasonMSOAAgegroup"] <- "ili_count" # Rename outcome
influenza_dat$ili_pred <- (influenza_dat$hospSARI2_seasonAgegroup / influenza_dat$totalpop_seasonAgegroup) * influenza_dat$totalpop_seasonMSOAAgegroup # Predicted level of ILI admissions (offset)

# Get a version of 2019 IMD for MSOAs
imd_msoa <- fread("./covariates/imd2019_msoa_level_data.csv") # Data downloaded from https://research.mysociety.org/sites/imd2019/about/
names(imd_msoa)[names(imd_msoa) == "MSOAC"] <- "MSOACode" # Rename variables
imd_msoa <- imd_msoa[, c("MSOACode", "MSOADECILE", "MSOAQUINTILE", "IMD19 SCORE", "POPMID15")] # Subset just variables need

# Get covariates for outliers analysis
covars <- st_read("./covariates/cm_msoa_covariates.geojson") # Load

# Join together data
msoas_cm_age <- merge(influenza_dat, imd_msoa, by = "MSOACode", all.x = TRUE)
rm(influenza_dat, imd_msoa)


## Descriptive statistics ##

# Table 1 - summary statistics
total_stat <- msoas_cm_age[, list(total_ili_count = sum(ili_count, na.rm = TRUE), vacc_seasonMSOAAgegroup = sum(vacc_seasonMSOAAgegroup, na.rm = TRUE), totalpop_seasonMSOAAgegroup = sum(totalpop_seasonMSOAAgegroup, na.rm = TRUE)), by = c("agegroup", "season")] # Aggregate data to total counts
total_stat$total_pc_vax <- (total_stat$vacc_seasonMSOAAgegroup / total_stat$totalpop_seasonMSOAAgegroup) * 100 # Calculate percentage vaccinated
total_stat <- total_stat[, c(1,2,3,6)] # Keep only columns needed
median_stat <- msoas_cm_age[, list(median_ili_count = median(ili_count, na.rm = TRUE), median_pc_vax = median(pc_vax, na.rm = TRUE)), by = c("agegroup", "season")] # Aggregate data to median values
table1 <- merge(total_stat, median_stat, by = c("agegroup", "season"))# Merge tables together
fwrite(table1, "./Outputs/table1.csv")# Save
rm(total_stat, median_stat, table1) # Tidy
gc()

# Plot of vaccination uptake by IMD quintile and age band
vax_imd <- msoas_cm_age[, list(vacc_seasonMSOAAgegroup = sum(vacc_seasonMSOAAgegroup, na.rm = TRUE), totalpop_seasonMSOAAgegroup = sum(totalpop_seasonMSOAAgegroup, na.rm = TRUE), total_ili_count = sum(ili_count, na.rm = TRUE), total_pred_ili = sum(ili_pred, na.rm = TRUE)), by = c("MSOAQUINTILE", "agegroup", "season")] # Aggregate data to quintiles
vax_imd$pc_vax <- (vax_imd$vacc_seasonMSOAAgegroup / vax_imd$totalpop_seasonMSOAAgegroup) * 100 # Calculate percentage vaccinated
vax_imd$ratio <- vax_imd$total_ili_count / vax_imd$total_pred_ili
plot1 <- ggplot(vax_imd, aes(x = MSOAQUINTILE, y = pc_vax, color = agegroup)) + # Plot
  geom_point() +
  facet_wrap(~season) +
  scale_colour_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
  ylim(0,100) +
  labs(color = "Age band",
       x = "Deprivation Quintile (1 = most deprived, 5 = least deprived)",
       y = "Percent vaccinated (%)") +
  theme(text = element_text(size = 16))
ggsave(plot = plot1, filename = "./Outputs/vaccination_update_imd.jpeg", dpi = 300) # Save
plot1a <- ggplot(vax_imd, aes(x = MSOAQUINTILE, y = total_ili_count, color = season)) + # Plot
  geom_point() +
  facet_wrap(~agegroup) +
  scale_colour_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
  #ylim(0,1) +
  labs(color = "Season",
       x = "Deprivation Quintile (1 = most deprived, 5 = least deprived)",
       y = "Total number of influenza-like-illness admissions") +
  theme(text = element_text(size = 16))
ggsave(plot = plot1a, filename = "./Outputs/ili_imd.jpeg", dpi = 300) # Save
rm(vax_imd, plot1, plot1a) # Tidy
gc()


### Regression analyses ###


## Predicting count of ILI admissions ##

## Unadjusted model
model0a <- glmer.nb(ili_count ~ pc_vax + offset(log(ili_pred)) + (1|MSOACode), data = msoas_cm_age[msoas_cm_age$agegroup == "50-54",]) # Run each model per age group
model0b <- glmer.nb(ili_count ~ pc_vax + offset(log(ili_pred)) + (1|MSOACode), data = msoas_cm_age[msoas_cm_age$agegroup == "55-64",])
model0c <- glmer.nb(ili_count ~ pc_vax + offset(log(ili_pred)) + (1|MSOACode), data = msoas_cm_age[msoas_cm_age$agegroup == "65-74",])
model0d <- glmer.nb(ili_count ~ pc_vax + offset(log(ili_pred)) + (1|MSOACode), data = msoas_cm_age[msoas_cm_age$agegroup == "75-84",])
model0e <- glmer.nb(ili_count ~ pc_vax + offset(log(ili_pred)) + (1|MSOACode), data = msoas_cm_age[msoas_cm_age$agegroup == "85+",])

# Tidy up results
tidy_model1 <- tidy(model0a, exponentiate = TRUE, conf.int = TRUE) # Tidy each model
tidy_model2 <- tidy(model0b, exponentiate = TRUE, conf.int = TRUE) 
tidy_model3 <- tidy(model0c, exponentiate = TRUE, conf.int = TRUE)
tidy_model4 <- tidy(model0d, exponentiate = TRUE, conf.int = TRUE)
tidy_model5 <- tidy(model0e, exponentiate = TRUE, conf.int = TRUE)

tidy_model1$model <- "50-54" # Rename models
tidy_model2$model <- "55-64"
tidy_model3$model <- "65-74"
tidy_model4$model <- "75-84"
tidy_model5$model <- "85+"

# Create summary table of results
combined_unadj_results <- rbind(tidy_model1, tidy_model2, tidy_model3, tidy_model4, tidy_model5) # Join models together
combined_unadj_results$model_type <- "Unadjusted" # Define so can seperate out later
rm(model0a, model0b, model0c, model0d, model0e, tidy_model1, tidy_model2, tidy_model3, tidy_model4, tidy_model5) # Tidy
gc()


## Adjusted for deprivation 

# Set reference group for models
msoas_cm_age$MSOAQUINTILE <- as.factor(msoas_cm_age$MSOAQUINTILE)
msoas_cm_age$MSOAQUINTILE <- relevel(msoas_cm_age$MSOAQUINTILE, ref = 5)

# Regression models
model1a <- glmer.nb(ili_count ~ pc_vax + factor(MSOAQUINTILE) + offset(log(ili_pred)) + (1|MSOACode), data = msoas_cm_age[msoas_cm_age$agegroup == "50-54",]) # Run each model per age group
model1b <- glmer.nb(ili_count ~ pc_vax + factor(MSOAQUINTILE) + offset(log(ili_pred)) + (1|MSOACode), data = msoas_cm_age[msoas_cm_age$agegroup == "55-64",])
model1c <- glmer.nb(ili_count ~ pc_vax + factor(MSOAQUINTILE) + offset(log(ili_pred)) + (1|MSOACode), data = msoas_cm_age[msoas_cm_age$agegroup == "65-74",])
model1d <- glmer.nb(ili_count ~ pc_vax + factor(MSOAQUINTILE) + offset(log(ili_pred)) + (1|MSOACode), data = msoas_cm_age[msoas_cm_age$agegroup == "75-84",])
model1e <- glmer.nb(ili_count ~ pc_vax + factor(MSOAQUINTILE) + offset(log(ili_pred)) + (1|MSOACode), data = msoas_cm_age[msoas_cm_age$agegroup == "85+",])

# Tidy up results
tidy_model1 <- tidy(model1a, exponentiate = TRUE, conf.int = TRUE) # Tidy each model
tidy_model2 <- tidy(model1b, exponentiate = TRUE, conf.int = TRUE)
tidy_model3 <- tidy(model1c, exponentiate = TRUE, conf.int = TRUE)
tidy_model4 <- tidy(model1d, exponentiate = TRUE, conf.int = TRUE)
tidy_model5 <- tidy(model1e, exponentiate = TRUE, conf.int = TRUE)

tidy_model1$model <- "50-54" # Rename models
tidy_model2$model <- "55-64"
tidy_model3$model <- "65-74"
tidy_model4$model <- "75-84"
tidy_model5$model <- "85+"

# Create summary table of results
combined_results <- rbind(tidy_model1, tidy_model2, tidy_model3, tidy_model4, tidy_model5) # Join models together
combined_results$model_type <- "Adjusted" # Define model type
combined_results <- rbind(combined_results, combined_unadj_results) # Join tables together
rm(combined_unadj_results, tidy_model1, tidy_model2, tidy_model3, tidy_model4, tidy_model5) # Tidy

# Tidy table
combined_results <- combined_results[combined_results$term != "(Intercept)",] # Drop intercept
combined_results$term[combined_results$term == "pc_vax"] <- "Vaccinated" # Rename rows
combined_results$term[combined_results$term == "factor(MSOAQUINTILE)1"] <- "Quintile 1" 
combined_results$term[combined_results$term == "factor(MSOAQUINTILE)2"] <- "Quintile 2" 
combined_results$term[combined_results$term == "factor(MSOAQUINTILE)3"] <- "Quintile 3" 
combined_results$term[combined_results$term == "factor(MSOAQUINTILE)4"] <- "Quintile 4" 

# Add in Quintile 1 as reference
new_row <- data.frame(
  effect = rep("fixed", 5),
  group = rep(NA, 5),
  term = rep("Quintile 5", 5),
  estimate = rep(1, 5),
  std.error = rep(1, 5),
  statistic = rep(1, 5),
  p.value = rep(1, 5),
  conf.low = rep(1, 5),
  conf.high = rep(1, 5),
  model = c("50-54", "55-64", "65-74", "75-84", "85+"),
  model_type = rep("Adjusted", 5),
  stringsAsFactors = FALSE  # Prevents character columns from being converted to factors
)
combined_results <- rbind(combined_results, new_row) # Join on
fwrite(combined_results, "./Outputs/regression_model1.csv") # Save
rm(new_row) # Tidy
                          
# Plot
combined_results$var <- "Deprivation - model adjusted" # Add labels for presentation purpose - model type and variable
combined_results$var[combined_results$term == "Vaccinated" & combined_results$model_type == "Adjusted"] <- "Vaccination uptake (%) - model adjusted"
combined_results$var[combined_results$term == "Vaccinated" & combined_results$model_type == "Unadjusted"] <- "Vaccination uptake (%) - model unadjusted"
combined_results$var <- factor(combined_results$var, # Set as factor
  levels = c("Vaccination uptake (%) - model unadjusted", "Vaccination uptake (%) - model adjusted", "Deprivation - model adjusted")  # Define order for plotting purposes
)
plot2 <- ggplot(combined_results[combined_results$effect == "fixed",], aes(x = estimate, y = term, color = model)) + # Plot
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), 
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  facet_wrap(~var, scales = "free", nrow = 2) +
  scale_colour_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
  labs(color = "Age band",
       x = "Incidence Rate Ratio",
       y = "Variable") +
  theme(text = element_text(size = 16))
ggsave(plot = plot2, filename = "./Outputs/regression_model1.jpeg", dpi = 300) # Save
rm(combined_results, plot2) # Tidy
gc()


## Predicting the outliers ##

# Assess correlation of measures
hold <- covars[, c("MSOADECILE", "distance_to_gp", "pop_density", "mean_data_use_gb", "houses_pre1945_pc", "total_resident_communal", "multigen_pc", "poor_health_pc", "Asian..Asian.British.or.Asian.Welsh", "Mixed.or.Multiple.ethnic.groups", "Black..Black.British..Black.Welsh..Caribbean.or.African", "White", "Other.ethnic.group")] # Subset variables
hold$geometry <- NULL # Remove
colnames(hold) <- c("Deprivation", "Distance to GP", "Population density", "Data use", "Houses pre-1945", "Communal residents", "Multigeneration houses", "Poor health", "Asian", "Mixed", "Black", "White", "Other")
# Rename columns
cormat <- cor(hold) # Estimate correlation matrix
corrplot(cormat, type = "upper", method = "square", tl.col = "black", tl.srt = 45, addCoef.col = T, col = COL1('YlGn')) # Plot
# To save - do this manually in the 'plots' viewer as easier
rm(cormat, hold) # Tidy

# Get residuals
# Use Pearsons residuals rather than response why - "Particularly valuable when analyzing data with non-normal distributions, as they account for the variability in the response variable based on the mode. Due to their standardized nature, large absolute values of Pearson residuals readily point to potential outliers. Without standardization, it can be difficult to compare the magnitude of errors across different data points."
res_m1a <- augment(model1a, data = msoas_cm_age[msoas_cm_age$agegroup == "50-54",], type.residuals = "pearson") # Get residuals for first model
res_m1b <- augment(model1b, data = msoas_cm_age[msoas_cm_age$agegroup == "55-64",], type.residuals = "pearson") # Get residuals for second model
residuals <- rbind(res_m1a, res_m1b) # Join tables together
res_m1c <- augment(model1c, data = msoas_cm_age[msoas_cm_age$agegroup == "65-74",], type.residuals = "pearson") # Repeat process...
residuals <- rbind(residuals, res_m1c) # Join tables together
res_m1d <- augment(model1d, data = msoas_cm_age[msoas_cm_age$agegroup == "75-84",], type.residuals = "pearson")
residuals <- rbind(residuals, res_m1d) # Join tables together
res_m1e <- augment(model1e, data = msoas_cm_age[msoas_cm_age$agegroup == "85+",], type.residuals = "pearson")
residuals <- rbind(residuals, res_m1e) # Join tables together
rm(res_m1a, res_m1b, res_m1c, res_m1d, res_m1e, model1a, model1b, model1c, model1d, model1e, msoas_cm_age) # Tidy
gc()

# Join on residuals to covariates
outliers <- merge(residuals, covars, by.x = "MSOACode", by.y = "MSOA11CD", all.x = TRUE) # Join together
rm(residuals, covars) # Tidy 
gc()

# Regression models
outliers$pop_density <- outliers$pop_density / 1000 # Adjust else coefficient does not show
outliers$communal <- (outliers$total_resident_communal / outliers$POPMID15) * 100 # Adjust for population size
model2a <- lm(.resid ~ distance_to_gp + pop_density + mean_data_use_gb + houses_pre1945_pc + communal + multigen_pc + White, data = outliers[outliers$agegroup == "50-54",])
model2b <- lm(.resid ~ distance_to_gp + pop_density + mean_data_use_gb + houses_pre1945_pc + communal + multigen_pc + White, data = outliers[outliers$agegroup == "55-64",])
model2c <- lm(.resid ~ distance_to_gp + pop_density + mean_data_use_gb + houses_pre1945_pc + communal + multigen_pc + White, data = outliers[outliers$agegroup == "65-74",])
model2d <- lm(.resid ~ distance_to_gp + pop_density + mean_data_use_gb + houses_pre1945_pc + communal + multigen_pc + White, data = outliers[outliers$agegroup == "75-84",])
model2e <- lm(.resid ~ distance_to_gp + pop_density + mean_data_use_gb + houses_pre1945_pc + communal + multigen_pc + White, data = outliers[outliers$agegroup == "85+",])

# Tidy up results
tidy_model2a <- tidy(model2a) # Tidy each model
tidy_model2b <- tidy(model2b)
tidy_model2c <- tidy(model2c)
tidy_model2d <- tidy(model2d)
tidy_model2e <- tidy(model2e)

tidy_model2a$model <- "50-54" # Rename models
tidy_model2b$model <- "55-64"
tidy_model2c$model <- "65-74"
tidy_model2d$model <- "75-84"
tidy_model2e$model <- "85+"

# Create tidy summary table
combined_results2 <- rbind(tidy_model2a, tidy_model2b, tidy_model2c, tidy_model2d, tidy_model2e) # Join models together
combined_results2 <- combined_results2[combined_results2$term != "(Intercept)",] # Drop intercept
combined_results2$term <- with(combined_results2, recode(term, # Rename the rows 
                                                         "distance_to_gp" = "Distance to GP",
                                                         "pop_density" = "Population density",
                                                         "mean_data_use_gb" = "Mean data use (GB)",
                                                         "houses_pre1945_pc" = "Pre-1945 houses (%)",
                                                         "communal" = "People in communal establishments (%)",
                                                         "multigen_pc" = "Multigenerational households (%)",
                                                         "White" = "White ethnicity (%)"
))
combined_results2$conf.low <- combined_results2$estimate - (1.96 * combined_results2$std.error) # Confidence intervals
combined_results2$conf.high <- combined_results2$estimate + (1.96 * combined_results2$std.error)
fwrite(combined_results2, "./Outputs/regression_model2.csv") # Save
rm(tidy_model2a, tidy_model2b, tidy_model2c, tidy_model2d, tidy_model2e, model2a, model2b, model2c, model2d, model2e) # Tidy
gc()

# Plot results
plot4 <- ggplot(combined_results2, aes(x = estimate, y = term, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = (estimate - 1.96 * std.error), xmax = (estimate + 1.96 * std.error)), 
                position = position_dodge(width = 0.5), width = 0.2) +
  scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
  labs(color = "Age band",
       x = "Estimate",
       y = "Variable") +
  theme(text = element_text(size = 16))
ggsave(plot = plot4, filename = "./Outputs/regression_model2.jpeg", dpi = 300) # Save
rm(combined_results2, plot4) # Tidy
gc()


# Map residuals
sf_data <- st_sf(outliers, crs = 27700) # Convert object back to sf object
sf_data$.resid_transformed <- sign(sf_data$.resid) * log1p(abs(sf_data$.resid)) # Transform residuals by logging them - reduce visual impact of outliers
map <- ggplot() + # Plot
  geom_sf(data = sf_data, mapping = aes(fill = .resid_transformed), lwd = 0, color = NA) +
  facet_wrap(~agegroup, nrow = 2, scales = "fixed") +
  scale_fill_viridis_c(option = "viridis") +
  labs(fill = "Log Residuals") +
  theme_void() +
  theme(text = element_text(size = 16))
ggsave(plot = map, filename = "./Outputs/residuals_map.jpeg", dpi = 300) # Save
rm(map, sf_data, outliers) # Tidy
gc()
