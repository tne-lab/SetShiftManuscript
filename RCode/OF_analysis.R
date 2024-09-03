# This script generates summary tables for the Open-field test
# Created 11/01/2023 by AER
# Last Updated 08/31/2024

# Load required libraries - not all of these are used in this script
libraries <- c("gdata", "lme4", "emmeans", "metamisc", "psych", "sjPlot", "sjstats", 
               "sjmisc", "sjlabelled", "plyr", "phia", "dplyr", "jtools", "interactions", 
               "lmerTest", "car", "nloptr", "optimx", "data.table", "glmmTMB", 
               "effectsize", "fitdistrplus", "tidyr", "grid", "gridExtra", "broom.mixed")
lapply(libraries, library, character.only = TRUE)


# IMPORTANT: Update this path to match your local file structure
# Load and preprocess data
dat_OF <- read.csv("OF_data.csv", header = TRUE)

# Data cleaning and preparation
dat_OF <- drop.levels(dat_OF)  # Remove unused factor levels
dat_OF <- dat_OF[complete.cases(dat_OF), ]  # Remove rows with NA values

# Convert variables to factors
dat_OF$Stim <- factor(dat_OF$Stim)
dat_OF$Session <- factor(dat_OF$Session)
dat_OF$Subject <- factor(dat_OF$Subject)

# Set up contrasts for factorial variables
contrasts(dat_OF$Stim) <- "contr.treatment"
contrasts(dat_OF$Session) <- "contr.treatment"
contrasts(dat_OF$Subject) <- "contr.treatment"

dat_OF <- drop.levels(dat_OF)  # Remove any newly created unused levels

# Transform response variables
dat_OF$log_travelled_distance_cm <- log(dat_OF$traveled_distance_cm)
dat_OF$log_speed_moving_cm_s <- log(dat_OF$speed_moving_cm.s)

# Model fitting

# 1. Immobility time model (using Gamma distribution)
m1.immobility <- glmmPQL(middle_body.time.stationary ~ 1 + Stim, 
                         random = ~ 1|Subject, 
                         family = Gamma(link=identity), 
                         data = dat_OF, 
                         control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

# 2. Distance traveled model (log-transformed, using Gaussian distribution)
m1.distance <- glmmPQL(log_travelled_distance_cm ~ 1 + Stim, 
                       random = ~ 1|Subject, 
                       family = gaussian, 
                       data = dat_OF)

# 3. Speed model (log-transformed, using Gaussian distribution)
m1.speed <- glmmPQL(log_speed_moving_cm_s ~ 1 + Stim, 
                    random = ~ 1|Subject, 
                    family = gaussian, 
                    data = dat_OF)



# Custom theme for tables
my_theme <- ttheme_default(
  core = list(fg_params = list(hjust = 0.5, x = 0.5),
              bg_params = list(fill = "transparent", col = "black", lwd = 1)),
  colhead = list(fg_params = list(hjust = 0.5, x = 0.5, fontsize = 11),
                 bg_params = list(fill = "transparent", col = "black", lwd = 1)),
  rowhead = list(fg_params = list(hjust = 0.5, x = 0.5, fontsize = 11),
                 bg_params = list(fill = "transparent", col = "black", lwd = 1)),
  rowhead_fringe = list(bg_params = list(fill = "transparent", col = NA, lwd = 1)),
  colhead_fringe = list(bg_params = list(fill = "transparent", col = NA, lwd = 1))
)

# Function to create summary tables
create_summary_table <- function(model, model_name) {
  model_summary <- tidy(model, effects = "fixed", conf.int = TRUE)
  model_summary$p.value <- as.numeric(model_summary$p.value)
  
  # Add Bonferroni-corrected p-values
  model_summary <- model_summary %>%
    mutate(
      Bonferroni.p = ifelse(row_number() == 2,
                            pmin(p.value * 3, 1, na.rm = TRUE),
                            NA_real_)
    )
  
  # Format p-values and round numeric columns
  model_summary <- model_summary %>%
    mutate(
      p.value = ifelse(p.value > 1e-3, 
                       format(round(p.value, 3), nsmall = 3, scientific = FALSE),
                       sprintf("%.3e", p.value)),
      estimate = round(estimate, 3),
      std.error = round(std.error, 3),
      statistic = round(statistic, 3),
      Bonferroni.p = ifelse(Bonferroni.p > 1e-3, 
                            format(round(Bonferroni.p, 3), nsmall = 3, scientific = FALSE),
                            sprintf("%.3e", Bonferroni.p))
    )
  
  model_summary$Bonferroni.p <- as.character(model_summary$Bonferroni.p)
  model_summary[is.na(model_summary)] <- ""
  
  final_table <- dplyr::select(model_summary, term, estimate, std.error, statistic, df, p.value, Bonferroni.p)
  final_table <- dplyr::rename(final_table,
                               Name = term, 
                               Coefficient = estimate, 
                               SE = std.error,
                               tstat = statistic, 
                               DF = df,
                               p.value = p.value,
                               Adjusted.p = Bonferroni.p)
  
  final_table <- final_table %>% mutate(Name = str_replace_all(Name, "Stim", "Stim_"))
  
  return(final_table)
}

# Create summary tables for each model
final_table_distance <- create_summary_table(m1.distance, "Distance")
final_table_immobility <- create_summary_table(m1.immobility, "Immobility")
final_table_speed <- create_summary_table(m1.speed, "Speed")
