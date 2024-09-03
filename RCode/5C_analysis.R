# This script generates summary tables for the 5-choice serial reaction time task
# Created 11/01/2023 by AER
# Last Updated 08/31/2024

# Load required libraries - not all of these are used in this script
libraries <- c("emmeans", "metamisc", "psych", "sjPlot", "sjstats", "broom.mixed", 
               "sjmisc", "sjlabelled", "plyr", "phia", "dplyr", "jtools", "interactions", 
               "gtsummary", "lmerTest", "car", "nloptr", "optimx", "data.table", "glmmTMB", 
               "effectsize", "fitdistrplus", "rcompanion", "huxtable", "grid", "gridExtra", 
               "gdata", "stringr", "tidyr", "officer", "flextable", "magrittr", "readxl")

# Load all required libraries
lapply(libraries, library, character.only = TRUE)

# Set global options
options("jtools-digits" = 3)

# IMPORTANT: Update this path to match your local file structure
dat <- read.csv("5C_data.csv", header = TRUE)

# Convert columns to factors
cols_factor <- c('Animal', 'Session.Number', "Outcome", "Protocol", "Stim")
dat[cols_factor] <- lapply(dat[cols_factor], factor)

# Process data for all the variables
process_data <- function(data, type) {
  if (type == "premature") {
    data %>% 
      filter(Omission_or_not == "Not_omission") %>% 
      mutate(Premature_binary = if_else(Premature_or_not == "Not_premature", 0, 1)) %>% 
      droplevels()
  } else if (type == "omission") {
    data %>% 
      mutate(Omission_binary = if_else(Omission_or_not == "Omission", 1, 0)) %>% 
      droplevels()
  } else if (type == "accuracy") {
    data %>%
      filter(Accuracy %in% c("Reward", "Error")) %>%
      mutate(Accuracy_binary = ifelse(Accuracy == "Reward", 1, 0)) %>%
      droplevels()
  } else if (type == "rt") {
    data %>%
      filter(Accuracy %in% c("Reward", "Error"),
             Reaction.Time < 5.5,
             Reaction.Time > 0.01) %>%
      mutate(Reaction.Time = as.numeric(Reaction.Time)) %>%
      droplevels() %>%
      {rownames(.) <- NULL; .}
  }
}

# Process data for different analyses
dat_premature <- process_data(dat, "premature")
dat <- process_data(dat, "omission")
dat_accuracy <- process_data(dat, "accuracy")
dat_rt <- process_data(dat, "rt")

# Set contrasts
cols_contrast <- c("Stim", "Protocol", "Animal")
for(col in cols_contrast) {
  contrasts(dat_premature[[col]]) <- "contr.treatment"
  contrasts(dat[[col]]) <- "contr.treatment"
  contrasts(dat_accuracy[[col]]) <- "contr.treatment"
  contrasts(dat_rt[[col]]) <- "contr.treatment"
}

# Fit models
fit_model <- function(formula, data, family = binomial) {
  glmmPQL(formula, random = ~ 1 | Animal, family = family, data = data)
}

m1.premature_glmmPQL <- fit_model(Premature_binary ~ 1 + Stim + Protocol, dat_premature)
m1.omissions_glmmPQL <- fit_model(Omission_binary ~ 1 + Stim + Protocol, dat)
m1.accuracy_glmmPQL <- fit_model(Accuracy_binary ~ 1 + Stim + Protocol, dat_accuracy)
m1.RT_glmmPQL <- fit_model(Reaction.Time ~ 1 + Stim + Protocol, dat_rt, family = Gamma(link = identity))

# Create a custom theme for the tables
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

# Function to process model summary
process_model_summary <- function(model, n_comparisons = 4) {
  model_summary <- tidy(model, effects = "fixed", conf.int = TRUE)
  model_summary$p.value <- as.numeric(model_summary$p.value)
  
  model_summary <- model_summary %>%
    mutate(
      Bonferroni.p = ifelse(row_number() == 2,
                            pmin(p.value * n_comparisons, 1, na.rm = TRUE),
                            NA_real_),
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
  
  return(model_summary)
}

# Function to create final table
create_final_table <- function(model_summary) {
  final_table <- dplyr::select(model_summary, term, estimate, std.error, statistic, df, p.value, Bonferroni.p)
  final_table <- dplyr::rename(final_table,
                               Name = term, 
                               Coefficient = estimate, 
                               SE = std.error,
                               tstat = statistic, 
                               DF = df,
                               p.value = p.value,
                               Adjusted.p = Bonferroni.p)
  
  final_table <- final_table %>%
    mutate(
      Name = str_replace_all(Name, "Stim", "Stim_"),
      Name = str_replace_all(Name, "on", "ON"),
      Name = str_replace_all(Name, "Protocol5CSRTT Day ", "Schedule_")
    )
  
  return(final_table)
}


# Process model summaries and create final tables
model_summaries <- list(
  premature = process_model_summary(m1.premature_glmmPQL),
  omissions = process_model_summary(m1.omissions_glmmPQL),
  accuracy = process_model_summary(m1.accuracy_glmmPQL),
  rt = process_model_summary(m1.RT_glmmPQL)
)

final_tables <- lapply(model_summaries, create_final_table)

# Function to format flextable
format_flextable <- function(table) {
  flextable(table) %>%
    set_table_properties(layout = "autofit") %>%
    theme_box() %>%
    theme_vanilla() %>%
    align(align = "center", part = "all") %>%
    border(border = fp_border(color = "black", width = 0.5), part = "all") %>%
    autofit()
}


final_tables$premature
final_tables$omissions
final_tables$accuracy
final_tables$rt

