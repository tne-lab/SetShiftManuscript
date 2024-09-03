# This script generates summary tables for the main set-shifting experiment
# Created 11/01/2023 by AER
# Last Updated 08/31/2024


# Load required libraries - not all of these are used in this script
libraries <- c("emmeans", "metamisc", "psych", "sjPlot", "sjstats", "broom.mixed", 
               "sjmisc", "sjlabelled", "plyr", "phia", "dplyr", "jtools", "interactions", 
               "gtsummary", "lmerTest", "car", "nloptr", "optimx", "data.table", "glmmTMB", 
               "effectsize", "fitdistrplus", "rcompanion", "huxtable", "grid", "gridExtra", 
               "gdata", "stringr", "tidyr", "stringr", "officer", "flextable", "magrittr")

# Load all libraries
lapply(libraries, library, character.only = TRUE)

# Set the number of digits for jtools package output
options("jtools-digits" = 3)


# IMPORTANT: Update this path to match your local file structure
dat <- read.csv("Main_SS_dataset.csv", header = TRUE)

# Data Preprocessing and Transformation
dat <- dat %>%
  # Handle missing values in Changing_Point column
  mutate(Correctness = ifelse(is.na(Changing_Point), Correctness, Changing_Point),
         # Convert variables to factors for categorical analysis
         Subject = factor(Subject),
         Stim = factor(Stim),
         Site = factor(Site),
         Protocol = factor(Protocol),
         Task_Type = factor(Task_Type),
         Session_factor = factor(Session)) %>%
  # Remove unused factor levels
  droplevels() %>%

# Set up contrasts for the model
model_contrasts <- c("Stim", "Site", "Protocol", "Task_Type", "Session_factor")
dat[model_contrasts] <- lapply(dat[model_contrasts], function(x) {contrasts(x) <- "contr.treatment"; x})

# Function to create summary statistics
create_summary <- function(data, grouping_vars, summary_fn) {
  data %>%
    group_by(!!!syms(grouping_vars)) %>%
    summarise(
      correct = summary_fn(Trial == "correct"),
      incorrect = summary_fn(Trial == "incorrect"),
      .groups = "drop") %>%
    droplevels() %>%
    within(Site <- relevel(Site, ref = "ventral"))
}

# Function to create summary of maximum values
create_summary_max <- function(data, grouping_vars, summary_var) {
  data %>%
    group_by(!!!syms(grouping_vars)) %>%
    summarise(max_time = max(!!sym(summary_var), na.rm = TRUE),
              .groups = "drop") %>%
    droplevels() %>%
    within(Site <- relevel(Site, ref = "ventral"))
}

# Create summaries for maximum time and correct/incorrect trials
summary_max_dat <- create_summary_max(dat, c("Subject", "Stim", "Site", "Protocol", "Session_factor"), "Time_in_sec_decision_made")
summary_incorrect_correct <- create_summary(dat, 
                                            c("Subject", "Stim", "Site", "Protocol", "Session_factor"), 
                                            sum)

# Set up control parameters for generalized linear mixed models
control <- glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))

# Function to build generalized linear mixed models
build_model <- function(formula, data, family, random = "~ 1 | Subject") {
  glmmPQL(formula, random = random, family = family, data = data)
}

# Function to summarize model results and add adjusted p-values
summarize_model <- function(model) {
  tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(p.value = as.numeric(p.value),
           p.value = ifelse(p.value < 1e-26, 1e-26, p.value),
           Bonferroni.p = ifelse(row_number() > (n() - 4),
                                 pmin(p.value * 4, 1, na.rm = TRUE),
                                 NA_real_),
           p.value = ifelse(p.value > 1e-3, 
                            format(round(p.value, 3), nsmall = 3, scientific = FALSE),
                            sprintf("%.3e", p.value)),
           estimate = round(estimate, 3),
           std.error = round(std.error, 3),
           statistic = round(statistic, 3),
           Bonferroni.p = ifelse(Bonferroni.p > 1e-3, 
                                 format(round(Bonferroni.p, 3), nsmall = 3, scientific = FALSE),
                                 sprintf("%.3e", Bonferroni.p))) %>%
    mutate(Bonferroni.p = as.character(Bonferroni.p)) %>%
    mutate_all(~ifelse(is.na(.), "", .))
}

# Define and fit the models

# Reaction Time (RT) Model
RT_model_glmmPQL <- build_model(RT ~ 1 + Stim:Site + Site + Task_Type + Protocol, random = ~ 1 | Subject, dat, Gamma(link = identity))

# Accuracy Model
Accuracy_model_glmmPQL <- build_model(Correctness ~ 1 + Stim:Site + Site + Task_Type + Protocol, random = ~ 1 | Subject, dat, binomial)

# Session Duration Model
Session_duration_model_glmmPQL <- build_model(max_time ~ 1 + Stim:Site + Site + Protocol, random = ~ 1 | Subject, summary_max_dat, gaussian)

# Error Model
Errors_model_glmmPQL <- build_model(incorrect ~ 1 + Stim:Site + Site + Protocol, random = ~ 1 | Subject, summary_incorrect_correct, negative.binomial(theta = 0.0001))

# Summarize the models
summary_RT <- summarize_model(RT_model_glmmPQL)
summary_Accuracy <- summarize_model(Accuracy_model_glmmPQL)
summary_Duration <- summarize_model(Session_duration_model_glmmPQL)
summary_Errors <- summarize_model(Errors_model_glmmPQL)

# Function to modify summary tables for better readability
modify_summary <- function(model_summary) {
  dplyr::select(model_summary, term, estimate, std.error, statistic, df, p.value, Bonferroni.p) %>%
    dplyr::rename(
      Name = term, 
      Coefficient = estimate, 
      SE = std.error,
      tstat = statistic, 
      DF = df,
      p.value = p.value,
      Bonferroni.p = Bonferroni.p
    ) %>%
    dplyr::mutate(
      Name = stringr::str_replace_all(Name, "Stim", "Stim_"),
      Name = stringr::str_replace_all(Name, "Site", "Site_"),
      Name = stringr::str_replace_all(Name, "Task_Type", "Rule_"),
      Name = stringr::str_replace_all(Name, "Protocol45", "Schedule_2"),
      Name = stringr::str_replace_all(Name, "Protocol46", "Schedule_3"),
      Name = stringr::str_replace_all(Name, "Protocol47", "Schedule_4"),
      DF = format(DF, big.mark = ""),
      tstat = format(tstat, big.mark = ""),
      Coefficient = format(Coefficient, big.mark = "")
    )
}

# Apply the modify_summary function to create final summary tables
final_table_RT <- modify_summary(summary_RT)
final_table_Accuracy <- modify_summary(summary_Accuracy)
final_table_Duration <- modify_summary(summary_Duration)
final_table_Errors <- modify_summary(summary_Errors)
