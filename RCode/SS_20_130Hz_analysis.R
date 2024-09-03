# This script generates summary tables for the 20 x 130 Hz set-shifting experiment
# Created 11/01/2023 by AER
# Last Updated 08/31/2024


# Load required libraries - not all of these are used in this script
libraries <- c("gdata", "lme4", "emmeans", "metamisc", "psych", "sjPlot", "sjstats", 
               "sjmisc", "sjlabelled", "plyr", "phia", "dplyr", "jtools", "interactions", 
               "lmerTest", "car", "nloptr", "optimx", "data.table", "glmmTMB", 
               "effectsize", "fitdistrplus", "tidyr", "grid", "gridExtra", "broom.mixed")
lapply(libraries, library, character.only = TRUE)


# TODO: Update the path to match your file structure
dat <- read.csv("SS_dataset_20_130Hz.csv", header = TRUE)


# Data preprocessing
dat <- dat %>%
  mutate(Subject = factor(Subject),
         Stim = factor(Stim),
         Protocol = factor(Protocol),
         Task_Type = factor(Task_Type),
         Session_factor = factor(Session),
         Site = factor(Site),
         Stim_frequency_group = factor(Stim_frequency_group),
         Stim_frequency = factor(Stim_frequency)) %>%
  droplevels() %>%
  within(Stim_frequency_group <- relevel(Stim_frequency_group, ref = "20")) %>%
  within(Stim <- relevel(Stim, ref = "ON")) %>%
  within(Stim_frequency <- relevel(Stim_frequency, ref = "0"))

# Set up model contrasts
model_contrasts <- c("Stim", "Protocol", "Task_Type", "Session_factor", "Stim_frequency_group", "Stim_frequency")
dat[model_contrasts] <- lapply(dat[model_contrasts], function(x) {contrasts(x) <- "contr.treatment"; x})

# Control settings for models
control <- glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))

# Model building function
build_model <- function(formula, data, family, random = "1|Stim_frequency_group/Subject") {
  glmmPQL(formula, random = random, family = family, data = data, control = control)
}

# Remove unused factor levels
dat <- drop.levels(dat)

# Preprocess data for errors model
dat.errors <- dat %>%
  group_by(Subject, Stim, Stim_frequency, Stim_frequency_group, Protocol, Session) %>%
  summarize(correct = sum(Trial == "correct"),
            incorrect = sum(Trial == "incorrect"),
            .groups = "drop") %>%
  mutate(incorrect = ifelse(is.na(incorrect), 0, incorrect),
         across(c(Stim, Protocol, Session, Stim_frequency, Stim_frequency_group), factor)) %>%
  droplevels()

# Set up contrasts for error model
model_contrasts2 <- c("Stim", "Protocol", "Session", "Stim_frequency_group", "Stim_frequency")
dat.errors[model_contrasts2] <- lapply(dat.errors[model_contrasts2], function(x) {contrasts(x) <- "contr.treatment"; x})

# Summarize max duration of the test
summary_max_dat <- dat %>%
  group_by(Subject, Stim, Stim_frequency_group, Stim_frequency, Protocol, Session) %>%
  summarise(max_duration = max(Time_in_sec, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(c(Stim, Stim_frequency_group, Stim_frequency, Protocol, Session), factor)) %>%
  droplevels()

# Set up contrasts for max duration model
summary_max_dat[model_contrasts2] <- lapply(summary_max_dat[model_contrasts2], function(x) {contrasts(x) <- "contr.treatment"; x})
summary_max_dat <- within(summary_max_dat, Stim_frequency <- relevel(Stim_frequency, ref = "0"))
summary_max_dat <- within(summary_max_dat, Stim <- relevel(Stim, ref = "OFF"))
summary_max_dat <- within(summary_max_dat, Stim_frequency_group <- relevel(Stim_frequency_group, ref = "130"))

# Build models using glmer
# Reaction Time (RT) Model
bothgroups.rt_glmer <- glmer(RT ~ 1 + Stim_frequency + Task_Type + Protocol + (1|Stim_frequency_group:Subject), 
                             family = Gamma(link=identity), data = dat, 
                             control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

# Accuracy Model
bothgroups.accuracy_glmer <- glmer(Correctness ~ 1 + Stim_frequency + Task_Type + Protocol + (1|Stim_frequency_group:Subject), 
                                   family = binomial, data = dat, 
                                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

# Errors Model
bothgroups.errors_glmer <- glmer(incorrect ~ 1 + Stim_frequency + Protocol + (1|Stim_frequency_group:Subject), 
                                 family = poisson, data = dat.errors, 
                                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

# Max Duration Model
bothgroups.max_duration_glmer <- glmer(max_duration ~ 1 + Stim_frequency + Protocol + (1|Stim_frequency_group:Subject), 
                                       family = Gamma(link=identity), data = summary_max_dat, 
                                       control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

# Function to summarize model results
summarize_model <- function(model) {
  tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      p.value = as.numeric(p.value),
      Bonferroni.p = ifelse(row_number() %in% 2:3, pmin(p.value * 4, 1, na.rm = TRUE), NA_real_),
      p.value = ifelse(p.value > 1e-3, 
                       format(round(p.value, 3), nsmall = 3, scientific = FALSE),
                       sprintf("%.3e", p.value)),
      estimate = round(estimate, 3),
      std.error = round(std.error, 3),
      statistic = round(statistic, 3),
      Bonferroni.p = ifelse(Bonferroni.p > 1e-3, 
                            format(round(Bonferroni.p, 3), nsmall = 3, scientific = FALSE),
                            sprintf("%.3e", Bonferroni.p))
    ) %>%
    mutate(Bonferroni.p = as.character(Bonferroni.p)) %>%
    mutate_all(~ifelse(is.na(.), "", .))
}

# Build models using glmmPQL (for degrees of freedom)
bothgroups.rt <- build_model(RT ~ 1 + Stim_frequency + Task_Type + Protocol, random = ~ 1|Stim_frequency_group/Subject, dat, Gamma(link = identity))
bothgroups.accuracy <- build_model(Correctness ~ 1 + Stim_frequency + Task_Type + Protocol, random = ~ 1|Stim_frequency_group/Subject, dat, binomial)
bothgroups.errors <- build_model(incorrect ~ 1 + Stim_frequency + Protocol, random = ~ 1 | Stim_frequency_group/Subject, dat.errors, poisson)
bothgroups.max_duration <- glmmPQL(max_duration ~ 1 + Stim_frequency + Protocol, random = ~1 | Stim_frequency_group/Subject, summary_max_dat, family = Gamma(link = identity))

# Summarize models
bothgroups.rt_summary <- summarize_model(bothgroups.rt)
bothgroups.accuracy_summary <- summarize_model(bothgroups.accuracy)
bothgroups.errors_summary <- summarize_model(bothgroups.errors)
bothgroups.max_duration_summary <- summarize_model(bothgroups.max_duration)

# Function to modify summary tables
modify_summary <- function(model_summary) {
  dplyr::select(model_summary, term, estimate, std.error, statistic, df, p.value, Bonferroni.p) %>%
    dplyr::rename(
      Name = term, 
      Coefficient = estimate, 
      SE = std.error,
      tstat = statistic, 
      DF = df,
      p.value = p.value,
      Adjusted.p = Bonferroni.p
    ) %>%
    dplyr::mutate(
      Name = str_replace_all(Name, c("Stim_frequency20" = "Stim_frequency_20",
                                     "Stim_frequency130" = "Stim_frequency_130",
                                     "Site" = "Site_",
                                     "Task_Type" = "Rule_",
                                     "frequency_group20" = "20Hz",
                                     "frequency_group130" = "130Hz",
                                     "Protocol45" = "Schedule_2",
                                     "Protocol46" = "Schedule_3",
                                     "Protocol47" = "Schedule_4")),
      DF = format(DF, big.mark = ""),
      tstat = format(tstat, big.mark = ""),
      Coefficient = format(Coefficient, big.mark = "")
    )
}

# Apply the function to the model summaries
final_table_lf_vs_hf_RT <- modify_summary(bothgroups.rt_summary)
final_table_lf_vs_hf_Accuracy <- modify_summary(bothgroups.accuracy_summary)
final_table_lf_vs_hf_Duration <- modify_summary(bothgroups.max_duration_summary)
final_table_lf_vs_hf_Errors <- modify_summary(bothgroups.errors_summary)

# Summarize glmer Models
bothgroups.rt_summary_glmer <- summarize_model(bothgroups.rt_glmer)
bothgroups.accuracy_summary_glmer <- summarize_model(bothgroups.accuracy_glmer)
bothgroups.errors_summary_glmer <- summarize_model(bothgroups.errors_glmer)
bothgroups.max_duration_summary_glmer <- summarize_model(bothgroups.max_duration_glmer)

# Modify glmer summary tables
modify_summary_glmer <- function(model_summary) {
  dplyr::select(model_summary, term, estimate, std.error, statistic, p.value, Bonferroni.p) %>%
    dplyr::rename(
      Name = term, 
      Coefficient = estimate, 
      SE = std.error,
      tstat = statistic, 
      p.value = p.value,
      Adjusted.p = Bonferroni.p
    ) %>%
    dplyr::mutate(
      Name = str_replace_all(Name, c("Stim_frequency20" = "Stim_frequency_20",
                                     "Stim_frequency130" = "Stim_frequency_130",
                                     "Site" = "Site_",
                                     "Task_Type" = "Rule_",
                                     "frequency_group20" = "20Hz",
                                     "frequency_group130" = "130Hz",
                                     "Protocol45" = "Schedule_2",
                                     "Protocol46" = "Schedule_3",
                                     "Protocol47" = "Schedule_4")),
      tstat = format(tstat, big.mark = ""),
      Coefficient = format(Coefficient, big.mark = "")
    )
}

# Apply the function to the glmer model summaries
final_table_lf_vs_hf_RT_glmer <- modify_summary_glmer(bothgroups.rt_summary_glmer)
final_table_lf_vs_hf_Accuracy_glmer <- modify_summary_glmer(bothgroups.accuracy_summary_glmer)
final_table_lf_vs_hf_Duration_glmer <- modify_summary_glmer(bothgroups.max_duration_summary_glmer)
final_table_lf_vs_hf_Errors_glmer <- modify_summary_glmer(bothgroups.errors_summary_glmer)

# Add DF column to glmer summary tables and reorder columns
add_df_and_reorder <- function(glmer_table, pql_table) {
  glmer_table %>% 
    dplyr::mutate(DF = pql_table$DF) %>%
    dplyr::select(Name, Coefficient, SE, tstat, DF, p.value, Adjusted.p)
}

final_table_lf_vs_hf_RT_glmer <- add_df_and_reorder(final_table_lf_vs_hf_RT_glmer, final_table_lf_vs_hf_RT)
final_table_lf_vs_hf_Accuracy_glmer <- add_df_and_reorder(final_table_lf_vs_hf_Accuracy_glmer, final_table_lf_vs_hf_Accuracy)
final_table_lf_vs_hf_Duration_glmer <- add_df_and_reorder(final_table_lf_vs_hf_Duration_glmer, final_table_lf_vs_hf_Duration)
final_table_lf_vs_hf_Errors_glmer <- add_df_and_reorder(final_table_lf_vs_hf_Errors_glmer, final_table_lf_vs_hf_Errors)
