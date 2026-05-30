install.packages(
  "simtrial",
  repos = c(
    "https://merck.r-universe.dev",
    "https://cloud.r-project.org"
  )
);

install.packages(c(
  "survival",
  "survminer",
  "tidyverse",
  "gsDesign",
  "rpact",
  "gtsummary",
  "gt",
  "flextable"
))



library(simtrial)
library(survival)
library(survminer)
library(tidyverse)
library(gsDesign)
library(rpact)



data(ex3_cure_with_ph)

head(ex3_cure_with_ph)


str(ex3_cure_with_ph)
summary(ex3_cure_with_ph)

#Missing values
colSums(is.na(ex3_cure_with_ph))
#No of patients
nrow(ex3_cure_with_ph)
#Treatment distribution
table(ex3_cure_with_ph$trt)
#Event distribution
table(ex3_cure_with_ph$evntd)

#Summary
summary(ex3_cure_with_ph)

#Treatment Group Counts
table(ex3_cure_with_ph$trt)

#Event V/s censoring count
table(ex3_cure_with_ph$evntd)


#Event proportion by treatment
table(ex3_cure_with_ph$trt,ex3_cure_with_ph$evntd)

#Histogram of survival times
ggplot(ex3_cure_with_ph,
       aes(x = month)) +
  geom_histogram(
    bins = 30
  ) +
  labs(
    title = "Distribution of Survival Times",
    x = "Months",
    y = "Frequency"
  )


#Survival Time by Treatment
ggplot(ex3_cure_with_ph,
       aes(x = month,
           fill = factor(trt))) +
  geom_histogram(
    bins = 30,
    alpha = 0.6,
    position = "identity"
  ) +
  labs(
    title = "Survival Time by Treatment Group",
    x = "Months",
    fill = "Treatment"
  )




#Boxplot of Survival Times
ggplot(ex3_cure_with_ph,
       aes(x = factor(trt),
           y = month,
           fill = factor(trt))) +
  geom_boxplot() +
  labs(
    title = "Survival Times by Treatment",
    x = "Treatment Group",
    y = "Months"
  )


#Censoring visualisation
#Event v/s sensored
ggplot(ex3_cure_with_ph,
       aes(x = factor(evntd),
           fill = factor(evntd))) +
  geom_bar() +
  labs(
    title = "Event vs Censoring",
    x = "Event Indicator",
    y = "Count"
  )



#Follow-Up Time by Event Status
ggplot(ex3_cure_with_ph,
       aes(x = factor(evntd),
           y = month,
           fill = factor(evntd))) +
  geom_boxplot() +
  labs(
    title = "Follow-Up Time by Event Status",
    x = "Event Status",
    y = "Months"
  )

##upto added in report


#Survival Object
surv_object <- Surv(
  time  = ex3_cure_with_ph$month,
  event = ex3_cure_with_ph$evntd
)



#Fit Overall KM Curve

km_overall <- survfit(surv_object ~ 1)

summary(km_overall)



#Plot Overall KM Curve
ggsurvplot(
  km_overall,
  data = ex3_cure_with_ph,
  
  conf.int = TRUE,
  risk.table = TRUE,
  censor = TRUE,
  
  xlab = "Months",
  ylab = "Survival Probability",
  
  title = "Overall Kaplan-Meier Survival Curve"
)





#KM CURVES BY TREATMENT GROUP
km_treatment <- survfit(
  Surv(month, evntd) ~ trt,
  data = ex3_cure_with_ph
)

summary(km_treatment)


#Plot Treatment vs Control KM Curves
ggsurvplot(
  km_treatment,
  data = ex3_cure_with_ph,
  
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  censor = TRUE,
  
  legend.title = "Treatment",
  legend.labs = c("Control", "Treatment"),
  
  xlab = "Months",
  ylab = "Survival Probability",
  
  title = "Kaplan-Meier Curves by Treatment Group"
)

#Median Survival Time
summary(km_treatment)$table


#Survival Probabilities at Fixed Times
summary(
  km_treatment,
  times = c(6, 12, 24)
)



#Advanced Publication-Style Plot
ggsurvplot(
  km_treatment,
  data = ex3_cure_with_ph,
  
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  
  pval = TRUE,
  
  surv.median.line = "hv",
  
  ggtheme = theme_bw(),
  
  palette = c("#E7B800", "#2E9FDF"),
  
  title = "Treatment vs Control Survival",
  
  xlab = "Time (Months)",
  ylab = "Overall Survival Probability"
)

#Perform Log-Rank Test
logrank_test <- survdiff(
  Surv(month, evntd) ~ trt,
  data = ex3_cure_with_ph
)

logrank_test

logrank_test$chisq

p_value <- 1 - pchisq(
  logrank_test$chisq,
  df = 1
)

p_value

#P-value to KM Plot
ggsurvplot(
  km_treatment,
  data = ex3_cure_with_ph,
  
  conf.int = TRUE,
  risk.table = TRUE,
  
  pval = TRUE,
  
  legend.title = "Treatment",
  legend.labs = c("Control", "Treatment"),
  
  title = "Kaplan-Meier Curves with Log-Rank Test"
)



#Fit cox model
cox_fit <- coxph(
  Surv(month, evntd) ~ trt,
  data = ex3_cure_with_ph
)

summary(cox_fit)
#Extract hazard ratio
exp(coef(cox_fit))
#Confidence interval
exp(confint(cox_fit))
#Forest Plot
ggforest(cox_fit)
#Survival Curves from Cox Model
ggsurvplot(
  survfit(cox_fit),
  data = ex3_cure_with_ph,
  conf.int = TRUE
)

#Schoenfeld Residual Test
ph_test <- cox.zph(cox_fit)

ph_test

#Plot PH Diagnostics
plot(ph_test)

#load package
library(gsDesign)

#Create O’Brien–Fleming Design
design_of <- gsDesign(
  k = 4,
  test.type = 2,
  alpha = 0.05,
  beta = 0.20,
  sfu = "OF"
)

design_of

#Examine Boundaries
design_of$upper$bound

#Nominal P-values
design_of$upper$prob

#Plot Sequential Boundaries
plot(design_of)

#Pocock Design
design_pocock <- gsDesign(
  k = 4,
  test.type = 2,
  alpha = 0.05,
  beta = 0.20,
  sfu = "Pocock"
)

design_pocock

#Compare Boundaries
design_of$upper$bound

design_pocock$upper$bound



library(survival)
library(gsDesign)
library(ggplot2)
library(dplyr)


design_of <- gsDesign(
  k = 4,
  test.type = 2,
  alpha = 0.05,
  beta = 0.20,
  sfu = "OF"
)

print(design_of)



boundaries <- design_of$upper$bound

boundaries



looks <- c(0.25, 0.50, 0.75, 1.00)


simulate_trial <- function(

    n_control = 150,
    n_treatment = 150,

    lambda_control = 0.10,
    lambda_treatment = 0.07,

    censor_max = 24,

    looks = c(0.25,0.50,0.75,1.00),

    boundaries
){


  time_control <- rexp(
    n_control,
    rate = lambda_control
  )

  time_treatment <- rexp(
    n_treatment,
    rate = lambda_treatment
  )


  censor_control <- runif(
    n_control,
    min = 0,
    max = censor_max
  )

  censor_treatment <- runif(
    n_treatment,
    min = 0,
    max = censor_max
  )


  obs_control <- pmin(
    time_control,
    censor_control
  )

  obs_treatment <- pmin(
    time_treatment,
    censor_treatment
  )


  event_control <- ifelse(
    time_control <= censor_control,
    1,
    0
  )

  event_treatment <- ifelse(
    time_treatment <= censor_treatment,
    1,
    0
  )


  trial_data <- data.frame(

    time = c(
      obs_control,
      obs_treatment
    ),

    event = c(
      event_control,
      event_treatment
    ),

    trt = c(
      rep(0, n_control),
      rep(1, n_treatment)
    )
  )


  trial_data <- trial_data[
    sample(1:nrow(trial_data)),
  ]


  stop_trial <- FALSE
  stop_look  <- length(looks)
  final_z    <- NA
  n_used     <- nrow(trial_data)



  for(i in 1:length(looks)){


    n_current <- floor(
      nrow(trial_data) * looks[i]
    )


    interim_data <- trial_data[
      1:n_current,
    ]

    cox_fit <- coxph(
      Surv(time, event) ~ trt,
      data = interim_data
    )

    z_stat <- summary(
      cox_fit
    )$coef[,"z"]

    final_z <- z_stat

    if(abs(z_stat) >= boundaries[i]){

      stop_trial <- TRUE

      stop_look <- i

      n_used <- n_current

      break
    }
  }


  return(

    data.frame(

      reject_null = stop_trial,

      stop_look = stop_look,

      final_z = final_z,

      sample_size = n_used
    )
  )
}

one_trial <- simulate_trial(
  boundaries = boundaries
)

print(one_trial)


set.seed(123)

nsim <- 1000

results <- replicate(

  nsim,

  simulate_trial(
    boundaries = boundaries
  ),

  simplify = FALSE
)


results_df <- do.call(
  rbind,
  results
)

head(results_df)


power_estimate <- mean(
  results_df$reject_null
)

cat("\nEstimated Power =", power_estimate, "\n")


early_stop_prob <- mean(
  results_df$stop_look < 4
)

cat(
  "\nEarly Stopping Probability =",
  early_stop_prob,
  "\n"
)

avg_sample_size <- mean(
  results_df$sample_size
)

cat(
  "\nAverage Sample Size =",
  avg_sample_size,
  "\n"
)

avg_z <- mean(
  results_df$final_z
)

cat(
  "\nAverage Final Z Statistic =",
  avg_z,
  "\n"
)


summary(results_df)


ggplot(
  results_df,
  aes(x = final_z)
) +

  geom_histogram(
    bins = 30
  ) +

  labs(
    title = "Distribution of Final Z Statistics",
    x = "Final Z Statistic",
    y = "Frequency"
  )


ggplot(
  results_df,
  aes(x = factor(stop_look))
) +

  geom_bar() +

  labs(
    title = "Distribution of Stopping Looks",
    x = "Interim Look",
    y = "Count"
  )


# To estimate Type I error:
# rerun simulation with
#
# lambda_control  = 0.10
# lambda_treatment = 0.10
#
# because HR = 1 under H0






# PHASE 8
# Dropout & Missing Data in Sequential Survival Trial


library(survival)
library(survminer)
library(gsDesign)
library(ggplot2)
library(dplyr)

# STEP 2 — Simulation Parameters

set.seed(123)

n_control   <- 150
n_treatment <- 150

lambda_control   <- 0.10
lambda_treatment <- 0.07

dropout_rate <- 0.20

censor_max <- 24


# STEP 3 — Generate True Survival Times


time_control <- rexp(
  n_control,
  rate = lambda_control
)

time_treatment <- rexp(
  n_treatment,
  rate = lambda_treatment
)
# STEP 4 — Generate Administrative Censoring


censor_control <- runif(
  n_control,
  min = 0,
  max = censor_max
)

censor_treatment <- runif(
  n_treatment,
  min = 0,
  max = censor_max
)

# STEP 5 — Generate Dropout Indicators

dropout_control <- rbinom(
  n_control,
  size = 1,
  prob = dropout_rate
)

dropout_treatment <- rbinom(
  n_treatment,
  size = 1,
  prob = dropout_rate
)


# STEP 6 — Generate Dropout Times

dropout_time_control <- ifelse(

  dropout_control == 1,

  runif(
    n_control,
    min = 0,
    max = 18
  ),

  Inf
)

dropout_time_treatment <- ifelse(

  dropout_treatment == 1,

  runif(
    n_treatment,
    min = 0,
    max = 18
  ),

  Inf
)


# STEP 7 — Observed Follow-Up Times


obs_control <- pmin(
  time_control,
  censor_control,
  dropout_time_control
)

obs_treatment <- pmin(
  time_treatment,
  censor_treatment,
  dropout_time_treatment
)

# STEP 8 — Event Indicators
event_control <- ifelse(

  time_control <= censor_control &
  time_control <= dropout_time_control,

  1,
  0
)

event_treatment <- ifelse(

  time_treatment <= censor_treatment &
  time_treatment <= dropout_time_treatment,

  1,
  0
)

# STEP 9 — Create Trial Dataset

trial_data <- data.frame(

  time = c(
    obs_control,
    obs_treatment
  ),

  event = c(
    event_control,
    event_treatment
  ),

  trt = c(
    rep(0, n_control),
    rep(1, n_treatment)
  ),

  dropout = c(
    dropout_control,
    dropout_treatment
  )
)
# STEP 10 — Summary of Dropout


cat("\n==============================\n")
cat("Dropout Summary\n")
cat("==============================\n")

table(trial_data$dropout)

cat("\nDropout Rate = ",
    mean(trial_data$dropout),
    "\n")


# STEP 11 — Event Summary

cat("\n==============================\n")
cat("Event Summary\n")
cat("==============================\n")

table(trial_data$event)

# STEP 12 — Kaplan-Meier Analysis


km_fit <- survfit(
  Surv(time, event) ~ trt,
  data = trial_data
)


# STEP 13 — KM Plot


ggsurvplot(
  km_fit,
  data = trial_data,

  risk.table = TRUE,
  conf.int = TRUE,
  pval = TRUE,

  legend.title = "Treatment",
  legend.labs = c(
    "Control",
    "Treatment"
  ),

  title = "Kaplan-Meier Curves with Dropout",

  xlab = "Time",
  ylab = "Survival Probability"
)


# STEP 14 — Cox Proportional Hazards Model


cox_fit <- coxph(
  Surv(time, event) ~ trt,
  data = trial_data
)

# STEP 15 — Cox Model Summary
summary(cox_fit)


# STEP 16 — Hazard Ratio

hr <- exp(
  coef(cox_fit)
)

cat("\n==============================\n")
cat("Hazard Ratio\n")
cat("==============================\n")

print(hr)


# STEP 17 — Confidence Interval


ci <- exp(
  confint(cox_fit)
)

cat("\n==============================\n")
cat("95% Confidence Interval\n")
cat("==============================\n")

print(ci)


# STEP 18 — PH Assumption Check


ph_test <- cox.zph(cox_fit)

cat("\n==============================\n")
cat("PH Assumption Test\n")
cat("==============================\n")

print(ph_test)


# STEP 19 — PH Diagnostic Plot


plot(ph_test)


# STEP 20 — Forest Plot


ggforest(cox_fit)


# STEP 21 — Sequential Design

design_of <- gsDesign(
  k = 4,
  test.type = 2,
  alpha = 0.05,
  beta = 0.20,
  sfu = "OF"
)


# STEP 22 — Sequential Boundaries


boundaries <- design_of$upper$bound

print(boundaries)


# STEP 23 — Interim Looks


looks <- c(
  0.25,
  0.50,
  0.75,
  1.00
)

# STEP 24 — Sequential Monitoring


trial_data <- trial_data[
  sample(1:nrow(trial_data)),
]

stop_trial <- FALSE

stop_look <- NA

final_z <- NA


# STEP 25 — Interim Analysis Loop


for(i in 1:length(looks)){

    # Current Data Size


  n_current <- floor(
    nrow(trial_data) * looks[i]
  )

  # Interim Dataset


  interim_data <- trial_data[
    1:n_current,
  ]


  # Cox Model at Interim


  cox_interim <- coxph(
    Surv(time, event) ~ trt,
    data = interim_data
  )

 
  # Z Statistic
 

  z_stat <- summary(
    cox_interim
  )$coef[,"z"]

  final_z <- z_stat

  # Print Interim Results


  cat("\n================================\n")
  cat("Look:", i, "\n")
  cat("Information Fraction:",
      looks[i], "\n")
  cat("Z-statistic:",
      round(z_stat,3), "\n")
  cat("Boundary:",
      round(boundaries[i],3), "\n")

  # Stopping Rule


  if(abs(z_stat) >= boundaries[i]){

    stop_trial <- TRUE

    stop_look <- i

    cat("Trial Stopped Early\n")

    break
  }
}


# STEP 26 — Final Decision

cat("\n================================\n")

if(stop_trial){

  cat("Final Decision: Early Stop\n")
  cat("Stopping Look:",
      stop_look, "\n")

} else {

  cat("Final Decision: Continue to Final Analysis\n")
}


# STEP 27 — Missing Covariate Example


age <- rnorm(
  nrow(trial_data),
  mean = 60,
  sd = 10
)


# Introduce Missing Age Values


missing_index <- sample(
  1:nrow(trial_data),
  size = 30
)

age[missing_index] <- NA

trial_data$age <- age

# STEP 28 — Missing Data Summary


cat("\n================================\n")
cat("Missing Age Values\n")
cat("================================\n")

sum(is.na(trial_data$age))


# STEP 29 — Complete Case Dataset


complete_case_data <- na.omit(
  trial_data
)


# STEP 30 — Cox Model with Age

cox_age <- coxph(
  Surv(time, event) ~ trt + age,
  data = complete_case_data
)


# STEP 31 — Final Model Summary


summary(cox_age)

