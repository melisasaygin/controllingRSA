#install the required packages
#install.packages("lme4")
#install.packages("ggplot2")
#install.packages("lmerTest")
#install.packages("performance")
#install.packages("dplyr")
#install.packages("HLMdiag")
#install.packages("jtools")
#install.packages("r2mlm")
#devtools::install_github("mkshaw/r2mlm")
#install.packages("tidyverse")

#importing the required packages
library(lme4)
library(ggplot2)
library(lmerTest)
library(performance)
library(dplyr)
library(readxl)
library(HLMdiag)
library(jtools)
library(r2mlm)
library(tidyverse)

getwd()
final_df_sitting <- read_excel("/Users/melisasaygin/Documents/Validation of RSA Paper Code (for GitHub)/final_df_sitting.xlsx")

#checking data classes (types) and correcting as needed
# _45 subscripts near the RSA metrics just refers to the fact that while making the relevant regressions for the approaches, the points 4.5SDs away from the mean (y-axis) were not included to create more robust regressions
sapply(final_df_sitting[, c("Average_IBI_msec", "RSA0_msec", "BMI", "IPAQ_SF", 
                            "Age", "Sex", "Current_Smoking", "participant_id",
                            "delta_RSA0_over_Ttot_laboratory_45", "delta_RSA0_over_Ttot_ambulatory_45", #
                            "RSA0_residual_45", "RSA0_over_Vt")], class) 

final_df_sitting$participant_id <- as.factor(final_df_sitting$participant_id)
final_df_sitting$Sex <- factor(final_df_sitting$Sex) #make sex into a factor variable
final_df_sitting$Current_Smoking <- factor(final_df_sitting$Current_Smoking, levels = c(0, 1), labels = c("Non-Smoker", "Smoker"))
#make numeric current_smoking and sex variables to later be able to use the r2mlm package
final_df_sitting$Current_Smoking_numeric <- ifelse(final_df_sitting$Current_Smoking == "Smoker", 1, 0)  # 0 becomes non-smoker, 1 becomes smoker
final_df_sitting$Sex_numeric <- ifelse(final_df_sitting$Sex == "Male", 1, 0) # 0 becomes female, 1 becomes male

#re-check the data types
sapply(final_df_sitting[, c("Average_IBI_msec", "RSA0_msec", "BMI", "IPAQ_SF", 
                            "Age", "Sex", "Current_Smoking", "participant_id",
                            "delta_RSA0_over_Ttot_laboratory_45", "delta_RSA0_over_Ttot_ambulatory_45",
                            "RSA0_residual_45", "RSA0_over_Vt", "Current_Smoking_numeric", "Sex_numeric")], class)

# remove outliers (values deviating more than 4.5 SDs for IBI and RSA metrics)
final_df_sitting <- final_df_sitting %>%
  group_by(participant_id) %>%
  filter(
    abs(RSA0_msec - mean(RSA0_msec)) / sd(RSA0_msec) <= 4.5,  # Remove RSA0_msec outliers
    abs(Average_IBI_msec - mean(Average_IBI_msec)) / sd(Average_IBI_msec) <= 4.5,  # Remove IBI outliers
    abs(delta_RSA0_over_Ttot_laboratory_45 - mean(delta_RSA0_over_Ttot_laboratory_45)) / sd(delta_RSA0_over_Ttot_laboratory_45) <= 4.5,  # Remove delta_RSA0_over_Ttot_laboratory_45 outliers
    abs(delta_RSA0_over_Ttot_ambulatory_45 - mean(delta_RSA0_over_Ttot_ambulatory_45)) / sd(delta_RSA0_over_Ttot_ambulatory_45) <= 4.5,  # Remove delta_RSA0_over_Ttot_ambulatory_45 outliers
    abs(RSA0_residual_45 - mean(RSA0_residual_45)) / sd(RSA0_residual_45) <= 4.5,  # Remove RSA0_residual outliers
    abs(RSA0_over_Vt - mean(RSA0_over_Vt)) / sd(RSA0_over_Vt) <= 4.5  # Remove RSA0_over_Vt outliers
  ) %>%
  ungroup()

#checking how many epochs does each participant have
epoch_counts <- final_df_sitting %>%
  count(participant_id, name = "n_epochs")

print(epoch_counts)

# get mean and standard deviation of epochs across participants
mean_epochs <- mean(epoch_counts$n_epochs)
sd_epochs <- sd(epoch_counts$n_epochs)

cat("Mean number of epochs per participant:", round(mean_epochs, 2), "\n")
cat("Standard deviation of epochs per participant:", round(sd_epochs, 2), "\n")

## MAKING THE WITHIN-PERSON CORRELATION MATRIX
library(tidyverse)
library(reshape2)

# defining a list of variables to correlate
vars_to_correlate <- c(
  "RSA0_msec", 
  "delta_RSA0_over_Ttot_laboratory_45", 
  "delta_RSA0_over_Ttot_ambulatory_45", 
  "RSA0_residual_45", 
  "RSA0_over_Vt", 
  "Average_IBI_msec", 
  "Respiration_Rate_bpm",             
  "Calibrated_Tidal_Volume"
)

# custom labels for each variables to map these labels
var_labels <- c(
  "RSA0_msec" = "RSA (Approach 0)",
  "delta_RSA0_over_Ttot_laboratory_45" = "(ΔRSA/Vt)_paced (Approach 1)",
  "delta_RSA0_over_Ttot_ambulatory_45" = "(ΔRSA/Vt)_ambulatory (Approach 2)",
  "RSA0_residual_45" = "RSA_residual (Approach 3)",
  "RSA0_over_Vt" = "RSA/Vt (Approach 4)",
  "Average_IBI_msec" = "Heart Period",
  "Respiration_Rate_bpm" = "Respiration Rate",
  "Calibrated_Tidal_Volume" = "Tidal Volume"
)

# extracts the list of unique participant ids
participants <- unique(final_df_sitting$participant_id)

# initialize a list to store the individual correlation matrix of each participant
corr_matrices <- list()

# start a loop that goes through each participant to compute their correlation matrix
for (p in participants) {
  df_sub <- final_df_sitting %>% filter(participant_id == p) #create dataframe for participant
  df_sub <- df_sub %>% select(all_of(vars_to_correlate)) %>% na.omit() #select the variables of interest for the participant
  corr <- cor(df_sub, method = "pearson")  # calculate correlation matrix for this participant's data
  corr_matrices[[as.character(p)]] <- corr #save the resulting correlation matrix into a list
  
}

# add all the correlation matrices and divides by the number of matrices
mean_corr_matrix <- Reduce("+", corr_matrices) / length(corr_matrices)

# convert square to long format dataframe
melted_corr <- melt(mean_corr_matrix)
melted_corr <- melted_corr %>%
  filter(Var1 != Var2) %>% #remove self-correlations
  filter(as.integer(factor(Var1, levels = vars_to_correlate)) > as.integer(factor(Var2, levels = vars_to_correlate)))
melted_corr 

# plot the lower-left triangle heatmap with custom variable names and box positions on the left and bottom
ggplot(melted_corr, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", value)), size = 3) +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limit = c(-1, 1),
    name = "Mean of within-person\nPearson's r"
  ) +
  scale_x_discrete(limits = vars_to_correlate, labels = var_labels) +  
  scale_y_discrete(limits = rev(vars_to_correlate), labels = var_labels) +  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(hjust = 1)  
  ) +
  coord_fixed()  


ggsave(
  filename = "corr_matrix_within.jpg",
  plot = last_plot(),            
  width = 8.5, height = 8.5,         
  dpi = 600,                     # High resolution
  units = "in",                 
  device = "jpeg"
)

#automatically make all the by-participant plots for each pair of variables in the correlation matrix (for Supplementary material)

library(tidyverse)

# Define your correlation variable list and final dataset
vars_to_correlate <- c(
  "RSA0_msec", 
  "delta_RSA0_over_Ttot_laboratory_45", 
  "delta_RSA0_over_Ttot_ambulatory_45", 
  "RSA0_residual_45", 
  "RSA0_over_Vt", 
  "Average_IBI_msec", 
  "Respiration_Rate_bpm",             
  "Calibrated_Tidal_Volume"
)

# Optional: custom labels for prettier plot titles
var_labels <- c(
  "RSA0_msec" = "RSA (Approach 0)",
  "delta_RSA0_over_Ttot_laboratory_45" = "(ΔRSA/Vt)_paced (Approach 1)",
  "delta_RSA0_over_Ttot_ambulatory_45" = "(ΔRSA/Vt)_ambulatory (Approach 2)",
  "RSA0_residual_45" = "RSA_residual (Approach 3)",
  "RSA0_over_Vt" = "RSA/Vt (Approach 4)",
  "Average_IBI_msec" = "Heart Period",
  "Respiration_Rate_bpm" = "Respiration Rate",
  "Calibrated_Tidal_Volume" = "Tidal Volume"
)

# Loop over all unique variable pairs (lower triangle, no self-pairs)
for (i in 1:(length(vars_to_correlate)-1)) {
  for (j in (i+1):length(vars_to_correlate)) {
    
    var_x <- vars_to_correlate[i]
    var_y <- vars_to_correlate[j]
    
    plot_data <- final_df_sitting %>%
      select(participant_id, all_of(c(var_x, var_y))) %>%
      drop_na()
    
    p <- ggplot(plot_data, aes_string(x = var_x, y = var_y)) +
      geom_point(alpha = 0.1, color = "darkblue") +
      geom_smooth(method = "lm", se = TRUE, color = "red", fill = "red", alpha = 0.3, linewidth = 0.8) +
      facet_wrap(~ participant_id, scales = "free") +
      labs(
        title = paste0(
          var_labels[[var_y]], " vs ", var_labels[[var_x]], " by Participant"
        ),
        x = var_labels[[var_x]],
        y = var_labels[[var_y]]
      ) +
      theme_classic(base_size = 12)
    
    filename <- paste0(
      "/Users/melisasaygin/Documents/Validation of RSA Paper Code (for GitHub)/",
      var_y, "_vs_", var_x, "_by_participant.png"
    )
    
    ggsave(
      filename = filename,
      plot = p,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
}



####################################################################################
##### RESEARCH QUESTION 2 (PREDICTING IBI FROM RSA METRICS) FOR SITTING POSTURE ####
####################################################################################

IBI_from_covariates_sitting_cwc <- lmer(  # just the covariates
  Average_IBI_msec ~ BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 | participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(IBI_from_covariates_sitting_cwc)
r2mlm(IBI_from_covariates_sitting_cwc) #no within variance explained as there is no level 1 variable yet

###predicting from uncontrolled RSA (Approach 0) with no respiratory controlling
final_df_sitting <- final_df_sitting %>%   #centering within clusters (individuals), also calculating the individual mean
  group_by(participant_id) %>%
  mutate(
    RSA0_msec_mean = mean(RSA0_msec),
    RSA0_msec_cwc = RSA0_msec - RSA0_msec_mean  # CWC = centered within cluster (i.e., participant)
  ) %>%
  ungroup()


IBI_from_unadjusted_RSA_sitting_cwc <- lmer(  #introduce the centered within cluster RSA metric (level1) along with the cluster means (level2)
  Average_IBI_msec ~ RSA0_msec_cwc + RSA0_msec_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 + RSA0_msec_cwc | participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(IBI_from_unadjusted_RSA_sitting_cwc)
summ(IBI_from_unadjusted_RSA_sitting_cwc)
r2mlm(IBI_from_unadjusted_RSA_sitting_cwc)

###predicting from the metric of Approach 1 - ΔRSA0/Vt via baseline regression made with paced breathing

final_df_sitting <- final_df_sitting %>%   #doing group mean centering
  group_by(participant_id) %>%
  mutate(
    delta_RSA0_over_Ttot_lab45_mean = mean(delta_RSA0_over_Ttot_laboratory_45),
    delta_RSA0_over_Ttot_lab45_cwc = delta_RSA0_over_Ttot_laboratory_45 - delta_RSA0_over_Ttot_lab45_mean
  ) %>%
  ungroup()

IBI_from_Approach1_RSA_sitting_cwc <- lmer(
  Average_IBI_msec ~ delta_RSA0_over_Ttot_lab45_cwc + delta_RSA0_over_Ttot_lab45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 + delta_RSA0_over_Ttot_lab45_cwc | participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(IBI_from_Approach1_RSA_sitting_cwc)
summ(IBI_from_Approach1_RSA_sitting_cwc)
r2mlm(IBI_from_Approach1_RSA_sitting_cwc)

# histogram of residuals
hist(res, breaks = 30, main = "Histogram of Residuals", xlab = "Residuals")

#predicting from the metric of Approach 2 - ΔRSA0/Vt via baseline regression made with ambulatory data
final_df_sitting <- final_df_sitting %>%    #group mean centering
  group_by(participant_id) %>%
  mutate(
    delta_RSA0_over_Ttot_amb45_mean = mean(delta_RSA0_over_Ttot_ambulatory_45),
    delta_RSA0_over_Ttot_amb45_cwc = delta_RSA0_over_Ttot_ambulatory_45 - delta_RSA0_over_Ttot_amb45_mean
  ) %>%
  ungroup()

IBI_from_Approach2_RSA_sitting_cwc <- lmer(
  Average_IBI_msec ~ delta_RSA0_over_Ttot_amb45_cwc + delta_RSA0_over_Ttot_amb45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 + delta_RSA0_over_Ttot_amb45_cwc | participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(IBI_from_Approach2_RSA_sitting_cwc)
summ(IBI_from_Approach2_RSA_sitting_cwc)
r2mlm(IBI_from_Approach2_RSA_sitting_cwc)

#predicting from the metric of Approach 3 - residual RSA

final_df_sitting <- final_df_sitting %>%
  group_by(participant_id) %>%
  mutate(
    RSA0_residual_45_mean = mean(RSA0_residual_45),
    RSA0_residual_45_cwc = RSA0_residual_45 - RSA0_residual_45_mean
  ) %>%
  ungroup()

IBI_from_Approach3_RSA_sitting_cwc <- lmer(
  Average_IBI_msec ~ RSA0_residual_45_cwc + RSA0_residual_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 + RSA0_residual_45_cwc| participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(IBI_from_Approach3_RSA_sitting_cwc)
summ(IBI_from_Approach3_RSA_sitting_cwc)
r2mlm(IBI_from_Approach3_RSA_sitting_cwc)

#predicting from the metric of Approach 4 - RSA/Vt

final_df_sitting <- final_df_sitting %>%   #group mean centering
  group_by(participant_id) %>%
  mutate(
    RSA0_over_Vt_mean = mean(RSA0_over_Vt),
    RSA0_over_Vt_cwc = RSA0_over_Vt - RSA0_over_Vt_mean
  ) %>%
  ungroup()


IBI_from_Approach4_RSA_sitting_cwc <- lmer(
  Average_IBI_msec ~ RSA0_over_Vt_cwc + RSA0_over_Vt_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 + RSA0_over_Vt_cwc| participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(IBI_from_Approach4_RSA_sitting_cwc)
summ(IBI_from_Approach4_RSA_sitting_cwc)
r2mlm(IBI_from_Approach4_RSA_sitting_cwc)

#predicting from the metric of Approach 5 - EXPLORATORY - RSA/Ttot

#make the exploratory (Approach "5" variable)
final_df_sitting <- final_df_sitting %>%
  mutate(RSA0_over_Ttot = RSA0_msec / Ttot)

final_df_sitting <- final_df_sitting %>%
  group_by(participant_id) %>%
  mutate(
    RSA0_over_Ttot_mean = mean(RSA0_over_Ttot, na.rm = TRUE),
    RSA0_over_Ttot_cwc = RSA0_over_Ttot - RSA0_over_Ttot_mean
  ) %>%
  ungroup()

IBI_from_Approach5_RSA_sitting_cwc <- lmer(
  Average_IBI_msec ~ RSA0_over_Ttot_cwc + RSA0_over_Ttot_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 + RSA0_over_Ttot_cwc | participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)


summary(IBI_from_Approach5_RSA_sitting_cwc)
summ(IBI_from_Approach5_RSA_sitting_cwc)
r2mlm(IBI_from_Approach5_RSA_sitting_cwc)

#Variance Decomposition into within and between cluster (individual) levels using the r2mlm package
par(mar = c(6.75, 10.5, 2.625, 10.5))
r2mlm(IBI_from_unadjusted_RSA_sitting_cwc)
r2mlm(IBI_from_Approach1_RSA_sitting_cwc)
r2mlm(IBI_from_Approach2_RSA_sitting_cwc)
r2mlm(IBI_from_Approach3_RSA_sitting_cwc)
r2mlm(IBI_from_Approach4_RSA_sitting_cwc)
r2mlm(IBI_from_Approach5_RSA_sitting_cwc)


####################################################################################
##### RESEARCH QUESTION 3 (PREDICTING PSYCHOLOGICAL STATES FROM RSA METRICS) ####
####################################################################################
#importing the required packages
library(lme4)
library(ggplot2)
library(lmerTest)
library(performance)
library(dplyr)
library(readxl)
library(HLMdiag)
library(jtools)

getwd()
df_ema <- read_excel("/Users/melisasaygin/Documents/Validation of RSA Paper Code (for GitHub)/df_ema_merged.xlsx")
#View(df_ema)

#checking data classes (types) and correcting as needed

sapply(df_ema[, c("Tot_Stress", "Pos_Aff", "Neg_Aff", "Safety", "BMI", 
                  "IPAQ_SF", "Age", "Sex", "Current_Smoking", "participant_id",
                  "RSA0_msec", "delta_RSA0_over_Ttot_ambulatory_45", 
                  "delta_RSA0_over_Ttot_laboratory_45", "RSA0_residual_45", 
                  "RSA0_over_Vt")], class)

df_ema$participant_id <- as.factor(df_ema$participant_id)
df_ema$Sex <- factor(df_ema$Sex) #make sex into a factor variable
df_ema$Current_Smoking <- factor(df_ema$Current_Smoking, levels = c(0, 1), labels = c("Non-Smoker", "Smoker"))
#make numeric variables for current_smoking and sex variables to later be able to use the r2mlm package
df_ema$Current_Smoking_numeric <- ifelse(df_ema$Current_Smoking == "Smoker", 1, 0)  #0 becomes non-smoker, 1 becomes smoker
df_ema$Sex_numeric <- ifelse(df_ema$Sex == "Male", 1, 0) # 0 becomes female, 1 becomes male

#re-check the data types
sapply(df_ema[, c("Tot_Stress", "Pos_Aff", "Neg_Aff", "Safety", "BMI", 
                            "IPAQ_SF", "Age", "Sex", "Current_Smoking", "participant_id",
                            "RSA0_msec", "delta_RSA0_over_Ttot_ambulatory_45", 
                            "delta_RSA0_over_Ttot_laboratory_45", "RSA0_residual_45", 
                            "RSA0_over_Vt", "Current_Smoking_numeric", "Sex_numeric")], class)
posture_counts <- table(df_ema$Majority_posture_5min)

# calculate the percentage for each posture
posture_percentage <- prop.table(posture_counts) * 100
posture_summary <- data.frame(
  Posture = names(posture_counts),
  Count = as.integer(posture_counts),
  Percentage = round(posture_percentage, 2)
)

print(posture_summary)

#filter to only keep the epochs with Sitting posture
df_ema_sitting <- df_ema[df_ema$Majority_posture_5min == "Sitting", ]

# counting responses per participant
participant_counts <- table(df_ema_sitting$participant_id)
participant_counts

# only keep participants with at least 8 responses
valid_participants <- names(participant_counts[participant_counts >= 8])
df_ema_filtered <- df_ema_sitting[df_ema_sitting$participant_id %in% valid_participants, ]

#Visualization of data ditribution for Safety

# add an ordered factor to both datasets based on SD
summary_safety <- df_ema_filtered %>%
  group_by(participant_id) %>%
  summarise(mean_safety = mean(Safety, na.rm = TRUE),
            sd_safety = sd(Safety, na.rm = TRUE)) %>%
  mutate(participant_id_ordered = reorder(participant_id, sd_safety))

df_ema_filtered <- df_ema_filtered %>%
  left_join(summary_safety %>% select(participant_id, sd_safety), by = "participant_id") %>%
  mutate(participant_id_ordered = reorder(participant_id, sd_safety))

ggplot() +
  # subtle individual points
  geom_jitter(data = df_ema_filtered, 
              aes(x = participant_id_ordered, y = Safety), 
              width = 0.2, alpha = 0.3, color = "grey40", size = 1) +
  
  # prominent mean points
  geom_point(data = summary_safety, 
             aes(x = participant_id_ordered, y = mean_safety), 
             color = "black", size = 2.5) +
  
  # prominent SD error bars
  geom_errorbar(data = summary_safety, 
                aes(x = participant_id_ordered, ymin = mean_safety - sd_safety, ymax = mean_safety + sd_safety), 
                width = 0.2, color = "black", linewidth = 1) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  labs(title = "", 
       x = "Participant", y = "Safety at the moment")

# save the plot as a .jpg file
ggsave("safety_distribution_byparticipant.jpg", 
       plot = last_plot(), 
       width = 10, 
       height = 6, 
       dpi = 300)


#visualization of data distribution for Stress
summary_stress <- df_ema_filtered %>%
  group_by(participant_id) %>%
  summarise(mean_stress = mean(Tot_Stress, na.rm = TRUE),
            sd_stress = sd(Tot_Stress, na.rm = TRUE)) %>%
  mutate(participant_id_ordered = reorder(participant_id, sd_stress))

df_ema_filtered <- df_ema_filtered %>%
  left_join(summary_stress %>% select(participant_id, sd_stress), by = "participant_id") %>%
  mutate(participant_id_ordered_stress = reorder(participant_id, sd_stress))

ggplot() +
  geom_jitter(data = df_ema_filtered, 
              aes(x = participant_id_ordered_stress, y = Tot_Stress), 
              width = 0.2, alpha = 0.3, color = "grey40", size = 1) +
  geom_point(data = summary_stress, 
             aes(x = participant_id_ordered, y = mean_stress), 
             color = "black", size = 2.5) +
  geom_errorbar(data = summary_stress, 
                aes(x = participant_id_ordered, ymin = mean_stress - sd_stress, ymax = mean_stress + sd_stress), 
                width = 0.2, color = "black", linewidth = 1) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  labs(title = "",
       x = "Participant", y = "Psychological Stress at the moment")
# save the plot as a .jpg file
ggsave("stress_distribution_byparticipant.jpg", 
       plot = last_plot(), 
       width = 10, 
       height = 6, 
       dpi = 300)


#visualization of data distribution for positive affect
summary_pos_aff <- df_ema_filtered %>%
  group_by(participant_id) %>%
  summarise(mean_pos_aff = mean(Pos_Aff, na.rm = TRUE),
            sd_pos_aff = sd(Pos_Aff, na.rm = TRUE)) %>%
  mutate(participant_id_ordered = reorder(participant_id, sd_pos_aff))

df_ema_filtered <- df_ema_filtered %>%
  left_join(summary_pos_aff %>% select(participant_id, sd_pos_aff), by = "participant_id") %>%
  mutate(participant_id_ordered_pos_aff = reorder(participant_id, sd_pos_aff))

ggplot() +
  geom_jitter(data = df_ema_filtered, 
              aes(x = participant_id_ordered_pos_aff, y = Pos_Aff), 
              width = 0.2, alpha = 0.3, color = "grey40", size = 1) +
  geom_point(data = summary_pos_aff, 
             aes(x = participant_id_ordered, y = mean_pos_aff), 
             color = "black", size = 2.5) +
  geom_errorbar(data = summary_pos_aff, 
                aes(x = participant_id_ordered, ymin = mean_pos_aff - sd_pos_aff, ymax = mean_pos_aff + sd_pos_aff), 
                width = 0.2, color = "black", linewidth = 1) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  labs(title = "",
       x = "Participant", y = "Positive Affect at the moment")

# save the plot as a .jpg file
ggsave("pos_aff_distribution_byparticipant.jpg", 
       plot = last_plot(), 
       width = 10, 
       height = 6, 
       dpi = 300)

#visualization for the data distribution of negative affect
summary_neg_aff <- df_ema_filtered %>%
  group_by(participant_id) %>%
  summarise(mean_neg_aff = mean(Neg_Aff, na.rm = TRUE),
            sd_neg_aff = sd(Neg_Aff, na.rm = TRUE)) %>%
  mutate(participant_id_ordered = reorder(participant_id, sd_neg_aff))

df_ema_filtered <- df_ema_filtered %>%
  left_join(summary_neg_aff %>% select(participant_id, sd_neg_aff), by = "participant_id") %>%
  mutate(participant_id_ordered_neg_aff = reorder(participant_id, sd_neg_aff))

ggplot() +
  geom_jitter(data = df_ema_filtered, 
              aes(x = participant_id_ordered_neg_aff, y = Neg_Aff), 
              width = 0.2, alpha = 0.3, color = "grey40", size = 1) +
  geom_point(data = summary_neg_aff, 
             aes(x = participant_id_ordered, y = mean_neg_aff), 
             color = "black", size = 2.5) +
  geom_errorbar(data = summary_neg_aff, 
                aes(x = participant_id_ordered, ymin = mean_neg_aff - sd_neg_aff, ymax = mean_neg_aff + sd_neg_aff), 
                width = 0.2, color = "black", linewidth = 1) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  labs(title = " ",
       x = "Participant", y = "Negative Affect at the moment")

# save the plot as a .jpg file
ggsave("neg_aff_distribution_byparticipant.jpg", 
       plot = last_plot(), 
       width = 10, 
       height = 6, 
       dpi = 300)

install.packages("jpeg")
library(gridExtra)
library(jpeg)
library(grid)
library(gridExtra)
library(grid)
library(jpeg)

# read images
image1 <- readJPEG("safety_distribution_byparticipant.jpg")
image2 <- readJPEG("stress_distribution_byparticipant.jpg")
image3 <- readJPEG("pos_aff_distribution_byparticipant.jpg")
image4 <- readJPEG("neg_aff_distribution_byparticipant.jpg")

# convert to grobs
grob1 <- rasterGrob(image1, interpolate = TRUE)
grob2 <- rasterGrob(image2, interpolate = TRUE)
grob3 <- rasterGrob(image3, interpolate = TRUE)
grob4 <- rasterGrob(image4, interpolate = TRUE)

# arrange into a 2x2 layout
combined_plot <- grid.arrange(grob1, grob2, grob3, grob4, ncol = 2)

# save combined psychological state distributions plot
ggsave("combined_distribution_byparticipant.jpg", 
       plot = combined_plot, 
       width = 12, height = 10, dpi = 300)



###STRESS###

###predicting from unadjusted RSA with no respiratory controlling (Approach 0)
df_ema_filtered <- df_ema_filtered %>%  
  group_by(participant_id) %>%
  mutate(
    RSA0_msec_mean = mean(RSA0_msec),       # person's mean
    RSA0_msec_cwc = RSA0_msec - RSA0_msec_mean     # variations of the individual
  ) %>%
  ungroup()

Stress_from_unadjusted_RSA_cwc <- lmer(
  Tot_Stress ~ RSA0_msec_cwc + RSA0_msec_mean +  BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1  | participant_id),  
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

#from Approach 1 - delta_RSA0_over_Ttot_laboratory_45
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    delta_RSA0_over_Ttot_laboratory_45_mean = mean(delta_RSA0_over_Ttot_laboratory_45),  # person's mean
    delta_RSA0_over_Ttot_laboratory_45_cwc = delta_RSA0_over_Ttot_laboratory_45 - delta_RSA0_over_Ttot_laboratory_45_mean  # centered within person
  ) %>%
  ungroup()

Stress_from_Approach1_RSA_cwc <- lmer(
  Tot_Stress ~ delta_RSA0_over_Ttot_laboratory_45_cwc + delta_RSA0_over_Ttot_laboratory_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1  | participant_id),
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

#the model with also random slopes yields a singular fit, while the fixed slopes only model converges successfully


summary(Stress_from_Approach1_RSA_cwc)
summ(Stress_from_Approach1_RSA_cwc)
r2mlm(Stress_from_Approach1_RSA_cwc)

#from Approach 2 - delta_RSA0_over_Ttot_ambulatory_45
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    delta_RSA0_over_Ttot_ambulatory_45_mean = mean(delta_RSA0_over_Ttot_ambulatory_45),  # person's mean
    delta_RSA0_over_Ttot_ambulatory_45_cwc = delta_RSA0_over_Ttot_ambulatory_45 - delta_RSA0_over_Ttot_ambulatory_45_mean  # centered within person
  ) %>%
  ungroup()

Stress_from_Approach2_RSA_cwc <- lmer(
  Tot_Stress ~ delta_RSA0_over_Ttot_ambulatory_45_cwc + delta_RSA0_over_Ttot_ambulatory_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1  | participant_id),
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

#the model with also the random slopes yielded a singular fit, while the fixed slopes only model converges successfully

summary(Stress_from_Approach2_RSA_cwc)
summ(Stress_from_Approach2_RSA_cwc)
r2mlm(Stress_from_Approach2_RSA_cwc)

#from Approach 3 - RSA0_residual_45
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    RSA0_residual_45_mean = mean(RSA0_residual_45),  # person's mean
    RSA0_residual_45_cwc = RSA0_residual_45 - RSA0_residual_45_mean  # centered within person
  ) %>%
  ungroup()

Stress_from_Approach3_RSA_cwc <- lmer(
  Tot_Stress ~ RSA0_residual_45_cwc + RSA0_residual_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1 | participant_id),  
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Stress_from_Approach3_RSA_cwc)
summ(Stress_from_Approach3_RSA_cwc)
r2mlm(Stress_from_Approach3_RSA_cwc)

#from Approach 4 - RSA0_over_Vt
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    RSA0_over_Vt_mean = mean(RSA0_over_Vt),  # person's mean
    RSA0_over_Vt_cwc = RSA0_over_Vt - RSA0_over_Vt_mean    # centered within person
  ) %>%
  ungroup()

Stress_from_Approach4_RSA_cwc <- lmer(
  Tot_Stress ~ RSA0_over_Vt_cwc + RSA0_over_Vt_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1 | participant_id),  # random intercept only
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Stress_from_Approach4_RSA_cwc)
summ(Stress_from_Approach4_RSA_cwc)
r2mlm(Stress_from_Approach4_RSA_cwc)

r2mlm(Stress_from_unadjusted_RSA_cwc)
r2mlm(Stress_from_Approach1_RSA_cwc)
r2mlm(Stress_from_Approach2_RSA_cwc)
r2mlm(Stress_from_Approach3_RSA_cwc)
r2mlm(Stress_from_Approach4_RSA_cwc)

#Predicting stress from respiratory rate and tidal volume
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    TV_mean = mean(Calibrated_Tidal_Volume),
    TV_cwc = Calibrated_Tidal_Volume - TV_mean,
    RR_mean = mean(Respiration_Rate_bpm),
    RR_cwc = Respiration_Rate_bpm - RR_mean
  ) %>%
  ungroup()
  
Stress_from_RR_TV_cwc <- lmer(
  Tot_Stress ~ TV_cwc + TV_mean + RR_cwc + RR_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 | participant_id),
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

# Summary and model metrics
summary(Stress_from_RR_TV_cwc)
summ(Stress_from_RR_TV_cwc)
r2(Stress_from_RR_TV_cwc)
r2mlm(Stress_from_RR_TV_cwc)

###SAFETY###
###predicting from unadjusted RSA with no respiratory controlling
Safety_from_unadjusted_RSA_cwc <- lmer(
  Safety ~ RSA0_msec_cwc + RSA0_msec_mean +  BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 | participant_id),  
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)
#including a random slopes results in a singular fit while the fixed slopes only model converges

summary(Safety_from_unadjusted_RSA_cwc)  
summ(Safety_from_unadjusted_RSA_cwc)
r2mlm(Safety_from_unadjusted_RSA_cwc)

#from Approach 1 - delta_RSA0_over_Ttot_laboratory_45
Safety_from_Approach1_RSA_cwc <- lmer(
  Safety ~ delta_RSA0_over_Ttot_laboratory_45_cwc + delta_RSA0_over_Ttot_laboratory_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1  | participant_id),
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)
#including a random slopes variable results in a singular fit while the fixed slopes only model does converge

summary(Safety_from_Approach1_RSA_cwc)
summ(Safety_from_Approach1_RSA_cwc)
r2(Safety_from_Approach1_RSA_cwc)
r2mlm(Safety_from_unadjusted_RSA_cwc)
r2mlm(Safety_from_Approach1_RSA_cwc)

#from Approach 2 - delta_RSA0_over_Ttot_ambulatory_45
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    delta_RSA0_over_Ttot_ambulatory_45_mean = mean(delta_RSA0_over_Ttot_ambulatory_45),  # person's mean
    delta_RSA0_over_Ttot_ambulatory_45_cwc = delta_RSA0_over_Ttot_ambulatory_45 - delta_RSA0_over_Ttot_ambulatory_45_mean  # centered within person
  ) %>%
  ungroup()

Safety_from_Approach2_RSA_cwc <- lmer(
  Safety ~ delta_RSA0_over_Ttot_ambulatory_45_cwc + delta_RSA0_over_Ttot_ambulatory_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1  | participant_id),
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Safety_from_Approach2_RSA_cwc)
summ(Safety_from_Approach2_RSA_cwc)
r2(Safety_from_Approach2_RSA_cwc)
plot(Safety_from_Approach2_RSA_cwc)

#from Approach 3 - RSA0_residual_45
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    RSA0_residual_45_mean = mean(RSA0_residual_45),  # person's mean
    RSA0_residual_45_cwc = RSA0_residual_45 - RSA0_residual_45_mean  # centered within person
  ) %>%
  ungroup()

Safety_from_Approach3_RSA_cwc <- lmer(
  Safety ~ RSA0_residual_45_cwc + RSA0_residual_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1 | participant_id),  
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Safety_from_Approach3_RSA_cwc)
summ(Safety_from_Approach3_RSA_cwc)
r2(Safety_from_Approach3_RSA_cwc)
plot(Safety_from_Approach3_RSA_cwc)

#from Approach 4 - RSA0_over_Vt
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    RSA0_over_Vt_mean = mean(RSA0_over_Vt),  # person's mean
    RSA0_over_Vt_cwc = RSA0_over_Vt - RSA0_over_Vt_mean    # centered within person
  ) %>%
  ungroup()

Safety_from_Approach4_RSA_cwc <- lmer(
  Safety ~ RSA0_over_Vt_cwc + RSA0_over_Vt_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1 | participant_id),  # random intercept only
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Safety_from_Approach4_RSA_cwc)
summ(Safety_from_Approach4_RSA_cwc)

r2mlm(Safety_from_unadjusted_RSA_cwc)
r2mlm(Safety_from_Approach1_RSA_cwc)
r2mlm(Safety_from_Approach2_RSA_cwc)
r2mlm(Safety_from_Approach3_RSA_cwc)
r2mlm(Safety_from_Approach4_RSA_cwc)

#safety from respiration rate and tidal volume
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    TV_mean = mean(Calibrated_Tidal_Volume),
    TV_cwc = Calibrated_Tidal_Volume - TV_mean,
    RR_mean = mean(Respiration_Rate_bpm),
    RR_cwc = Respiration_Rate_bpm - RR_mean
  ) %>%
  ungroup()


Safety_from_RR_TV_cwc <- lmer(
  Safety ~ TV_cwc + TV_mean + RR_cwc + RR_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1  | participant_id),
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Safety_from_RR_TV_cwc)
r2(Safety_from_RR_TV_cwc)
r2mlm(Safety_from_RR_TV_cwc)

# predict from uncontrolled RSA alongside respiratory metrics
Safety_from_RSA_RR_TV <- lmer(
  Safety ~ RSA0_msec_cwc + RSA0_msec_mean +
    RR_cwc + RR_mean +
    TV_cwc + TV_mean +
    BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1 | participant_id),
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

# Model summary and R²
summary(Safety_from_RSA_RR_TV)
r2mlm(Safety_from_RSA_RR_TV)  


#predicting from the metric of Approach 5 - EXPLORATORY - RSA/Ttot
df_ema_filtered <- df_ema_filtered %>%
  mutate(RSA0_over_Ttot = RSA0_msec / Ttot)


df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    RSA0_over_Ttot_mean = mean(RSA0_over_Ttot, na.rm = TRUE),
    RSA0_over_Ttot_cwc = RSA0_over_Ttot - RSA0_over_Ttot_mean
  ) %>%
  ungroup()

Safety_from_Approach5_RSA_sitting_cwc <- lmer(
  Safety ~ RSA0_over_Ttot_cwc + RSA0_over_Ttot_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 | participant_id),
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Safety_from_Approach5_RSA_sitting_cwc)
r2mlm(Safety_from_Approach5_RSA_sitting_cwc)



###POSITIVE AFFECT###
###predicting from unadjusted RSA with no respiratory controlling
Pos_Aff_from_unadjusted_RSA_cwc <- lmer(
  Pos_Aff ~ RSA0_msec_cwc + RSA0_msec_mean +  BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1| participant_id),  
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Pos_Aff_from_unadjusted_RSA_cwc)  
summ(Pos_Aff_from_unadjusted_RSA_cwc)
r2mlm(Pos_Aff_from_unadjusted_RSA_cwc)


#from Approach 1 - delta_RSA0_over_Ttot_laboratory_45
Pos_Aff_from_Approach1_RSA_cwc <- lmer(
  Pos_Aff ~ delta_RSA0_over_Ttot_laboratory_45_cwc + delta_RSA0_over_Ttot_laboratory_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1  | participant_id),
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Pos_Aff_from_Approach1_RSA_cwc)
summ(Pos_Aff_from_Approach1_RSA_cwc)
r2mlm(Pos_Aff_from_Approach1_RSA_cwc)


#from Approach 2 - delta_RSA0_over_Ttot_ambulatory_45
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    delta_RSA0_over_Ttot_ambulatory_45_mean = mean(delta_RSA0_over_Ttot_ambulatory_45),  # person's mean
    delta_RSA0_over_Ttot_ambulatory_45_cwc = delta_RSA0_over_Ttot_ambulatory_45 - delta_RSA0_over_Ttot_ambulatory_45_mean  # centered within person
  ) %>%
  ungroup()

Pos_Aff_from_Approach2_RSA_cwc <- lmer(
  Pos_Aff ~ delta_RSA0_over_Ttot_ambulatory_45_cwc + delta_RSA0_over_Ttot_ambulatory_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1  | participant_id),
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Pos_Aff_from_Approach2_RSA_cwc)
summ(Pos_Aff_from_Approach2_RSA_cwc)
r2mlm(Pos_Aff_from_Approach2_RSA_cwc)

#from Approach 3 - RSA0_residual_45
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    RSA0_residual_45_mean = mean(RSA0_residual_45),  # person's mean
    RSA0_residual_45_cwc = RSA0_residual_45 - RSA0_residual_45_mean  # centered within person
  ) %>%
  ungroup()

Pos_Aff_from_Approach3_RSA_cwc <- lmer(
  Pos_Aff ~ RSA0_residual_45_cwc + RSA0_residual_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1 | participant_id),  
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Pos_Aff_from_Approach3_RSA_cwc)
summ(Pos_Aff_from_Approach3_RSA_cwc)
r2mlm(Pos_Aff_from_Approach3_RSA_cwc)


#from Approach 4 - RSA0_over_Vt
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    RSA0_over_Vt_mean = mean(RSA0_over_Vt),  # person's mean
    RSA0_over_Vt_cwc = RSA0_over_Vt - RSA0_over_Vt_mean    # centered within person
  ) %>%
  ungroup()

Pos_Aff_from_Approach4_RSA_cwc <- lmer(
  Pos_Aff ~ RSA0_over_Vt_cwc + RSA0_over_Vt_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1 | participant_id),  # random intercept only
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Pos_Aff_from_Approach4_RSA_cwc)
summ(Pos_Aff_from_Approach4_RSA_cwc)
r2mlm(Pos_Aff_from_Approach4_RSA_cwc)

r2mlm(Pos_Aff_from_unadjusted_RSA_cwc)
r2mlm(Pos_Aff_from_Approach1_RSA_cwc)
r2mlm(Pos_Aff_from_Approach2_RSA_cwc)
r2mlm(Pos_Aff_from_Approach3_RSA_cwc)
r2mlm(Pos_Aff_from_Approach4_RSA_cwc)

###NEGATIVE AFFECT###
###predicting from unadjusted RSA with no respiratory controlling
Neg_Aff_from_unadjusted_RSA_cwc <- lmer(
  Neg_Aff ~ RSA0_msec_cwc + RSA0_msec_mean +  BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1| participant_id),  
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Neg_Aff_from_unadjusted_RSA_cwc)  
summ(Neg_Aff_from_unadjusted_RSA_cwc)
r2mlm(Neg_Aff_from_unadjusted_RSA_cwc)


#from Approach 1 - delta_RSA0_over_Ttot_laboratory_45
Neg_Aff_from_Approach1_RSA_cwc <- lmer(
  Neg_Aff ~ delta_RSA0_over_Ttot_laboratory_45_cwc + delta_RSA0_over_Ttot_laboratory_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1  | participant_id),
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Neg_Aff_from_Approach1_RSA_cwc)
summ(Neg_Aff_from_Approach1_RSA_cwc)
r2mlm(Neg_Aff_from_Approach1_RSA_cwc)

r2mlm(Neg_Aff_from_unadjusted_RSA_cwc)
r2mlm(Neg_Aff_from_Approach1_RSA_cwc)


#from Approach 2 - delta_RSA0_over_Ttot_ambulatory_45
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    delta_RSA0_over_Ttot_ambulatory_45_mean = mean(delta_RSA0_over_Ttot_ambulatory_45),  # person's mean
    delta_RSA0_over_Ttot_ambulatory_45_cwc = delta_RSA0_over_Ttot_ambulatory_45 - delta_RSA0_over_Ttot_ambulatory_45_mean  # centered within person
  ) %>%
  ungroup()

Neg_Aff_from_Approach2_RSA_cwc <- lmer(
  Neg_Aff ~ delta_RSA0_over_Ttot_ambulatory_45_cwc + delta_RSA0_over_Ttot_ambulatory_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1  | participant_id),
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Neg_Aff_from_Approach2_RSA_cwc)
summ(Neg_Aff_from_Approach2_RSA_cwc)
r2mlm(Neg_Aff_from_Approach2_RSA_cwc)

#from Approach 3 - RSA0_residual_45
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    RSA0_residual_45_mean = mean(RSA0_residual_45),  # person's mean
    RSA0_residual_45_cwc = RSA0_residual_45 - RSA0_residual_45_mean  # centered within person
  ) %>%
  ungroup()

Neg_Aff_from_Approach3_RSA_cwc <- lmer(
  Neg_Aff ~ RSA0_residual_45_cwc + RSA0_residual_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1 | participant_id),  
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Neg_Aff_from_Approach3_RSA_cwc)
summ(Neg_Aff_from_Approach3_RSA_cwc)
r2mlm(Neg_Aff_from_Approach3_RSA_cwc)


#from Approach 4 - RSA0_over_Vt
df_ema_filtered <- df_ema_filtered %>%
  group_by(participant_id) %>%
  mutate(
    RSA0_over_Vt_mean = mean(RSA0_over_Vt),  # person's mean
    RSA0_over_Vt_cwc = RSA0_over_Vt - RSA0_over_Vt_mean    # centered within person
  ) %>%
  ungroup()

Neg_Aff_from_Approach4_RSA_cwc <- lmer(
  Neg_Aff ~ RSA0_over_Vt_cwc + RSA0_over_Vt_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1 | participant_id),  # random intercept only
  data = df_ema_filtered,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(Neg_Aff_from_Approach4_RSA_cwc)
summ(Neg_Aff_from_Approach4_RSA_cwc)
r2mlm(Neg_Aff_from_Approach4_RSA_cwc)

r2mlm(Neg_Aff_from_unadjusted_RSA_cwc)
r2mlm(Neg_Aff_from_Approach1_RSA_cwc)
r2mlm(Neg_Aff_from_Approach2_RSA_cwc)
r2mlm(Neg_Aff_from_Approach3_RSA_cwc)
r2mlm(Neg_Aff_from_Approach4_RSA_cwc)


##############################################################
####CALCULATION OF THE SUMMARY STATS FOR THE EMA VARIABLES####
##############################################################

library(dplyr)

# calculate per-participant mean and SD
psychological_stats <- df_ema_filtered %>%
  group_by(participant_id) %>% #calculating the sd and mean per participant
  summarise(
    mean_Tot_Stress = mean(Tot_Stress, na.rm = TRUE),
    sd_Tot_Stress = sd(Tot_Stress, na.rm = TRUE),
    mean_Safety = mean(Safety, na.rm = TRUE),
    sd_Safety = sd(Safety, na.rm = TRUE),
    mean_Pos_Aff = mean(Pos_Aff, na.rm = TRUE),
    sd_Pos_Aff = sd(Pos_Aff, na.rm = TRUE),
    mean_Neg_Aff = mean(Neg_Aff, na.rm = TRUE),
    sd_Neg_Aff = sd(Neg_Aff, na.rm = TRUE)
  )

# view participant-level stats 
print(psychological_stats, n=37)
View(psychological_stats)
# Compute average of participant-level means and SDs
overall_psychological_stats <- psychological_stats %>%
  summarise(
    avg_mean_Tot_Stress = mean(mean_Tot_Stress, na.rm = TRUE),
    avg_sd_Tot_Stress = mean(sd_Tot_Stress, na.rm = TRUE),
    avg_mean_Safety = mean(mean_Safety, na.rm = TRUE),
    avg_sd_Safety = mean(sd_Safety, na.rm = TRUE),
    avg_mean_Pos_Aff = mean(mean_Pos_Aff, na.rm = TRUE),
    avg_sd_Pos_Aff = mean(sd_Pos_Aff, na.rm = TRUE),
    avg_mean_Neg_Aff = mean(mean_Neg_Aff, na.rm = TRUE),
    avg_sd_Neg_Aff = mean(sd_Neg_Aff, na.rm = TRUE)
  )

# View overall means and SDs
print(overall_psychological_stats)

##############################################################
##CALCULATION OF THE SUMMARY STATS FOR THE PHYSIO VARIABLES###
##############################################################

#calculating the summary statistics

physio_stats <- final_df_sitting %>%
  group_by(participant_id) %>%
  summarise(
    mean_RSA0 = mean(RSA0_msec, na.rm = TRUE),
    sd_RSA0 = sd(RSA0_msec, na.rm = TRUE),
    
    mean_Avg_IBI = mean(Average_IBI_msec, na.rm = TRUE),
    sd_Avg_IBI = sd(Average_IBI_msec, na.rm = TRUE),
    
    mean_Tidal_Volume = mean(Calibrated_Tidal_Volume, na.rm = TRUE),
    sd_Tidal_Volume = sd(Calibrated_Tidal_Volume, na.rm = TRUE),
    
    mean_Resp_Rate = mean(Respiration_Rate_bpm, na.rm = TRUE),
    sd_Resp_Rate = sd(Respiration_Rate_bpm, na.rm = TRUE)
  ) %>%
  ungroup()

print(physio_stats, n=42)


# 2. Compute overall averages of those participant-level stats
overall_physio_stats <- physio_stats %>%
  summarise(
    avg_mean_RSA0 = mean(mean_RSA0, na.rm = TRUE),
    avg_sd_RSA0 = mean(sd_RSA0, na.rm = TRUE),
    
    avg_mean_Avg_IBI = mean(mean_Avg_IBI, na.rm = TRUE),
    avg_sd_Avg_IBI = mean(sd_Avg_IBI, na.rm = TRUE),
    
    avg_mean_Tidal_Volume = mean(mean_Tidal_Volume, na.rm = TRUE),
    avg_sd_Tidal_Volume = mean(sd_Tidal_Volume, na.rm = TRUE),
    
    avg_mean_Resp_Rate = mean(mean_Resp_Rate, na.rm = TRUE),
    avg_sd_Resp_Rate = mean(sd_Resp_Rate, na.rm = TRUE)
  )

# View overall means and SDs
print(overall_physio_stats)
View(overall_physio_stats)

#QUADRATIC MODELS FOR THE PREDICTION OF IBI: these do not converge, and so this added complexity is not needed/cannot be accommodated

# Center variables and add quadratic terms
final_df_sitting <- final_df_sitting %>%
  group_by(participant_id) %>%
  mutate(
    # Approach 0: RSA0_msec
    RSA0_msec_mean = mean(RSA0_msec),
    RSA0_msec_cwc = RSA0_msec - RSA0_msec_mean,
    RSA0_msec_cwc2 = RSA0_msec_cwc^2,
    
    # Approach 1: ΔRSA0/Vt - lab
    delta_RSA0_over_Ttot_lab45_mean = mean(delta_RSA0_over_Ttot_laboratory_45),
    delta_RSA0_over_Ttot_lab45_cwc = delta_RSA0_over_Ttot_laboratory_45 - delta_RSA0_over_Ttot_lab45_mean,
    delta_RSA0_over_Ttot_lab45_cwc2 = delta_RSA0_over_Ttot_lab45_cwc^2,
    
    # Approach 2: ΔRSA0/Vt - ambulatory
    delta_RSA0_over_Ttot_amb45_mean = mean(delta_RSA0_over_Ttot_ambulatory_45),
    delta_RSA0_over_Ttot_amb45_cwc = delta_RSA0_over_Ttot_ambulatory_45 - delta_RSA0_over_Ttot_amb45_mean,
    delta_RSA0_over_Ttot_amb45_cwc2 = delta_RSA0_over_Ttot_amb45_cwc^2,
    
    # Approach 3: residual RSA
    RSA0_residual_45_mean = mean(RSA0_residual_45),
    RSA0_residual_45_cwc = RSA0_residual_45 - RSA0_residual_45_mean,
    RSA0_residual_45_cwc2 = RSA0_residual_45_cwc^2,
    
    # Approach 4: RSA/Vt
    RSA0_over_Vt_mean = mean(RSA0_over_Vt),
    RSA0_over_Vt_cwc = RSA0_over_Vt - RSA0_over_Vt_mean,
    RSA0_over_Vt_cwc2 = RSA0_over_Vt_cwc^2
  ) %>%
  ungroup()

IBI_from_unadjusted_RSA_sitting_cwc_quad <- lmer(
  Average_IBI_msec ~ RSA0_msec_cwc + RSA0_msec_cwc2 + RSA0_msec_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1 + RSA0_msec_cwc + RSA0_msec_cwc2 | participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(IBI_from_unadjusted_RSA_sitting_cwc_quad)
r2mlm(IBI_from_unadjusted_RSA_sitting_cwc_quad)


IBI_from_Approach1_RSA_sitting_cwc_quad <- lmer(
  Average_IBI_msec ~ delta_RSA0_over_Ttot_lab45_cwc + delta_RSA0_over_Ttot_lab45_cwc2 + delta_RSA0_over_Ttot_lab45_mean +
    BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric +
    (1 + delta_RSA0_over_Ttot_lab45_cwc + delta_RSA0_over_Ttot_lab45_cwc2 | participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

r2mlm(IBI_from_unadjusted_RSA_sitting_cwc_quad)



###PREDICTING IBI USING A VERY CONSERVATIVE PHYSICAL ACTIVITY THRESHOLD###
#the dataframe name is the same, make sure you clear the global environment before running
#the dataframe that should be imported is: 

#importing the required packages
library(lme4)
library(ggplot2)
library(lmerTest)
library(performance)
library(dplyr)
library(readxl)
library(HLMdiag)
library(jtools)
library(r2mlm)
library(tidyverse)

final_df_sitting <- read_excel("/Users/melisasaygin/Documents/Validation of RSA Paper Code (for GitHub)/final_df_sitting_conservative.xlsx")

#checking data classes (types) and correcting as needed
sapply(final_df_sitting[, c("Average_IBI_msec", "RSA0_msec", "BMI", "IPAQ_SF", 
                            "Age", "Sex", "Current_Smoking", "participant_id",
                            "delta_RSA0_over_Ttot_laboratory_45", "delta_RSA0_over_Ttot_ambulatory_45",
                            "RSA0_residual", "RSA0_over_Vt")], class)

final_df_sitting$participant_id <- as.factor(final_df_sitting$participant_id)
final_df_sitting$Sex <- factor(final_df_sitting$Sex) #make sex into a factor variable
final_df_sitting$Current_Smoking <- factor(final_df_sitting$Current_Smoking, levels = c(0, 1), labels = c("Non-Smoker", "Smoker"))
#make the current_smoking and sex variables numeric to later be able to use the r2mlm package
final_df_sitting$Current_Smoking_numeric <- ifelse(final_df_sitting$Current_Smoking == "Smoker", 1, 0)  #0 becomes non-smoker, 1 becomes smoker
final_df_sitting$Sex_numeric <- ifelse(final_df_sitting$Sex == "Male", 1, 0) # 0 becomes female, 1 becomes male

#re-check the data types
sapply(final_df_sitting[, c("Average_IBI_msec", "RSA0_msec", "BMI", "IPAQ_SF", 
                            "Age", "Sex", "Current_Smoking", "participant_id",
                            "delta_RSA0_over_Ttot_laboratory_45", "delta_RSA0_over_Ttot_ambulatory_45",
                            "RSA0_residual", "RSA0_over_Vt", "Current_Smoking_numeric", "Sex_numeric")], class)

# remove outliers (values deviating more than 4.5 SDs for IBI and RSA metrics)
final_df_sitting <- final_df_sitting %>%
  group_by(participant_id) %>%
  filter(
    abs(RSA0_msec - mean(RSA0_msec)) / sd(RSA0_msec) <= 4.5,  # Remove RSA0_msec outliers
    abs(Average_IBI_msec - mean(Average_IBI_msec)) / sd(Average_IBI_msec) <= 4.5,  # Remove IBI outliers
    abs(delta_RSA0_over_Ttot_laboratory_45 - mean(delta_RSA0_over_Ttot_laboratory_45)) / sd(delta_RSA0_over_Ttot_laboratory_45) <= 4.5,  # Remove delta_RSA0_over_Ttot_laboratory_45 outliers
    abs(delta_RSA0_over_Ttot_ambulatory_45 - mean(delta_RSA0_over_Ttot_ambulatory_45)) / sd(delta_RSA0_over_Ttot_ambulatory_45) <= 4.5,  # Remove delta_RSA0_over_Ttot_ambulatory_45 outliers
    abs(RSA0_residual - mean(RSA0_residual)) / sd(RSA0_residual) <= 4.5,  # Remove RSA0_residual outliers
    abs(RSA0_over_Vt - mean(RSA0_over_Vt)) / sd(RSA0_over_Vt) <= 4.5  # Remove RSA0_over_Vt outliers
  ) %>%
  ungroup()



####################################################################################
##### RESEARCH QUESTION 2 (PREDICTING IBI FROM RSA METRICS) FOR SITTING POSTURE ####
####################################################################################
IBI_from_covariates_sitting_cwc <- lmer(  #seeing how much between-persons variance is explained by just the covariates - exploratory
  Average_IBI_msec ~ BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 | participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(IBI_from_covariates_sitting_cwc)
r2mlm(IBI_from_covariates_sitting_cwc)

###predicting from unadjusted RSA with no respiratory controlling

final_df_sitting <- final_df_sitting %>%   #centering within clusters (individuals), also calculating the individual mean
  group_by(participant_id) %>%
  mutate(
    RSA0_msec_mean = mean(RSA0_msec),
    RSA0_msec_cwc = RSA0_msec - RSA0_msec_mean  # CWC = centered within cluster
  ) %>%
  ungroup()


IBI_from_unadjusted_RSA_sitting_cwc <- lmer(  #introduce the centered within cluster RSA metric (level1) along with the cluster means (level2)
  Average_IBI_msec ~ RSA0_msec_cwc + RSA0_msec_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 + RSA0_msec_cwc | participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")   #when fitted with REML set as TRUE, no warning messages but the variances are essentially the same
)
plot(IBI_from_unadjusted_RSA_sitting_cwc)
summary(IBI_from_unadjusted_RSA_sitting_cwc)
summ(IBI_from_unadjusted_RSA_sitting_cwc)
r2(IBI_from_unadjusted_RSA_sitting_cwc)
r2mlm(IBI_from_unadjusted_RSA_sitting_cwc)

###predicting from the metric of Approach 1 - ΔRSA0/Vt via baseline regression made with paced breathing

final_df_sitting <- final_df_sitting %>%   #doing group mean centering
  group_by(participant_id) %>%
  mutate(
    delta_RSA0_over_Ttot_lab45_mean = mean(delta_RSA0_over_Ttot_laboratory_45),
    delta_RSA0_over_Ttot_lab45_cwc = delta_RSA0_over_Ttot_laboratory_45 - delta_RSA0_over_Ttot_lab45_mean
  ) %>%
  ungroup()

IBI_from_Approach1_RSA_sitting_cwc <- lmer(
  Average_IBI_msec ~ delta_RSA0_over_Ttot_lab45_cwc + delta_RSA0_over_Ttot_lab45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 + delta_RSA0_over_Ttot_lab45_cwc | participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(IBI_from_Approach1_RSA_sitting_cwc)
summ(IBI_from_Approach1_RSA_sitting_cwc)
r2(IBI_from_Approach1_RSA_sitting_cwc)
r2mlm(IBI_from_Approach1_RSA_sitting_cwc)



#predicting from the metric of Approach 2 - ΔRSA0/Vt via baseline regression made with ambulatory data
final_df_sitting <- final_df_sitting %>%    #group mean centering
  group_by(participant_id) %>%
  mutate(
    delta_RSA0_over_Ttot_amb45_mean = mean(delta_RSA0_over_Ttot_ambulatory_45),
    delta_RSA0_over_Ttot_amb45_cwc = delta_RSA0_over_Ttot_ambulatory_45 - delta_RSA0_over_Ttot_amb45_mean
  ) %>%
  ungroup()

IBI_from_Approach2_RSA_sitting_cwc <- lmer(
  Average_IBI_msec ~ delta_RSA0_over_Ttot_amb45_cwc + delta_RSA0_over_Ttot_amb45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 + delta_RSA0_over_Ttot_amb45_cwc | participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(IBI_from_Approach2_RSA_sitting_cwc)
summ(IBI_from_Approach2_RSA_sitting_cwc)
r2mlm(IBI_from_Approach2_RSA_sitting_cwc)

#predicting from the metric of Approach 3 - residual RSA

final_df_sitting <- final_df_sitting %>%
  group_by(participant_id) %>%
  mutate(
    RSA0_residual_45_mean = mean(RSA0_residual_45),
    RSA0_residual_45_cwc = RSA0_residual_45 - RSA0_residual_45_mean
  ) %>%
  ungroup()

IBI_from_Approach3_RSA_sitting_cwc <- lmer(
  Average_IBI_msec ~ RSA0_residual_45_cwc + RSA0_residual_45_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 + RSA0_residual_45_cwc| participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(IBI_from_Approach3_RSA_sitting_cwc)
summ(IBI_from_Approach3_RSA_sitting_cwc)
r2(IBI_from_Approach3_RSA_sitting_cwc)
r2mlm(IBI_from_Approach3_RSA_sitting_cwc)

#predicting from the metric of Approach 4 - RSA/Vt

final_df_sitting <- final_df_sitting %>%   #group mean centering
  group_by(participant_id) %>%
  mutate(
    RSA0_over_Vt_mean = mean(RSA0_over_Vt),
    RSA0_over_Vt_cwc = RSA0_over_Vt - RSA0_over_Vt_mean
  ) %>%
  ungroup()


IBI_from_Approach4_RSA_sitting_cwc <- lmer(
  Average_IBI_msec ~ RSA0_over_Vt_cwc + RSA0_over_Vt_mean + BMI + IPAQ_SF + Age + Sex_numeric + Current_Smoking_numeric + 
    (1 + RSA0_over_Vt_cwc| participant_id),
  data = final_df_sitting,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

r2mlm(IBI_from_Approach4_RSA_sitting_cwc)
summary(IBI_from_Approach4_RSA_sitting_cwc)
summ(IBI_from_Approach4_RSA_sitting_cwc)

par(mar = c(6.75, 10.5, 2.625, 10.5))
r2mlm(IBI_from_unadjusted_RSA_sitting_cwc)
r2mlm(IBI_from_Approach1_RSA_sitting_cwc)
r2mlm(IBI_from_Approach2_RSA_sitting_cwc)
r2mlm(IBI_from_Approach3_RSA_sitting_cwc)
r2mlm(IBI_from_Approach4_RSA_sitting_cwc)




