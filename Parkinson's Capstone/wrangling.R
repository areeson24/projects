# load packages
library(tidyverse)
### PATIENT META DATA ###
# load patient meta data
pat_data <- read.csv("~/Desktop/pads-parkinsons-disease-smartwatch-dataset-1.0.0/preprocessed/file_list.csv")
# format demographic info and patient characteristics
pat_data <- pat_data |>
  select(id, age, height, weight, gender, handedness, label) |>
  mutate(gender = factor(gender) |> relevel(ref = "male"),
         handedness = factor(handedness) |> relevel(ref = "right"),
         label = factor(label, labels = c("Healthy", "Other", "PD")) |> 
           relevel(ref = "Healthy"))
### MOVEMENT SIGNAL DATA ###
# function to load binary movement file
read_pads_move <- function(file_path, sample_rate = 100, cut_initial_sec = 0.5) {
  # connection
  con <- file(file_path, "rb")
  # read in raw data
  raw_data <- readBin(con, what = "numeric", n = 1e7, size = 4, endian = "little")
  # close connection
  close(con)
  # define knowns
  n_channels <- 132  
  n_samples <- length(raw_data) / n_channels
  # reshape into matrix: samples x channels
  data_mat <- matrix(raw_data, nrow = n_samples, ncol = n_channels, byrow = TRUE)
  # generate channel names
  tasks <- c("Relaxed1", "Relaxed2", "RelaxedTask1", "RelaxedTask2", "StretchHold", 
             "HoldWeight", "DrinkGlas", "CrossArms", "TouchNose", "Entrainment1", "Entrainment2")
  wrists <- c("Left", "Right")
  sensors <- c("Acceleration", "Rotation")
  axes <- c("X", "Y", "Z")
  channel_names <- c()
  for (task in tasks) {
    for (wrist in wrists) {
      for (sensor in sensors) {
        for (axis in axes) {
          channel_names <- c(channel_names, paste(task, wrist, sensor, axis, sep = "_"))
        }
      }
    }
  }
  # assign column names
  colnames(data_mat) <- channel_names
  # cut off first 0.5 seconds (assuming 100Hz -> 50 samples)
  cut_samples <- sample_rate * cut_initial_sec
  if (n_samples > cut_samples) {
    data_mat <- data_mat[-(1:cut_samples), ]
  } else {
    warning("Not enough samples to cut the initial 0.5 seconds.")
  }
  return(data_mat)
}
# list all subject file paths
subject_files <- list.files(path = "~/Desktop/pads-parkinsons-disease-smartwatch-dataset-1.0.0/preprocessed/movement/", 
                            pattern = "*.bin", full.names = TRUE)
# load all subject movement data
all_subject_data <- lapply(subject_files, read_pads_move)
# summarize time series for each subject
summarize_subject <- function(data_mat) {
  # compute means and standard deviations for each channel
  feature_means <- colMeans(data_mat, na.rm = TRUE)
  feature_sds <- apply(data_mat, 2, sd, na.rm = TRUE)
  # label the mean and sd features
  mean_labels <- paste(colnames(data_mat), "_mean", sep = "")
  sd_labels <- paste(colnames(data_mat), "_sd", sep = "")
  # combine means and standard deviations into a single feature vector
  feature_vector <- c(feature_means, feature_sds)
  # assign the proper column names for both mean and sd features
  names(feature_vector) <- c(mean_labels, sd_labels)
  return(feature_vector)
}
# apply to all subjects
movement_data <- lapply(all_subject_data, summarize_subject)
movement_data <- do.call(rbind, movement_data)
movement_data <- as.data.frame(movement_data)
# fit PCA on movement data
fit_pca <- prcomp(movement_data, center = TRUE, scale. = TRUE)
summary(fit_pca)
fit_pca$rotation[,1:10] |> round(3)
# store PCs 1-10 (>90% of variation)
movement_PCA_data <- fit_pca$x[,1:10] |> as.data.frame()
colnames(movement_PCA_data) <- paste("move_", colnames(movement_PCA_data), sep = "")
### QUESTIONAIRE DATA ###
# function to load individual's binary questionnaire file
read_pads_quest <- function(file_path) {
  # connection
  con <- file(file_path, "rb")
  # read in raw data
  raw_data <- readBin(con, what = "numeric", n = 1e7, size = 4, endian = "little")
  # close connection
  close(con)
  # generate question names
  q_names <- paste("nmsymp_Q", 1:30, sep = "")
  # assign question names
  names(raw_data) <- q_names
  return(raw_data)
}
# list all subject questionnaire file paths
subject_q_files <- list.files(path = "~/Desktop/pads-parkinsons-disease-smartwatch-dataset-1.0.0/preprocessed/questionnaire/", 
                            pattern = "*.bin", full.names = TRUE)
# load all subject questionnaire data
questionaire_data <- lapply(subject_q_files, read_pads_quest)
questionaire_data <- do.call(rbind, questionaire_data)
questionaire_data <- as.data.frame(questionaire_data) |>
  mutate(across(everything(), ~factor(ifelse(. == 1, "Yes", "No"))))
### FINAL DATASET ###
pads_clean <- cbind(pat_data, movement_PCA_data, questionaire_data)
saveRDS(pads_clean, "~/Desktop/2024-2025/School/bst263/finalproj/pads_clean.rds")
