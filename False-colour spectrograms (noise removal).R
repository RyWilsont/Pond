library(tuneR)
library(seewave)
library(signal)
library(entropy)

# Define the directory containing the .wav files
directory <- "C:/Users/Administrador/Downloads/Ordesa sound files/Daily Spectrograms/J016/Summer"

# Print directory path to confirm
print(paste("Checking directory:", directory))

# Check if the directory exists
if (!dir.exists(directory)) {
  stop("The specified directory does not exist.")
}

# List all .WAV files in the directory (considering case sensitivity)
file_list <- list.files(directory, pattern = "\\.WAV$", full.names = TRUE)

# Print file list to debug
print("List of .WAV files found:")
print(file_list)

# Check if there are any .WAV files
if (length(file_list) == 0) {
  stop("No .WAV files found in the directory.")
}

# Initialize lists to store results
entropy_list <- list()
aci_list <- list()
background_noise_list <- list()

# Function definitions for calculating metrics
calculate_entropy <- function(freq_bin) {
  total_energy <- sum(freq_bin)
  if (total_energy == 0) {
    return(0)  # Avoid division by zero
  }
  prob_distribution <- freq_bin / total_energy
  shannon_entropy <- -sum(prob_distribution * log(prob_distribution + 1e-10))
  return(shannon_entropy)
}

calculate_aci <- function(freq_bin) {
  local_maxima <- sum(diff(sign(diff(freq_bin))) == -2)
  aci_value <- local_maxima / length(freq_bin)
  return(aci_value)
}

calculate_background_noise <- function(freq_bin) {
  median_noise <- median(freq_bin)
  return(median_noise)
}

# Loop through each .wav file
for (audio_file in file_list) {
  # Load the audio file
  wave <- readWave(audio_file)
  
  # Extract the audio signal and the sampling rate
  y <- wave@left
  sr <- wave@samp.rate
  
  # Parameters for STFT
  n_fft <- 16384  # Number of FFT points
  hop_length <- n_fft / 2  # Hop length (50% overlap)
  
  # Perform Short-Time Fourier Transform (STFT) using specgram
  specgram_result <- specgram(y, n = n_fft, Fs = sr, overlap = n_fft - hop_length)
  
  # Get the magnitude spectrogram
  S <- abs(specgram_result$S)
  
  # Normalize the spectrogram
  S <- S / max(S)
  
  # Convert to dB scale
  S_db <- 10 * log10(S + 1e-10)
  
  # Calculate metrics
  entropy_values <- apply(S_db, 1, calculate_entropy)
  aci_values <- apply(S_db, 1, calculate_aci)
  background_noise_values <- apply(S_db, 1, calculate_background_noise)
  
  # Compute frequency bins
  freq_bins <- seq(0, sr / 2, length.out = n_fft / 2 + 1)[1:nrow(S_db)]
  
  # Ensure the length of freq_bins matches the length of metric values
  if (length(freq_bins) != length(entropy_values) ||
      length(freq_bins) != length(aci_values) ||
      length(freq_bins) != length(background_noise_values)) {
    warning(paste("Length mismatch in file:", basename(audio_file)))
    next
  }
  
  # Create data frames for each metric
  entropy_df <- data.frame(Frequency = freq_bins, Entropy = entropy_values, File = basename(audio_file))
  aci_df <- data.frame(Frequency = freq_bins, ACI = aci_values, File = basename(audio_file))
  background_noise_df <- data.frame(Frequency = freq_bins, Background_Noise = background_noise_values, File = basename(audio_file))
  
  # Append results to lists
  entropy_list[[basename(audio_file)]] <- entropy_df
  aci_list[[basename(audio_file)]] <- aci_df
  background_noise_list[[basename(audio_file)]] <- background_noise_df
}

# Combine all results into single data frames
entropy_all <- do.call(rbind, entropy_list)
aci_all <- do.call(rbind, aci_list)
background_noise_all <- do.call(rbind, background_noise_list)

# Print resulting data frames
print(head(entropy_all))
print(head(aci_all))
print(head(background_noise_all))

# Load necessary library
library(dplyr)

# Extract time and create the 'Time' column for entropy_all, and remove the 'File' column
entropy_all <- entropy_all %>%
  mutate(
    # Extracting the time part from the File column
    TempTime = substr(File, 10, 15),
    # Converting extracted time to hh:mm format
    Time = paste0(substr(TempTime, 1, 2), ":", substr(TempTime, 3, 4))
  ) %>%
  select(-TempTime, -File)  # Remove the intermediate 'TempTime' and 'File' columns

# View the updated data frame
head(entropy_all)

# Extract time and create the 'Time' column for aci_all, and remove the 'File' column
aci_all <- aci_all %>%
  mutate(
    # Extracting the time part from the File column
    TempTime = substr(File, 10, 15),
    # Converting extracted time to hh:mm format
    Time = paste0(substr(TempTime, 1, 2), ":", substr(TempTime, 3, 4))
  ) %>%
  select(-TempTime, -File)  # Remove the intermediate 'TempTime' and 'File' columns

# View the updated data frame
head(aci_all)

# Extract time and create the 'Time' column for background_noise_all, and remove the 'File' column
background_noise_all <- background_noise_all %>%
  mutate(
    # Extracting the time part from the File column
    TempTime = substr(File, 10, 15),
    # Converting extracted time to hh:mm format
    Time = paste0(substr(TempTime, 1, 2), ":", substr(TempTime, 3, 4))
  ) %>%
  select(-TempTime, -File)  # Remove the intermediate 'TempTime' and 'File' columns

# View the updated data frame
head(background_noise_all)

#### Scale acoustic indices ####

# Entropy

# Find min and max of Entropy
min_entropy <- min(entropy_all$Entropy)
max_entropy <- max(entropy_all$Entropy)

# Scale the Entropy column
entropy_all$Scaled_Entropy <- (entropy_all$Entropy - min_entropy) / (max_entropy - min_entropy)

# Print the result
print(data)

# Acoustic complexity

# Find min and max of aci
min_aci <- min(aci_all$ACI)
max_aci <- max(aci_all$ACI)

# Scale the aci column
aci_all$Scaled_ACI <- (aci_all$ACI - min_aci) / (max_aci - min_aci)

# Background noise

# Find min and max of background_noise
min_noise <- min(background_noise_all$Background_Noise)
max_noise <- max(background_noise_all$Background_Noise)

# Scale the background_noise column
background_noise_all$Scaled_Background_Noise <- (background_noise_all$Background_Noise - min_noise) / (max_noise - min_noise)

#### Mapping to RBG ####

# Assuming that the scaled columns are already present in the respective data frames
# Combine the relevant columns into a single data frame

combined_data <- data.frame(
  Scaled_Entropy = entropy_all$Scaled_Entropy,
  Scaled_ACI = aci_all$Scaled_ACI,
  Scaled_Background_Noise = background_noise_all$Scaled_Background_Noise
)

# Initialize RGB_Data with the same number of rows
RGB_Data <- data.frame(matrix(ncol = 1, nrow = nrow(combined_data)))
colnames(RGB_Data) <- "Color"

# Assuming combined_data is your data frame with the scaled variables
library(grDevices) # For the rgb function

# Function to generate RGB values based on different variable combinations
generate_rgb_combinations <- function(data, var3, var2, var1) {
  rgb(data[[var2]], data[[var1]], data[[var3]])
}

# All possible combinations of the variables
combinations <- list(
  c("Scaled_Entropy", "Scaled_ACI", "Scaled_Background_Noise"),
  c("Scaled_Entropy", "Scaled_Background_Noise", "Scaled_ACI"),
  c("Scaled_ACI", "Scaled_Entropy", "Scaled_Background_Noise"),
  c("Scaled_ACI", "Scaled_Background_Noise", "Scaled_Entropy"),
  c("Scaled_Background_Noise", "Scaled_Entropy", "Scaled_ACI"),
  c("Scaled_Background_Noise", "Scaled_ACI", "Scaled_Entropy")
)

# Create an empty list to store the results
RGB_Combinations <- list()

# Loop through the combinations and generate RGB colors
for (i in 1:length(combinations)) {
  combo <- combinations[[i]]
  RGB_Combinations[[i]] <- generate_rgb_combinations(combined_data, combo[1], combo[2], combo[3])
}

# Example: Storing the results in a data frame
RGB_Data <- data.frame(
  Color1 = RGB_Combinations[[1]],
  Color2 = RGB_Combinations[[2]],
  Color3 = RGB_Combinations[[3]],
  Color4 = RGB_Combinations[[4]],
  Color5 = RGB_Combinations[[5]],
  Color6 = RGB_Combinations[[6]]
)

# Print the RGB data frame to check the results
print(RGB_Data)

# Print the result
head(RGB_Data)

head(entropy_all)

##### Fit HEX codes to freq and time #####

# Extract the required columns from entropy_all and RGB_Data
frequency <- entropy_all$Frequency
time <- entropy_all$Time
color <- RGB_Data$Color1

# Combine these columns into a new data frame
combined_df <- data.frame(Frequency = frequency, Time = time, Color = color)

# Print the result
print(head(combined_df))

#### Plot the spectrogram ####

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Ensure 'Time' is in a time format (for ggplot2 to handle it correctly)
combined_df <- combined_df %>%
  mutate(Time = as.POSIXct(Time, format = "%H:%M", tz = "UTC"))

#### Removing background noise ####

# Step 1: Count the occurrences of each HEX code
hex_table <- table(unlist(RGB_Data))

# Convert the table to a data frame for easier manipulation
hex_df <- as.data.frame(hex_table)

# Rename columns for clarity
colnames(hex_df) <- c("HEX", "Count")

# Step 2: Determine the threshold for the dominant 10%
threshold_count <- quantile(hex_df$Count, 0.9)  # Get the count for the 90th percentile

# Identify dominant HEX codes (those that have counts greater than or equal to the threshold)
dominant_hexes <- hex_df$HEX[hex_df$Count >= threshold_count]

# Print dominant HEX codes for reference
print("Dominant HEX Codes:")
print(dominant_hexes)

# Step 3: Replace dominant HEX codes with black in RGB_Data
# Create a function to replace dominant HEX codes
replace_dominant_hex <- function(color) {
  if (color %in% dominant_hexes) {
    return("#000000")  # Replace with black
  } else {
    return(color)  # Keep the original color
  }
}

# Apply the function to replace colors in the first column of RGB_Data
RGB_Data$Adjusted_Color <- sapply(RGB_Data$Color1, replace_dominant_hex)

# Create a new combined_df with the adjusted colors
combined_df_adjusted <- data.frame(
  Frequency = frequency,
  Time = time,
  Color = RGB_Data$Adjusted_Color
)

# Ensure 'Time' is in POSIXct format
combined_df_adjusted <- combined_df_adjusted %>%
  mutate(Time = as.POSIXct(Time, format = "%H:%M", tz = "UTC"))

# Step 4: Plot the adjusted spectrogram
p <- ggplot(combined_df_adjusted, aes(x = Time, y = Frequency)) +
  geom_tile(aes(fill = Color)) +
  scale_fill_identity() +
  scale_x_datetime(date_breaks = "2 hours", date_labels = "%H:%M") +
  labs(x = "Time", y = "Frequency (Hz)") +
  theme_bw()

ggsave("J016 Summer (321) .jpeg", plot = p, width = 8, height = 6, dpi = 300)

