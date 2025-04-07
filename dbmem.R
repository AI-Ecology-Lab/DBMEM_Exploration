# SHR DBMEM Script----
# Southern Hydrate Ridge Marine Data
# Date: April , 4, 2025


# ---------------------------
# 1. Read and Clean the Data
# ---------------------------

#
# RUN THE SCIPEN!!!!
#
  
options(scipen=999) 

rm(list=ls()) 
gc()

setwd("C:/Users/Attca/OneDrive/Desktop/Edpsy538/DBMEM")
getwd()  ## confirm that the working directory is set correctly

# Load required packages
library(adespatial)   # for dbmem() and forward.sel()
library(vegan)        # like rda() and varpart()
library(lubridate)    # for handling dates/times

data <- read.csv("SHR_2022-2023_fauna_ctd_calcO2_pressure.csv", header = TRUE, stringsAsFactors = FALSE)

# Remove the "Bubble" column entirely
data$Bubble <- NULL

# ---------------------------
# 2. Change datetime
# ---------------------------
data$datetime <- as.POSIXct(data$Timestamp, format = "%m/%d/%Y %H:%M", tz = "UTC")

# Alternatively:
# data$datetime <- as.POSIXct(paste(data$Date, data$Time),
#                             format = "%d-%b-%y %H:%M:%S", tz = "UTC")

data <- data[order(data$datetime), ]

# ---------------------------
# 3. Create a Temporal Coordinate and DBMEM
# ---------------------------
# Convert datetime to numeric (seconds since epoch)
time_numeric <- as.numeric(data$datetime)
# DBMEM expects a matrix of coordinates, here a 1D time vector is used, we may switch to 4D if we have enough df
time_matrix <- as.matrix(time_numeric)

# Generate DBMEM eigenfunctions 
dbmem_vars <- dbmem(time_matrix, silent = FALSE)

# ---------------------------
# 4. Define Community and Environmental Data
# ---------------------------

species_columns <- c("Sablefish", "Sea.Stars", "Crabs", "Hagfish", "Euphausia",
                     "Liponema", "Flatfish", "Rockfish", "Eelpout")

comm <- data[, species_columns]

# Define environmental variables to be used.
env_columns <- c("Temperature", "Conductivity", "Pressure", "Salinity",
                 "Oxygen..ml.l.", "PressurePSI")

env <- data[, env_columns]

# ---------------------------
# 5. Forward Selection of DBMEM Eigenfunctions
# ---------------------------
# Select significant DBMEM eigenfunctions using forward selection against community data.
sel <- forward.sel(comm, as.matrix(dbmem_vars), alpha = 0.05, nperm = 999)
if (!is.null(sel$order) && length(sel$order) > 0) {
  selected_dbmem <- dbmem_vars[, sel$order]
} else {
  warning("No significant DBMEM eigenfunctions were selected. Using all DBMEM eigenfunctions.")
  selected_dbmem <- dbmem_vars
}

# ---------------------------
# 6. Analysis 1: Community ~ Temporal DBMEM Only
# ---------------------------
# Build an RDA model using only temporal DBMEM eigenfunctions as predictors.
rda_dbmem <- rda(comm ~ ., data = as.data.frame(selected_dbmem))
summary(rda_dbmem)

# Plot the RDA ordination
plot(rda_dbmem, main = "RDA: Community ~ Temporal DBMEM Eigenfunctions")
# Add a title, and use default biplot settings for interpretation

# ---------------------------
# 7. Analysis 2: Community ~ Environmental Data Only
# ---------------------------
# Build an RDA model using only the environmental data.
rda_env <- rda(comm ~ ., data = env)
summary(rda_env)

# Plot the RDA ordination
plot(rda_env, main = "RDA: Community ~ Environmental Data")
# Default biplot for environmental RDA

# ---------------------------
# 8. Analysis 3: Combined Model and Variation Partitioning
# ---------------------------
# Combine the DBMEM eigenfunctions and environmental data in a single RDA model.
rda_combined <- rda(comm ~ ., data = cbind(as.data.frame(selected_dbmem), env))
summary(rda_combined)

# Plot the combined RDA
plot(rda_combined, main = "RDA: Community ~ Temporal DBMEM + Environmental Data")

# Variation Partitioning to assess the unique and shared contributions, maybe we need to zscore them later?:
varpart_model <- varpart(comm, selected_dbmem, env)
plot(varpart_model, main = "Variation Partitioning: DBMEM vs. Environmental Data")

# ---------------------------
# 9. Additional Plots (Optional)
# ---------------------------
# Scree plot for DBMEM eigenvalues to inspect the distribution of variance. Should be normal, skewed left(?)
eigenvalues <- attributes(dbmem_vars)$values
plot(eigenvalues, type = "b", xlab = "Eigenfunction Number", 
     ylab = "Eigenvalue", main = "DBMEM Eigenvalues Scree Plot")
