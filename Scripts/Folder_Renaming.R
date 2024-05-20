setwd("C:/Users/axi313/Documents/TARA_Entry/Raw_Data")

# Load necessary library

library(dplyr)
old_names = c("DeArmas-17585-001", "DeArmas-17585-002", "DeArmas-17585-003", "DeArmas-17585-004",
              "DeArmas-17586-001", "DeArmas-17586-002", "DeArmas-17587-001", "DeArmas-17587-002",
              "DeArmas-17587-003", "DeArmas-17587-004", "DeArmas-17588-002")

# Create a dataframe with the existing folder names and corresponding new names
data <- data.frame(
  old_names = c("DeArmas-17585-001", "DeArmas-17585-002", "DeArmas-17585-003", "DeArmas-17585-004",
                "DeArmas-17586-001", "DeArmas-17586-002", "DeArmas-17587-001", "DeArmas-17587-002",
                "DeArmas-17587-003", "DeArmas-17587-004", "DeArmas-17588-002"),
  new_names = c("CP013_entry", "CP016_entry", "CP022_entry", "SA-AH-29_entry", 
                "CP011_entry", "SA-TY-021_entry", "CP017_entry", "CP042_entry", 
                "CS005_entry", "CS015_entry", "CP006_12m")
)

# Function to rename folders
rename_folders <- function(old_name, new_name) {
  if (file.exists(old_name)) {
    file.rename(old_name, new_name)
  } else {
    message(paste("Folder", old_name, "does not exist."))
  }
}

# Apply the renaming function to each folder
apply(data, 1, function(row) rename_folders(row[1], row[2]))

# Print the updated data frame to confirm renaming
print(data)

# Function to rename internal folders within each "multi" folder
rename_internal_folders <- function(main_folder) {
  multi_path <- file.path(main_folder, "multi")
  count_path <- file.path(multi_path, "count")
  vdj_b_path <- file.path(multi_path, "vdj_b")
  vdj_t_path <- file.path(multi_path, "vdj_T")
  
  # Rename subfolders if they exist
  if (file.exists(count_path)) {
    file.rename(count_path, file.path(multi_path, "GEX_FB"))
  } else {
    message(paste("Folder", count_path, "does not exist."))
  }
  
  if (file.exists(vdj_b_path)) {
    file.rename(vdj_b_path, file.path(multi_path, "BCR"))
  } else {
    message(paste("Folder", vdj_b_path, "does not exist."))
  }
  
  if (file.exists(vdj_t_path)) {
    file.rename(vdj_t_path, file.path(multi_path, "TCR"))
  } else {
    message(paste("Folder", vdj_t_path, "does not exist."))
  }
}

# List of main folders
main_folders <- c("CP013_entry", "CP016_entry", "CP022_entry", "SA-AH-29_entry", 
                  "CP011_entry", "SA-TY-021_entry", "CP017_entry", "CP042_entry", 
                  "CS005_entry", "CS015_entry", "CP006_12m")

# Apply the renaming function to each main folder
sapply(main_folders, rename_internal_folders)

#######
# Load necessary library
library(dplyr)

# Function to move subfolders up one level and remove the old_names folder
move_and_remove_old_names <- function(new_name, old_name) {
  old_names_path <- file.path(new_name, "per_sample_outs", old_name)
  
  if (file.exists(old_names_path)) {
    # List subfolders in the old_names folder
    subfolders <- list.files(old_names_path, full.names = TRUE)
    
    # Move each subfolder up one level
    for (subfolder in subfolders) {
      subfolder_name <- basename(subfolder)
      new_path <- file.path(new_name, "per_sample_outs", subfolder_name)
      file.rename(subfolder, new_path)
    }
    
    # Remove the old_names folder
    unlink(old_names_path, recursive = TRUE)
  } else {
    message(paste("Folder", old_names_path, "does not exist."))
  }
}

# Apply the function to each row in the data frame
apply(data, 1, function(row) move_and_remove_old_names(row["new_names"], row["old_names"]))


# Create the data frame
data <- data.frame(
  new_names = c("CP013_entry", "CP016_entry", "CP022_entry", "SA-AH-29_entry", 
                "CP011_entry", "SA-TY-021_entry", "CP017_entry", "CP042_entry", 
                "CS005_entry", "CS015_entry", "CP006_12m")
)

# Function to rename subfolders within per_sample_outs
rename_subfolders <- function(main_folder) {
  per_sample_outs_path <- file.path(main_folder, "per_sample_outs")
  
  count_path <- file.path(per_sample_outs_path, "count")
  vdj_b_path <- file.path(per_sample_outs_path, "vdj_b")
  vdj_t_path <- file.path(per_sample_outs_path, "vdj_t")
  
  # Rename subfolders if they exist
  if (file.exists(count_path)) {
    file.rename(count_path, file.path(per_sample_outs_path, "GEX_FB"))
  } else {
    message(paste("Folder", count_path, "does not exist."))
  }
  
  if (file.exists(vdj_b_path)) {
    file.rename(vdj_b_path, file.path(per_sample_outs_path, "BCR"))
  } else {
    message(paste("Folder", vdj_b_path, "does not exist."))
  }
  
  if (file.exists(vdj_t_path)) {
    file.rename(vdj_t_path, file.path(per_sample_outs_path, "TCR"))
  } else {
    message(paste("Folder", vdj_t_path, "does not exist."))
  }
}

# Apply the function to each new_name in the data frame
apply(data, 1, function(row) rename_subfolders(row["new_names"]))

# Function to append _entry to every subfolder within OLD
append_entry_to_subfolders <- function(base_folder) {
  old_path <- file.path(base_folder, "OLD")
  
  if (file.exists(old_path)) {
    # List subfolders in the OLD folder
    subfolders <- list.dirs(old_path, full.names = TRUE, recursive = FALSE)
    
    # Append _entry to each subfolder name
    for (subfolder in subfolders) {
      subfolder_name <- basename(subfolder)
      new_subfolder_name <- paste0(subfolder_name, "_entry")
      new_path <- file.path(old_path, new_subfolder_name)
      file.rename(subfolder, new_path)
    }
  } else {
    message(paste("Folder", old_path, "does not exist."))
  }
}

# Set your base directory
base_directory <- "C:/Users/axi313/Documents/TARA_Entry/Raw_Data"
append_entry_to_subfolders(base_directory)

###
# Function to create per_sample_outs subfolder and move all other subfolders into it
move_subfolders_to_per_sample_outs <- function(base_folder) {
  old_path <- file.path(base_folder, "OLD")
  
  if (file.exists(old_path)) {
    # List all subfolders in the OLD folder
    main_folders <- list.dirs(old_path, full.names = TRUE, recursive = FALSE)
    
    for (main_folder in main_folders) {
      # Create the per_sample_outs subfolder
      per_sample_outs_path <- file.path(main_folder, "per_sample_outs")
      dir.create(per_sample_outs_path, showWarnings = FALSE)
      
      # List all subfolders within the main folder
      subfolders <- list.dirs(main_folder, full.names = TRUE, recursive = FALSE)
      subfolders <- subfolders[subfolders != per_sample_outs_path]  # Exclude the new per_sample_outs folder
      
      # Move each subfolder into the per_sample_outs subfolder
      for (subfolder in subfolders) {
        subfolder_name <- basename(subfolder)
        new_path <- file.path(per_sample_outs_path, subfolder_name)
        file.rename(subfolder, new_path)
      }
    }
  } else {
    message(paste("Folder", old_path, "does not exist."))
  }
}

# Set your base directory
base_directory <- "C:/Users/axi313/Documents/TARA_Entry/Raw_Data"
# Call the function to move subfolders into per_sample_outs within each folder in OLD
move_subfolders_to_per_sample_outs(base_directory)

###
# Function to create the folder structure and move the file
move_raw_feature_bc_matrix <- function(base_folder) {
  old_path <- file.path(base_folder, "OLD")
  
  if (file.exists(old_path)) {
    # List all subfolders in the OLD folder
    main_folders <- list.dirs(old_path, full.names = TRUE, recursive = FALSE)
    
    for (main_folder in main_folders) {
      # Create the multi and GEX_FB subfolders
      multi_path <- file.path(main_folder, "multi")
      gex_fb_path <- file.path(multi_path, "GEX_FB")
      dir.create(gex_fb_path, recursive = TRUE, showWarnings = FALSE)
      
      # Define the source and destination file paths
      source_file <- file.path(main_folder, "per_sample_outs", "GEX_FB", "raw_feature_bc_matrix.h5")
      destination_file <- file.path(gex_fb_path, "raw_feature_bc_matrix.h5")
      
      # Move the file if it exists
      if (file.exists(source_file)) {
        file.rename(source_file, destination_file)
      } else {
        message(paste("File", source_file, "does not exist."))
      }
    }
  } else {
    message(paste("Folder", old_path, "does not exist."))
  }
}


# Call the function to create the folder structure and move the file
move_raw_feature_bc_matrix(base_directory)

