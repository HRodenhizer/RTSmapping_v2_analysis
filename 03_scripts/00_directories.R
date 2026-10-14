# List of directories to be created
directories <- c(
  "00_raw_data",
  "01_processed_data",
  "03_scripts",
  "03_graphs",
  "04_tables"
)

# Create each directory
for (dir in directories) {
  dir.create(dir)
}

# List of script directories to be created
script_directories <- c(
  "00_functions",
  "01_data_prep",
  "03_analysis"
)

# Create each directory
for (dir in script_directories) {
  dir.create(paste0("03_scripts/", dir))
}
