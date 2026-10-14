# List of directories to be created
directories <- c(
  "00_raw_data",
  "01_processed_data",
  "03_scripts",
  "04_graphs",
  "05_tables"
)

# Create each directory
for (dir in directories) {
  dir.create(dir)
}
