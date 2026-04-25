# Build script with dependency checking
setwd("D:/QClaw")

# Step 1: Generate documentation
message("=== Step 1: Generating documentation ===")
if (requireNamespace("roxygen2", quietly = TRUE)) {
  roxygen2::roxygenise("OmniSciKit")
  message("✅ Documentation generated successfully")
} else {
  message("❌ roxygen2 not available. Install with: install.packages('roxygen2')")
  quit(status = 1)
}

# Step 2: Build package
message("\n=== Step 2: Building package ===")
result <- system("R CMD build OmniSciKit", intern = TRUE)
cat(paste(result, collapse = "\n"), "\n")

# Step 3: Check if build succeeded
files <- list.files(pattern = "OmniSciKit_.*\\.tar\\.gz")
if (length(files) > 0) {
  file_info <- file.info(files[1])
  message("\n✅ SUCCESS! Package built: ", files[1])
  message("📦 File size: ", round(file_info$size / 1024, 2), " KB")
  message("📅 Built at: ", file_info$mtime)
  
  # Step 4: Optional - check the package
  message("\n=== Step 3: Running R CMD check (optional) ===")
  message("Run the following to check the package:")
  message("R CMD check ", files[1])
} else {
  message("\n❌ Build failed - no .tar.gz file found")
}
