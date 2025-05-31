library(httr)

args <- commandArgs(trailingOnly = TRUE)
uuid_file <- args[1]
output_dir <- args[2]

dir.create(output_dir, showWarnings = FALSE)

uuids <- readLines(uuid_file)

for (uuid in uuids) {
  out_file <- file.path(output_dir, paste0(uuid, ".tsv"))
  url <- paste0("https://api.gdc.cancer.gov/data/", uuid)

  tryCatch({
    res <- GET(url)
    if (status_code(res) == 200) {
      writeBin(content(res, "raw"), out_file)
      cat("✅ Downloaded", uuid, "\n")
    } else {
      cat("❌ Failed", uuid, "Status:", status_code(res), "\n")
    }
  }, error = function(e) {
    cat("❗ Error downloading", uuid, ":", e$message, "\n")
  })
}
