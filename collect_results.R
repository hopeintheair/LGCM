library(dplyr)

files <- list.files("results", pattern = "\\.RData$", full.names = TRUE)

all_summary <- lapply(files, function(f) {
  env <- new.env()
  load(f, envir = env)
  tbl <- env$result[["summary_table"]]
  tbl$model      <- env$row$model
  tbl$rho        <- env$row$rho_val
  tbl$source_file <- basename(f)
  tbl
}) |> bind_rows()

print(all_summary)

# CSV로도 저장
write.csv(all_summary, "results/all_summary.csv", row.names = FALSE)
cat("Saved: results/all_summary.csv\n")
