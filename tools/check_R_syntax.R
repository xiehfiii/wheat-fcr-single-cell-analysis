args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript check_R_syntax.R <directory>")
files <- list.files(args[[1]], pattern = "[.]R$", recursive = TRUE,
                    full.names = TRUE)
failed <- character()
for (file in files) {
  result <- tryCatch({parse(file = file); NULL}, error = function(e) e)
  if (is.null(result)) {
    cat("OK", file, "\n")
  } else {
    cat("FAIL", file, conditionMessage(result), "\n")
    failed <- c(failed, file)
  }
}
cat("Checked", length(files), "R scripts; failures:", length(failed), "\n")
if (length(failed)) quit(status = 1L)
