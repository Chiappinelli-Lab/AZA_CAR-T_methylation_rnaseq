# Build the canonical sample sheet for the Haydar mouse methylation array project
# from "Methyl Array Sample Manifest EV.xlsx" (donor/dose/CAR-arm metadata, Sheet1)
# and its IDAT-basename cross-reference (Sheet2). 

suppressMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
})

source(file.path(Sys.getenv("REPO_ROOT"), "R", "paths.R"))
manifest_file <- METH_MANIFEST
idat_dir <- IDAT_DIR

sheet1 <- read_excel(manifest_file, sheet = 1, n_max = 24)
design <- sheet1 %>%
  select(Code, ID, Donor, CAR, Dose) %>%
  mutate(Drug = ifelse(Donor == "EV03", "GSK-3484862", "5-Azacytidine"))

sheet2 <- read_excel(manifest_file, sheet = 2)
basenames <- sheet2 %>%
  filter(!is.na(external_sample_id)) %>%
  mutate(Basename_rel = str_remove(file_name, "_(Grn|Red)\\.idat$")) %>%
  distinct(external_sample_id, Basename_rel) %>%
  rename(Code = external_sample_id)

design <- left_join(design, basenames, by = "Code") %>%
  mutate(
    Sentrix_ID = str_extract(Basename_rel, "^[0-9]+"),
    Basename = file.path(idat_dir, Sentrix_ID, Basename_rel)
  )

stopifnot(
  nrow(design) == 24,
  !any(is.na(design$Basename)),
  all(file.exists(paste0(design$Basename, "_Grn.idat"))),
  all(file.exists(paste0(design$Basename, "_Red.idat")))
)

out_file <- METH_SHEET
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
write.csv(design, out_file, row.names = FALSE)
cat("Wrote", nrow(design), "samples to", out_file, "\n")
