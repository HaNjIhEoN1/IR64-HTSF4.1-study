library(dplyr)
library(tidyr)
library(readr)

# Usage:
#   Rscript make_IR64HTSF4_origin_bed.R \
#     SJ_10kb.f.coords N22_10kb.f.coords IR64HTSF4_origin.bed
#
# If arguments are omitted, the three filenames above are used in the
# current working directory.
args <- commandArgs(trailingOnly = TRUE)

ir64_file <- if (length(args) >= 1) args[1] else "SJ_10kb.f.coords"
n22_file  <- if (length(args) >= 2) args[2] else "N22_10kb.f.coords"
output_bed <- if (length(args) >= 3) args[3] else "IR64HTSF4_origin.bed"

coords_columns <- paste0("S", 1:11)

read_coords <- function(path) {
  read_table(
    path,
    col_names = coords_columns,
    show_col_types = FALSE,
    progress = FALSE
  )
}

# Retain the highest alignment identity (S7) for each 10-kb query window (S11).
ir64_best <- read_coords(ir64_file) %>%
  group_by(S11) %>%
  summarise(IDY_ir64 = max(S7, na.rm = TRUE), .groups = "drop")

n22_best <- read_coords(n22_file) %>%
  group_by(S11) %>%
  summarise(IDY_n22 = max(S7, na.rm = TRUE), .groups = "drop")

# Assign each window to the reference with the higher alignment identity.
# Ties are assigned to IR64, following the original script.
origin_bed <- full_join(ir64_best, n22_best, by = "S11") %>%
  mutate(
    origin = case_when(
      is.na(IDY_ir64) ~ "N22",
      is.na(IDY_n22) ~ "IR64",
      IDY_ir64 >= IDY_n22 ~ "IR64",
      TRUE ~ "N22"
    )
  ) %>%
  extract(
    S11,
    into = c("chrom", "start", "end"),
    regex = "^(IR64-htsf4_Chr\\d+):(\\d+)-(\\d+)$",
    convert = TRUE,
    remove = TRUE
  )

if (anyNA(origin_bed[, c("chrom", "start", "end")])) {
  stop("Some S11 window names do not match IR64-htsf4_ChrNN:start-end.")
}

origin_bed <- origin_bed %>%
  mutate(chr_number = as.integer(sub("^IR64-htsf4_Chr", "", chrom))) %>%
  arrange(chr_number, start, end) %>%
  select(chrom, start, end, origin)

# Four-column, headerless BED-like output:
# chromosome, start, end, inferred origin
write_tsv(origin_bed, output_bed, col_names = FALSE, na = "")

message("Wrote ", nrow(origin_bed), " windows to ", output_bed)
