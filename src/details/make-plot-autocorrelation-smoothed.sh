#!/bin/bash
# author : Mitchell Vollger and sjn
# date : July.2025

set -exuo pipefail

if [[ $# != 4 ]]; then
  printf "Expect $0 <sample-name> <input-file> <output-pdf> <output-stat.txt>\n"
  exit 1
fi

samplenm=$1
inp=$2
outpdf=$3
outstat=$4

if [ ! -s ${inp} ]; then
  printf "Problem finding 1 file: %s\n" ${inp}
  exit 1
fi

ftype=m6a.autocor
tmpd=${TMPDIR}/$(whoami)/$$
stat_name=$(basename "$outstat" | cut -f2 -d'.')
rm -rf $tmpd
mkdir -p $tmpd
mkdir -p $(dirname "${outpdf}")
mkdir -p $(dirname "${outstat}")

nreads=5000
(set +eo pipefail && samtools view -h --subsample-seed 41 --subsample 0.01 $inp \
  | awk -v nrds=$nreads 'BEGIN {i=0} ; { if($1 ~ /^@/) print; else if (++i <= nrds) print; else { exit 0; } }' \
  | samtools view -b \
  | ft extract --all - \
  > $tmpd/$samplenm.$ftype)

R --no-save --quiet <<__R__
  library(pracma)
  library(ggplot2)
  library(tidyverse)
  library(data.table)
  library(zoo)

  read_m6a = function(file, my_tag = "", min_ml = 200, nrows=Inf, ref=FALSE){
      tmp = fread(glue::glue(file), nrows=nrows) %>%
          filter(en - st > 0.5 * fiber_length | (en == 0 & st == 0)) %>%
          filter(ec > 3.9)
      if ("sam_flag" %in% colnames(tmp)) {
          tmp = tmp %>% filter(sam_flag <= 16)
      }
      tmp = tmp %>%
          mutate(fake_end = ifelse(en == 0, 11, as.numeric(en)),
                 st = as.numeric(st),
                 index = row_number())
      m6a = tmp %>%
          select(m6a, fiber, fiber_length, m6a_qual, st, en) %>%
          filter(m6a != ".") %>%
          separate_rows(m6a, m6a_qual, sep=",") %>%
          filter(as.numeric(m6a_qual) > min_ml) %>%
          mutate(type = "m6A", start = as.numeric(m6a))
      m6a
  }

  # Read and process methylation calls
  dists_df = read_m6a("$tmpd/$samplenm.$ftype") %>%
    group_by(fiber) %>%
    arrange(start) %>%
    mutate(dist = start - lag(start)) %>%
    ungroup()

  calls_per_fiber = dists_df %>% group_by(fiber) %>% summarise(n = n())
  keep_fibers = calls_per_fiber %>% filter(n > 400)
  dists_df = dists_df %>% filter(fiber %in% keep_fibers[["fiber"]])

  # Compute ACF
  LAG = 2000
  acf_list = dists_df %>%
    group_by(fiber) %>%
    summarise(index = list(seq(0, max(start)) %in% start)) %>%
    mutate(auto = map(index, ~ acf(as.numeric(.x), lag.max=LAG, plot=FALSE)[["acf"]])) %>%
    pull(auto)

  min_len <- min(sapply(acf_list, length))
  acf_list <- lapply(acf_list, function(x) x[1:min_len])
  auto <- Reduce("+", acf_list) / length(acf_list)

  lag = 0:(length(auto)-1)
  tdf = data.table(lag=lag, auto=auto)

  # Pad to lag 2000
  if (max(tdf[["lag"]]) < LAG) {
    padding = data.table(
      lag = (max(tdf[["lag"]])+1):LAG,
      auto = 0
    )
    tdf = rbind(tdf, padding)
  }

  # Smooth ACF
  tdf[["smooth"]] = zoo::rollmean(tdf[["auto"]], k=7, fill=0, align="center")

  # === Peak Detection ===
  lag_cut <- 50
  smoothed_adj <- tdf[["smooth"]]
  smoothed_adj[tdf[["lag"]] < lag_cut] <- 0
  ### smooth a bit more still to reduce spurious peaks; only used for peak-calling
  smoothed_adj <- as.numeric(stats::filter(smoothed_adj, rep(1/3, 3), sides = 2))

  peaks_mat <- pracma::findpeaks(smoothed_adj, minpeakdistance = 10)

  if (!is.null(peaks_mat)) {
    min_width <- 40
    peak_widths <- peaks_mat[, 4] - peaks_mat[, 3]
    wide_enough <- peak_widths >= min_width
    peaks_mat <- peaks_mat[wide_enough, , drop = FALSE]

    # Remove any peak < lag 75
    peak_lags_raw <- tdf[["lag"]][peaks_mat[, 2]]
    keep_idx <- which(peak_lags_raw >= 75)
    peaks_mat <- peaks_mat[keep_idx, , drop=FALSE]
  }

  num_peaks <- if (is.null(peaks_mat)) 0L else nrow(peaks_mat)
  peak_lags <- if (num_peaks == 0) numeric(0) else tdf[["lag"]][peaks_mat[, 2]]

  # === Area Under Curve ===
  asymptote = mean(tdf[lag > 1800 & lag <= 2000][["auto"]], na.rm = TRUE)
  first_crossing = min(tdf[lag > 0 & auto < asymptote][["lag"]])
  tdf[, region := ifelse(lag >= first_crossing, "included", "excluded")]
  area_between = trapz(tdf[region == "included"][["lag"]], abs(tdf[region == "included"][["auto"]] - asymptote))

  # === Peak Area ===
  tdf[["peak_area"]] <- FALSE
  area_peaks <- 0
  if (!is.null(peaks_mat)) {
    for (i in seq_len(nrow(peaks_mat))) {
      idxs <- peaks_mat[i, 3]:peaks_mat[i, 4]
      tdf[["peak_area"]][idxs] <- TRUE
      area_peaks <- area_peaks + trapz(tdf[lag %in% tdf[["lag"]][idxs]][["lag"]], smoothed_adj[idxs])
    }
  }

  # === Total Difference ===
  valid_idx <- which(tdf[["lag"]] > 75)
  total_diff <- sum(abs(tdf[["auto"]][valid_idx] - tdf[["smooth"]][valid_idx]), na.rm = TRUE)

  # === Output Stats ===
  stats_file = "$outstat"
  cat("# Note: ***Autocorrelation stats***\n", file=stats_file, append=FALSE)
  cat("# Stats:", "$stat_name", "\n", file=stats_file, sep="", append=TRUE)
  cat(sprintf("Asymptote(Curve)=%.3f\n", asymptote), file=stats_file, append=TRUE)
  cat(sprintf("Area(Curve)=%.3f\n", area_between), file=stats_file, append=TRUE)
  cat(sprintf("Count(Peaks)=%d\n", num_peaks), file=stats_file, append=TRUE)
  #cat(sprintf("Area(Peaks)=%.3f\n", area_peaks), file=stats_file, append=TRUE)
  cat(sprintf("Noise(Curve)=%.3f\n", total_diff), file=stats_file, append=TRUE)
  cat(sprintf("Area(Curve)/Noise(Curve)=%.3f\n", area_between/(total_diff+0.01)), file=stats_file, append=TRUE)
  cat("\n", file=stats_file, append=TRUE)

  # === Plot ===
  pdf("$outpdf")
  p <- ggplot(tdf, aes(x = lag)) +
    geom_line(aes(y = auto), color = "gray40") +
    geom_line(aes(y = smooth), color = "blue") +
    geom_hline(yintercept = asymptote, color = "red", linetype = "dashed") +
    geom_vline(xintercept = peak_lags, color = "green", linetype = "dotted") +
    geom_ribbon(data = subset(tdf, peak_area), aes(ymin = 0, ymax = smooth), fill = "orange", alpha = 0.4) +
    coord_cartesian(xlim = c(0, 2000)) +
    theme_minimal() +
    labs(title = "Smoothed autocorrelation with peaks", y = "Autocorrelation", x = "Lag") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2.5,
      label = paste0(
        "Peaks=", sprintf("%d", num_peaks), "\n",
        #"PeakArea=", sprintf("%.3f", area_peaks), "\n",
        "Area=", sprintf("%.3f", area_between), "\n",
        "Noise=", sprintf("%.3f", total_diff), "\n",
        "Area/Noise=", sprintf("%.3f", area_between / (total_diff + 0.01))
      ),
      size = 3.5,
      color="darkgreen"
    )
  print(p)
  dev.off()
__R__

rm -rf ${tmpd}

exit 0
