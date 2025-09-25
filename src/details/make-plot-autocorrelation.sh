#!/bin/bash
# author : Mitchell Vollger and sjn

set -exuo pipefail

if [[ $# != 5 ]]; then
  printf "Expect $0 <sample-name> <input-file> <ec-filter-boolean> <output-pdf> <output-stat.txt>\n"
  exit 1
fi

samplenm=$1
inp=$2 # fiberseq.bam
ecfilter=$3 # set to 0 (or neg) to remove filtering by ec
outpdf=$4
outstat=$5

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

# Mitchell suggests using 10k reads or fewer
nreads=5000
(set +eo pipefail && samtools view -h --subsample-seed 41 --subsample 0.01 $inp \
  | awk -v nrds=$nreads 'BEGIN {i=0} ; { if($1 ~ /^@/) print; else if (++i <= nrds) print; else { exit 0; } }' \
  | samtools view -b \
  | ft extract --all - \
  > $tmpd/$samplenm.$ftype)

R --no-save --quiet <<__R__
library(ggplot2)
library(tidyverse)
library(data.table)
library(scales)
library(ggforce)
library(cowplot)
library(dplyr)
library(splitstackshape)
library(ggridges)
library(IRanges)
library(ggrepel)
library(ggnewscale)
library(ggside)
library(glue)
library("tidylog", warn.conflicts = FALSE)
library(patchwork)

read_m6a = function(file, my_tag = "", min_ml = 200, nrows = Inf, ref = TRUE, ecfilt = 3.9){
    tmp = fread(glue(file), nrows = nrows) %>%
      filter(en - st > 0.5 * fiber_length | (en == 0 & st == 0)) %>%
      { if(ecfilt > 0) filter(., ec > ecfilt) else . }

    if ("sam_flag" %in% colnames(tmp)){
        print("filtering by sam flag")
        tmp = tmp %>% filter(sam_flag <= 16)
    }

    tmp = tmp %>% 
        mutate(
            fake_end = dplyr::case_when(
                en == 0 ~ 11,
                TRUE ~ as.numeric(en)
            ),
            st=as.numeric(st),
            index = row_number(),
            bin = IRanges::disjointBins(
                IRanges(st+1, fake_end) + 150
            ),
        )

    if(ref){
        m6a = tmp %>%
            select(ref_m6a, fiber, fiber_length, m6a_qual, bin, st, en, strand) %>%
            filter(ref_m6a!=".") %>%
            cSplit(c("m6a_qual", "ref_m6a"), direction = "long") %>%
            filter(m6a_qual > min_ml) %>%
            mutate(type = "m6A", start = ref_m6a, end = ref_m6a + 1, alpha = 1.0, size = 2)

        nuc = tmp %>%
            select(fiber, fiber_length, ref_nuc_starts, ref_nuc_lengths, bin, st, en, strand) %>%
            filter(ref_nuc_starts!=".") %>%
            cSplit(c("ref_nuc_starts", "ref_nuc_lengths"), direction = "long") %>%
            mutate(type = "Nucleosome", start = ref_nuc_starts, end = ref_nuc_starts + ref_nuc_lengths, alpha=0.8, size=1)

        msp = tmp %>%
            select(fiber, fiber_length, ref_msp_starts, ref_msp_lengths, bin, st, en, strand) %>%
            filter(ref_msp_starts!=".") %>%
            cSplit(c("ref_msp_starts", "ref_msp_lengths"), direction = "long") %>%
            mutate(type="MSP", start=ref_msp_starts, end=ref_msp_starts + ref_msp_lengths, alpha=0.8, size=1)
    } else {
        m6a = tmp %>%
            select(m6a, fiber, fiber_length, m6a_qual, bin, st, en, strand) %>%
            filter(m6a!=".") %>%
            cSplit(c("m6a_qual", "m6a"), direction = "long") %>%
            filter(m6a_qual > min_ml) %>%
            mutate(type="m6A", start=m6a, end=m6a + 1, alpha=1.0, size=2)

        nuc = tmp %>%
            select(fiber, fiber_length, nuc_starts, nuc_lengths, bin, st, en, strand) %>%
            filter(nuc_starts!=".") %>%
            cSplit(c("nuc_starts", "nuc_lengths"), direction = "long") %>%
            mutate(type="Nucleosome", start=nuc_starts, end=nuc_starts+nuc_lengths, alpha=0.8, size=1)

        msp = tmp %>%
            select(fiber, fiber_length, msp_starts, msp_lengths, bin, st, en, strand) %>%
            filter(msp_starts!=".") %>%
            cSplit(c("msp_starts", "msp_lengths"), direction = "long") %>%
            mutate(type="MSP", start=msp_starts, end=msp_starts + msp_lengths, alpha=0.8, size=1)
    }

    print(my_tag)
    print(dim(tmp))
    print(dim(m6a))
    bind_rows(list(nuc, m6a, msp)) %>%
        mutate(tag = my_tag) %>%
        filter(start != -1) %>%
        group_by(type, fiber) %>%
        arrange(start) %>%
        mutate(dist = start - lag(start)) %>%
        data.table()
}

ecf <- if (${ecfilter} > 0) 3.9 else 0
dists_df = read_m6a("$tmpd/$samplenm.$ftype", my_tag="$samplenm", ref=FALSE, ecfilt=ecf)
dists_df[
  order(start),
  c("dist", "count_per_fiber"):= list(start - lag(start), .N),
  by=list(tag, type, fiber)
]
dists_df = dists_df[count_per_fiber/fiber_length > 0.02 & type == "m6A"]
dists_df = dists_df[order(tag, type, fiber, start)]

if(length(unique(dists_df[["fiber"]])) < 50) {
    pdf("$outpdf")
    plot.new()
    text(0.5, 0.5, "Not enough reads with features for autocorrelation", cex=1.5)
    dev.off()

    # write stats
    stats_file <- "$outstat"
    cat("# Note: ***Autocorrelation stats***\n", file=stats_file, append=FALSE)
    cat("# Stats:", "$stat_name", "\n", file=stats_file, sep="", append=TRUE)
    cat("X[Y==0]=NaN", "\n", file=stats_file, sep="", append=TRUE)
    cat("\n", file=stats_file, append=TRUE)
} else {
    LAG = 220

    tdf <- dists_df %>%
      filter(type == "m6A") %>%
      group_by(tag, fiber) %>%
      reframe(
        count = n(),
        index = if(n() > 0 && all(is.finite(start))) list(seq(0, max(start, na.rm = TRUE)) %in% start) else list(integer(0))
      ) %>%
      group_by(tag) %>%
      summarise(
        auto = acf(as.numeric(unlist(index)), lag.max = LAG, plot = FALSE)[["acf"]],
        .groups = "drop"
      ) %>%
      ungroup() %>%
      mutate(lag = seq_len(n())) %>%
      ungroup() %>%
      as.data.table()  # fully convert to data.table

    # safe in-place assignment
    tdf[, change := c(0, (auto[1:(.N-1)] * auto[2:.N]) <= 0)]
    xcross <- tdf[change == 1, lag]

    # write stats
    stats_file <- "$outstat"
    cat("# Note: ***Autocorrelation stats***\n", file=stats_file, append=FALSE)
    cat("# Stats:", "$stat_name", "\n", file=stats_file, sep="", append=TRUE)
    cat(paste("X[Y==0]=", paste(xcross, collapse=","), "\n", sep=""), file=stats_file, append=TRUE)
    cat("\n", file=stats_file, append=TRUE)

    # make the autocorrelation plot
    pdf("$outpdf")
    p <- tdf %>%
      filter(lag > 25) %>%
      ggplot(aes(x=lag, y=auto, color=tag)) +
      geom_hline(aes(yintercept=0), color="darkblue", size=1, linetype="dashed") +
      geom_vline(data=NULL, aes(xintercept=147), linetype="dashed", color="black", alpha=0.5) +
      geom_line() +
      geom_label_repel(
          data = . %>% filter(change == 1),
          aes(y=0, x = lag, label=lag),
          min.segment.length = 0,
          nudge_x=5,
          nudge_y=0.01,
          show.legend = FALSE
      ) +
      theme_minimal_grid() +
      scale_y_continuous("Autocorrelation between m6A events") +
      scale_x_continuous("Lag between m6A events") +
      guides(color = guide_legend(override.aes = list(size = 2, shape="") ) ) +
      theme(
          legend.position = "top",
          legend.text=element_text(size=8)
      )
    print(p)
    dev.off()
}

__R__

rm -rf ${tmpd}
exit 0
