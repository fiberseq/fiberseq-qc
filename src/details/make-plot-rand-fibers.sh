#!/bin/bash
# author : sjn
# date : Aug.2023

set -exuo pipefail

if [[ $# != 4 ]]; then
  printf "Expect $0 <sample-name> <input-file> <output-pdf> <output-stat.txt>\n"
  exit 1
fi

samplenm=$1
inp=$2 # fiberseq.bam
outpdf=$3
outstat=$4

ftype=nucs.sample
tmpd=${TMPDIR}/$(whoami)/$$
rm -rf ${tmpd}
mkdir -p ${tmpd}
mkdir -p $(dirname "${outpdf}")
mkdir -p $(dirname "${outstat}")

nreads=100
(set +eo pipefail && samtools view -h --subsample-seed 43 --subsample 0.01 $inp \
  | awk -v nrds=$nreads 'BEGIN {i=0} ; { if($1 ~ /^@/) print; else if (++i <= nrds) print; else { exit 0; } }' \
  | samtools view -b \
  > ${tmpd}/sample.bam)

# plot end of fiber as 10 nt rectangle
cat ${tmpd}/sample.bam \
  | ft extract --all - \
  | cutnm m6a,5mC,msp_starts,msp_lengths,nuc_starts,nuc_lengths,fiber_length,fiber \
  | awk 'NR > 1' \
  | awk 'BEGIN {OFS="\t"; print "Fiber", "Feature", "Start", "End"} ; { \
          n_m6a=split($1, m6a, ","); \
          n_cpg=split($2, cpg, ","); \
          n_msp=split($3, msp, ","); \
          split($4, lmsp, ","); \
          n_nuc=split($5, nuc, ","); \
          split($6, lnuc, ","); \
          f_length=$(NF-1); \
          split($NF, fnm, "/"); \
          $NF=fnm[2]; \
          for(i=1;i<n_m6a;++i) { print $NF, "m6a", m6a[i], m6a[i]+1; } \
          for(i=1;i<n_cpg;++i) { print $NF, "5mc", cpg[i], cpg[i]+1; } \
          for(i=1;i<n_msp;++i) { print $NF, "msp", msp[i], msp[i]+lmsp[i]; } \
          for(i=1;i<n_nuc;++i) { print $NF, "nuc", nuc[i], nuc[i]+lnuc[i]; } \
          print $NF, "xfiber-end", f_length, f_length+10; \
        }' \
  > ${tmpd}/${ftype}.samples

R --no-save --quiet <<__R__
  library(ggplot2)
  library(cowplot)
  library(ggforce)
  library(tidyverse)
  library(RColorBrewer)
  library(ggpubr)

  df <- read.table("${tmpd}/${ftype}.samples", header=TRUE, row.names=NULL, sep="\t")
  custom_colors <- scale_fill_manual(name="Feature", values = brewer.pal(5, "Spectral"))

  stats_file <- "${outstat}"
  cat("# Note: ***Random fiber stats***\n", file=stats_file, append=FALSE)
  cat("Number(Fibers)=", length(unique(df[["Fiber"]])), "\n", file=stats_file, sep="", append=TRUE)

  df <- df %>% 
    mutate(
      y=case_when(
        Feature=="5mc" ~ 1.25,
        Feature=="m6a" ~ 1,
        Feature=="msp" ~ 0.5,
        Feature=="nuc" ~ 0.25,
        Feature=="xfiber-end" ~ 1.25,
        TRUE ~ 0
      ),
    )

  pdf("${outpdf}", height=1+length(unique(df[["Fiber"]])), width=15)

  df %>%
    ggplot(
      aes(
        xmin=Start, xmax=End,
        color=NULL,
        fill=Feature,
        ymin=-y, ymax=y
      ),
    ) +
    facet_col(~Fiber, scales="free_y", strip.position="left") + 
    geom_rect(alpha=1) +
    geom_vline(xintercept=seq(0, max(df[["End"]]), 5000), color="white") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    theme(legend.position="top") +
    custom_colors

  dev.off()
  #ggsave("${outpdf}", height=1+length(unique(df[["Fiber"]])), width=15, limitsize=FALSE)
__R__

rm -rf ${tmpd}

exit 0
