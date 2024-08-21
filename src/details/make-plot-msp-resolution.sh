#!/bin/bash
# author : sjn
# date : Aug 22, 2022

set -euox pipefail

if [[ $# != 4 ]]; then
  printf "Expect $0 <sample-name> <input-table> <output-pdf> <output-stat.txt>\n"
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

ftype=msp
tmpd=${TMPDIR}/$(whoami)/$$
rm -rf ${tmpd}
mkdir -p ${tmpd}
mkdir -p $(dirname "${outpdf}")
mkdir -p $(dirname "${outstat}")

# bins set by Andrew
bins="20,40,60,80,100,125,150,175,200,225,250,300,350,400,500"

# Mitchell suggested just using a sample for this qc plot
#  $inp is from unaligned bam and is in arbitrary order
# subshell to get around awk's exit early return, pipe-related problems
#  without exit 0, the whole input is read through
nreads=1000000
BASEDIR=$(dirname "$0")
(${BASEDIR}/cutnm msp_starts,msp_lengths,m6a ${inp} |
  awk 'NR > 1' |
  awk '$1 != "."' |
  awk -v n=${nreads} 'NR <= n {print} NR > n {exit 0}' || true) |
  awk 'BEGIN {OFS="\t"} ; { \
        nmsp=split($1, a, ","); \
        split($2, b, ","); \
        nm6a=split($3, c, ","); \
        lst=1; \
        for ( i=1; i<nmsp; ++i ) { \
          m6a_count = 0; \
          start=a[i]; end=a[i]+b[i]; \
          for ( j=lst; j<nm6a; ++j ) { \
            if ( c[j]>=start ) { \
              if ( c[j]<end ) { \
                m6a_count += 1; \
                lst = j; \
              } else { \
                break; \
              } \
            } \
          } \
          print m6a_count/b[i], b[i]; \
        } \
      }' |
  awk -v bins=${bins} \
      'BEGIN {OFS="\t"; split(bins,b,",")} ; { \
        if($2>=b[1]) { print } \
      }' |
  awk -v bins=${bins} \
      'BEGIN {OFS="\t"; nbins=split(bins,b,",")} ; { \
        for ( i=1; i<=nbins; ++i ) { \
          if($2>=b[i]) { \
            if(i==nbins || $2<b[i+1]) { \
              print $1, i; \
              break; \
            } \
          } \
        } \
      }' \
   >${tmpd}/${samplenm}.${ftype}

R --no-save --quiet <<__R__
  df <- read.table("$tmpd/$samplenm.$ftype", header=FALSE, sep="\t", row.names=NULL)

  # bin labels
  bins_list <- strsplit("$bins", ",")[[1]]
  bin_labels <- c(sapply(seq_along(bins_list)[-length(bins_list)], function(i) {
    paste(bins_list[i], bins_list[i+1], sep="-")
  }), "500+")

  data <- data.frame(
    m6a_count = df[["V1"]],
    bin = factor(df[["V2"]], levels=1:length(bin_labels), labels=bin_labels)
  )

  pdf("$outpdf")
  mycol <- "darkgreen"
  boxplot(m6a_count ~ bin, data = data,
          xlab = "MSP Size", 
          ylab = "Avg. Resolution for Bin", 
          main = "Resolving Power",
          ylim = c(0,0.5),
          col  = mycol,
          outline = FALSE,
          las = 2,
          cex.axis = 0.7,
          pch = 16,
          cex = 0.5)
  dev.off()

  stats_file <- "$outstat"
  target_bin <- "175-200"
  median_value <- median(data[["m6a_count"]][data[["bin"]] == target_bin])
  cat("***MSP Resolution stats***\n", file=stats_file, sep="", append=FALSE)
  cat("Median(", target_bin, ")=", median_value, "\n", file=stats_file, sep="", append=TRUE)
  cat("\n", file=stats_file, append=TRUE)
__R__

rm -rf ${tmpd}

exit 0
