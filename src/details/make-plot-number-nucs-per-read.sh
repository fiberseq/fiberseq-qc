#!/bin/bash -efx
# author : sjn
# date : Aug 22, 2022

# histogram of per read #nucleosomes
# add a line and values for the median, 10%ile and 90%ile

set -euo pipefail

if [ $# != 4 ]; then
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

ftype=nuc
tmpd=${TMPDIR}/$(whoami)/$$
rm -rf ${tmpd}
mkdir -p ${tmpd}
mkdir -p $(dirname "${outpdf}")
mkdir -p $(dirname "${outstat}")

BASEDIR=$(dirname "$0")
${BASEDIR}/cutnm nuc_lengths ${inp} |
  awk 'NR > 1' |
  awk '$1 != "."' |
  rev |
  sed 's;,;;' |
  rev |
  awk '{ print split($1, a, ",") }' |
  sort -gk1,1 |
  uniq -c |
  awk '{ print $2"\t"$1 }' \
 >${tmpd}/${samplenm}.${ftype}

R --no-save --quiet <<__R__
  # 0.0 <= quantile <= 1.0
  fast_kth <- function(array_2d, lower_b, upper_b, quantile) {
    scores <- subset(array_2d, V1>=lower_b & V1<=upper_b)
    sn <- sum(scores[,2])
    marker <- round(sn * quantile, 0)
    i <- 1
    cntr <- 0
    repeat {
      cntr <- cntr + scores[i,2]
      if ( cntr >= marker ) {
        break
      }
      i <- i + 1
    }
    marker <- scores[i,1]
    return(marker)
  }

  s <- read.table("$tmpd/$samplenm.$ftype", header=FALSE, sep="\t", row.names=NULL)
  mxh <- 150

  f <- subset(s, V1>mxh)
  g <- subset(s, V1<mxh)
  h <- subset(s, V1==mxh)
  h[,2] <- h[,2] + sum(f[,2])
  p <- round(100*sum(f[,2])/sum(s[,2]), 2)
  s <- rbind(g, h)
  mxc <- max(s[,2])

  scores_10 <- fast_kth(s, 0, mxh, 0.1)
  scores_50 <- fast_kth(s, 0, mxh, 0.5)
  scores_90 <- fast_kth(s, 0, mxh, 0.9)

  mycol <- "darkgreen"
  pdf("$outpdf")
  plot(s, axes=F, xlim=c(0,mxh), type="h", main="$samplenm", xlab="# nucleosomes per read", ylab="Count")
  abline(v=scores_10, col=mycol, lty=1)
  abline(v=scores_50, col=mycol, lty=1)
  abline(v=scores_90, col=mycol, lty=1)

  rtoff <- 4
  div <- 2
  msg1 <- paste(p, "% > ", mxh, sep="")
  msg2 <- paste(scores_10)
  msg3 <- paste(scores_50)
  msg4 <- paste(scores_90)

  text(mxh-20, mxc/div, msg1, col=mycol)
  text(scores_10-rtoff, mxc, msg2, col=mycol)
  text(scores_50+rtoff, mxc, msg3, col=mycol)
  text(scores_90+rtoff, mxc, msg4, col=mycol)

  axis(1)
  axis(2)

  dev.off()

  stats_file <- "$outstat"
  cat("# Note: ***per read number of nucleosomes***\n", file=stats_file, append=FALSE)
  cat("Percent(NucsPerRead)", ">", mxh, "=", p, "%\n", file=stats_file, sep="", append=TRUE)
  cat("Quantile10%(NucsPerRead)=", scores_10, "\n", file=stats_file, sep="", append=TRUE)
  cat("Median(NucsPerRead)=", scores_50, "\n", file=stats_file, sep="", append=TRUE)
  cat("Quantile90%(NucsPerRead)=", scores_90, "\n", file=stats_file, sep="", append=TRUE)
__R__

rm -rf ${tmpd}

exit 0
