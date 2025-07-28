#!/bin/bash
# author : sjn
# date : Aug 22, 2022

set -euo pipefail

if [[ $# != 5 ]]; then
  printf "Expect $0 <sample-name> <input-table> <max-x-scale> <output-pdf> <output-stat.txt>\n"
  exit 1
fi

samplenm=$1
inp=$2 # fiber-all-table.tbl.gz
max_x_scale=$3 # should be an integer
outpdf=$4
outstat=$5

if [ ! -s ${inp} ]; then
  printf "Problem finding 1 file: %s\n" ${inp}
  exit 1
fi

ftype=readlengths
tmpd=${TMPDIR}/$(whoami)/$$
stat_name=$(basename "$outstat" | cut -f2 -d'.')
rm -rf ${tmpd}
mkdir -p ${tmpd}
mkdir -p $(dirname "${outpdf}")
mkdir -p $(dirname "${outstat}")

# putting things in bins of size 10
BASEDIR=$(dirname "$0")
${BASEDIR}/cutnm fiber_length ${inp} |
  awk 'NR > 1' |
  awk '$1 != "."' |
  rev |
  sed 's;,;;' |
  rev |
  awk '{ $1=int($1/10); $1=$1*10; print $1 }' |
  sort -gk1,1 |
  uniq -c |
  awk '{ print $2"\t"$1 }' \
 >${tmpd}/${samplenm}.${ftype}

R --no-save --quiet <<__R__
  # 0.0 <= quantile <= 1.0
  fast_kth <- function(array_2d, lower_b, upper_b, quantile) {
    reads <- subset(array_2d, V1>=lower_b & V1<=upper_b)
    sn <- sum(reads[,2])
    marker <- round(sn * quantile, 0)
    i <- 1
    cntr <- 0
    repeat {
      cntr <- cntr + reads[i,2]
      if ( cntr >= marker ) {
        break
      }
      i <- i + 1
    }
    marker <- reads[i,1]
    return(marker)
  }

  s <- read.table("$tmpd/$samplenm.$ftype", header=FALSE, sep="\t", row.names=NULL)
  mxh <- as.integer("$max_x_scale")

  # those reads with size > mxh, are put in the mxh bin
  f <- subset(s, V1>mxh)
  g <- subset(s, V1<mxh)
  h <- subset(s, V1==mxh)
  h[,2] <- h[,2] + sum(f[,2])
  p <- round(100*dim(f)[1]/dim(s)[1], 2)
  s <- rbind(g, h)

  reads_10 <- fast_kth(s, 0, mxh, 0.1)
  reads_50 <- fast_kth(s, 0, mxh, 0.5)
  reads_90 <- fast_kth(s, 0, mxh, 0.9)

  mycol <- "darkgreen"
  pdf("$outpdf")
  plot(s[,1], s[,2], xlim=c(0, mxh), axes=F, type="h", main=paste("$samplenm", "$ftype", sep=":"), xlab="Unaligned read length (bp)", ylab="Count")
  abline(v=reads_10, col=mycol, lty=1)
  msg1 <- paste(round(reads_10/1000,1), "k", sep="")
  abline(v=reads_50, col=mycol, lty=1)
  msg2 <- paste(round(reads_50/1000,1), "k", sep="")
  abline(v=reads_90, col=mycol, lty=1)
  msg3 <- paste(round(reads_90/1000,1), "k", sep="")
  msg4 <- paste(p, "% > ", paste(mxh/1000, "kb", sep=""), sep="")

  mxc <- max(s[,2])
  xv <- 75
  offxv <- 8
  rtoff <- 2000
  div <- 3

  lines(c(reads_10-rtoff, reads_10), c(mxc/div-strheight(msg1, units="user"), 0), col=mycol)
  text(reads_10-rtoff-3*offxv, mxc/div, msg1, col=mycol, cex=0.9)

  lines(c(reads_50+rtoff, reads_50), c(mxc/div-strheight(msg2, units="user"), 0), col=mycol)
  text(reads_50+rtoff+3*offxv, mxc/div, msg2, col=mycol, cex=0.9)

  lines(c(reads_90+rtoff, reads_90), c(mxc/div-strheight(msg3, units="user"), 0), col=mycol)
  text(reads_90+rtoff+3*offxv, mxc/div, msg3, col=mycol, cex=0.9)

  text(mxh-mxh/16, mxc/div, msg4, col=mycol)

  axis(1, seq(0, mxh, 4000), labels=paste(seq(0, 50, 4), "k", sep=""))
  axis(2)
  gg <- c("Quantile10%", "Median", "Quantile90%")
  gg <- format(gg, justify="left")
  nnum <- prettyNum(c(reads_10, reads_50, reads_90), big.mark=",", scientific=FALSE)
  nnum <- format(nnum, justify="right")
  nnum <- paste(gg, nnum)

  par(cex=0.8)
  legend("topright",
         legend=nnum,
         col=mycol,
         box.lty=0,
        )

  dev.off()

  stats_file <- "$outstat"
  cat("# Note: ***Read length stats***\n", file=stats_file, append=FALSE)
  cat("# Stats:", "$stat_name", "\n", file=stats_file, sep="", append=TRUE)
  cat("Quantile10%(ReadLength)=", reads_10, "\n", file=stats_file, sep="", append=TRUE)
  cat("Median(ReadLength)=", reads_50, "\n", file=stats_file, sep="", append=TRUE)
  cat("Quantile90%(ReadLength)=", reads_90, "\n", file=stats_file, sep="", append=TRUE)
  cat(paste("Percent(ReadLength>", paste(mxh/1000, "kb)", sep=""), "=", p, "%\n", sep=""), file=stats_file, sep="", append=TRUE)
  cat("\n", file=stats_file, append=TRUE)
__R__

rm -rf ${tmpd}

exit 0
