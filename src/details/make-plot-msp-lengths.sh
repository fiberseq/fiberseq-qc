#!/bin/bash
# author : sjn
# date : Aug 22, 2022

set -euox pipefail

if [[ $# != 4 ]]; then
  printf "Expect $0 <sample-name> <input-file> <output-pdf> <output-stat.txt>\n"
  exit 1
fi

samplenm=$1
inp=$2
outpdf=$3
outstat=$4

ftype=msp
tmpd=${TMPDIR}/$(whoami)/$$

if [ ! -s $inp ]; then
  printf "Problem finding 1 file: %s\n" $inp
  exit 1
fi

rm -rf $tmpd
mkdir -p $tmpd
mkdir -p $(dirname "${outpdf}")
mkdir -p $(dirname "${outstat}")

BASEDIR=$(dirname "$0")
zcat $inp |
    ${BASEDIR}/cutnm msp_lengths |
    awk 'NR > 1' |
    awk '$1 != "."' |
    rev |
    sed 's;,;;' |
    rev |
    awk 'BEGIN {OFS="\t"; vals[-1]=-1; delete vals[-1]} ; { \
          lng=split($0, a, ","); \
          for ( i=1; i<=lng; ++i ) { \
            if ( a[i] in vals ) { vals[a[i]] += 1 } \
            else { vals[a[i]] = 1 } \
          } \
        } END { \
          for ( v in vals ) { print v, vals[v] } \
        }' |
    sort -gk1,1 \
        >$tmpd/$samplenm.$ftype

R --no-save --quiet <<__R__
  fast_median <- function(array_2d, lower_b, upper_b) {
    msps <- subset(array_2d, V1>=lower_b & V1<=upper_b)
    sn <- sum(msps[,2])
    if ( sn %% 2 == 1 ) {
      mds <- (sn+1)/2
    } else {
      mds <- (sn/2 + sn/2 + 1)/2
    }
    i <- 1
    cntr <- 0
    repeat {
      cntr <- cntr + msps[i,2]
      if ( cntr >= mds ) {
        break
      }
      i <- i + 1
    }
    mds <- msps[i,1]
    return(mds)
  }

  s <- read.table("$tmpd/$samplenm.$ftype", header=FALSE, sep="\t", row.names=NULL)
  mxh = 400
  mnv <- min(s[,1])
  if ( mnv > 1 ) {
    v <- 1:(mnv-1)
    v <- cbind(v, rep(0, mnv-1))
    colnames(v) <- c("V1", "V2")
    colnames(s) <- c("V1", "V2")
    s <- rbind(v, s)
  }

  length_ones <- subset(s, V1==1)[,2]
  length_all <- sum(s[,2])
  p <- round(100*length_ones/length_all, 1)
  all_greater <- sum(subset(s, V1>150)[,2])
  p2 <- round(100*all_greater/length_all)

  # those msps with size > mxh, are put in the mxh bin
  f <- subset(s, V1>mxh)
  g <- subset(s, V1<mxh)
  h <- subset(s, V1==mxh)
  h[,2] <- h[,2] + sum(f[,2])
  pl <- round(100*sum(f[,2])/sum(s[,2]),1)
  s <- rbind(g, h)

  # want to exclude MSP regions of size <= 1
  ss <- subset(s, V1>1)
  msp_med <- fast_median(ss, 0, mxh)
  stats_file <- "$outstat"
  cat("# Note: ***MSP stats***\n", file=stats_file, append=FALSE)
  cat("MSPs>150=", p2, "%\n", file=stats_file, sep="", append=TRUE)
  cat("Median(MSPs!=1bp)=", msp_med, "\n", file=stats_file, sep="", append=TRUE)

  s <- s[,2]

  mycol <- "darkgreen"
  pdf("$outpdf")
  plot(s, xlim=c(0, mxh), axes=F, type="h", main=paste("$samplenm", "$ftype", sep=":"), xlab="Methyltransferase-sensitive patch (MSP) length (bp)", ylab="Count")
  abline(v=msp_med, col=mycol, lty=1)
  msg <- paste(msp_med)
  msg2 <- paste(p, "% of MSPs are 1 bp", sep="")
  msg3 <- paste(pl, "% > ", mxh, "bp", sep="")

  mxc <- max(s)
  xv <- 75
  offxv <- 8
  rtoff <- 50
  div <- 5
  lines(c(msp_med+rtoff, msp_med), c(mxc/div, 0), col=mycol)
  lines(c(msp_med+rtoff, msp_med+rtoff+offxv), c(mxc/div, mxc/div), col=mycol)
  text(msp_med+rtoff+3*offxv, mxc/div, msg, col=mycol)
  text(mxh/2, mxc/2, msg2, col=mycol)
  text(5/6*mxh, mxc/(div+1), msg3, col=mycol)

  axis(1, seq(0, mxh, 50))
  axis(2)

  dev.off()
__R__

rm -rf $tmpd

exit 0
