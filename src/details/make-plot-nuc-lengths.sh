#!/bin/bash
# author : sjn
# date : Aug 22, 2022

# histogram of nucleosome lengths
# report stats: median length of 80-220 bp nucs
#               (#80-220bp)/(#250-450bp) mono vs. di-nucleosome ratio

set -euox pipefail

if [[ $# != 4 ]]; then
    printf "Expect $0 <sample-name> <input-file> <output-pdf> <output-stat.txt>\n"
    exit 1
fi

samplenm=$1
inp=$2
outpdf=$3
outstat=$4

if [ ! -s $inp ]; then
    printf "Problem finding 1 file: %s\n" $inp
    exit 1
fi

ftype=nuc
tmpd=${TMPDIR}/$(whoami)/$$
rm -rf $tmpd
mkdir -p $tmpd
mkdir -p $(dirname "${outpdf}")

# in some fiberseq-smk runs, there are 'nucleosomes' of size 1 that are the first and last listed always.  I am throwing out size=1 nucleosomes.
BASEDIR=$(dirname "$0")
zcat $inp |
    ${BASEDIR}/cutnm nuc_lengths |
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
    awk '$1 > 1' |
    sort -gk1,1 \
        >$tmpd/$samplenm.$ftype

R --no-save --quiet <<__R__
  library(scales)

  fast_median <- function(array_2d, lower_b, upper_b) {
    nucs <- subset(array_2d, V1>=lower_b & V1<=upper_b)
    sn <- sum(nucs[,2])
    if ( sn %% 2 == 1 ) {
      mds <- (sn+1)/2
    } else {
      mds <- (sn/2 + sn/2 + 1)/2
    }
    i <- 1
    cntr <- 0
    repeat {
      cntr <- cntr + nucs[i,2]
      if ( cntr >= mds ) {
        break
      }
      i <- i + 1
    }
    mds <- nucs[i,1]
    return(mds)
  }

  s <- read.table("$tmpd/$samplenm.$ftype", header=FALSE, sep="\t", row.names=NULL)
  mxh = 500
  mnv <- min(s[,1])
  if ( mnv > 1 ) {
    v <- 1:(mnv-1)
    v <- cbind(v, rep(0, mnv-1))
    colnames(v) <- c("V1", "V2")
    colnames(s) <- c("V1", "V2")
    s <- rbind(v, s)
  }

  # lower/upper bounds on mono and di-nucleosome sizes
  mono_lb <- 80
  mono_ub <- 220
  di_lb   <- 250
  di_ub   <- 450
  mono_nuc_med <- fast_median(s, mono_lb, mono_ub)
  di_nuc_med <- fast_median(s, di_lb, di_ub)

  # those nucs with size > mxh, are put into the mxh bin
  f <- subset(s, V1>mxh)
  g <- subset(s, V1<mxh)
  h <- subset(s, V1==mxh)
  h[,2] <- h[,2] + sum(f[,2])
  numf <- sum(f[,2])
  s <- rbind(g, h)
  mxc <- max(s[,2])

  stats_file <- "$outstat"
  cat("# Note: ***Nucleosome stats***\n", file=stats_file, sep="", append=FALSE)
  cat("# Note: ", mono_lb, "bp <= mono_nuc <= ", mono_ub, "bp\n", file=stats_file, sep="", append=TRUE)
  cat("# Note: ", di_lb, "bp <= di_nuc <= ", di_ub, "bp\n", file=stats_file, sep="", append=TRUE)
  cat("Median(mono_nuc)=", mono_nuc_med, "\n", file=stats_file, sep="", append=TRUE)
  cat("Median(di_nuc)=", di_nuc_med, "\n", file=stats_file, sep="", append=TRUE)
  ratio <- round(sum(subset(s, V1>=mono_lb & V1<=mono_ub)[,2])/sum(subset(s, V1>=di_lb & V1<=di_ub)[,2]), 4)
  cat("Ratio(#mono_nuc/#di_nuc)=", ratio, "\n", file=stats_file, sep="", append=TRUE)
  r1 <- round(sum(subset(s, V1>=mono_lb & V1<=mono_ub)[,2])/sum(s[,2]), 4)
  cat("Percent(mono/all)=", round(100*r1,1), "%\n", file=stats_file, sep="", append=TRUE)
  r2 <- round(sum(subset(s, V1>=di_lb & V1<=di_ub)[,2])/sum(s[,2]), 4)
  cat("Percent(di/all)=", round(100*r2,1), "%\n", file=stats_file, sep="", append=TRUE)
  cat("Percent(NucLength)>", mxh, "bp=", round(numf/sum(s[,2])*100, 1), "%\n", file=stats_file, sep="", append=TRUE)

  mycol <- "darkgreen"
  pdf("$outpdf")
  pp <- plot(s, xlim=c(0, mxh), axes=F, type="h", 
             main=paste("$samplenm", "$ftype", sep=":"), xlab="Length (bp)", ylab="Count",
             panel.first=c(rect(mono_lb, 0, mono_ub, mxc*r1, col=alpha(mycol, 0.25), border=NA),
                           rect(di_lb, 0, di_ub, mxc*r2, col=alpha(mycol, 0.25), border=NA)
                          )
            )
  lines(c(mono_nuc_med, mono_nuc_med), c(0, mxc*0.95), col=mycol, lty=1)
  lines(c(di_nuc_med, di_nuc_med), c(0, mxc*0.95), col=mycol, lty=1)
  msg1 <- paste(mono_nuc_med)
  msg2 <- paste(di_nuc_med)
  msg3 <- paste(round(numf/sum(s[,2])*100, 1), "%>", mxh, "bp", sep="")

  xv <- 75
  offxv <- 8
  rtoff <- 50
  div <- 5
  lines(c(xv, mono_nuc_med), c(mxc/div, 0), col=mycol)
  lines(c(xv-offxv, xv), c(mxc/div, mxc/div), col=mycol)
  text(xv-3*offxv, mxc/div, msg1, col=mycol)
  lines(c(di_nuc_med+rtoff, di_nuc_med), c(mxc/div, 0), col=mycol)
  lines(c(di_nuc_med+rtoff, di_nuc_med+rtoff+offxv), c(mxc/div, mxc/div), col=mycol)
  text(di_nuc_med+rtoff+3*offxv, mxc/div, msg2, col=mycol)
  text((mono_lb+mono_ub)/2, mxc, paste("%(mono/all)=", round(100*r1,1), "%", sep=""))
  text((di_lb+di_ub)/2, mxc, paste("%(di/all)=", round(100*r2,1), "%", sep=""))
  text(mxh*0.9, mxc*0.5, msg3, col=mycol)

  axis(1, seq(0, mxh, 50))
  axis(2)

  dev.off()
__R__

rm -rf $tmpd

exit 0
