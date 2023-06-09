#!/bin/tcsh -ef
# author : sjn
# date : May.2023

if ( $#argv != 4 ) then
  printf "%s\n" $0:t
  printf "  <base-output-dir>\n"
  printf "  <sample-name>\n"
  printf "  <fiberseq.bam>\n"
  printf "  <fiberseq.all.tbl.gz>\n"
  exit -1
endif

set echo

set baseoutd = $1
set samplenm = $2
set baminp   = $3 # <samplenm>.fiberseq.bam with corresponding bai file
set tableinp = $4 # <samplenm>.fiberseq.all.tbl.gz

mkdir -p $baseoutd

set src_dir = `readlink -f $0`
set src_dir = $src_dir:h

set statsfs = ()
set pdfs = ()

# qc_msp
($src_dir/details/make-plot-msp-lengths.sh \
  $samplenm \
  $tableinp \
  $baseoutd/$samplenm.qc_msp_lengths.pdf \
  $baseoutd/$samplenm.qc_msp_lengths.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_msp_lengths.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_msp_lengths.pdf)

# qc_nuc
($src_dir/details/make-plot-nuc-lengths.sh \
  $samplenm \
  $tableinp \
  $baseoutd/$samplenm.qc_nuc_lengths.pdf \
  $baseoutd/$samplenm.qc_nuc.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_nuc.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_nuc_lengths.pdf)

# qc_m6a
($src_dir/details/make-plot-number-m6a-per-read.sh \
  $samplenm \
  $tableinp \
  $baseoutd/$samplenm.qc_m6a_per_read.pdf \
  $baseoutd/$samplenm.qc_ccs_passes.pdf \
  $baseoutd/$samplenm.qc_m6a.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_m6a.intermediate.stat.txt $baseoutd/$samplenm.qc_m6a.intermediate.stat.txt) # must duplicate
set pdfs = ($pdfs $baseoutd/$samplenm.qc_m6a_per_read.pdf $baseoutd/$samplenm.qc_ccs_passes.pdf)

# qc_nucs_per_read
($src_dir/details/make-plot-number-nucs-per-read.sh \
  $samplenm \
  $tableinp \
  $baseoutd/$samplenm.qc_number_nucs_per_read.pdf \
  $baseoutd/$samplenm.qc_number_nucs_per_read.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_number_nucs_per_read.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_number_nucs_per_read.pdf)

# qc_readlength_per_nuc
($src_dir/details/make-plot-readlength-per-nuc.sh \
  $samplenm \
  $tableinp \
  $baseoutd/$samplenm.qc_readlength_per_nuc.pdf \
  $baseoutd/$samplenm.qc_readlength_per_nuc.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_readlength_per_nuc.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_readlength_per_nuc.pdf)

# qc_readlengths
($src_dir/details/make-plot-readlengths.sh \
  $samplenm \
  $tableinp \
  $baseoutd/$samplenm.qc_readlengths.pdf \
  $baseoutd/$samplenm.qc_readlengths.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_readlengths.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_readlengths.pdf)

# qc_rq
($src_dir/details/make-plot-rq.sh \
  $samplenm \
  $tableinp \
  $baseoutd/$samplenm.qc_readquality.pdf \
  $baseoutd/$samplenm.qc_readquality.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_readquality.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_readquality.pdf)

# qc_autocorrelation
($src_dir/details/make-plot-autocorrelation.sh \
  $samplenm \
  $baminp \
  $baseoutd/$samplenm.qc_autocorrelation.pdf \
  $baseoutd/$samplenm.qc_autocorrelation.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_autocorrelation.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_autocorrelation.pdf)

wait

# qc_combine_stats
cat $statsfs \
 >! $baseoutd/$samplenm.qc_stats.txt

# create html outputs
set fs = ($pdfs $statsfs)
$src_dir/details/make-html.tcsh \
  $samplenm \
  $baseoutd/$samplenm.overview.html \
  $baseoutd/$samplenm.qc.html \
  $fs

foreach s ($statsfs)
  rm -f $s
end

exit 0
