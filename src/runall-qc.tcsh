#!/bin/tcsh -ef
# author : sjn
# date : May.2023

if ( $#argv != 3 ) then
  printf "%s\n" $0:t
  printf "  <base-output-dir>\n"
  printf "  <sample-name>\n"
  printf "  <fiberseq.bam>\n"
  exit -1
endif

set echo

set baseoutd = $1
set samplenm = $2
set baminp   = $3 # <samplenm>.fiberseq.bam with corresponding bai file

if (! $?TMPDIR) then # this should typically be defined already
  setenv TMPDIR "/tmp"
  if ( ! -d $TMPDIR ) then # try again
    setenv TMPDIR $baseoutd/tmp
  endif
endif

set tmpd = $TMPDIR/`whoami`/$$
rm -rf $tmpd
mkdir -p $tmpd

mkdir -p $baseoutd

set src_dir = `readlink -f $0`
set src_dir = $src_dir:h

set statsfs = ()
set pdfs = ()

# create table; cut out what is needed to save disk space
set table = $tmpd/$samplenm.fiberseq.all.tbl
ft extract -t 8 $baminp --all - \
  | $src_dir/details/cutnm 5mC,ec,fiber,fiber_length,m6a,msp_lengths,msp_starts,nuc_lengths,nuc_starts,rq,total_5mC_bp,total_AT_bp,total_m6a_bp \
 >! $table

# get #reads
(samtools view -c $baminp >! $tmpd/nreads.txt) &

# qc_msp
($src_dir/details/make-plot-msp-lengths.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.qc_msp_lengths.pdf \
  $baseoutd/$samplenm.qc_msp_lengths.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_msp_lengths.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_msp_lengths.pdf)

# qc_nuc
($src_dir/details/make-plot-nuc-lengths.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.qc_nuc_lengths.pdf \
  $baseoutd/$samplenm.qc_nuc.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_nuc.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_nuc_lengths.pdf)

# qc_m6a
($src_dir/details/make-plot-number-m6a-per-read.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.qc_m6a_per_read.pdf \
  $baseoutd/$samplenm.qc_ccs_passes.pdf \
  $baseoutd/$samplenm.qc_m6a.intermediate.stat.txt \
  $baseoutd/$samplenm.qc_m6a.intermediate.stat_ccs_passes.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_m6a.intermediate.stat.txt $baseoutd/$samplenm.qc_m6a.intermediate.stat_ccs_passes.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_m6a_per_read.pdf $baseoutd/$samplenm.qc_ccs_passes.pdf)

# qc_nucs_per_read
($src_dir/details/make-plot-number-nucs-per-read.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.qc_number_nucs_per_read.pdf \
  $baseoutd/$samplenm.qc_number_nucs_per_read.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_number_nucs_per_read.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_number_nucs_per_read.pdf)

# qc_5mCs_per_read
($src_dir/details/make-plot-number-mcpgs-per-read.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.qc_number_cpgs_per_read.pdf \
  $baseoutd/$samplenm.qc_number_cpgs_per_read.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_number_cpgs_per_read.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_number_cpgs_per_read.pdf)

# qc_readlength_per_nuc
($src_dir/details/make-plot-readlength-per-nuc.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.qc_readlength_per_nuc.pdf \
  $baseoutd/$samplenm.qc_readlength_per_nuc.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_readlength_per_nuc.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_readlength_per_nuc.pdf)

# qc_readlengths
($src_dir/details/make-plot-readlengths.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.qc_readlengths.pdf \
  $baseoutd/$samplenm.qc_readlengths.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_readlengths.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_readlengths.pdf)

# qc_rq
($src_dir/details/make-plot-rq.sh \
  $samplenm \
  $table \
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

# qc_randfibers
set maxz = 20000
($src_dir/details/make-plot-rand-fibers.sh \
  $samplenm \
  $baminp \
  0 \
  $maxz \
  $baseoutd/$samplenm.qc_randfibers.pdf \
  $baseoutd/$samplenm.qc_randfibers.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.qc_randfibers.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.qc_randfibers.pdf)

# zoomed qc_randfibers : 2k-4k to start
set lastz = 2000
set zooms = (4000)
foreach z ($zooms)
  set znm = `echo $z | numfmt --to-unit=1000 --suffix=K`
  set lastznm = `echo $lastz | numfmt --to-unit=1000 --suffix=K`
  ($src_dir/details/make-plot-rand-fibers.sh \
    $samplenm \
    $baminp \
    $lastz \
    $z \
    $baseoutd/$samplenm.qc_randfibers.$lastznm"-"$znm.pdf \
    $baseoutd/$samplenm.qc_randfibers.$lastznm"-"$znm.intermediate.stat.txt) &

  set statsfs = ($statsfs $baseoutd/$samplenm.qc_randfibers.$lastznm"-"$znm.intermediate.stat.txt)
  set pdfs = ($pdfs $baseoutd/$samplenm.qc_randfibers.$lastznm"-"$znm.pdf)
  set lastz = $z
end

wait

# qc_combine_stats
cat $tmpd/nreads.txt \
  | awk '{ printf "# Note: ***Unaligned Reads***\nNumber(Reads)=%s\n\n", $1; }' \
  | cat - $statsfs \
 >! $baseoutd/$samplenm.qc_stats.txt

# create html outputs
if ( ! -s $tmpd/nreads.txt ) then
  printf "Problem with counting reads in %s\n" $baminp
  printf "This could be a truncation problem\n"
  printf "Failed command: samtools view -c <bam>\n"
  exit -1
endif
set nreads = `cat $tmpd/nreads.txt | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
set samplenm_html = $samplenm":nreads="$nreads
set fs = ($pdfs $statsfs)
$src_dir/details/make-html.tcsh \
  $samplenm_html \
  $baseoutd/$samplenm.overview.html \
  $baseoutd/$samplenm.qc.html \
  $fs

foreach s ($statsfs)
  rm -f $s
end

rm -rf $tmpd

exit 0
