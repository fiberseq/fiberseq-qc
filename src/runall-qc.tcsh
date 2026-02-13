#!/bin/tcsh -ef
# author : sjn

if ( $#argv != 3 ) then
  printf "%s\n" $0:t
  printf "  <base-output-dir>\n"
  printf "  <sample-name>\n"
  printf "  <fiberseq.bam>\n"
  exit -1
endif

# check if correct conda environment is activated
set req_env = "fiberseq-qc"
set env_check = `env | grep CONDA_DEFAULT_ENV | awk -F"=" '{ print $2 }'`
if ( "$env_check" != "$req_env" ) then
  printf "Require conda environment activated: %s\n" $req_env
  exit -1
endif

set echo

set baseoutd = $1
set samplenm = $2
set baminp   = $3 # <samplenm>.fiberseq.bam with corresponding bai file

if (! -s $baminp ) then
  printf "Bad <fiberseq.bam> file given:\n  %s\n" $baminp
  exit -1
endif

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

# PacBio or ONT?
set tech = "ignore"
set tcntr = `samtools view -H $baminp | grep -m1 "PL:" |& tr '\t' '\n' | awk '$1 ~ /^PL:/' | cut -f2 -d':' | tr '[[:upper:]]' '[[:lower:]]'`
if ( "$tcntr" != "" ) then
  set tech = "$tcntr"
endif

if ( "$tech" != "pacbio" && "$tech" != "ont" && "$tech" != "ignore" ) then
  printf "Unknown tech found: %s.\nExpect value PacBio or ONT\n" $tech
  exit -1
endif

@ methrate_div = 1
@ ecfilter = 1
if ( "$tech" != "pacbio" ) then
  # assume ONT -> no PL: flag
  # Mitchell:
  #  ONT sequences one strand, so we will have either As or Ts on a given molecule coming out of ft extract (depending on forward or reverse alignment).
  #   dividing the total AT count by 2 across the experiment will give a very good estimate. But it may be misleading in super small datasets.
  @ methrate_div = 2
  @ ecfilter = 0
endif


# create table; cut out what is needed to save disk space
set table = $tmpd/$samplenm.fiberseq.all.tbl
ft extract -t 8 $baminp --all - \
  | $src_dir/details/cutnm 5mC,ec,fiber,fiber_length,m6a,msp_lengths,msp_starts,nuc_lengths,nuc_starts,rq,total_5mC_bp,total_AT_bp,total_m6a_bp \
 >! $table

# get #reads
(samtools view -c $baminp >! $tmpd/nreads.txt) &

# msp
set stat_name = 'msp_lengths'
($src_dir/details/make-plot-msp-lengths.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.$stat_name.pdf \
  $baseoutd/$samplenm.$stat_name.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.$stat_name.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.$stat_name.pdf)

# nuc
set stat_name = 'nuc_lengths'
($src_dir/details/make-plot-nuc-lengths.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.$stat_name.pdf \
  $baseoutd/$samplenm.$stat_name.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.$stat_name.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.$stat_name.pdf)

# m6a
set stat1_name = 'm6a_per_read'
set stat2_name = 'ccs_passes'
($src_dir/details/make-plot-number-m6a-per-read.sh \
  $samplenm \
  $table \
  $methrate_div \
  $baseoutd/$samplenm.$stat1_name.pdf \
  $baseoutd/$samplenm.$stat2_name.pdf \
  $baseoutd/$samplenm.$stat1_name.intermediate.stat.txt \
  $baseoutd/$samplenm.$stat2_name.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.$stat1_name.intermediate.stat.txt $baseoutd/$samplenm.$stat2_name.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.$stat1_name.pdf $baseoutd/$samplenm.$stat2_name.pdf)

# nucs_per_read
set stat_name = 'number_nucs_per_read'
($src_dir/details/make-plot-number-nucs-per-read.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.$stat_name.pdf \
  $baseoutd/$samplenm.$stat_name.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.$stat_name.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.$stat_name.pdf)

# 5mCs_per_read
set stat_name = 'number_cpgs_per_read'
($src_dir/details/make-plot-number-mcpgs-per-read.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.$stat_name.pdf \
  $baseoutd/$samplenm.$stat_name.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.$stat_name.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.$stat_name.pdf)

# readlength_per_nuc
set stat_name = 'readlength_per_nuc'
($src_dir/details/make-plot-readlength-per-nuc.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.$stat_name.pdf \
  $baseoutd/$samplenm.$stat_name.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.$stat_name.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.$stat_name.pdf)

# readlengths
if ( "$tech" == "pacbio" ) then
  set max_scale = 50000
else # ONT
  set max_scale = 250000
endif
set stat_name = 'readlengths'
($src_dir/details/make-plot-readlengths.sh \
  $samplenm \
  $table \
  $max_scale \
  $baseoutd/$samplenm.$stat_name.pdf \
  $baseoutd/$samplenm.$stat_name.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.$stat_name.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.$stat_name.pdf)

# resolution
set stat_name = 'msp_resolution'
($src_dir/details/make-plot-msp-resolution.sh \
  $samplenm \
  $table \
  $baseoutd/$samplenm.$stat_name.pdf \
  $baseoutd/$samplenm.$stat_name.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.$stat_name.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.$stat_name.pdf)

# rq is for pacbio data only
if ( "$tech" == "pacbio" ) then
  set stat_name = 'readquality'
  ($src_dir/details/make-plot-rq.sh \
    $samplenm \
    $table \
    $baseoutd/$samplenm.$stat_name.pdf \
    $baseoutd/$samplenm.$stat_name.intermediate.stat.txt) &

  set statsfs = ($statsfs $baseoutd/$samplenm.$stat_name.intermediate.stat.txt)
  set pdfs = ($pdfs $baseoutd/$samplenm.$stat_name.pdf)
endif

# autocorrelation
set stat_name = 'autocorrelation'
($src_dir/details/make-plot-autocorrelation.sh \
  $samplenm \
  $baminp \
  $ecfilter \
  $baseoutd/$samplenm.$stat_name.pdf \
  $baseoutd/$samplenm.$stat_name.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.$stat_name.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.$stat_name.pdf)

# autocorrelation smoothed stats
set stat_name = 'autocorrelation_smoothed'
($src_dir/details/make-plot-autocorrelation-smoothed.sh \
  $samplenm \
  $baminp \
  $ecfilter \
  $baseoutd/$samplenm.$stat_name.pdf \
  $baseoutd/$samplenm.$stat_name.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.$stat_name.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.$stat_name.pdf)

# randfibers
set stat_name = 'randfibers'
set maxz = 20000
($src_dir/details/make-plot-rand-fibers.sh \
  $samplenm \
  $baminp \
  0 \
  $maxz \
  $baseoutd/$samplenm.$stat_name.pdf \
  $baseoutd/$samplenm.$stat_name.intermediate.stat.txt) &

set statsfs = ($statsfs $baseoutd/$samplenm.$stat_name.intermediate.stat.txt)
set pdfs = ($pdfs $baseoutd/$samplenm.$stat_name.pdf)

# zoomed randfibers : 2k-4k to start
set lastz = 2000
set zooms = (4000)
foreach z ($zooms)
  set znm = `echo $z | numfmt --to-unit=1000 --suffix=K`
  set lastznm = `echo $lastz | numfmt --to-unit=1000 --suffix=K`
  set stat_name = "randfibers."$lastznm"-"$znm
  ($src_dir/details/make-plot-rand-fibers.sh \
    $samplenm \
    $baminp \
    $lastz \
    $z \
    $baseoutd/$samplenm.$stat_name.pdf \
    $baseoutd/$samplenm.$stat_name.intermediate.stat.txt) &

  set statsfs = ($statsfs $baseoutd/$samplenm.$stat_name.intermediate.stat.txt)
  set pdfs = ($pdfs $baseoutd/$samplenm.$stat_name.pdf)
  set lastz = $z
end

wait

if ( ! -s $tmpd/nreads.txt ) then
  printf "Problem with counting reads in %s\n" $baminp
  printf "This could be a truncation problem\n"
  printf "Failed command: samtools view -c <bam>\n"
  exit -1
endif

# combine stats
cat $tmpd/nreads.txt \
  | awk '{ printf "# Note: ***Unaligned Reads***\n# Stats:nreads\nNumber(Reads)=%s\n\n", $1; }' \
  | cat - $statsfs \
 >! $baseoutd/$samplenm.qc_stats.txt

# create html outputs
set nreads = `cat $tmpd/nreads.txt | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
set n50 = `grep -F -e "N50(ReadLength)" $baseoutd/$samplenm.qc_stats.txt | cut -f2 -d'=' | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
set totGB = `grep -F -e "Total_GB(Bases)" $baseoutd/$samplenm.qc_stats.txt | cut -f2 -d'=' | awk '{printf "%.2f\n", $1}'`
set medm6A = `grep -F -e "Median(#m6A/#ATs)" $baseoutd/$samplenm.qc_stats.txt | cut -f2 -d'=' | awk '{printf "%.2f%%\n", $1*100}'`
set s2n = `grep -F -e "Area(Curve)/Noise(Curve)" $baseoutd/$samplenm.qc_stats.txt | cut -f2 -d'=' | awk '{printf "%.2f\n", $1}'`
set blockm6A = `grep -F -e "Percent(NucLength>500bp)" $baseoutd/$samplenm.qc_stats.txt | cut -f2 -d'='`
set samplenm_html = $samplenm":Number_of_reads="$nreads":N50_length_(bp)="$n50":Total_sequencing_yield_(Gbp)="$totGB":Methylation_rate_(m6A/AT)="$medm6A":Chromatin_signal_to_noise_(CSN)="$s2n":Untreated_fiber_rate_(MSP>500bp)="$blockm6A
set fs = ($pdfs $statsfs)
$src_dir/details/make-html.tcsh \
  $samplenm_html \
  $baseoutd/$samplenm.overview.html \
  $baseoutd/$samplenm.qc.html \
  $fs

rm -f $statsfs
rm -rf $tmpd

exit 0
