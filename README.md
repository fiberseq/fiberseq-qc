# fiberseq-qc
QC plots and statistics

<hr/>
<h4>Software Setup</h4>
# download the software:<br/>
<code>git clone https://github.com/fiberseq/fiberseq-qc.git</code><br/>
<code>cd fiberseq-qc</code><br/><br/>
# install external dependencies using:<br/>
<code>conda create -n fiberseq-qc</code><br/>
<code>mamba env update -n fiberseq-qc --file env/qc.yaml</code><br/><br/>
Then, activate your environment<br/>
<code>conda activate fiberseq-qc</code>

<hr/>
<h4>Example Run</h4>
<code>src/runall-qc.tcsh my-output-dir my-sample-name fiberseq.ft.bam</code>
<br/><br/>where fiberseq.ft.bam was produced using the <code>ft predict-m6a</code> or <code>ft add-nucleosome</code> command<br/>

<hr/>
<h4>Output</h4>
The main output is a page that ends in overview.html that you can open in your browser.<br/>
On the left side is a series of pertinent histogram and other analysis images with statistics.<br/>
On the right side is a blow-up of any image that you mouse over on the left.<br/>
<br/><br/>
<img src="./share/ss-full.png">
<br/><br/>

<hr/>
<h4>Stats</h4>
You can look at your sample's statistics relative to the distribution of all of our lab's human samples' <a href="https://s3.kopah.uw.edu/prod-reporter/index.html">statistics</a>.<br/>
When loaded, you will be under the 'QC Reports' menu.  Near the top, click on the 'Meta' tab, and then press the 'Upload-QC' button.<br/>
Upload the file that ends with qc_stats.txt.  Your sample's values are shown in red.<br/>
<br/><br/>
<img src="./share/meta.png">
<br/><br/>

Some salient statistics appear at the top of your sample's QC output, which may look like:<br/>
Number of reads=6,081,165<br/>
N50 length (bp)=17,770<br/>
Total sequencing yield (Gbp)=101.16<br/>
Methylation rate (m6A/AT)=7.19%<br/>
Chromatin signal to noise (CSN)=105.30<br/>
Untreated fiber rate (MSP>500bp)=0.8%<br/>
<br/>
In addition to the Number of reads, the final 3 shown are important.<br/>
The Methylation rate shows the proportion of A/T's that were methylated.<br/>
The Untreated fiber rate is related and shows the proportion of fibers with very little methylation (a high value indicates undermethylation).<br/>
The Chromatin signal to noise (CSN) uses autocorrelation and is a metric related to the Nuclei extract step.  You want CSN>10.<br/>
<br/>
Additionally, we like to look closely at the 2 autocorrelation plots and the 2 randfibers plots.  Here are plots from a successful sample.<br/>
<br/><br/>
<img src="./share/autocorr.png">
<br/><br/>
<img src="./share/randfibers.png">
<br/><br/>
The randfibers plots show a bird's eye view of sampled fibers.  The plots are from 0-20 kb and (zoomed in) 2-4 kb.<br/>
The autocorrelation for a successful sample shows a signature curve profile, often crossing at y=0.<br/>
When it does, it often marks the typical size of a nucleosome.<br/>
The second autocorrelation plot is the same but extends out to a lag of 2kb.  We use to this to derive the<br/>
Chromatin signal to noise (CSN) metric.<br/>
