# fiberseq-qc
QC plots and statistics

<hr/>
<h4>External Software Requirements</h4>
# install external dependencies using:<br/>
<code>conda create -n fiberseq-qc</code><br/>
<code>mamba env update -n fiberseq-qc --file env/qc.yaml</code><br/><br/>
Then, activate your environment<br/>
<code>conda activate fiberseq-qc</code>

<hr/>
<h4>Example Run</h4>
<code>runall-qc.tcsh my-output-dir PS00246 fiberseq.bam fiberseq.all.tbl.gz</code>
<br/><br/>where fiberseq.bam and fiberseq.all.tbl.gz were produced using the fiberseq-smk or fibertools-rs pipeline<br/>

<hr/>
<h4>Output</h4>
The main output is a page that ends in overview.html that you can open in your browser.<br/>
On the left side is a series of pertinent histogram images with statistics.<br/>
On the right side is a blow-up of any image that you mouse over on the left.<br/>
<br/><br/>
<img src="./share/ss-full.png">

