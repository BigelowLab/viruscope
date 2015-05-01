# viralscan

# Install
```
git clone git@github.com:BigelowLab/viralscan.git
cd viralscan
python setup.py install
```

# Requires
+ [bedtools](https://github.com/arq5x/bedtools2)
+ blastp
+ [diamond](http://ab.inf.uni-tuebingen.de/software/diamond/)
+ [prodigal](https://github.com/hyattpd/Prodigal)
+ [tRNAscan-SE](http://selab.janelia.org/tRNAscan-SE/)
    + [Fixing](http://happykhan.com/getting-trnascan-to-work-on-linux.html) a common issue when installing

# Usage
```
$ viralscan -n AAA015-L03 -t 40 --verbose AAA015-L03_2kb_pass_contigs.fasta /dev/shm/vscan_test /mnt/scgc_nfs/ref/viral_dbs/*.fasta.gz
Running Prodigal on AAA015-L03_2kb_pass_contigs.fasta
Running tRNAscan on AAA015-L03_2kb_pass_contigs.fasta
Calculating GC content of AAA015-L03_2kb_pass_contigs.fasta
Running blastp on /dev/shm/vscan_test/prodigal/AAA015-L03_proteins.fasta using:
    -db nr
    -num_alignments 10
    -evalue 0.001000
    -num_threads 40
```


[setup]
input_path=/mnt/scgc_nfs/lab/jlabonte/Test_FR_methods/Test_pipeline_diamond/
name=AAA164A21
output_path=/mnt/scgc_nfs/lab/jlabonte/Test_FR_methods/Test_pipeline_diamond/Summary
[inputs]
fasta_file=Sequence/AAA164A21.fasta #original
gc_content_file=Output/AAA164A21/AAA164A21_gc_content.tsv.gz
proteins_file=Output/AAA164A21/prodigal/AAA164A21_proteins.fasta
blastp_file=Output/AAA164A21/AAA164A21_blastp.tsv
trna_file=Output/AAA164A21/tRNAscan/AAA164A21-tRNAscan.txt
tetramer_file=Output/AAA164A21/tetramerPCA/AAA164A21-tetramer-PC.csv
[genes]
name=genes
tRNA=blue
hypothetical=black
viral=red
viral2=darkred
[Coverage]
log=false
color=black
[GC]
percent=blue
skew=red
log=false
[Tetramer PC]
PC2=blue
PC1=red
log=false
[Similarity_1]
relative_path=true
name1=POV
name2=LineP
dir1=Output/AAA164A21/diamond/
dir2=Output/AAA164A21/diamond/
flag1=POV.tsv.gz
flag2=LineP-all.tsv.gz
log=false
color1=blue
color2=red
reads1=5922080
reads2=8279226
[Pileup_1]
relative_path=true
name=POV
dir=Output/AAA164A21/coverage/
flag=POV.tsv.gz
log=false
color=blue
[Pileup_2]
relative_path=true
name=LineP
dir=Output/AAA164A21/coverage/
flag=LineP-all.tsv.gz
log=false
color=red
