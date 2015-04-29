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
