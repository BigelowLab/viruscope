# viralscan
Jessica should be the one writing this up. Simple, 3-sentence description.

# Install
```
git clone git@github.com:BigelowLab/viralscan.git
cd viralscan
python setup.py install
```

# Requires
+ [bedtools](https://github.com/arq5x/bedtools2)
+ [blastp](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST)
+ [diamond](http://ab.inf.uni-tuebingen.de/software/diamond/)
+ [prodigal](https://github.com/hyattpd/Prodigal)
+ [R](http://cran.r-project.org/)
+ [samtools](https://github.com/samtools/samtools)
+ [tRNAscan-SE](http://selab.janelia.org/tRNAscan-SE/)
    + [Fixing](http://happykhan.com/getting-trnascan-to-work-on-linux.html) a common issue when installing
+ tetramer (R packages)
    + [Biostrings](http://www.bioconductor.org/packages/release/bioc/html/Biostrings.html)
    + [argparser](https://bitbucket.org/djhshih/argparser) be sure you have the version from this site
+ graphsignals (R packages)
    + [futile.logger](https://github.com/zatonovo/futile.logger)
    + [data.table](https://github.com/Rdatatable/data.table)

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
Creating DIAMOND database for /dev/shm/vscan_test/prodigal/AAA015-L03_proteins.fasta
Running DIAMOND BLASTX on GoM_DeepSeaSediment_023.fna.gz across AAA015-L03_proteins
Converting DIAMOND database GoM_DeepSeaSediment_023.daa to tabular (GoM_DeepSeaSediment_023.tsv.gz)
Finding coverages per contig based coverage on hits above 50 percent in GoM_DeepSeaSediment_023
Running DIAMOND BLASTX on GoM_DeepSeaSediment_278.fna.gz across AAA015-L03_proteins
Converting DIAMOND database GoM_DeepSeaSediment_278.daa to tabular (GoM_DeepSeaSediment_278.tsv.gz)
Finding coverages per contig based coverage on hits above 50 percent in GoM_DeepSeaSediment_278
Running DIAMOND BLASTX on GoM_DeepSeaSediment_315.fna.gz across AAA015-L03_proteins
Converting DIAMOND database GoM_DeepSeaSediment_315.daa to tabular (GoM_DeepSeaSediment_315.tsv.gz)
Finding coverages per contig based coverage on hits above 50 percent in GoM_DeepSeaSediment_315
Running Tetramer PCA on AAA015-L03_2kb_pass_contigs.fasta
```
