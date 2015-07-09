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
+ tetramer (R packages, install as described at each link)
    + [Biostrings](http://www.bioconductor.org/packages/release/bioc/html/Biostrings.html)
    + [argparser](https://bitbucket.org/djhshih/argparser)
+ graphsignals (R packages, install as described at each link)
    + [futile.logger](https://github.com/zatonovo/futile.logger)
    + [data.table](https://github.com/Rdatatable/data.table)

# Usage
```
$ viral-scan -h
usage: viral-scan.py [-h] [-v] [-n NAME] [-t THREADS] [-i IDENTITY] [--verbose]
                    [--db DB] [--num-alignments NUM_ALIGNMENTS]
                    [--evalue EVALUE] [--script-path SCRIPT_PATH]
                    [--window-size WINDOW_SIZE] [--step-size STEP_SIZE]
                    fasta output query [query ...]

To be determined.

positional arguments:
  fasta                 Fasta file of sequences to analyze.
  output                Location to write output files.
  query                 Reference fasta(s). Specify as many as you would like
                        or use a file pattern.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -n NAME, --name NAME  Name of sample being processed. (default: None)
  -t THREADS, --threads THREADS
                        Number of threads to use. (default: 12)
  -i IDENTITY, --identity IDENTITY
                        passing blastx percent identity per hit (default: 50)
  --verbose

BLAST options:
  --db DB               BLAST database (default: nr)
  --num-alignments NUM_ALIGNMENTS
                        Number of database sequences to show alignments for.
                        (default: 10)
  --evalue EVALUE       Expectation value threshold for saving hits. (default:
                        0.001)

TetramerPCA options:
  --script-path SCRIPT_PATH
                        tetramer PCA R script path; testing ignore if this
                        option is skipped (default: None)
  --window-size WINDOW_SIZE
  --step-size STEP_SIZE
```

Example:
```
$ viral-scan -n AAA015-L03 \
    -t 40 \
    --verbose \
    AAA015-L03.fasta \
    AAA015-L03_output \
    /mnt/scgc_nfs/ref/viral_dbs/*.fasta.gz
```
