# ViruSCope
A tool specifically developed to identify viral sequences in single amplified genomes. ViruSCope uses a
combination of BLAST annotations, genomic anomalies (GC content and tetramer frequencies), and contrasting fragment
recruitment of viral and bacterial metagenomic reads to identify viral sequences in the generally fragmented single
cell genomes. ViruSCope is a two-step program that was written in python for the data processing (`viruscope`) and R
(`graphsignals.Rscript`) to produce the final graphical and tabular output identifying the viral contigs.

# Install
```
git clone git@github.com:BigelowLab/viruscope.git
cd viruscope
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
    + [class](https://cran.r-project.org/web/packages/class/index.html)

# Usage
## Single or batch  
There are two ways to run viruscope.You can run viruSCope on a single genome, or run in batch mode, which cuts the run time per genome considerably when running many genomes through the workflow. Either workflow identifies similar statistics on potential viral contigs, but the structure of output directories differs between the two.  Batch mode conducts virus sequence prediction and returns statistics on likely viral sequences, while single genome mode returns statistics on likely viral sequences as well as a graphical output.  

## Single Genome Mode
### viruscope

```
$ viruscope -h

positional arguments:
  fasta                 FASTA file of sequences to analyze.
  output                location to write output files
  query                 reference FASTA(s) -- specify as many as you would
                        like or use a file pattern

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -n NAME, --name NAME  name of sample being processed (default: None)
  -t THREADS, --threads THREADS
                        number of threads to use. (default: 12)
  -i IDENTITY, --identity IDENTITY
                        passing blastx percent identity per hit (default: 50)
  --verbose

BLAST options:
  --db DB               BLAST database (default: nr)
  --num-alignments NUM_ALIGNMENTS
                        number of database sequences for which to show
                        alignments (default: 10)
  --evalue EVALUE       expectation value threshold for saving hits (default:
                        0.001)

TetramerPCA options:
  --script-path SCRIPT_PATH
                        tetramer PCA R script path; testing ignored if this
                        option is skipped (default: None)
  --window-size WINDOW_SIZE
                        sequence window size (default: 1600)
  --step-size STEP_SIZE
                        step size for window (default: 200)

Classifier options:
  --training-file TRAINING_FILE
                        classifier training file (default: None)
  --knn-k KNN_K         nearest neighbors for classifier (default: 27)
```

Example:
```
$ viruscope -n AAA015-L03 \
    -t 40 \
    --verbose \
    AAA015-L03.fasta \
    AAA015-L03_output \
    /mnt/scgc_nfs/ref/viral_dbs/*.fasta.gz
```

### graphsignals

graphsignals is a standalone script that draws upon a configuration file generated by `viruscope`.  See the [example](https://github.com/BigelowLab/viruscope/blob/master/example.cfg)

```
$ Rscript --vanilla graphsignals.Rscript -h
useage: Rscript --vanilla graphsignals.Rscript [-h] config_file

optional arguments:
  -h, --help            show this help message and exit

positional arguments:
  config_file           configuration file
```

## Batch mode

Batch mode involves running several different scripts.  This workflow is in development and currently requires manually triggering several different components.

#### ORF setup
command:
```batch-viruscope orf-setup --threads {threads} {path to directory containing input genomes} {viruscope_outdir}```

Identifies and clusters ORFs from all input genomes.  Clustered ORFs are placed in a subdirectory of the viruscope_outdir as such: viruscope_outdir/clustering/for_mica/

#### Comparison of ORFs to nr:
Clustered ORFs are separated into files containing 1000 ORFs each, to be compared to NCBI's NR database.

This step takes clustered ORFs, and compares them to NCBI's nr database using either BLAST or mica, a BLAST accellorator.  This is designed to be used with a job scheduler, so that each sub-set of clustered ORFS may be submitted separately to speed up the process step.  The command for each ORF file is:

```batch-viruscope blast viruscope_outdir/clustering/for_mica/subset_XX.fasta --mica --threads {threads}```

#### Recruitment of microbial and viral metagenomic reads to each input genome using Diamond
Viral and bacterial metagenomic reads are recruited to each genome invidivually.  This step is also designed to be used with a job scheduler for parallellization.  The command for each individual genome is:

```batch-viruscope recruit-single --threads {# threads} --output {viruscope_dir}/diamond/ --sag-contigs {path_to_genomic_contigs} {path_to_prokka_generated_proteins} {path_to_virus_db} {path_to_bac_db}```

#### Summarise Results
Takes results of comparison to NR and metagenomic read recruitments and compares to a training set of viral sequences via kmeans clustering, returns ID as viral or not viral with a summary statistic.

command:
```batch-viruscope summarize {directory_containing_input_genomes} {directory_containing_viruscope_results}```
