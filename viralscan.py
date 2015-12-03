"""
viral-SCan, a tool specifically developed to identify viral sequences in single amplified genomes. Viral-SCan uses a
combination of BLAST annotations, genomic anomalies (GC content and tetramer frequencies), and contrasting fragment
recruitment of viral and bacterial metagenomic reads to identify viral sequences in the generally fragmented single
cell genomes.
"""

from __future__ import print_function

import argparse
import contextlib
import gzip
import itertools
import os
import shutil
import subprocess
import sys
import tempfile
from distutils.spawn import find_executable

__version__ = "0.1.7"
REQUIRES = ["bedtools", "samtools", "prodigal", "tRNAscan-SE", "blastp",
            "diamond", "gzip", "gunzip", "Rscript"]


@contextlib.contextmanager
def file_transaction(*rollback_files):
    """
    Wrap file generation in a transaction, moving to output if finishes.
    """
    safe_names, orig_names = _flatten_plus_safe(rollback_files)
    # remove any half-finished transactions
    remove_files(safe_names)
    try:
        if len(safe_names) == 1:
            yield safe_names[0]
        else:
            yield tuple(safe_names)
    # failure -- delete any temporary files
    except:
        remove_files(safe_names)
        remove_tmpdirs(safe_names)
        raise
    # worked -- move the temporary files to permanent location
    else:
        for safe, orig in zip(safe_names, orig_names):
            if os.path.exists(safe):
                shutil.move(safe, orig)
        remove_tmpdirs(safe_names)


def remove_tmpdirs(fnames):
    for x in fnames:
        xdir = os.path.dirname(os.path.abspath(x))
        if xdir and os.path.exists(xdir):
            shutil.rmtree(xdir, ignore_errors=True)


def remove_files(fnames):
    for x in fnames:
        if x and os.path.exists(x):
            if os.path.isfile(x):
                os.remove(x)
            elif os.path.isdir(x):
                shutil.rmtree(x, ignore_errors=True)


def _flatten_plus_safe(rollback_files):
    """
    Flatten names of files and create temporary file names.
    """
    tx_files, orig_files = [], []
    for fnames in rollback_files:
        if isinstance(fnames, basestring):
            fnames = [fnames]
        for fname in fnames:
            basedir = safe_makedir(os.path.dirname(fname))
            tmpdir = safe_makedir(tempfile.mkdtemp(dir=basedir))
            tx_file = os.path.join(tmpdir, os.path.basename(fname))
            tx_files.append(tx_file)
            orig_files.append(fname)
    return tx_files, orig_files


def safe_makedir(dname):
    """
    Make a directory if it doesn't exist, handling concurrent race conditions.
    """
    if not dname:
        return dname
    num_tries = 0
    max_tries = 5
    while not os.path.exists(dname):
        try:
            os.makedirs(dname)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dname


def file_exists(fnames):
    """
    Check if a file or files exist and are non-empty.

    parameters
        fnames : file path as string or paths as list; if list, all must exist

    returns
        boolean
    """
    if isinstance(fnames, basestring):
        fnames = [fnames]
    for f in fnames:
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            return False
    return True


def check_dependencies(executables):
    for exe in executables:
        if not find_executable(exe):
            sys.exit("`%s` not found in PATH." % exe)


def name_from_path(path):
    file, ext = os.path.splitext(os.path.basename(path))
    if ext == ".gz":
        file, ext = os.path.splitext(file)
    return file


def prodigal(fasta, out_files, verbose=False):
    """Expected order of 4 out_files is proteins, genes, genbank, and score."""
    if file_exists(out_files):
        return out_files

    if verbose:
        print("Running Prodigal on %s" % fasta, file=sys.stderr)

    with file_transaction(out_files) as tx_out_files:
        cmd = ("prodigal -a {proteins} -d {genes} "
               "-i {fasta} -o {genbank} -p meta -s {score}"
              ).format(proteins=tx_out_files[0],
                       genes=tx_out_files[1],
                       fasta=fasta,
                       genbank=tx_out_files[2],
                       score=tx_out_files[3])
        subprocess.check_call(cmd, shell=True)
    return out_files


def create_fasta_index(fasta, verbose=False):
    out_file = fasta + ".fai"
    if file_exists(out_file):
        return out_file

    if verbose:
        print("Creating an index of", fasta, file=sys.stderr)

    cmd = "samtools faidx {fasta}".format(fasta=fasta)
    subprocess.check_call(cmd, shell=True)
    return out_file


def trnascan(fasta, out_file, verbose=False):
    if file_exists(out_file):
        return out_file

    if verbose:
        print("Running tRNAscan on %s" % fasta, file=sys.stderr)

    with file_transaction(out_file) as tx_out_file:
        # TODO parameterize
        cmd = "tRNAscan-SE -Q -B -o {out} {fasta}".format(out=tx_out_file,
                                                          fasta=fasta)
        subprocess.check_call(cmd, shell=True)
    return out_file


def readfa(fh):
    for header, group in itertools.groupby(fh, lambda line: line[0] == '>'):
        if header:
            line = group.next()
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq


def gc_skew_and_content(seq, window_size=500):
    """sliding window implementation of gc_skew and gc_content.

    >>> seq = 'TTATCCTATTCAGCCACCCGATTTGGAAACCCGGATTGCAATTCTCCAAA'
    >>> window_size = 40
    >>> vals = [(mid, skew, content) for mid, skew, content in gc_skew_and_content(seq, window_size)]
    >>> vals[0][0]
    20
    >>> vals[0][1]
    -0.263...
    >>> vals[0][2]
    0.47...
    """
    half_window = window_size / 2
    seq = seq.upper()
    seq_len = len(seq)
    assert seq_len >= window_size

    for start in xrange(0, seq_len - window_size + 1):
        stop = start + window_size
        mid = stop - half_window
        s = seq[start:stop]

        g = s.count('G')
        c = s.count('C')
        gc = g + c

        content = float(gc) / window_size
        skew = 0 if gc == 0 else (g - c) / float(g + c)
        yield mid, skew, content


def gc_content(fasta, out_file, verbose=False):
    if file_exists(out_file):
        return out_file

    if verbose:
        print("Calculating GC content of %s" % os.path.basename(fasta))

    header = ["SEQUENCE_NAME", "POSITION", "SKEW", "CONTENT"]
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file.rpartition(".")[0],
                  'w') as wfh, open(fasta) as rfh:
            print(*header, sep="\t", file=wfh)
            for name, seq in readfa(rfh):
                for point, skew, content in gc_skew_and_content(seq, 500):
                    print("%s\t%i\t%0.3f\t%0.3f" %
                          (name, point, skew, content),
                          file=wfh)
        cmd = "gzip {tsv}".format(tsv=tx_out_file.rpartition(".")[0])
        subprocess.check_call(cmd, shell=True)
    return out_file


def blastp(fasta, out_file, db,
           num_alignments=10,
           evalue=0.001,
           threads=1,
           verbose=False):
    if file_exists(out_file):
        return out_file

    if verbose:
        print("Running blastp on %s using:" % fasta, "    -db %s" % db,
              "    -num_alignments %d" % num_alignments, "    -evalue %f" %
              evalue, "    -num_threads %d" % threads,
              sep="\n",
              file=sys.stderr)

    with file_transaction(out_file) as tx_out_file:
        nongz = tx_out_file.rpartition(".")[0]

        fields = ["qseqid", "sseqid", "pident", "length", "mismatch",
                  "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                  "bitscore", "sallseqid", "score", "nident", "positive",
                  "gaps", "ppos", "qframe", "sframe", "qseq", "sseq", "qlen",
                  "slen", "salltitles"]

        with open(nongz, 'w') as fo:
            print(*fields, sep="\t", file=fo)

        cmd = ("blastp -db {db} -query {query} -outfmt "
               "'6 {fields}' "
               "-num_threads {threads} "
               "-num_alignments {alignments} "
               "-evalue {evalue} >> {out}").format(db=db,
                                                   query=fasta,
                                                   fields=" ".join(fields),
                                                   threads=threads,
                                                   alignments=num_alignments,
                                                   evalue=evalue,
                                                   out=nongz)
        subprocess.check_call(cmd, shell=True)
        subprocess.check_call("gzip {tsv}".format(tsv=nongz), shell=True)
    return out_file


def make_diamond_db(fasta, db, threads=1, verbose=False):
    out_file = db + ".dmnd"
    if file_exists(out_file):
        return db

    if verbose:
        print("Creating DIAMOND database for", fasta, file=sys.stderr)

    with file_transaction(out_file) as tx_out_file:
        cmd = ("diamond makedb --in {fasta} -d {db} "
               "-p {threads}").format(fasta=fasta,
                                      db=tx_out_file,
                                      threads=threads)
        subprocess.check_call(cmd, shell=True)
    return db


def diamond_blastx(fasta, out_file, db, threads=1, verbose=False):
    if file_exists(out_file):
        return out_file

    if verbose:
        print("Running DIAMOND BLASTX on %s across %s" %
              (os.path.basename(fasta), os.path.basename(db)),
              file=sys.stderr)

    with file_transaction(out_file) as tx_out_file:
        cmd = ("diamond blastx -d {db} -q {fasta} "
               "-a {out} -p {threads} -t {tmpdir}").format(db=db,
                                                           fasta=fasta,
                                                           out=tx_out_file,
                                                           threads=threads,
                                                           tmpdir=tempfile.gettempdir())
        subprocess.check_call(cmd, shell=True)
    return out_file


def diamond_view(daa, out_file, verbose=False):
    if file_exists(out_file):
        return out_file

    if verbose:
        print("Converting DIAMOND database %s to tabular (%s)" %
              (os.path.basename(daa), os.path.basename(out_file)),
              file=sys.stderr)

    with file_transaction(out_file) as tx_out_file:
        nongz = tx_out_file.rpartition(".")[0]
        cmd = ("diamond view -a {daa} -o {out}").format(daa=daa, out=nongz)
        subprocess.check_call(cmd, shell=True)
        cmd = "gzip {tsv}".format(tsv=nongz)
        subprocess.check_call(cmd, shell=True)
    return out_file


def alignment_coverage(daa, tsv, fasta_index, out_file,
                       identity=50.0,
                       verbose=False):
    """
    daa is diamond alignment archive
    tsv is gzipped diamond tab view file
    out_file is a tsv of chrome, 1-based position, count
    """

    def _hit_set_from_tsv(tsv, identity=50.0):
        hits = set()
        with gzip.open(tsv) as fh:
            for l in fh:
                l = l.strip().split("\t")
                if float(l[2]) > identity:
                    hits.add("%s:%s" % (l[0], l[1]))
        return hits

    def _daa_to_sam(daa, sam):
        cmd = "diamond view -f sam -a {daa} -o {sam}".format(daa=daa, sam=sam)
        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError:
            print("Conversion from DAA to SAM has failed.", file=sys.stderr)
            raise
        return sam

    def _filter_sam(sam, hits, out_file):
        with open(out_file, 'w') as filtered_sam, open(sam) as unfiltered_sam:
            for l in unfiltered_sam:
                l = l.strip()
                if l.startswith("@"):
                    print(l, file=filtered_sam)
                    continue
                l = l.split("\t")
                if "%s:%s" % (l[0], l[2]) in hits:
                    print(*l, sep="\t", file=filtered_sam)
        return out_file

    def _sam_to_bam(sam, idx, out_file):
        cmd = ("samtools view -Sbht {index} {sam} 2> /dev/null "
               "| samtools sort -m 8G - {bam} 2> /dev/null"
              ).format(index=idx,
                       sam=sam,
                       bam=out_file.rpartition(".")[0])
        subprocess.check_call(cmd, shell=True)
        return out_file

    if file_exists(out_file):
        return out_file

    name = name_from_path(daa)
    if verbose:
        print("Finding coverages per contig based coverage on hits above",
              identity, "percent in", name,
              file=sys.stderr)

    with file_transaction(out_file) as tx_out_file:
        tmpdir = os.path.dirname(tx_out_file)
        filtered_sam = os.path.join(tmpdir, "filtered.sam")
        # create dictionary of passing hits
        passing_hits = _hit_set_from_tsv(tsv, identity)
        if not passing_hits:
            print("No hits passed filter for", name, file=sys.stderr)
            return

        # create a sam from daa
        unfiltered_sam = _daa_to_sam(daa, os.path.join(tmpdir,
                                                       "unfiltered.sam"))
        # create new filtered sam
        filtered_sam = _filter_sam(unfiltered_sam, passing_hits,
                                   os.path.join(tmpdir, "filtered.sam"))
        # convert sam to bam
        bam = _sam_to_bam(filtered_sam, fasta_index,
                          os.path.join(tmpdir, "filtered.bam"))
        cmd = "bedtools genomecov -d -ibam {bam} | gzip > {out}".format(
            bam=bam,
            out=tx_out_file)
        subprocess.check_call(cmd, shell=True)
    return out_file


def read_count(fasta):
    """Count the number of entries within a fasta file."""
    total = 0
    if fasta.endswith(".gz"):
        count_file = fasta.rsplit(".gz", 1)[0] + ".count"
        cat = "gunzip -c"
    else:
        count_file = fasta + ".count"
        cat = "cat"
    if not os.path.exists(count_file):
        cmd = '{cat} {fasta} | grep -c "^>" > {out}'.format(cat=cat,
                                                            fasta=fasta,
                                                            out=count_file)
        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError:
            # TODO: have this be a warning and grab the count without
            # having to write the file
            print(fasta, "needs to be in a writable folder!")
            raise
    for count in open(count_file):
        total = int(count.strip())
        break
    return total


def tetramer_pca(fasta, out_file, script_path,
                 window_size=1600,
                 step_size=200,
                 verbose=False):
    if file_exists(out_file):
        return out_file

    if verbose:
        print("Running Tetramer PCA on", os.path.basename(fasta),
              file=sys.stderr)

    # using this to create the directory structure only
    with file_transaction(out_file) as tx_out_file:
        cmd = ("Rscript --vanilla {script} --input {fasta} "
               "--outputPC {out} --window {window} "
               "--step {step}").format(script=script_path,
                                       fasta=fasta,
                                       out=out_file,
                                       window=window_size,
                                       step=step_size)
        if verbose:
            print("$> ", cmd, file=sys.stderr)
        subprocess.check_call(cmd, shell=True)
    return out_file


def write_signals_conf(cfg_file, output, name, input_fasta, gc_output,
                       protein_fasta, blastp_tsv, trna_output, pca_output,
                       query_results, training_file, knn_k, verbose=False):
    if verbose:
        print("Building viral signals plot configuration file", cfg_file,
              file=sys.stderr)
    pallette = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                "#ffff33", "#a65628", "#f781bf", "#999999"]
    num_cols = len(pallette)
    names = sorted(query_results.keys())
    # blast tabular output
    comps = []
    for i, (j, k) in enumerate(itertools.combinations(names, 2), start=1):
        if i <= 2:
            comps.append("[Similarity_%d]" % i)
            comps.append("name1=%s" % j)
            comps.append("name2=%s" % k)
            comps.append("dir1=%s" % os.path.dirname(query_results[j]['blast']))
            comps.append("dir2=%s" % os.path.dirname(query_results[k]['blast']))
            comps.append("flag1=%s" %
                         os.path.basename(query_results[j]['blast']))
            comps.append("flag2=%s" %
                         os.path.basename(query_results[k]['blast']))
            comps.append("log=false")
            comps.append('color1=%s' % pallette[names.index(j) % num_cols])
            comps.append('color2=%s' % pallette[names.index(k) % num_cols])
            comps.append("reads1=%d" % read_count(query_results[j]['fasta']))
            comps.append("reads2=%d" % read_count(query_results[k]['fasta']))
        else:
            comps.append("#[Similarity_%d]" % i)
            comps.append("#name1=%s" % j)
            comps.append("#name2=%s" % k)
            comps.append("#dir1=%s" %
                         os.path.dirname(query_results[j]['blast']))
            comps.append("#dir2=%s" %
                         os.path.dirname(query_results[k]['blast']))
            comps.append("#flag1=%s" %
                         os.path.basename(query_results[j]['blast']))
            comps.append("#flag2=%s" %
                         os.path.basename(query_results[k]['blast']))
            comps.append("#log=false")
            comps.append('#color1=%s' % pallette[names.index(j) % num_cols])
            comps.append('#color2=%s' % pallette[names.index(k) % num_cols])
            comps.append("#reads1=%d" % read_count(query_results[j]['fasta']))
            comps.append("#reads2=%d" % read_count(query_results[k]['fasta']))
    # genomic coverage output
    pileups = []
    for i, j in enumerate(names, start=1):
        if i <= 5:
            pileups.append("[Pileup_%d]" % i)
            pileups.append("name=%s" % j)
            pileups.append("dir=%s" %
                           os.path.dirname(query_results[j]['coverage']))
            pileups.append("flag=%s" %
                           os.path.basename(query_results[j]['coverage']))
            pileups.append("log=false")
            pileups.append('color=%s' % pallette[names.index(j) % num_cols])
        else:
            pileups.append("#[Pileup_%d]" % i)
            pileups.append("#name=%s" % j)
            pileups.append("#dir=%s" %
                           os.path.dirname(query_results[j]['coverage']))
            pileups.append("#flag=%s" %
                           os.path.basename(query_results[j]['coverage']))
            pileups.append("#log=false")
            pileups.append('#color=%s' % pallette[names.index(j) % num_cols])

    conf = """[Weights]
VIRAL_WEIGHTS=0:-1e6,2:0.2,3:0.3,5:0.5
POV_FR_WEIGHTS=0:-1e6,1:0.01,3:0.03,5:0.05
FR_RATIO_WEIGHTS=0:-1e6,2:2.25,3:5,5:25
[keywords]
VIRAL_KEYWORDS=phage,virus,prophage,terminase,t4-like,lambda-like,mu-like,capsid,tail,fiber,lambdoid,portal,tail,virion,lysis,podovirus,podo-like,head,baseplate,myovirus,siphovirus,structural
VIRAL2_KEYWORDS=integrase,transposase
[setup]
input_path={output}
name={name}
output_path={output}/summary
[inputs]
fasta_file={fasta}
gc_content_file={gc_output}
proteins_file={p_proteins}
blastp_file={blastp_tsv}
trna_file={trna_output}
tetramer_file={pca_output}
[classifier]
training_file={training_file}
knn_k={knn_k}
[genes]
name=genes
tRNA=#377eb8
hypothetical=black
viral=#e41a1c
viral2=#e41a1c
[Coverage]
log=false
color=black
[GC]
percent=#377eb8
skew=#e41a1c
log=false
[Tetramer PC]
PC2=#377eb8
PC1=#e41a1c
log=false
{comparisons}
{pileups}
""".format(output=output,
           name=name,
           fasta=input_fasta,
           gc_output=gc_output,
           p_proteins=protein_fasta,
           blastp_tsv=blastp_tsv,
           trna_output=trna_output,
           pca_output=pca_output,
           training_file=training_file or os.path.join(output, "YOUR_TRAINING_FILE.csv"),
           knn_k=knn_k,
           comparisons="\n".join(comps),
           pileups="\n".join(pileups))
    with open(cfg_file, 'w') as ofh:
        print(conf, file=ofh)
    return cfg_file


def viralscan(fasta, output, query, name, threads, identity, verbose, db,
              num_alignments, evalue, script_path, window_size, step_size,
              training_file, knn_k):
    check_dependencies(REQUIRES)
    if name is None:
        name = name_from_path(fasta)

    # prodigal
    p_proteins, p_genes, p_genbank, p_score = prodigal(
        fasta, [os.path.join(output, "prodigal", name + "_proteins.fasta"),
                os.path.join(output, "prodigal", name + "_genes.fasta"),
                os.path.join(output, "prodigal", name + ".gbk"),
                os.path.join(output, "prodigal", name + ".scores")], verbose)
    p_proteins_index = create_fasta_index(p_proteins)
    # tRNAscan-SE
    trna_output = trnascan(fasta, os.path.join(output, "tRNAscan",
                                               name + "-tRNAscan.txt"),
                           verbose)
    # gc content and skew
    gc_output = gc_content(fasta,
                           os.path.join(output, name + "_gc_content.tsv.gz"),
                           verbose)
    # blastp
    blastp_tsv = blastp(p_proteins,
                        os.path.join(output, name + "_blastp.tsv.gz"), db,
                        num_alignments, evalue, threads, verbose)
    # create database
    protein_db = make_diamond_db(p_proteins,
                                 os.path.join(output, "diamond",
                                              os.path.basename(p_proteins).partition(".")[0]),
                                 threads, verbose)
    # diamond
    query_results = {}
    for q in query:
        qname = name_from_path(q)
        diamond_daa = diamond_blastx(q, os.path.join(output, "diamond",
                                                     qname + ".daa"),
                                     protein_db, threads, verbose)
        diamond_tsv = diamond_view(diamond_daa,
                                   os.path.join(output, "diamond",
                                                qname + ".tsv.gz"), verbose)
        coverage_tsv = alignment_coverage(diamond_daa, diamond_tsv,
                                          p_proteins_index,
                                          os.path.join(output, "coverage",
                                                       qname + ".tsv.gz"),
                                          identity, verbose)
        if coverage_tsv:
            query_results[qname] = dict(coverage=coverage_tsv,
                                        fasta=q,
                                        blast=diamond_tsv)
    # tetramerPCA
    if script_path:
        pca_output = tetramer_pca(fasta,
                                  os.path.join(output, "tetramerPCA",
                                               name + "-tetramer-PC.csv"),
                                  script_path, window_size, step_size, verbose)
    else:
        pca_output = os.path.join(output, "tetramerPCA",
                                  "YOUR_SAMPLE-tetramer-PC.csv")
    # write viralsignals config file
    signals_cfg = write_signals_conf(os.path.join(output, "signals.cfg"),
                                     output, name, fasta, gc_output,
                                     p_proteins, blastp_tsv, trna_output,
                                     pca_output, query_results,
                                     training_file, knn_k)


def main():
    def _file_exists(parser, arg):
        if isinstance(arg, basestring):
            tmparg = [arg]
        for a in tmparg:
            if not os.path.exists(a):
                parser.error("The file %s does not exist." % a)
            if os.path.isdir(a):
                parser.error("Expecting file, but got folder %s." % a)
        return arg

    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        version=__version__)
    p.add_argument('fasta',
                   type=lambda x: _file_exists(p, x),
                   help="Fasta file of sequences to analyze.")
    p.add_argument('output', help="Location to write output files.")
    p.add_argument('query',
                   type=lambda x: _file_exists(p, x),
                   nargs="+",
                   help=("Reference fasta(s). Specify as many as "
                         "you would like or use a file pattern."))

    p.add_argument('-n', '--name', help="Name of sample being processed.")
    p.add_argument('-t', '--threads',
                   default=12,
                   type=int,
                   help="Number of threads to use.")
    p.add_argument('-i', '--identity',
                   default=50,
                   type=float,
                   help="passing blastx percent identity per hit")
    p.add_argument('--verbose', action='store_true')

    blasto = p.add_argument_group('BLAST options')
    blasto.add_argument('--db', default='nr', help="BLAST database")
    blasto.add_argument('--num-alignments',
                        type=int,
                        default=10,
                        help=("Number of database sequences to "
                              "show alignments for."))
    blasto.add_argument('--evalue',
                        type=float,
                        default=0.001,
                        help="Expectation value threshold for saving hits.")

    tetramero = p.add_argument_group('TetramerPCA options')
    tetramero.add_argument('--script-path',
                           type=lambda x: _file_exists(p, x),
                           help=("tetramer PCA R script path; testing "
                                 "ignore if this option is skipped"))
    tetramero.add_argument('--window-size', default=1600, type=int)
    tetramero.add_argument('--step-size', default=200, type=int)

    classifiero = p.add_argument_group('Classifier options')
    # having these options here implies viralscan runs the classifier
    # rather than simply filling the config file
    # TODO: write config file separately or include classifier as part of viral-scan
    classifiero.add.argument('--training-file',
                             type=lambda x: _file_exists(p, x),
                             help="classifier training file")
    classifiero.addargument("--knn-k",
                            type=int,
                            default=3,
                            help="nearest neighbors for classifier")

    args = vars(p.parse_args())
    viralscan(**args)


if __name__ == '__main__':
    main()
