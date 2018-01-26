from __future__ import print_function

import argparse
import contextlib
import fileinput
import gzip
import itertools
import os
import shutil
import subprocess
import six
import sys
import tempfile
import time
from collections import defaultdict
from distutils.spawn import find_executable
import pandas as pd

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
        if isinstance(fnames, six.string_types):
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
    if isinstance(fnames, six.string_types):
        fnames = [fnames]
    for f in fnames:
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            return False
    return True


def check_dependencies(executables):
    exes = []
    for exe in executables:
        if not find_executable(exe):
            exes.append(exe)
    if len(exes) > 0:
        for exe in exes:
            print("`%s` not found in PATH." % exe)
        sys.exit(1)


def name_from_path(path):
    file, ext = os.path.splitext(os.path.basename(path))
    if ext == ".gz":
        file, ext = os.path.splitext(file)
    return file

def readfa(fh):
    for header, group in itertools.groupby(fh, lambda line: line[0] == '>'):
        if header:
            line = next(group)
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq
            
def format_fasta_record(name, seq, wrap=100):
    """Fasta __str__ method.

    Convert fasta name and sequence into wrapped fasta format.

    Args:
        name (str): name of the record
        seq (str): sequence of the record
        wrap (int): length of sequence per line

    Returns:
        tuple: name, sequence

    >>> format_fasta_record("seq1", "ACTG")
    ">seq1\nACTG"
    """
    record = ">{name}\n".format(name=name)
    if wrap:
        for i in range(0, len(seq), wrap):
            record += seq[i:i+wrap] + "\n"
    else:
        record += seq + "\n"
    return record.strip()

def run_cd_hit(input_fa, output_fa, c=0.9, G=1, b=20, M=800,
    T=1, n=5, l=10, t=2, d=20, s=0.0, S=999999, aL=0.0, AL=99999999, aS=0.0,
    AS=99999999, A=0, uL=1.0, uS=1.0, U=99999999, g=1, sc=0, sf=0):
    """Run CD-HIT to cluster input FASTA.

    Args:
        input_fa (str): File path to fasta.
        output_fa (str): File path of output fasta.
        c (Optional[float]): sequence identity threshold, default 0.9
 	        this is the default cd-hit's "global sequence identity" calculated as:
 	        number of identical amino acids in alignment
            divided by the full length of the shorter sequence
        G (Optional[int]): use global sequence identity, default 1
 	        if set to 0, then use local sequence identity, calculated as :
            number of identical amino acids in alignment
 	        divided by the length of the alignment
 	        NOTE!!! don't use -G 0 unless you use alignment coverage controls
 	        see options -aL, -AL, -aS, -AS
        b (Optional[int]): band_width of alignment, default 20
        M (Optional[int]): memory limit (in MB) for the program, default 800; 0 for unlimited
        T (Optional[int]): number of threads, default 1; with 0, all CPUs will be used
        n (Optional[int]): word_length, default 5, see user's guide for choosing it
        l (Optional[int]): length of throw_away_sequences, default 10
        t (Optional[int]): tolerance for redundance, default 2
        d (Optional[int]): length of description in .clstr file, default 20
 	        if set to 0, it takes the fasta defline and stops at first space
        s (Optional[float]): length difference cutoff, default 0.0
 	        if set to 0.9, the shorter sequences need to be
            at least 90% length of the representative of the cluster
        S (Optional[int]): length difference cutoff in amino acid, default 999999
 	        if set to 60, the length difference between the shorter sequences
 	        and the representative of the cluster can not be bigger than 60
        aL (Optional[float]): alignment coverage for the longer sequence, default 0.0
 	        if set to 0.9, the alignment must covers 90% of the sequence
        AL (Optional[int]): alignment coverage control for the longer sequence, default 99999999
 	        if set to 60, and the length of the sequence is 400,
 	        then the alignment must be >= 340 (400-60) residues
        aS (Optional[float]): alignment coverage for the shorter sequence, default 0.0
        	if set to 0.9, the alignment must covers 90% of the sequence
        AS (Optional[int]): alignment coverage control for the shorter sequence, default 99999999
        	if set to 60, and the length of the sequence is 400,
        	then the alignment must be >= 340 (400-60) residues
        A (Optional[int]): minimal alignment coverage control for the both sequences, default 0
        	alignment must cover >= this value for both sequences
        uL (Optional[float]): maximum unmatched percentage for the longer sequence, default 1.0
        	if set to 0.1, the unmatched region (excluding leading and tailing gaps)
        	must not be more than 10% of the sequence
        uS (Optional[float]): maximum unmatched percentage for the shorter sequence, default 1.0
        	if set to 0.1, the unmatched region (excluding leading and tailing gaps)
        	must not be more than 10% of the sequence
        U (Optional[int]): maximum unmatched length, default 99999999
        	if set to 10, the unmatched region (excluding leading and tailing gaps)
        	must not be more than 10 bases
        g (Optional[int]): 1 or 0, default 1
        	when 0 a sequence is clustered to the first
        	cluster that meet the threshold (fast cluster). If set to 1, the program
        	will cluster it into the most similar cluster that meet the threshold
        	(accurate but slow mode)
        	but either 1 or 0 won't change the representatives of final clusters
        sc (Optional[int]): sort clusters by size (number of sequences), default 0, output clusters by decreasing length
        	if set to 1, output clusters by decreasing size
        sf (Optional[int]): sort fasta/fastq by cluster size (number of sequences), default 0, no sorting
        	if set to 1, output sequences by decreasing cluster size

    Returns:
        list, [file path of output fasta, file path of output cluster definitions]

    """
    output_clstr = "{fa}.clstr".format(fa=output_fa)
    output_files = [output_fa, output_clstr]
    if file_exists(output_files):
        return output_files

    print("Running CD-HIT on {fa}".format(fa=input_fa), file=sys.stderr)

    contig_name_map = {}
    tmp_fasta = "{fa}.rename.tmp".format(fa=input_fa)
    with open(input_fa) as f_in, open(tmp_fasta, "w") as f_out:
        for i, (name, seq) in enumerate(readfa(f_in), start=1):
            contig_name_map["%d" % i] = name
            print(format_fasta_record(i, seq), file=f_out)

    with file_transaction(output_files) as tx_out_files:
        cmd = ("cd-hit -i {input_fasta} -o {output_fasta} -c {c} "
                "-G {G} -b {b} -M {M} -T {T} -n {n} -l {l} -t {t} "
                "-d {d} -s {s} -S {S} -aL {aL} -AL {AL} -aS {aS} "
                "-AS {AS} -A {A} -uL {uL} -uS {uS} -U {U} "
                "-p 1 -g {g} -sc {sc} -sf {sf}").format(input_fasta=tmp_fasta,
                                                        output_fasta=tx_out_files[0],
                                                        **locals())
        subprocess.check_call(cmd, shell=True)
        # copy the clstr output to its temp file location; let file_transaction clean up
        shutil.copyfile("{fa}.clstr".format(fa=tx_out_files[0]), tx_out_files[1])

        # edit the output files in place back to their original names
        # changes the format of the cluster file

        # update change the contig names in the cluster file back to original
        for line in fileinput.input(tx_out_files[0], inplace=True):
            line = line.strip()
            if line.startswith(">"):
                name = contig_name_map[line.strip(">")]
                print(">{name}".format(name=name))
            else:
                print(line)

        # change the contig names in the cluster file
        for line in fileinput.input(tx_out_files[1], inplace=True):
            line = line.strip()
            if not line.startswith(">"):
                # changes:
                # 1	382aa, >6... at 1:382:1:382/100.00%
                # to just the original contig name
                if "*" in line:
                    print('{}*'.format(contig_name_map[line.partition(">")[-1].partition("...")[0]]))
                else:
                    print(contig_name_map[line.partition(">")[-1].partition("...")[0]])
            else:
                # this is the cluster ID
                print(line)

    if file_exists(tmp_fasta):
        os.remove(tmp_fasta)

    return output_files


def cluster_map(clstr, singles=False):
    ''' Create dict from cd-hit function output
    
    Args:
        clstr (string): path to cd-hit cluster file
    Returns:
        cluster_map (defaultdict): clutster_map[seed_cluster] = member list
    '''
    
    cluster_map = defaultdict(list)
    with open(clstr) as fh:
        for cluster_start, group in itertools.groupby(fh, lambda l: l[0] == '>'):
            members = []
            if not cluster_start: 
                for line in group:
                    if "*" in line: 
                        rep_seq = line.strip().replace("*", "").split(' ')[0]
                    else:
                        members.append(line.strip().split(" ")[0])
            if singles:
                cluster_map[rep_seq] = members
            elif len(members) > 0:
                cluster_map[rep_seq] = members
            else:
                continue
    return cluster_map

def swap_cluster_map(cm):
    cm_swap = {}

    for c in cm: 
        for k in cm[c]: cm_swap[k] = c
    return cm_swap


def write_fa_record(name, seq, oh, line_len=60):
    print(">{}".format(name), file=oh)
    for i in range(0, len(seq), line_len):
        print(seq[i:i+line_len], file=oh)
        
        
def run_mica(fasta, out, db='/mnt/scgc/simon/databases/mica/nr-20150620-mica', num_alignments=10,
           evalue=0.001,
           threads=20, fields = ["qseqid", "sseqid", "pident", "length", "mismatch",
                  "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                  "bitscore", "sallseqid", "score", "nident", "positive",
                  "gaps", "ppos", "qframe", "sframe", "qseq", "sseq", "qlen",
                  "slen", "salltitles"]):
    cmd = ("mica-search --p='{threads}' --blastp 'blastp' {db} {query} "
                   "--blast-args -outfmt '6 {fields}' "
                   "-num_alignments {alignments} -evalue {evalue} -out {out}").format(db=db,
                                                      query=fasta,
                                                      fields=" ".join(fields),
                                                      threads=threads,
                                                      alignments=num_alignments,
                                                      evalue=evalue,
                                                      out=out)
    print(cmd)
    return cmd

def run_blast(fasta, out, db='nr', num_alignments=10,
           evalue=0.001,
           threads=20, fields = ["qseqid", "sseqid", "pident", "length", "mismatch",
                  "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                  "bitscore", "sallseqid", "score", "nident", "positive",
                  "gaps", "ppos", "qframe", "sframe", "qseq", "sseq", "qlen",
                  "slen", "salltitles"]):
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
                                                      out=out)
    print(cmd)
    return cmd

def prodigal(fasta, out_files, verbose=False):
    """Expected order of 4 out_files is proteins, genes, genbank, and score."""
    if file_exists(out_files):
        return out_files

    if verbose:
        print("Running Prodigal on %s" % fasta, file=sys.stderr)

    with file_transaction(out_files) as tx_out_files:
        cmd = ("prodigal -f gff -a {proteins} -d {genes} "
               "-i {fasta} -o {genbank} -p meta -s {score}"
              ).format(proteins=tx_out_files[0],
                       genes=tx_out_files[1],
                       fasta=fasta,
                       genbank=tx_out_files[2],
                       score=tx_out_files[3])
        subprocess.check_call(cmd, shell=True)
    return out_files        


## from original viruscope

def diamond_view(daa, out_file, verbose=False):
    ''' coverts diamond result to a tabular output '''
    if file_exists(out_file):
        return out_file

    if verbose:
        print("Converting DIAMOND database %s to tabular (%s)" %
              (os.path.basename(daa), os.path.basename(out_file)),
              file=sys.stderr)

    with file_transaction(out_file) as tx_out_file:
        nongz = tx_out_file.rpartition(".")[0]
        subprocess.check_call(["diamond", "view", "-a", daa, "-o", nongz])
        subprocess.check_call(["gzip", nongz])
    return out_file


def create_fasta_index(fasta, verbose=False):
    out_file = fasta + ".fai"
    if file_exists(out_file):
        return out_file

    if verbose:
        print("Creating an index of", fasta, file=sys.stderr)

    cmd = "samtools faidx {fasta}".format(fasta=fasta)
    subprocess.check_call(cmd, shell=True)
    return out_file


def alignment_coverage(daa, fasta_index, out_file,
                       identity=50.0, threads=8,
                       verbose=False):
    """
    daa is diamond alignment archive
    tsv is gzipped diamond tab view file
    out_file is a tsv of chrom, 1-based position, count
    """

    def _daa_to_tsv(daa, tsv):
        cmd = "diamond view -a {daa} -o {tsv}".format(daa=daa, tsv=tsv)
        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError:
            print("Conversion from DAA to tsv has failed.", file=sys.stderr)
            raise
        return tsv
    
    def _hit_set_from_tsv(tsv, identity=50.0):
        hits = set()
        with open(tsv) as fh:
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

    def _sam_to_bam(sam, idx, out_file, threads=8):
        cmd = ("samtools view -@ {t} -bSht {index} {sam} 2> /dev/null "
               "| samtools sort -@ {t} -m 2G - > {bam} 2> /dev/null"
              ).format(t=threads,
                       index=idx,
                       sam=sam,
                       bam=out_file.rpartition(".")[0] + ".bam")
        subprocess.check_call(cmd, shell=True)
        return out_file

    if file_exists(out_file):
        return out_file

    name = name_from_path(daa)
    if verbose:
        print("Finding coverages per contig based coverage on hits above",
              identity, "percent in", name, file=sys.stderr)

    with file_transaction(out_file) as tx_out_file:
        tmpdir = os.path.dirname(tx_out_file)
        tsv = os.path.join(tmpdir, "daa.tsv")
        filtered_sam = os.path.join(tmpdir, "filtered.sam")
        
        # create a tsv from daa
        tsv = _daa_to_tsv(daa, tsv)
        
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
                          os.path.join(tmpdir, "filtered.bam"), threads=threads)
        cmd = "bedtools genomecov -d -ibam {bam} | gzip > {out}".format(
            bam=bam,
            out=tx_out_file)
        subprocess.check_call(cmd, shell=True)
    return out_file


def mean_coverage(daa, fasta_index, out_cov):
    if out_cov.endswith(".gz") == False:
        out_cov = "{}.gz".format(out_cov)
    
    if op.exists(out_cov):
        cov = out_cov
    else:
        cov = alignment_coverage(daa, fasta_index, out_cov)
    
    df = pd.read_csv(cov, sep="\t", names=['orf','position','coverage'])
    ov_per_orf = df.groupby('orf')['coverage'].mean()
    return pd.DataFrame(cov_per_orf).reset_index()


def orf_map(gff):
    gdf = pd.read_csv(gff, comment='#', sep="\t", names = ['contig','app','type','start','stop','dot','strand','val','notes']).dropna()
    gdf['id'] = [i.split(";")[0].replace("ID=",'') for i in gdf['notes']]
    gdf['len'] = gdf['stop'] - gdf['start']
    return gdf[['contig','id','len']]

# adapted from graphsignals by Ben:
def weighted_score(
    x = [0.901, 0.317, 0.653, 0.423, 0.000, 0.419, 0.299, 0.917, 0.195], 
    lut = {'0': -1e6, '2': 0.2, '3': 0.3, '5': 0.5}):
    
    
    """
    Score a set of values according to a weighted look-up-table
    
    Final score is cumulative in the sense that a score is the sum of 
    of a values score and all of the possible lower scores
    
    @param x numeric vector, of values in the range of 0-1
    @param lut dict, look up tables where items are cut-offs between
      weights and keys (when converted to numeric) are the weights
    """
    
    #weights
    kys = sorted(lut.keys())
    w = [float(k) for k in kys]
    
    #values
    v = [lut[k] for k in kys]
    
    # indices
    ix = [bisect_left(v, xi) - 1 for xi in x]
    
    # we find the cumulative weights
    cumweights = np.cumsum(np.asarray(w)).tolist()
    
    # and then assign weighted scores using those
    ws = [cumweights[i] for i in ix]
    
    return ws


def id_virus_orfs(micaout, 
                  keywords1 = 'phage,virus,prophage,terminase,t4-like,lambda-like,mu-like,capsid,tail,fiber,lambdoid,portal,tail,virion,lysis,podovirus,podo-like,head,baseplate,myovirus,siphovirus,structural', 
                  keywords2 = 'integrase,transposase'):
    ''' takes mica output and returns dataframe with 1 if hits had viral signal, 0 if not
    Args:
        micaout (path): path to mica output file
        keywords1 (string): comma separated list of keywords to search for to confirm that orf is viral
        keywords2 (string): comma separated list of secondary keywords to search for to confirm that orf is viral
    Returns:
        pandas DataFrame of ORF classifications, 1 if viral, 0 if not (p1 for keywords1, p2 for kewords2)
    
    '''
    keywords1 = keywords1.split(",")
    keywords2 = keywords2.split(",")
    
    def search_phrases(string, keylist):
        for i, k in enumerate(keylist):
            if k in string:
                return True
            else:
                if i == len(keylist)-1:
                    return False
                else:
                    continue
    cnames = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles".split()
    df = pd.read_csv(micaout, sep="\t", names=cnames)

    p1_true = df[[search_phrases(i, keywords1) for i in df['salltitles']]]['qseqid'].unique()
    p2_true = df[[search_phrases(i, keywords2) for i in df['salltitles']]]['qseqid'].unique()
    qseqs = df['qseqid'].unique()

    newdf = pd.DataFrame(data = {'orf':qseqs,
                                'p1':[1 if q in p1_true else 0 for q in list(qseqs)],
                                'p2':[1 if q in p2_true else 0 for q in list(qseqs)]})
    return newdf

def orf_map_fa(fasta):
    names = []
    contigs = []
    lens = []

    for name, seq in readfa(open(fasta)):
        i = name.split(" ")[0]
        names.append(i)
        contig = "_".join(i.split("_")[:-1])
        contigs.append(contig)
        lens.append(len(seq) * 3)

    return pd.DataFrame(data={'contig':contigs, 'id':names, 'len':lens})
