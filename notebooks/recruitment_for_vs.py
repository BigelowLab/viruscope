from __future__ import print_function

import argparse
import contextlib
import fileinput
import gzip
import itertools
import os
import os.path as op
import shutil
import subprocess
import six
import sys
import tempfile
import time
from collections import defaultdict
from distutils.spawn import find_executable
import click
import pandas as pd
from scgc import readfx



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


def read_count(fname):
    """Count the number of reads and write metadata .count file.

    Args:
        fname (str): fastq or fasta file path

    Returns:
        read_count (int): integer of number of reads within fasta/fastq file
    """
    total_reads = 0
    
    if op.exists(fname) == False:
        logger.error("could not find file: %s" % fname)
        return 0
    count_file = '{}.count'.format(fname)

    if op.exists(count_file):
        total_reads = open(count_file).read().split("\n")[0]
        return total_reads

    for name, seq, qual in readfx(fname):
        total_reads += 1
    try:
        with open(count_file, "w") as oh:
            print(total_reads, file=oh)
        return total_reads

    except:
        return total_reads


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


def diamond_view(daa, out_file, threads, verbose=False):
    ''' coverts diamond result to a tabular output '''
    if file_exists(out_file):
        return out_file

    if verbose:
        print("Converting DIAMOND database %s to tabular (%s)" %
              (os.path.basename(daa), os.path.basename(out_file)),
              file=sys.stderr)

    with file_transaction(out_file) as tx_out_file:
        nongz = tx_out_file.rpartition(".")[0]
        subprocess.check_call(["diamond", "view","-a", daa, "-o", nongz])
        subprocess.check_call(["gzip", nongz])
    return out_file

## stats table construction
def import_diamond_tsv(tsv, pctid=50.0, best_hit=True):
    df = pd.read_csv(tsv, 
                     names="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(), 
                     sep="\t")
    df = df[df['pident'] >= pctid]
    if best_hit:
        df = df.sort_values(by=['qseqid', 'length','bitscore'], ascending=False).drop_duplicates(subset='qseqid', keep='first')
    return df


def summarize_by_contig(df, hitscol):
    return pd.Series(df.groupby('contig')[hitscol].sum(), name=hitscol)


def contig_lengths(infile):
    ''' create a dict with contig names as keys and lengths as values from gff file'''
    outdict = {}
    if "g" in infile.split(".")[-1]:
        filetype = 'gff'
        print("looks like input contig file is in gff format.")
    elif "f" in infile.split(".")[-1]:
        filetype = 'fasta'
        print("looks like input config fiel is in fasta format.")
    else:
        raise IOError("can't figure out what kind of file contig file is.  Make sure it's either in fasta or gff format.")
    if filetype == 'gff':
        with open(infile) as ih:
            for l in ih:
                if l.startswith("##sequence-region"):
                    vec = l.strip().split()
                    outdict[vec[1]] = vec[-1]
    
    elif filetype == 'fasta':
        for name, seq, qual in readfx(infile):
            outdict[name] = len(seq)  
    return outdict


def compute_fr(tbl, clens, mult=1e6):
    '''
    Args:
        tbl: output stats table with mg hit and read counts from diamond recruitment
        clens: dict of contig lengths
        mult: factor to multiply fraction by to make readiable
        
    Outputs:
        pandas DataFrame with mg_fr values calculated
    '''

    tbl['contig_length'] = [float(clens[i]) for i in tbl['contig']]
    
    hits_cols = [i for i in tbl.columns if 'hit' in i]
    count_cols = ["_".join(["reads",i.split("_")[1]]) for i in hits_cols]
    
    for h, c in zip(hits_cols, count_cols):
        fr = tbl[h]/(tbl[c] * tbl['contig_length']) * mult
        tbl[h.replace("hit_","fr_")] = fr
    return tbl


def orf_map(gff):
    gdf = pd.read_csv(gff, comment='#', sep="\t", names = ['contig','app','type','start','stop','dot','strand','val','notes']).dropna()
    gdf['orf'] = [i.split(";")[0].replace("ID=",'') for i in gdf['notes']]
    gdf['len'] = gdf['stop'] - gdf['start']
    return gdf[['contig','orf']]


def map_orfs_to_contigs(df, contig_file):
    if "g" in contig_file.split(".")[-1]:
        gff = True
        print("looks like input contig file is in gff format.  Will map ORFs to contigs using that.")
    else:
        print("doesn't look like input contig file is in gff format.  Will assume that contig name is embedded in the ORF name.")
        gff = False
    
    if gff:
        gdf = orf_map(contig_file)
        return pd.merge(df, gdf, on='orf', how='outer')
    else:
        df['contig'] = ["_".join(i.split("_")[:-1]) for i in df['orf']]
        return df
    

def construct_recruit_tbl(vir_tsv, bac_tsv, read_count_dict, contig_file):
    '''
    Args:
        vir_tsv: diamond recruitment converted to tsv for vir metagenome
        bac_tsv: diamond recruitment converted to tsv for bac metagenome
        read_count_dict: dict of mg read counts with two keys -- 'vir_reads' and 'bac_reads'
        contig_file: path to a file with sag contigs in it; either in fasta or gff format
    Returns:
        pandas dataframe with mg fraction calculated
    '''
    def _count_hits(df, name):
        outser = pd.Series(df.groupby('sseqid')['qseqid'].count(), name=name)
        if len(outser) == 0:
            outdf = pd.DataFrame(columns=['orf',name])
        else:
            outdf = pd.DataFrame(outser).reset_index().rename(columns={'sseqid':'orf'})
        return outdf
    
    cnames = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split()
    bac_df = import_diamond_tsv(bac_tsv)
    vir_df = import_diamond_tsv(vir_tsv)

    bac_sum = _count_hits(bac_df, 'hit_mg-bac')
    vir_sum = _count_hits(vir_df, 'hit_mg-vir')
    
    orfhits = pd.merge(bac_sum, vir_sum, on='orf', how='outer')[['orf','hit_mg-bac','hit_mg-vir']]
    orfhits = map_orfs_to_contigs(orfhits, contig_file)

    chits = pd.concat([summarize_by_contig(orfhits, 'hit_mg-bac'), summarize_by_contig(orfhits, 'hit_mg-vir')], axis=1)
    chits['reads_mg-vir'] = float(read_count_dict['vir_reads'])
    chits['reads_mg-bac'] = float(read_count_dict['bac_reads'])
    
    clens = contig_lengths(contig_file)
    
    out_tbl = compute_fr(chits.reset_index(), clens, mult=1e6)
    out_tbl['ratio_virus_bacteria'] = out_tbl['fr_mg-vir'] / out_tbl['fr_mg-bac']
    out_tbl['ratio_virus_bacteria'] = [1000 if i == float('inf') else i for i in out_tbl['ratio_virus_bacteria']]
    
    return out_tbl



@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('prot_fasta', nargs=1)
@click.argument('vir_mg', nargs=1)
@click.argument('bac_mg', nargs=1)
@click.option('--sag-contigs',
             help='location of sag contigs in either fasta or gff format, no output created table created if default of None used',
             default=None,
             show_default=True)
@click.option('--output',
              default=None,
              help='directory location to place output files',
             show_default=True)
@click.option('--threads',
              default=8,
              show_default=True,
              help='number of cores to run on')
@click.option('--verbose',
              default=True,
             show_default=True)
def main(prot_fasta, vir_mg, bac_mg, sag_contigs, output, threads, verbose):
    '''
python recruitment_for_vs.py --threads 10 --output /mnt/scgc/simon/simonsproject/bats248_vs/diamond/pergenome/ --sag-contigs /mnt/scgc/simon/simonsproject/bats248_contigs/coassemblies/AG-920/AG-920-P22_contigs.fasta  /mnt/scgc/simon/simonsproject/bats248_vs/prodigal/AG-920-P22_proteins.fasta /mnt/scgc_nfs/ref/viral_dbs/POV.fasta.gz /mnt/scgc_nfs/ref/viral_dbs/LineP-all.fasta.gz
    '''
    fa_name = op.basename(prot_fasta).split(".")[0]
    
    protein_db = make_diamond_db(prot_fasta,
                                 os.path.join(output,
                                              os.path.basename(prot_fasta).partition(".")[0]),
                                 threads, verbose)
    
    tsv_list = []
    if op.exists(op.join(output, "{}_mg_diamond_recruitment_tbl.csv".format(fa_name))):
        print("output table already exists for {}!".format(fa_name))
        return
                 
    
    
    for q in [vir_mg, bac_mg]:
        qname = op.basename(q).split(".")[0]
        
        diamond_daa = diamond_blastx(q, os.path.join(output, '{fa}_vs_{qname}.daa'.format(fa=fa_name, qname=qname)),
                                     protein_db, threads, verbose)
        diamond_tsv = diamond_view(diamond_daa,
                                   os.path.join(output, '{fa}_vs_{qname}.tsv.gz'.format(fa=fa_name, qname=qname)), 
                                   verbose)
        tsv_list.append(diamond_tsv)
    print(",".join(tsv_list))
    if sag_contigs is not None:

        # dict of metagenome read counts
        read_count_dict = {}
        read_count_dict['vir_reads'] = read_count(vir_mg)
        read_count_dict['bac_reads'] = read_count(bac_mg)
        out_tbl = construct_recruit_tbl(tsv_list[0], tsv_list[1], read_count_dict, sag_contigs)
        out_tbl.to_csv(op.join(output, "{}_mg_diamond_recruitment_tbl.csv".format(fa_name)), sep=",", index=False)

    os.remove('{}.dmnd'.format(protein_db))

        
if __name__ =='__main__':
    main()