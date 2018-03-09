import pandas as pd
from collections import defaultdict
import itertools
from .tools import readfa
from .recruit import summarize_by_contig


def map_clstr_raw(clstr, singles=False):
    ''' Map clusters from raw cd-hit cluster output
    Args:
        clstr (str): path to raw cd-hit output
        singles (bool): if true, return singles within output dict with empty list as values
    Returns:
        cluster_map (dict)
    '''
    cluster_map = defaultdict(list)
    with open(clstr) as fh:
        for cluster_start, group in itertools.groupby(fh, lambda l: l[0] == '>'):
            members = []
            rep_seq = ''
            if not cluster_start:
                for line in group:
                    if "*" in line:
                        rep_seq = line.split(",")[1].split("...")[0].replace(">",'').replace(" ","")
                    else:
                        members.append(line.split(",")[1].split("...")[0].replace(">",'').replace(" ",""))
            if len(rep_seq) == 0:
                continue

            if singles:
                cluster_map[rep_seq] = members
            elif len(members) > 0:
                cluster_map[rep_seq] = members
            else:
                continue
    return cluster_map


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


def phage_contig_table(clstr_map, prot_fasta, phage_hits_df, outfile=None):
    ''' create a summary of phage hits to contig orfs, mapping back cd-hit cluster seeds to orfs from individual genomes
    Args:
        clstr_map (dict): cluster map with structure {orf:seed}
        prot_fasta (path): path to protein fasta file for mapping orfs back to contigs
        phage_hits_df (pandas.DataFrame): phage hits dataframe from summary of all Mica results
        outfile (path): where to write output table, if None, none written
    Returns:
        pandas dataframe
    '''
    cm = clstr_map
    om = orf_map_fa(prot_fasta)
    om['lookup'] = [cm.get(i, i) for i in om['id']]
    omp = pd.merge(phage_hits_df, om, left_on='orf', right_on='lookup', how='right').fillna(0)[['contig','id','p1','p2','len']]
    orf_count = pd.DataFrame(om.groupby('contig')['id'].count()).rename(columns={'id':'total_orfs'})
    csum = pd.concat([summarize_by_contig(omp, 'p1'), summarize_by_contig(omp, 'p2'), orf_count], axis=1)
    csum['viral_phage_gene_fraction'] = csum['p1'] / csum['total_orfs']
    csum['viral2_phage_gene_fraction'] = csum['p2'] / csum['total_orfs']

    if outfile is not None:
        csum.to_csv(outfile)

    return csum
