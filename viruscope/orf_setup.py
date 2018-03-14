import subprocess
import sys
import glob
import os
import os.path as op
from collections import defaultdict
import itertools
from pyfaidx import Fasta
from .tools import (file_transaction, file_exists, safe_makedir, write_fa_record,
                    tmp_dir, readfa)


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

def run_prodigal(f, prod_dir):

    name = os.path.basename(f).split(".")[0].split("_")[0]
    outfiles = [os.path.join(prod_dir, name + "_proteins.fasta"),
                os.path.join(prod_dir, name + "_genes.fasta"),
                os.path.join(prod_dir, name + ".gbk"),
                os.path.join(prod_dir, name + ".scores")]
    try:
        p_proteins, p_genes, p_genbank, p_score = prodigal(f, outfiles, verbose=True)
    except Exception as inst:
        print(inst)
        print('Unable to run prodigal for {}'.format(name))
    return p_proteins


def run_batch_prodigal(falist, workingdir = "./"):
    prod_dir = os.path.join(workingdir, 'prodigal')
    safe_makedir(prod_dir)
    out_prots = []
    for f in falist:
        try:
            p_proteins = run_prodigal(f, prod_dir)
        except Exception as inst:
            print(inst)
            continue
        out_prots.append(p_proteins)
    return out_prots


def concat_orfs(orf_dir, other_fastas=None, outfile='all_orfs.fasta'):
    with open(outfile, "w") as oh:
        for f in glob.glob(os.path.join(orf_dir, '*proteins.fasta')):
            with open(f) as ih:
                for l in ih: print(l.strip(), file=oh)
        if other_fastas is not None:
            for f in other_fastas:
                with open(f) as ih:
                    for l in ih: print(l.strip(), file=oh)
    return outfile


def run_cd_hit(input_fa, output_fa, c=0.9, G=1, b=20, M=5000,
    T=1, n=5, l=10, t=2, d=0, s=0.0, S=999999, aL=0.0, AL=99999999, aS=0.0,
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
    # TODO: for some reason, *.clstr file not produced when using file_transaction function. 
    # it would be better to run it this way if we can figure out what the issue is.
    #with file_transaction(output_files) as tx_out_files:
    #    cmd = ("cd-hit -i {input_fasta} -o {output_fasta} -c {c} "
    #            "-G {G} -b {b} -M {M} -T {T} -n {n} -l {l} -t {t} "
    #            "-d {d} -s {s} -S {S} -aL {aL} -AL {AL} -aS {aS} "
    #            "-AS {AS} -A {A} -uL {uL} -uS {uS} -U {U} "
    #            "-p 1 -g {g} -sc {sc} -sf {sf}").format(input_fasta=input_fa,
    #                                                    output_fasta=tx_out_files[0],
    #                                                    **locals())
    #    subprocess.check_call(cmd, shell=True)
    cmd = ("cd-hit -i {input_fasta} -o {output_fasta} -c {c} "
                "-G {G} -b {b} -M {M} -T {T} -n {n} -l {l} -t {t} "
                "-d {d} -s {s} -S {S} -aL {aL} -AL {AL} -aS {aS} "
                "-AS {AS} -A {A} -uL {uL} -uS {uS} -U {U} "
                "-p 1 -g {g} -sc {sc} -sf {sf}").format(input_fasta=input_fa,
                                                        output_fasta=output_fa,
                                                        **locals())
    subprocess.check_call(cmd, shell=True)
    return output_files


def id_added_seeds(clstr, original_seed_fasta = None):
    seeds = set()
    if original_seed_fasta is not None:
        for name, seq in readfa(open(original_seed_fasta)):
            seeds.add(name.split()[0])

    cluster_map = defaultdict(list)
    new_seeds = []
    member_seeds = 0
    old_seeds = 0

    with open(clstr) as fh:
        for cluster_start, group in itertools.groupby(fh, lambda l: l[0] == '>'):
            members = []
            rep_seq = ''
            if not cluster_start:
                for line in group:
                    if "*" in line:
                        rep_seq = line.split(",")[1].split("...")[0].replace(">",'').replace(" ","")

                    members.append(line.split(",")[1].split("...")[0].replace(">",'').replace(" ",""))

            if len(rep_seq) == 0:
                continue

            seed_seq = None

            if rep_seq in seeds:
                seed_seq = rep_seq
                members.remove(rep_seq)
                old_seeds += 1
            else:
                for m in members:
                    if m in seeds:
                        seed_seq = m
                        member_seeds += 1
                        members.remove(m)
                        break

            if seed_seq is None:
                seed_seq = rep_seq
                members.remove(rep_seq)
                new_seeds.append(rep_seq)
            cluster_map[seed_seq] = members
    print('there are', len(new_seeds), 'new seeds', file=sys.stderr)
    print('there are', old_seeds, 'old seeds', file=sys.stderr)
    print('there are', member_seeds, 'old seeds included in sequence clusters that will serve as the cluster seed', file=sys.stderr)
    return cluster_map, new_seeds


def write_cluster_map(cmap, out_map):
    count = 0
    with open(out_map, "w") as oh:
        for k in cmap:
            print(k, k, sep = "\t", file=oh)
            for val in cmap[k]:
                print(val, k, sep="\t", file=oh)
    return out_map

def read_cluster_map(map_file):
    assert op.exists(map_file), "Could not find find map file for ORF clusters within filesystem"
    
    cmap = {}
    with open(map_file) as ih:
        for l in ih:
            vec = l.strip().split("\t")
            cmap[vec[0]] = vec[1]
    return cmap

def swap_cluster_map(cm):
    cm_swap = {}
    for c in cm:
        for k in cm[c]: cm_swap[k] = c
    return cm_swap


def write_new_seeds(new_seed_fa, clstr_fasta, new_seed_list):
    if os.path.exists(new_seed_fa):
        print("new seed fasta file already exists, skipping.", file=sys.stderr)
        return new_seed_fa
    with open(new_seed_fa, "w") as oh:
        for name, seq in readfa(open(clstr_fasta)):
            if name.split()[0] in new_seed_list:
                print(">{}".format(name), file=oh)
                for i in range(0, len(seq), 60):
                    print(seq[i:i+60], file=oh)
    return new_seed_fa


def cluster_split_fa(fasta, outdir, seq_num=1000, pctid=70):
    with tmp_dir() as tmp:
        outfasta = os.path.join(tmp, fasta.replace(".f","_{}.f".format(pctid)))
        outfiles = run_cd_hit(fasta, outfasta, c=0.7)
        count = 0
        fa = Fasta(fasta)
        number = 0
        missing = []
        seedout = id_added_seeds(outfiles[1], None)
        cluster_map = seedout[0]

        for k in cluster_map:
            seqs = [k] + cluster_map[k]
            count += len(seqs)
            with open(os.path.join(outdir, "subset_{}.fasta".format(number)), "a") as oh:
                for s in seqs:
                    try:
                        rec = fa[s]
                    except:
                        missing.append(s)
                        continue

                    write_fa_record(rec.long_name, str(rec), oh)

            if count > 1000:
                count = 0
                number += 1
        print(number, "files created", file=sys.stderr)
    return outdir


def prep_contigs(falist, wd='./viruscope', old_seeds=None):
    '''sends sequences through the orf_setup workflow
    Args:
        falist (list): list of fasta files to process (expected to be SAG contigs)
        wd (str): path to the output directory location
        old_seeds (list): list of fasta files to use as the old seed sequences to be included in contig clustering
    Returns:
        str: path to output directory containing clustering results (seed maps, seed fasta files, and fasta files ready for mica).
    '''
    wd = safe_makedir(wd)
    clustdir = safe_makedir(op.join(wd, 'clustering'))
    for_mica = safe_makedir(op.join(clustdir, "for_mica"))
    
    cmap_file = op.join(clustdir,"seed_map90.tsv")
    new_seeds = op.join(clustdir, "new_seeds.fasta")
    
    prods = run_batch_prodigal(falist, wd)
    all_orfs = concat_orfs(op.join(wd, 'prodigal'), other_fastas=[old_seeds], outfile = op.join(clustdir, "all_orfs.fasta"))
    cd_hit_out = run_cd_hit(all_orfs, all_orfs.replace(".fasta","_clust90.fasta"), c=0.9, d=0)
    seedout = id_added_seeds(cd_hit_out[1], old_seeds)
    out = write_cluster_map(seedout[0], cmap_file)
    cmap = swap_cluster_map(seedout[0])

    new_seeds = write_new_seeds(new_seeds, cd_hit_out[0], seedout[1])
    for_mica = cluster_split_fa(new_seeds, for_mica)
    return clustdir

