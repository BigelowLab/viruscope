import pandas as pd
from collections import defaultdict
import itertools
from nb_tools import readfa, swap_cluster_map, orf_map


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



def phage_contig_table(clstr_map, gff, phage_hits_df, outfile=None):
    ''' create a summary of phage hits to contig orfs, mapping back cd-hit cluster seeds to orfs from individual genomes
    Args:
        clstr_map (dict): cluster map
        gff (path): path to gff output file for mapping orfs back to contigs
        phage_hits_df (pandas.DataFrame): phage hits dataframe from summary of all Mica results
        outfile (path): where to write output table, if None, none written
    Returns:
        pandas dataframe 
    '''
    cm = swap_cluster_map(clstr_map)
    om = orf_map(gff)
    om['lookup'] = [cm.get(i, i) for i in om['id']]
    omp = pd.merge(phage_hits_df, om, left_on='orf', right_on='lookup', how='right').fillna(0)[['contig','id','p1','p2','len']]
    orf_count = pd.DataFrame(om.groupby('contig')['id'].count()).rename(columns={'id':'total_orfs'})
    csum = pd.concat([summarize_by_contig(omp, 'p1'), summarize_by_contig(omp, 'p2'), orf_count], axis=1)
    csum['viral_phage_gene_fraction'] = csum['p1'] / csum['total_orfs']
    csum['viral2_phage_gene_fraction'] = csum['p2'] / csum['total_orfs']

    if outfile is not None:
        csum.to_csv(outfile)

    return csum

