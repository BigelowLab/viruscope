import click
import sys
import os.path as op
import os
import glob
from viruscope import recruit, orf_setup, phage_count, summary, tools, pbs


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.0.0')
@click.pass_context
def cli(obj):
    """batch viruscope packages."""
    pass


@cli.command('orf-setup', short_help='set up orf clusters for viruSCope')
@click.argument('fadir', nargs=1)
@click.argument('wd', nargs=1)
@click.option('--old-seeds', 
              help='list of fasta files containing previously used seeds',
              default = None,
              show_default=True)
@click.option('--use-pbs',
              is_flag=True,
              help='output mica and diamond commands as a pbs sub (Bigelow only)')
@click.option('--threads',
              help='threads to use if writing a pbs sub',
              default=20,
              show_default=True)
@click.option('--queue',
              help='queue to use if writing pbs sub',
              default='scgc-route',
              show_default=True)
@click.option('--sub-out',
              help='directory to write pbs sub outfile to',
              default=None)
@click.option('--mem',
              help='amount of memory to request for pbs subs',
              default=5,
              show_default=True)
@click.option('--bac-mg',
              help='location of bacterial metagenome for recruitment step',
              default='/mnt/scgc_nfs/ref/viral_dbs/LineP-all.fasta.gz',
             show_default=True)
@click.option('--vir-mg',
              help='location of bacterial metagenome for recruitment step',
              default='/mnt/scgc_nfs/ref/viral_dbs/POV.fasta.gz',
             show_default=True)
def set_up_clusters(fadir, wd, old_seeds, use_pbs, threads, mem, queue, sub_out, bac_mg, vir_mg):
    wd = tools.safe_makedir(wd)
    falist = glob.glob(op.join(fadir, "*_*.fasta"))
    clust_dir = orf_setup.prep_contigs(falist, wd, old_seeds)
    
    if use_pbs:
        micadir = op.join(clust_dir, 'for_mica')
        assert op.exists(micadir), 'could not find the directory containing new seeds for mica'
        
        outdir = tools.safe_makedir(op.join(wd, 'blast'))
        mica_outfile = op.join(wd, 'mica_pbs_sub.sh')
        diamond_outfile = op.join(wd, 'diamond_pbs_sub.sh')
        
        if sub_out is None:
            mica_sub_out = op.join(wd, 'mica_pbs_sub.out')
            diamond_sub_out = op.join(wd, 'diamond_pbs_sub.out')
        else:
            mica_sub_out = op.join(sub_out, 'mica_pbs_sub.out')
            diamond_sub_out = op.join(sub_out, 'diamond_pbs_sub.out')
        
        mica_sub = pbs.write_mica_array_sub(micadir, 
                                            outdir, 
                                            outfile=mica_outfile, 
                                            subout=mica_sub_out,
                                            queue=queue, 
                                            threads=threads)
        diamond_sub = pbs.write_diamond_array_sub(fadir, wd, subout=diamond_sub_out, queue=queue, threads=threads, 
                                                  outfile=diamond_outfile, bac_mg=bac_mg, vir_mg=vir_mg)
        print('PBS array submission script for mica written to {}'.format(mica_sub), file=sys.stderr)
                                
        print('PBS array submission script for diamond written to {}'.format(diamond_sub), file=sys.stderr)
        


@cli.command('recruit-single', short_help='run diamond recruitment of vir and bac mgs against SAG')
@click.argument('prot-fasta', nargs=1)
@click.argument('vir-mg', nargs=1)
@click.argument('bac-mg', nargs=1)
@click.option('--sag-contigs',
             help='location of sag contigs in either fasta or gff format',
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
def recruit_single(prot_fasta, vir_mg, bac_mg, sag_contigs, output, threads, verbose):
    recruit.run_recruitment(prot_fasta, vir_mg, bac_mg, sag_contigs, output, threads, verbose)
    

@cli.command('summarize', short_help='summarize results of diamond and mica, create viruscope summary statistic')
@click.argument('contig-dir', nargs=1)
@click.argument('wd', nargs=1)
@click.option('--seed-classifications', default=None, nargs=1, show_default=True)
@click.option('--training-file', help='location of training file for knn classifier',
              default='/mnt/scgc_nfs/opt/viruscope/virus-training.csv',
              show_default=True)
def summarize_results(contig_dir, wd, seed_classifications, training_file):
    phage_count.write_blast_summaries(wd, seed_classifications)
    summary.write_batch_summaries(contig_dir, wd, training_file=training_file)
    

if __name__=='__main__':
    cli()