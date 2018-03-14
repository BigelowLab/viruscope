import glob
import os.path as op
from .phage_count import run_mica
from .tools import safe_makedir


def write_mica_array_sub(micadir, outdir, outfile=None, subout="./mica.out", queue='scgc-route', threads=20, mem=5):
    fastas = glob.glob(op.join(micadir, "*.fasta"))
    
    mica_cmd = run_mica(op.join(micadir, "subset_${num}.fasta"), op.join(outdir, "subset_${num}.mica"), threads=threads)
    sub = '''#!/bin/bash                                                                                                             
## set name of PBS job                                                                                                  
#PBS -N newseeds_mica

## set the job array variable                                                                                           
#PBS -J 1-{total}                                                                                        

## send the environment variables with job                                                                              
#PBS -V

## set the queue                                                                                                        
#PBS -q {queue}
#PBS -l walltime=24:00:00                                                                                               

#PBS -l ncpus={threads},mem={mem}G                                                                           

#PBS -j oe                                                                                              
#PBS -o {subout}

module load mica

num=$(($PBS_ARRAY_INDEX-1))

{mica_cmd}'''.format(total = len(fastas), queue=queue, subout=subout, mem=mem, threads=threads, mica_cmd=mica_cmd)
    if outfile:
        with open(outfile, "w") as oh:
            print(sub, file=oh)
        
    return outfile


def write_diamond_array_sub(contigdir, wd, subout, queue, threads, mem, outfile, bac_mg, vir_mg):
    assert op.exists(op.join(wd, 'prodigal')), "could not find prodigal output directory"
    outdiamond = safe_makedir(op.join(wd, 'diamond'))
    fasta_glob = op.join(contigdir,'*_contigs.fasta')
    total = len(glob.glob(fasta_glob))
    contigpath = '"${fastalist[$num]}"'
    proteinpath=op.join(wd, 'prodigal', '${name}_proteins.fasta')
    

    sub = '''#!/bin/bash

## set name of PBS job
#PBS -N diamond_array

## set the job array variable
#PBS -J 1-{total}

## send the environment variables with job
#PBS -V

## set the queue
#PBS -q {queue}
#PBS -l walltime=24:00:00

#PBS -l ncpus={threads},mem={mem}G

#PBS -j oe
#PBS -o {subout}

module purge
module load anaconda3
module load diamond

num=$(($PBS_ARRAY_INDEX-1))

fastalist=({fasta_glob})

fastapath={contigpath}

base=$(basename $fastapath)

name=$(echo $base | cut -f 1 -d '_')

proteins={proteinpath}
bac_db={bac_mg}
vir_db={vir_mg}

batch-viruscope recruit-single \
--threads 20 --output {outdiamond} \
--sag-contigs $fastapath $proteins $vir_db $bac_db
'''.format(total=total, subout=subout, threads=threads, mem=mem, queue=queue, fasta_glob=fasta_glob, 
           contigpath=contigpath, proteinpath=proteinpath, outdiamond=outdiamond, 
           bac_mg=bac_mg, vir_mg=vir_mg)
    if outfile:
        with open(outfile, "w") as oh:
            print(sub, file=oh)
        
    return outfile