import glob
from .phage_count.py import run_mica


def write_mica_array_sub(micadir, outdir, out_file=None, subout="./mica.out", queue='scgc-route', threads=20):
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

#PBS -l ncpus={threads}                                                                           

#PBS -j oe                                                                                              
#PBS -o {subout}

module load mica

num=$(($PBS_ARRAY_INDEX-1))

{mica_cmd}'''.format(total = len(fastas), queue=queue, subout=subout, threads=threads, mica_cmd=mica_cmd)
    if outfile:
        with open(outfile, "w"):
            print(sub)
        
    return sub


#def write_diamond_array_sub