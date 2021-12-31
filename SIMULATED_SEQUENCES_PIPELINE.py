#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""

BUILD_SIMULATED_SEQUENCES

Daniel Ariad (daniel@ariad.org)
Dec 30, 2020

"""

    
import time, sys, random, os, operator, collections

from random import sample, choices, seed, choice
from multiprocessing import Process

### import pypy3_popcounts.build

if sys.platform == "linux" or sys.platform == "linux2":
    print('Detected OS: linux.')
    HOME='home'
elif sys.platform == "darwin":
    print('Detected OS: macOS.')
    HOME='Users'
elif sys.platform == "win32":
    print('Detected OS: windows.')
    HOME='????????'
    
###sys.path.append('/home/ariad/Dropbox/postdoc_JHU/Project1_LD-PGTA/LD-PGTA_ecosystem/LD-PGTA_V2')
###os.chdir('/home/ariad/Dropbox/postdoc_JHU/Project1_LD-PGTA/LD-PGTA_ecosystem/LD-PGTA_V2')

from MIX_HAPLOIDS import MixHaploids_wrapper
from EXTRACT_GENOTYPES import extract

leg_tuple = collections.namedtuple('leg_tuple', ('chr_id', 'pos', 'ref', 'alt')) #Encodes the rows of the legend table
sam_tuple = collections.namedtuple('sam_tuple', ('sample_id', 'group1', 'group2', 'sex')) #Encodes the rows of the samples table
obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table

ref_path = f'/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/reference_panels/'

def read_ref(filename):
    with open(filename, 'r') as data_in:
        tab = tuple(str(line.replace('\n','')) for line in data_in)
    return tab

def runInParallel(*fns,**kwargs):
    proc = []
    for fn in fns:
        try:
            p = Process(target=fn,args=kwargs.get('args',tuple()))
            p.start()
            proc.append(p)
            time.sleep(5)
        except Exception as error:
            print('caution: a process failed!')
            print(error)
    for p in proc:
        try:
          p.join()
        except:
            None

def contrast_crossovers_wrapper(disomy_obs_filename,monosomy_obs_filename,ancestral_makeup,chr_id,model,min_reads,max_reads,output_dir,ref_dir,output_filename):
    from DETECT_CROSSOVERS import contrast_crossovers
    sp = '_'.join(sorted(ancestral_makeup))
    args = dict(disomy_obs_filename = disomy_obs_filename,
                monosomy_obs_filename = monosomy_obs_filename,
                hap_filename = ref_path + f'{sp:s}_panel/{chr_id:s}_{sp:s}_panel.hap.gz',
                leg_filename = ref_path + f'{sp:s}_panel/{chr_id:s}_{sp:s}_panel.legend.gz',
                samp_filename = ref_path + f'{sp:s}_panel/{sp:s}_panel.samples.gz',
                window_size = 0,
                subsamples = 300,
                offset = 0,
                min_reads = min_reads, #3,
                max_reads = max_reads, #8,
                min_HF = 0.05,
                min_score = 2,
                output_dir = output_dir, #f'results_{sp:s}/',
                output_filename = output_filename,
                compress = 'bz2',
                ancestral_makeup = ancestral_makeup,
                seed=0)
                #model = model)
    LLR_dict, info = contrast_crossovers(**args)

    return LLR_dict, info

def simulate_haploids(sample_id,sp,chr_id,genotypes,output_dir):
    """ Wraps the function 'extract'. """
    
    leg_filename = ref_path + f'{sp:s}_panel/{chr_id:s}_{sp:s}_panel.legend.gz'
    hap_filename = ref_path + f'{sp:s}_panel/{chr_id:s}_{sp:s}_panel.hap.gz'
    sam_filename = ref_path + f'{sp:s}_panel/{sp:s}_panel.samples.gz'
    return extract(leg_filename,hap_filename,sam_filename,chr_id,sample_id,genotypes=genotypes,output_dir=output_dir)

def sample_indv(ancestral_makeup):    
    seed(None, version=2)
    
    if len(ancestral_makeup)==1:
        sp = next(iter(ancestral_makeup))
        INDIVIDUALS = read_ref(ref_path + f"samples_per_panel/{sp:s}_panel.txt") 
        A = sample(INDIVIDUALS,k=2)
        B = choices(['A','B'],k=2)
    
    elif len(ancestral_makeup)==2 and type(ancestral_makeup)==set:
        B = choices(['A','B'],k=2)
        A = []
        for i,sp in enumerate(random.sample(tuple(ancestral_makeup), len(ancestral_makeup)),start=1): #Sampling without replacement.
            INDIVIDUALS = read_ref(ref_path + f"samples_per_panel/{sp:s}_panel.txt") 
            A.extend(sample(INDIVIDUALS,k=i))
    
    else:
        raise Exception('error: unsupported sp value.')

    C = [i+j for i,j in zip(A,B)]
    return A,B,C 


def main(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir):
    print(ancestral_makeup)
    
    work_dir = work_dir.rstrip('/') + '/' if len(work_dir)!=0 else ''
    
    A,B,C = sample_indv(ancestral_makeup)
    
    print('Simulating effective haploids:')    
    
    sp = '_'.join(sorted(ancestral_makeup))
    
    for indv,genotype in zip(A,B): 
        print(indv+genotype)
        simulate_haploids(indv, sp, chr_id, genotype, work_dir)
    
    sim_obs_tabs = [f'{work_dir:s}{c:s}.{chr_id:s}.hg38.obs.p' for c in C]
    
    fns = MixHaploids_wrapper(*sim_obs_tabs, read_length=read_length, depth=depth, scenarios=('disomy','monosomy'),
                              compress='bz2', output_dir=work_dir)
    
    fns += MixHaploids_wrapper(sim_obs_tabs[2], read_length=read_length, depth=depth, scenarios=('monosomy',),
                               compress='bz2', output_dir=work_dir)
    
    disomy_obs_filename, monosomy0_obs_filename, monosomy1_obs_filename = fns
            
    print(disomy_obs_filename)
    print(monosomy0_obs_filename)
    print(monosomy1_obs_filename)
    
    output_filename0 = f'simulated.SPH.{C[0]:s}.{C[1]:s}.{chr_id:s}.LLR.p'
    output_filename1 = f'simulated.BPH.{C[0]:s}.{C[1]:s}.{C[2]:s}.{chr_id:s}.LLR.p'    
    
    contrast_crossovers_wrapper(disomy_obs_filename,monosomy0_obs_filename,ancestral_makeup,chr_id,None,min_reads,max_reads,work_dir,ref_path,output_filename0)
    contrast_crossovers_wrapper(disomy_obs_filename,monosomy1_obs_filename,ancestral_makeup,chr_id,None,min_reads,max_reads,work_dir,ref_path,output_filename1)
    
    for c in C: os.remove(f'{work_dir:s}{c:s}.{chr_id:s}.hg38.obs.p')

    return 0


if __name__ == "__main__":
    random.seed(a=None, version=2)
    #admixture_proportions=[] #[0.5,0.5]
    depth=0.05
    
    mismatch = False #Partial mismatch of reference panel for admixtures.
    chr_id='chr16'
    read_length = 75
    min_reads,max_reads = 12,6
    ancestral_makeup={'EUR'} #'EUR_EAS'
    
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_{'_'.join(sorted(ancestral_makeup)):s}_{chr_id:s}"

    ###main(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir)
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    
    
    ancestral_makeup={'EAS'}
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_{'_'.join(sorted(ancestral_makeup)):s}_{chr_id:s}"
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    
    ancestral_makeup={'SAS'}
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_{'_'.join(sorted(ancestral_makeup)):s}_{chr_id:s}"
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    
    ancestral_makeup={'AMR'}
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_{'_'.join(sorted(ancestral_makeup)):s}_{chr_id:s}"
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    
    ancestral_makeup={'AFR'}
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_{'_'.join(sorted(ancestral_makeup)):s}_{chr_id:s}"
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    
    ancestral_makeup={'EAS_EUR'}
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_{'_'.join(sorted(ancestral_makeup)):s}_{chr_id:s}"
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    
    ancestral_makeup={'SAS_EUR'}
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_{'_'.join(sorted(ancestral_makeup)):s}_{chr_id:s}"
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    
    ancestral_makeup={'AFR_EUR'}
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_{'_'.join(sorted(ancestral_makeup)):s}_{chr_id:s}"
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    
    ancestral_makeup={'EAS_SAS'}
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_{'_'.join(sorted(ancestral_makeup)):s}_{chr_id:s}"
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    runInParallel(*([main]*32),args=(ancestral_makeup,depth,chr_id,read_length,min_reads,max_reads,work_dir) )
    
    #for n in [*range(22,0,-1),'X']:
        #chr_id = 'chr' + str(n)
    #    runInParallel(*([main]*12),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir) )

        #main(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,distant_admixture)
        #for sp in ('AFR_EUR','EAS_SAS','SAS_EUR','EAS_EUR'):
            #work_dir = f"/mybox/F1-simulations/results_mixed_{sp:s}" #'../results' #'results_EAS'
        #runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,distant_admixture) )
    print('DONE.')
    pass
else:
    print("The module BUILD_SIMULATED_SEQUENCES was imported.")
