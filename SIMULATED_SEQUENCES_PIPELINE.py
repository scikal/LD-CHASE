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

def transitions(chr_id):
    """ Generates transitions between BPH to SPH regions for trisomy of meiosis II origin. """
    x = int(chr_id[3:]) if chr_id[3:].isnumeric() else chr_id[3:]
    if type(x) is int and 1<=x<=6:
        #BPH-SPH-BPH-SPH
        result = ('BPH',random.uniform(0,.25),random.uniform(.5,.75),random.uniform(.75,1))
    elif type(x) is int and 7<=x<=12:
        #BPH-SPH-BPH
        result = ('BPH',random.uniform(0,.333),random.uniform(.666,1))
    elif x=='X' or (type(x) is int and 13<=x<=22):
        #SPH-BPH
        result = ('SPH',random.uniform(.5,1))
    else:
        result = ('SPH',1)
    return (result,)

def contrast_crossovers_wrapper(disomy_obs_filename,monosomy_obs_filename,chr_id,sp,ancestral_makeup,model,min_reads,max_reads,output_dir,ref_dir,output_filename):
    from DETECT_CROSSOVERS import contrast_crossovers
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

def sample_indv(sp):    
    seed(None, version=2)
    list_SP = sp.split('_')
    if len(list_SP)==1:
        INDIVIDUALS = read_ref(ref_path + f"samples_per_panel/{sp:s}_panel.txt") 
        A = sample(INDIVIDUALS,k=3)
        B = choices(['A','B'],k=3)
    
    elif len(list_SP)==2 and admixture_proportions==[]:
        B = choices(['A','B'],k=3)
        A = []
        for i,p in enumerate(random.sample(list_SP, len(list_SP)),start=1): #Sampling without replacement.
            INDIVIDUALS = read_ref(ref_path + f"samples_per_panel/{p:s}_panel.txt") 
            A.extend(sample(INDIVIDUALS,k=i))
    
    elif len(list_SP)==2 and len(admixture_proportions)==2:
        B = choices(['A','B'],k=6)
        A = []
        for p in list_SP:
            INDIVIDUALS = read_ref(ref_path + f"samples_per_panel/{p:s}_panel.txt") 
            A.extend(sample(INDIVIDUALS,k=3))
        A = operator.itemgetter(0,3,1,4,2,5)(A)
    
    else:
        print('error: unsupported sp value.')

    C = [i+j for i,j in zip(A,B)]
    return A,B,C 


def main(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions=[],mismatch=False):
    work_dir = work_dir.rstrip('/') + '/' if len(work_dir)!=0 else ''
    SPsorted = {('EUR_EAS'): 'EAS_EUR', ('EAS_EUR'): 'EAS_EUR', ('EUR_SAS'): 'SAS_EUR', ('SAS_EUR'): 'SAS_EUR', ('EAS_SAS'): 'EAS_SAS', ('SAS_EAS'): 'EAS_SAS', ('AFR_EUR'): 'AFR_EUR', ('EUR_AFR'): 'AFR_EUR', 'EUR': 'EUR', 'EAS': 'EAS', 'SAS': 'SAS', 'AMR': 'AMR', 'AFR': 'AFR'}
    
    A,B,C = sample_indv(sp)
    
    print('Simulating effective haploids:')    
    for a,b in zip(A,B): 
        print(a+b)
        simulate_haploids(a, SPsorted[sp], chr_id, b, work_dir)
    sim_obs_tabs = [f'{work_dir:s}{c:s}.{chr_id:s}.hg38.obs.p' for c in C]
    fns = MixHaploids_wrapper(*sim_obs_tabs[:2], read_length=read_length, depth=depth, scenarios=('disomy','monosomy'),
                                    transitions=transitions(chr_id), compress='bz2', output_dir=work_dir, 
                                    distant_admixture=admixture_proportions)
    
    fns += MixHaploids_wrapper(sim_obs_tabs[2], read_length=read_length, depth=depth, scenarios=('monosomy',),
                                    transitions=transitions(chr_id), compress='bz2', output_dir=work_dir, 
                                    distant_admixture=admixture_proportions)
    
    
    disomy_obs_filename, monosomy0_obs_filename, monosomy1_obs_filename = fns
        
    print(disomy_obs_filename)
    print(monosomy0_obs_filename)
    print(monosomy1_obs_filename)
    
    if mismatch:
        effective_SP = choice(sp.split('_'))
        ancestral_makeup = {}
    else:
        effective_SP = SPsorted[sp] 
        ancestral_makeup = dict(zip(sp.split('_'),admixture_proportions))
        
    #print(effective_SP)
    print(ancestral_makeup)

    output_filename0 = f'simulated.SPH.{C[0]:s}.{C[1]:s}.{chr_id:s}.LLR.p'
    output_filename1 = f'simulated.BPH.{C[0]:s}.{C[1]:s}.{C[2]:s}.{chr_id:s}.LLR.p'    
    
    contrast_crossovers_wrapper(disomy_obs_filename,monosomy0_obs_filename,chr_id,effective_SP,ancestral_makeup,None,min_reads,max_reads,work_dir,ref_path,output_filename0)
    contrast_crossovers_wrapper(disomy_obs_filename,monosomy1_obs_filename,chr_id,effective_SP,ancestral_makeup,None,min_reads,max_reads,work_dir,ref_path,output_filename1)
    for c in C: os.remove(f'{work_dir:s}{c:s}.{chr_id:s}.hg38.obs.p')
    return 0


if __name__ == "__main__":
    random.seed(a=None, version=2)
    admixture_proportions=[] #[0.5,0.5]
    depth=0.05
    sp='EUR' #'EUR_EAS'
    chr_id='chr16'
    read_length = 75
    min_reads,max_reads = 12,6
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_EUR_chr16"
    mismatch = False #Partial mismatch of reference panel for admixtures.
    ###main(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch)
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    
    
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_EAS_chr16"
    sp='EAS'
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_SAS_chr16"
    sp='SAS'
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_AMR_chr16"
    sp='AMR'
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    
    work_dir = f"/{HOME:s}/ariad/Dropbox/postdoc_JHU/Project2_Trace_Crossovers/CONTRAST-CROSSOVERS/results/nonadmixed_AFR_chr16"
    sp='AFR'
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
    runInParallel(*([main]*32),args=(depth,sp,chr_id,read_length,min_reads,max_reads,work_dir,admixture_proportions, mismatch) )
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
