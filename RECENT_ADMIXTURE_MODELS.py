#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RECENT_ADMIXTURE_MODELS

Given reads that originated form disomy and a monosomy of a sibling, the 
likelihood of observed reads within a genomic window under the SPH and BPH
scenarios is calculated. The calculation is done using a reference panel
of two populations. Moreover, the statistical models are based on the
assumption of a recent-admixtures, where the parents are associated with
different ancestral populations.

BPH (Both Parental Homologs) correspond to the presence of three unmatched
haplotypes, while SPH (Single Parental Homolog) correspond to chromosome gains
involving identical homologs.

Daniel Ariad (daniel@ariad.org)
Dec 21, 2020
"""

import pickle, os, sys, bz2, collections, gzip, platform

from functools import reduce
from operator import and_, itemgetter
from itertools import combinations

leg_tuple = collections.namedtuple('leg_tuple', ('chr_id', 'pos', 'ref', 'alt')) #Encodes the rows of the legend table
sam_tuple = collections.namedtuple('sam_tuple', ('sample_id', 'group1', 'group2', 'sex')) #Encodes the rows of the samples table
obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table

if platform.python_implementation()=='PyPy':
    class PopCount:
        def __init__(self):
            self.A = bytes((bin(i).count('1') for i in range(1<<20)))
    
        def __call__(self,x):
            result = 0
            while(x): result += self.A[x & 1048575]; x >>= 20
            return result
    popcount = PopCount()
else:
    try: 
        from gmpy2 import popcount
    except ModuleNotFoundError: 
        print('caution: the module gmpy2 is missing.')
        def popcount(x):
            """ Counts non-zero bits in positive integer. """
            return bin(x).count('1')

class recent_admixture:
    """ Based on the statisitcal models (models_dict) and the reference panel
    (leg_tab, hap_tab and sam_tab), it allows to calculate the likelihoods of
    observed alleles under various statistical models (monosomy, disomy, SPH
    and BPH). """


    def __init__(self, obs_tab, leg_tab, hap_tab, sam_tab, models_dict, total_number_of_haplotypes):
        """ Initialize the attributes of the class. """

        if len(leg_tab)!=len(hap_tab):
            raise Exception('Error: the number of SNPs in the LEGEND file differ from the number of SNPs in the HAP file.')

        if total_number_of_haplotypes!=2*len(sam_tab):
            raise Exception('Error: the number of diploid samples in the SAMPLE file differ from the number of haplotypes in the HAP file.')

        self.total_number_of_haplotypes_in_reference_panel = total_number_of_haplotypes
        self.models_dict = models_dict
        self.hap_dict, self.fraction_of_matches = self.build_hap_dict(obs_tab, leg_tab, hap_tab)
        self.flags, self.num_of_hap_in_ref_subpanel, self.name2id = self.subpanels(sam_tab)

    def subpanels(self, sam_tab):
        """ Differentiates between the two groups that compose the reference
        panel. Then, all the haplotypes that are associated with each group are
        flagged using a binary representation marks and counted. """

        differentiate = [row.group2 == sam_tab[0].group2 for row in sam_tab for i in (1,2)]
        flag0 = sum(v<<i for i, v in enumerate(reversed(differentiate)))
        flag1 = flag0 ^ ((1 << self.total_number_of_haplotypes_in_reference_panel) - 1)
        flags = (flag0, flag1)
        name2id = {sam_tab[0].group2:0,
                   sam_tab[differentiate[::2].index(False,1)].group2:1}
        N0 = differentiate.count(True)
        N1 = self.total_number_of_haplotypes_in_reference_panel - N0
        number_of_haplotypes_in_reference_subpanel = (N0,N1)
        return flags, number_of_haplotypes_in_reference_subpanel, name2id

    def build_hap_dict(self, obs_tab, leg_tab, hap_tab):
        """ Returns a dictionary that lists SNP alleles and gives their
        relevent row from haplotypes table. The row is stored as bits, where
        True means that the haplotype contains the allele. We denote the
        returned dictionary as the reference panel. """

        hap_dict = dict()
        mismatches = 0
        combined = {pos: (ref,alt,hap) for (chr_id,pos,ref,alt),hap in zip(leg_tab, hap_tab)}
        missing = 3*(None,)


        b = (1 << self.total_number_of_haplotypes_in_reference_panel) - 1 #### equivalent to int('1'*number_of_haplotypes,2)

        for (pos, read_id, base) in obs_tab:
            ref, alt, hap = combined.get(pos, missing)
            if base==alt:
                hap_dict[(pos,base)] = hap
            elif base==ref:
                hap_dict[(pos,base)] = hap ^ b ### ^b flips all bits of the binary number, hap_tab[ind] using bitwise xor operator.
            else:
                mismatches += 1

        fraction_of_matches = 1-mismatches/len(obs_tab)

        print('Algorithm for recent-admixtures: %.2f%% of the observed alleles matched the reference panel.' % (100*fraction_of_matches))

        return hap_dict, fraction_of_matches

    def intrenal_hap_dict(self, *alleles):
        """ This function allows treatment of alleles and haplotypes on an
        equal footing. This is done in three steps: (1) All the alleles and
        haplotypes are enumerated. (2) For each given haplotype, tuples in the
        reference panel, associated with the haplotype's alleles, are
        intersected. (3) A dictionary that lists all the alleles and haplotypes
        by their index is returned. The dictionary gives for each allele and
        haplotype their associated tuple and intersected tuple, respectively. """
        hap = dict()

        for i, X in enumerate(alleles):
            if type(X[0])==tuple: #Checks if X is a tuple/list of alleles.
                n = len(X)
                if n==1:
                    hap[1 << i] = self.hap_dict[X[0]]
                elif n==2:
                    hap[1 << i] = self.hap_dict[X[0]] & self.hap_dict[X[1]]
                else:
                    hap[1 << i] = reduce(and_,itemgetter(*X)(self.hap_dict))

            elif type(X[0])==int: #Checks if X is a single allele.
                hap[1 << i] = self.hap_dict[X]
            else:
                raise Exception('error: joint_frequencies only accepts alleles and tuple/list of alleles.')

        return hap

    def joint_frequencies_combo(self, *alleles, group2_id, normalize):
        """ Based on the reference panel, it calculates joint frequencies of
            observed alleles. The function arguments are alleles, that is,
            tuples of position and base, e.g., (100,'T'), (123, 'A') and
            (386, 'C'). Each allele is enumerated according to the order it
            was received by the function. The function returns a dictionary that
            lists all the possible subgroups of the given alleles. Each key in
            the dictionary is a tuple of intergers that are lexicographically
            sorted. Moreover, each integer within the keys corresponds to an
            enumerated allele. For each subgroup of alleles the dictionary
            gives the joint frequencies in a given population. The function
            arguments can also include haplotypes, that is, tuples of alleles;
            Haplotypes are treated in the same manner as alleles. """

        flag = self.flags[group2_id]
        hap = self.intrenal_hap_dict(*alleles)

        result = {c: popcount(A&flag) for c,A in hap.items()}

        for C in combinations(hap, 2):
            result[C[0]|C[1]] = popcount(hap[C[0]]&hap[C[1]]&flag)

        for C in combinations(hap, 3):
            result[C[0]|C[1]|C[2]] = popcount(hap[C[0]]&hap[C[1]]&hap[C[2]]&flag)

        for r in range(4,len(alleles)):
            for C in combinations(hap, r):
                result[sum(C)] = popcount(reduce(and_,itemgetter(*C)(hap))&flag)

        if len(alleles)>=4:
            result[sum(hap.keys())] = popcount(reduce(and_,hap.values())&flag)

        if normalize:
            N = self.num_of_hap_in_ref_subpanel[group2_id]
            result = {k: v/N for k,v in result.items()}

        return result

    def likelihoods(self, *alleles):
        """ Calculates the likelihood to observe a set with alleles
        and haplotypes under four scenarios, namely, monosomy, disomy, SPH
        and BPH. """

        model = self.models_dict[len(alleles) - 1] #The number of reads from the disomy is smaller by one from the total number of given reads.
        F = self.joint_frequencies_combo(*alleles, group2_id=0, normalize=False)
        M = self.num_of_hap_in_ref_subpanel[0] #Divide values by M to normalize the joint frequencies, F.
        G = self.joint_frequencies_combo(*alleles, group2_id=1, normalize=False)
        N = self.num_of_hap_in_ref_subpanel[1] #Divide values by N to normalize the joint frequencies, G.

        MONOSOMY = 1 << len(alleles) - 1 ### The last read is the tuple is from the monosomy.

        ### BPH ###
        (((A0, A1),((B0,),)),) = model[1].items()
        DISOMY = ( F[B0] / M + G[B0] / N ) * A0 / ( 2 * A1 )

        DISOMY += sum( sum( (F[B0] * G[B1] + G[B0] * F[B1]) for (B0, B1) in C) * A0 / A1
                   for (A0, A1), C in model[2].items()) / (2 * M * N)

        BPH = (F[MONOSOMY] / M + G[MONOSOMY] / N) * DISOMY / 2
        
        ### SPH ###

        (((A0, A1),((B0,),)),) = model[1].items()
        SPH = ( F[B0|MONOSOMY] / M + G[B0] * F[MONOSOMY] / N / M ) * A0 / ( 4 * A1 )
        
        SPH += ( F[B0] * G[MONOSOMY] / M / N + G[B0|MONOSOMY] / N ) * A0 / ( 4 * A1 )


        SPH += sum( sum( (F[B0|MONOSOMY] * G[B1] + G[B0] * F[B1|MONOSOMY]) for (B0, B1) in C) * A0 / A1
                   for (A0, A1), C in model[2].items()) / (4 * M * N)
        
        SPH += sum( sum( (F[B0] * G[B1|MONOSOMY] + G[B0|MONOSOMY] * F[B1]) for (B0, B1) in C) * A0 / A1
                   for (A0, A1), C in model[2].items()) / (4 * M * N)
        
        return BPH, SPH

    def likelihoods2(self, *alleles):
        """ Calculates the likelihood to observe two alleles/haplotypes
        under SPH and BPH scenarios. """

        F = self.joint_frequencies_combo(*alleles, group2_id=0, normalize=True)
        G = self.joint_frequencies_combo(*alleles, group2_id=1, normalize=True)
        a, b, ab = F[1], F[2], F[3]
        A, B, AB = G[1], G[2], G[3]
        BPH = (A+a)*(B+b)/4 
        SPH = (ab+A*b+a*B+AB)/4
        
        return BPH, SPH
    
    def likelihoods3(self, *alleles):
        """ Calculates the likelihood to observe three alleles/haplotypes
        under SPH and BPH scenarios. """

        F = self.joint_frequencies_combo(*alleles, group2_id=0, normalize=True)
        G = self.joint_frequencies_combo(*alleles, group2_id=1, normalize=True)
        
        a, b, ab, c, ac, bc, abc = F[1], F[2], F[3], F[4], F[5], F[6], F[7]
        A, B, AB, C, AC, BC, ABC = G[1], G[2], G[3], G[4], G[5], G[6], G[7]
        
        BPH = (c+C)*(ab+A*b+AB+a*B)/8
        
        SPH = (abc+A*bc+AB*c+ac*B)/8
        SPH += (ab*C+AC*b+ABC+a*BC)/8

        return BPH, SPH

    def likelihoods4(self, *alleles):
        """ Calculates the likelihood to observe four alleles/haplotypes
        under SPH and BPH scenarios. """

        F = self.joint_frequencies_combo(*alleles, group2_id=0, normalize=True)
        G = self.joint_frequencies_combo(*alleles, group2_id=1, normalize=True)
        
        a, b, c, d = F[1], F[2], F[4], F[8],
        ab, ac, ad, bc, bd, cd = F[3], F[5], F[9], F[6], F[10], F[12]
        abc, abd, acd, bcd = F[7], F[11], F[13], F[14]
        abcd = F[15]

        A, B, C, D = G[1], G[2], G[4], G[8],
        AB, AC, AD, BC, BD, CD = G[3], G[5], G[9], G[6], G[10], G[12]
        ABC, ABD, ACD, BCD = G[7], G[11], G[13], G[14]
        ABCD = G[15]

        BPH = (d+D)*(abc+ab*C+ac*B+bc*A+ABC+AB*c+AC*b+BC*a)/16 
        
        SPH = (abcd+abd*C+acd*B+bcd*A+ABC*d+AB*cd+AC*bd+BC*ad)/16 
        SPH += (abc*D+ab*CD+ac*BD+bc*AD+ABCD+ABD*c+ACD*b+BCD*a)/16
        
        return BPH, SPH

    def likelihoods5(self, *alleles):
        """ Calculates the likelihood to observe five alleles/haplotypes
        under SPH and BPH scenarios. """

        F = self.joint_frequencies_combo(*alleles, group2_id=0, normalize=True)
        G = self.joint_frequencies_combo(*alleles, group2_id=1, normalize=True)
        a, b, c, d, e = F[1], F[2], F[4], F[8], F[16]
        ab, ac, ad, ae, bc, bd, be, cd, ce, de = F[3], F[5], F[9], F[17], F[6], F[10], F[18], F[12], F[20], F[24]
        abc, abd, abe, acd, ace, ade, bcd, bce, bde, cde = F[7], F[11], F[19], F[13], F[21], F[25], F[14], F[22], F[26], F[28]
        abcd, abce, abde, acde, bcde = F[15], F[23], F[27], F[29], F[30]
        abcde = F[31]

        A, B, C, D, E = G[1], G[2], G[4], G[8], G[16]
        AB, AC, AD, AE, BC, BD, BE, CD, CE, DE = G[3], G[5], G[9], G[17], G[6], G[10], G[18], G[12], G[20], G[24]
        ABC, ABD, ABE, ACD, ACE, ADE, BCD, BCE, BDE, CDE = G[7], G[11], G[19], G[13], G[21], G[25], G[14], G[22], G[26], G[28]
        ABCD, ABCE, ABDE, ACDE, BCDE = G[15], G[23], G[27], G[29], G[30]
        ABCDE = G[31]

        BPH = (e+E)*(abcd+abc*D+bcd*A+acd*B+abd*C+ab*CD+ad*BC+ac*BD+ABCD+ABC*d+BCD*a+ACD*b+ABD*c+AB*cd+AD*bc+AC*bd)/32 
        
        SPH = (abcde+abce*D+bcde*A+acde*B+abde*C+abe*CD+ade*BC+ace*BD+ABCD*e+ABC*de+BCD*ae+ACD*be+ABD*ce+AB*cde+AD*bce+AC*bde)/32 
        SPH += (abcd*E+abc*DE+bcd*AE+acd*BE+abd*CE+ab*CDE+ad*BCE+ac*BDE+ABCDE+ABCE*d+BCDE*a+ACDE*b+ABDE*c+ABE*cd+ADE*bc+ACE*bd)/32 

        return BPH, SPH

    def get_likelihoods(self, *x):
        """ Uses the optimal function to calculate the likelihoods.
        In general, self.likelihoods can get less than five alleles but the
        dedicated functions are optimized to a certain number of alleles. """

        l = len(x)
        if l==2:
            result = self.likelihoods2(*x)
        elif l==3:
            result = self.likelihoods3(*x)
        elif l==4:
            result = self.likelihoods4(*x)
        elif l==5:
            result = self.likelihoods5(*x)
        else:
            result = self.likelihoods(*x)
        return result

def wrapper_of_recent_admixture_for_debugging(obs_filename,leg_filename,hap_filename,sample_filename,models_filename):
    """ Wrapper function of the class 'recent_admixture'. It receives an observations
    file, legend file, haplotypes file, samples file and a file with the
    statistical models. Based on the given data it creates and returns an
    instance of the class. """

    if not os.path.isfile(obs_filename): raise Exception('Error: OBS file does not exist.')
    if not os.path.isfile(leg_filename): raise Exception('Error: LEGEND file does not exist.')
    if not os.path.isfile(hap_filename): raise Exception('Error: HAP file does not exist.')
    if not os.path.isfile(sample_filename): raise Exception('Error: SAMPLE file does not exist.')
    if not os.path.isfile(models_filename): raise Exception('Error: MODELS file does not exist.')

    load = lambda filename: {'bz2': bz2.open, 'gz': gzip.open}.get(filename.rsplit('.',1)[1], open)  #Adjusts the opening method according to the file extension.

    open_hap = load(hap_filename)
    with open_hap(hap_filename,'rb') as hap_in:
        hap_tab, total_number_of_haplotypes = pickle.load(hap_in)

    open_leg = load(leg_filename)
    with open_leg(leg_filename,'rb') as leg_in:
        leg_tab = pickle.load(leg_in)


    open_samp = load(sample_filename)
    with open_samp(sample_filename,'rb') as samp_in:
        sam_tab = pickle.load(samp_in)

    open_obs = load(obs_filename)
    with open_obs(obs_filename, 'rb') as obs_in:
        obs_tab = pickle.load(obs_in)
        #info = pickle.load(obs_in)

    open_model = load(models_filename)
    with open_model(models_filename, 'rb') as model_in:
        models_dict = pickle.load(model_in)

    return recent_admixture(obs_tab, leg_tab, hap_tab, sam_tab, models_dict, total_number_of_haplotypes)

if __name__ != "__main__":
    print('The module RECENT_ADMIXTURE_MODELS was imported.')
else:
    print('The module RECENT_ADMIXTURE_MODELS was invoked directly.')
    sys.exit(0)

###############################   END OF FILE   ###############################

"""

if __name__ != "__main__":
    print("The module RECENT_ADMIXTURE_MODELS was imported.")
else:
    print("The module RECENT_ADMIXTURE_MODELS was invoked directly.")
    #sys.exit(0)
    import time, random
    t0 = time.time()
    obs_filename = 'test/test.obs.p'
    hap_filename = 'test/test.hap.p'
    leg_filename = 'test/test.leg.p'
    sam_filename = 'test/test.sam.p'
    models_filename = 'MODELS/MODELS18.p'

    A = wrapper_of_recent_admixture_for_debugging(obs_filename,leg_filename,hap_filename,sam_filename,models_filename)

    alleles = tuple(A.hap_dict.keys())

    frequencies0 = lambda *x: {bin(a)[2:]:b for a,b in A.joint_frequencies_combo(*x,group2_id=0,normalize=True).items()}
    frequencies1 = lambda *x: {bin(a)[2:]:b for a,b in A.joint_frequencies_combo(*x,group2_id=1,normalize=True).items()}

    likelihoods = A.likelihoods
    likelihoods2 = A.likelihoods2
    likelihoods3 = A.likelihoods3
    likelihoods4 = A.likelihoods4
    likelihoods5 = A.likelihoods5


    random.seed(a=0, version=2)
    x = random.randrange(len(alleles)-16)
    haplotypes = (alleles[x:x+4],alleles[x+4:x+8],alleles[x+8:x+12],alleles[x+12:x+16])

    print('-----joint_frequencies_combo-----')
    print(frequencies0(alleles[x+0]))
    print(frequencies0(*alleles[x:x+4]))
    print(frequencies1(alleles[x+0]))
    print(frequencies1(*alleles[x:x+4]))
    print('-----likelihoods4-haplotypes-----')
    print(haplotypes)
    print(frequencies0(*haplotypes))
    print(frequencies1(*haplotypes))
    print(likelihoods(*haplotypes))
    print(likelihoods4(*haplotypes))
    print('-----likelihoods2-----')
    print(alleles[x:x+2])
    print(frequencies0(*alleles[x:x+2]))
    print(frequencies1(*alleles[x:x+2]))
    print(likelihoods(*alleles[x:x+2]))
    print(likelihoods2(*alleles[x:x+2]))
    print('-----likelihoods3-----')
    print(alleles[x:x+3])
    print(frequencies0(*alleles[x:x+3]))
    print(frequencies1(*alleles[x:x+3]))
    print(likelihoods(*alleles[x:x+3]))
    print(likelihoods3(*alleles[x:x+3]))
    print('-----likelihoods4-----')
    print(alleles[x:x+4])
    print(frequencies0(*alleles[x:x+4]))
    print(frequencies1(*alleles[x:x+4]))
    print(likelihoods(*alleles[x:x+4]))
    print(likelihoods4(*alleles[x:x+4]))
    print('-----likelihoods5-----')
    print(alleles[x:x+5])
    print(frequencies0(*alleles[x:x+5]))
    print(frequencies1(*alleles[x:x+5]))
    print(likelihoods(*alleles[x:x+5]))
    print(likelihoods5(*alleles[x:x+5]))
    t1 = time.time()

    print('Done in %.3f sec.' % ((t1-t0)))
"""

