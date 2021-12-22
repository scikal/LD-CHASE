#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAKE_STATISTICAL_MODEL

Builds statistical models of disomy for non-admixtures.

Daniel Ariad (daniel@ariad.org)
AUG 31, 2020
"""

from itertools import product
from collections import defaultdict
from math import gcd as greatest_common_divisor
from pickle import dump
from time import time
from bz2 import BZ2File
import argparse, sys

def ENGINE(number_of_reads,degeneracies):
    """ Generates polysomy statistical models for n-reads, based of a list of
    the degeneracy of each homolog. """

    degeneracies_dict = {i:w for i,w in enumerate(degeneracies) if w>0}
    model = defaultdict(int)
    for sequence in product(degeneracies_dict, repeat=number_of_reads):
        haplotypes = defaultdict(list)
        weight = 1
        for read_ind,hap in enumerate(sequence):
            haplotypes[hap].append(read_ind)
            weight *=  degeneracies_dict[hap]
        key = tuple(tuple(indices) for indices in haplotypes.values())
        model[key] += weight
    return model

def COMPACT(model,number_of_reads,degeneracies):
    """ The partitioning of reads can be encoded efficiently using the
    occupation basis. In this representation all the reads are enumerated.
    Each subset of reads is represented by a binary sequence, where the
    $i^\mathrm{th}$ element is one when the $i^\mathrm{th}$ read is included
    in the subset and zero otherwise. In other words, bits in a binary sequence
    map to whether a read is included in the subset. """

    compact = {i+1: defaultdict(list) for i in range(len(degeneracies))}
    T = sum(degeneracies)**number_of_reads
    while(len(model)!=0):
        haplotypes, weight = model.popitem()
        HAPLOTYPES = tuple(sum(1 << x for x in h) for h in sorted(haplotypes, key=len))
        gcd = greatest_common_divisor(weight,T)
        compact[len(HAPLOTYPES)][weight//gcd,T//gcd].append(HAPLOTYPES)
    for partition_size in compact:
        compact[partition_size] = {normalized_weight:tuple(hap) for normalized_weight,hap in compact[partition_size].items()}
    return compact

def DISOMY(number_of_reads):
    """ Builds a statistical model for n-reads under the disomy scenario. """

    degeneracies = (1, 1)
    model = ENGINE(number_of_reads,degeneracies)
    return COMPACT(model,number_of_reads,degeneracies)

def representationA(model):
    """ Represents the model that is returned from the function ENGINE. """
    result = ''
    for partition,weight in model.items():
        if result!='':
            result += '+'
        l = [''.join((chr(read+65) for read in hap)) for hap in partition if len(hap)]
        sandwitch = ''.join('f('+hap+')' for hap in l)
        if weight!=1:
            result += f'{weight:d}'
        result += f'{sandwitch:s}'
    return result

def representationB(model):
    """ An alternative representation of the model. """
    result=''
    for partition,weight in model.items():
        if result!='':
            result += '+'
        l = [''.join((chr(65+read_ind) for read_ind in hap)) for hap in partition if len(hap)]
        if weight!=1:
            result += f'{weight:d}*'
        result += f"{'*'.join(l):s}"
    return result

def BUILD(x):
    """ Build and store a dictionary with statistical models of disomy for
        various number of reads. """

    models = dict()
    for i in range(1,x+1):
        print('Makes the statistical model for %d reads.' % i)
        a = time()
        models[i]= DISOMY(i)
        b = time()
        print('Done building in %.3f sec.' % ((b-a)))
    with open( f'MODELS{x:d}.p', 'wb') as f:
        dump(models, f, protocol=4)
    with BZ2File( f'MODELS{x:d}.p.bz2', 'wb') as f:
        dump(models, f, protocol=4)
    return models

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Makes the statistical model for disomy.')
    parser.add_argument('maximal_number_of_reads', metavar='n', type=int, default=16,
                        help='Maximal number of supported reads. The default is 16 reads.')

    n = vars(parser.parse_args())['maximal_number_of_reads']
    models = BUILD(n)
    sys.exit(0)

else:
    print('The module MAKE_STATISTICAL_MODEL was imported.')
