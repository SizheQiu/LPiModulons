import numpy as np
import pandas as pd
import pickle



def read_regprecise(file):
    regulon = {}
    ffile = open(file, "rt")
    lines = ffile.readlines()
    ffile.close()
    for line in lines:
        line = line.strip()
        if len(line)<1:
            continue
        elif '#' in line:
            RF = ((line.split('-')[1]).split(':')[0]).strip()
            regulon[ RF ] = []
        else:
            g = str(line.split('\t')[2]).strip()
            regulon[RF].append(g)
    return regulon 

def load_pickle(filename):
    temp = None
    with open(filename,'rb') as f:
        temp = pickle.load(f)
    return temp

def dump_pickle(file, filename):
    with open(filename, 'wb') as f:
        pickle.dump( file , f)
        
        
def get_rpkm( numReads, gene_length, totalNumReads ):
    return  numReads / ( gene_length/1000 * totalNumReads/1000000 )