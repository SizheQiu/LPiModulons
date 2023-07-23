import numpy as np
import pandas as pd
import pickle
from scipy import special, stats
from statsmodels.stats.multitest import fdrcorrection



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
    '''
    Convert read counts to RPKM.
    '''
    return  numReads / ( gene_length/1000 * totalNumReads/1000000 )


def revstrand(inseq):
    outseq=inseq[::-1]
    return outseq

def complement(inseq):
    '''
    complement for both complete and incomplete nucleotides.
    '''
    inseq.upper()
    complement_dict = {'A':'T','T':'A','C':'G','G':'C',
                       'B':'V','D':'H','H':'D','K':'M','M':'K',
                       'S':'S','V':'B','W':'W','N':'N','R':'Y','Y':'R' }
    clist=[]
    for ntide in inseq:
        clist.append( complement_dict[ ntide ] )
        
    complement=''.join(clist)
    return complement



def compute_threshold(S,k,cutoff=550):
    """Computes kurtosis-based threshold for a component of an S matrix
        S: Component matrix with gene weights
        k: Component name
        cutoff: Minimum test statistic value to determine threshold (550 is default from sensitivity analysis)
    """
    i = 0 
    # Sort genes based on absolute value
    ordered_genes = abs(S[k]).sort_values()
    K,p = stats.normaltest(S.loc[:,k])
    while K > cutoff:
        i -= 1
        # Check if K statistic is below cutoff
        K,p = stats.normaltest(S.loc[ordered_genes.index[:i],k])
    comp_genes = ordered_genes.iloc[i:]
    if len(comp_genes) == len(S.index):
        return max(comp_genes)+.05
    else:
        return np.mean([ordered_genes.iloc[i],ordered_genes.iloc[i-1]])




