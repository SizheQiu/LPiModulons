import numpy as np
import pandas as pd
import pickle
from scipy import special, stats
from statsmodels.stats.multitest import fdrcorrection
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
            RF = ((line.split(' - ')[1]).split(':')[0]).strip()
            regulon[ RF ] = []
        else:
            g = str(line.split('\t')[1]).strip()
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




def get_gene_region(gene_id, gb_features, ref_features):
    
    left_gb,left_ref = np.inf,np.inf
    right_gb,right_ref = -np.inf, -np.inf
    if gene_id in list(gb_features['locus_tag']):
        left_gb, right_gb, strand = gb_features[gb_features['locus_tag']==gene_id].values[0][0:3]
    if gene_id in list(ref_features['locus_tag']):
        left_ref, right_ref, strand = ref_features[ref_features['locus_tag']==gene_id].values[0][0:3]
        
    return (min(left_gb,left_ref), max(right_gb,right_ref), strand)
    
    

def get_upstream(seq_path, site, N_up, strand):
    sequence = SeqIO.read(seq_path, "fasta").seq
    site = site - 1
    if strand == '+':
        left = site - N_up
        right = site
        s = str(sequence)[left:right]
    else:
        left = site+1
        right = site + N_up+1
        s = complement(revstrand( str(sequence)[left:right] ))
    return s
        
def check_overlap(left,right, gb_features, ref_features):
    '''
    Check genes that overlap with the region( left, right ).
    gb_features, ref_features: feature tables form genbank and refseq.
    '''
    output_genes = []
    temp_gb = gb_features[ ( gb_features['start'] < right ) &  ( gb_features['end'] > left )].reset_index()
    temp_ref = ref_features[ ( ref_features['start'] < right ) &  ( ref_features['end'] > left )].reset_index()
    for i in range(len(temp_gb.index)):
        if 'pWCFS' in str(temp_gb['locus_tag'][i]):
            continue
        output_genes.append( temp_gb['locus_tag'][i] )
    for i in range(len(temp_ref.index)):
        if 'pWCFS' in str(temp_ref['locus_tag'][i]):
            continue
        output_genes.append( temp_ref['locus_tag'][i] )
    output_genes = list(set(output_genes))
    return output_genes
        

# def find_operon_start( site, strand, gb_features, ref_features ):
#     '''
#     find the first gene in the operon,
#     assuming that intergenetic space is < 200bp in the same operon.
#     '''
#     N_up = 200
#     if strand == '+':
#         overlap = check_overlap(site-N_up, site, gb_features, ref_features)
#     else:
#         overlap = check_overlap(site, site+N_up, gb_features, ref_features)
#     while len(overlap) > 0:
#         gene_id = overlap[0]
#         left, right, strand = get_gene_region(gene_id, gb_features, ref_features)
#         if strand == '+':
#             overlap = check_overlap(left-N_up, left, gb_features, ref_features)
#         else:
#             overlap = check_overlap(right, right+N_up, gb_features, ref_features)
            
#     return (gene_id, left, right, strand)
    
    
    
def get_im_promoters(im, seq_path, gb_features, ref_features):
    '''
    Get promoter sequences (200bp upstream to translation start site) 
    for genes in the I-modulon.
    '''
    promoters = {}
    N_up = 200
    for gene in im:
        left, right, strand = get_gene_region(gene, gb_features, ref_features)
        
        if strand == '+':
            overlap = check_overlap(left-N_up, left, gb_features, ref_features)
        else:
            overlap = check_overlap(right, right+N_up, gb_features, ref_features)
        
        if len(overlap) < 1:
            if strand == '+':
                promoters[gene] = get_upstream(seq_path, left, N_up, strand)
            else:
                promoters[gene] = get_upstream(seq_path, right, N_up, strand)
        elif overlap[0] in im:
            continue
#         else:
#             if strand == '+':
#                 gene_id, start_left, start_right, strand = find_operon_start( left, strand, gb_features, ref_features )
#                 promoters[gene_id] = get_upstream(seq_path, start_left, N_up, strand)
#             else:
#                 gene_id, start_left, start_right, strand = find_operon_start( right, strand, gb_features, ref_features )
#                 promoters[gene_id] = get_upstream(seq_path, start_right, N_up, strand)
                
    return promoters
             




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
    
    
def get_seqfeature( gene, table ):
    temp_pd = table[table['locus_tag']==gene].reset_index()
    output = {'start':int(temp_pd['start'][0]),'end':int(temp_pd['end'][0]),
            'strand':str(temp_pd['strand'][0]),'label': str(temp_pd['symbol'][0])}
    if str(temp_pd['strand'][0]) == '+':
        output['strand'] = +1
    else:
        output['strand'] = -1
    if str(temp_pd['symbol'][0]) == gene:
        output['label'] = gene
    else:
        output['label'] = str(temp_pd['symbol'][0])+'('+gene+')'
    return output




