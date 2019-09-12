from pyfaidx import Fasta

def read_fasta_file(filename, duplicate_action='stop'):
    """
    Function that reads in a FASTA file and returns a dictionary. Wrapper around
    pyfaidx FASTA reading code (https://pypi.org/project/pyfaidx/)

    Parameters
    ----------
    filename : string
        FASTA filename that is to be read in

    duplicate_action : string, optional
        selector that defines how to deal with duplicate entries. Valid options
        are "stop", "first", "last", "longest", "shortest". The default value
        is "stop"
           

    Returns
    -------

    dict
        Returns a dictionary where the keys are the FASTA headers for each sequence
        and the values are the sequences themselves.
    
    """
    

    IN = Fasta(filename, duplicate_action=duplicate_action, read_long_names=True)
    ID2seq={}
    KEYZ = IN.keys()

    for k in KEYZ:
        ID2seq[k] = str(IN[k][:])

    return ID2seq


## ------------------------------------
##
##        Script starts here
## 
## ------------------------------------

# read in human proteome as defined from swissprot
human_proteome = read_fasta_file('swissprot_human_proteome.fasta')

# define parameters for our search:
#
# - fragsize is the size of the region exmained
# - RS_threshold is the fraction of the fragsize which must be R+S residues
# - ratio_threshold defines the max that R:S or S:R can be 

fragsize = 60 
RS_threshold = 0.6
ratio_threshold = 2

# initialize an empty gene set that we're going to populate with entries that contain at
# least 1 RS fragment
geneset = set([])

for i in human_proteome:

    # get the sequence, but if the sequence is shorter than the fragsize then skip
    s = human_proteome[i]
    if len(s) < fragsize:
        continue

    # for each fragment inside the sequence extract out the fragment, calculate the 
    # fraction of R and S, and if it conforms to the definitions we've defined print
    # to the screen and add this gene to the geneset
    for p in range(0, len(s)-fragsize):
        fragment = s[p:p+fragsize]

        fraction_R =  fragment.count('R')/fragsize
        fraction_S =  fragment.count('S')/fragsize

        if fraction_R + fraction_S > RS_threshold:

            try:
                ratio1 = fraction_R/fraction_S
                ratio2 = fraction_S/fraction_R
            except ZeroDivisionError:
                continue

            if ratio1 < ratio_threshold and ratio2 < ratio_threshold:
                print("%s, %i, %i, %s" %(i, p, p+fragsize, fragment))
                geneset.add(i)

# Finally having scanned all the sequences write out the full fasta header
# for the proteins that matched the criteria
with open('SR_proteins.txt','w') as fh:
    gene_list = list(geneset)
    for i in gene_list:
        fh.write('%s\n'%(i))
        

# finally write out JUST the uniprot IDs (not this assumes the '|' character
# separates fields in the FASTA header!
with open('SR_proteins_IDONLY.txt','w') as fh:
    gene_list = list(geneset)
    for i in gene_list:
        fh.write('%s\n'%(i.split('|')[1]))
        
        
