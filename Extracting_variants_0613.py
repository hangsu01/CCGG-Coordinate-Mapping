#!/usr/bin/env python
# coding: utf-8

# In[1]:


import json
import pandas as pd
import numpy
import gzip
import argparse


def loadFasta(filename):
    """ Parses a classically formatted and possibly 
        compressed FASTA file into a list of headers 
        and fragment sequences for each sequence contained.
        The resulting sequences are 0-indexed! """
    if (filename.endswith(".gz")):
        fp = gzip.open(filename, 'rb')
    else:
        fp = open(filename, 'rb')
    # split at headers
    data = fp.read().split(b'>')
    fp.close()
    # ignore whatever appears before the 1st header
    data.pop(0)     
    headers = []
    sequences = []
    for sequence in data:
        sequence = sequence.decode()
        lines = sequence.split('\n')
        headers.append(lines.pop(0))
        sequences.append(''.join(lines))
    return (headers, sequences)

def processCigar(cigar):
    '''expand cigar string
    def processCigar(cigar:<str>):
    Parameters:
        cigar:<str> - Samtool Compatible Cigar String. e.g. 2=
    Output:
        out:<str> - expanded cigar string. e.g.:==
    '''

    out = ''
    N = 0
    for symbol in cigar:
        if symbol in '0123456789':
            N = 10*N + int(symbol)
        else:
            #if (symbol != 'D'):
            if (N == 0):
                out += symbol
            else:
                out += N*symbol
            N = 0
    return out

def combineCigar(cigar):
    '''compress expanded cigar string
    def combineCigar(cigar:<str>):
    Parameters:
        cigar:<str> - Expanded cigar strings. e.g. ==
    Output:
        out:<str> - Samtool Compatible Cigar String. e.g.:2=
    '''
    cigar = cigar +'$'
    out = ''
    N = 0
    start = 0
    for i in range(1,len(cigar)):
        if cigar[i-1] != cigar[i]:
            out += str(i-start) + cigar[i-1]
            start = i
    return out 



def report_variants(cigar, ref_i, alt_i, int_i):
    # ref:019, alt:037
    '''Given a compressed cigar string between two sequences, return the variants between the two sequences.
    def report_variants(cigar <str>, ref_i <int>, alt_i <int>):
    Parameters:
        cigar: <str> - Samtool Compatible Cigar String. e.g.:2=
        ref_i: <int> - integer, linear coordinates of the first base of the proximal anchor on Genome A (as ref)
        alt_i: <int> - integer, linear coordinates of the first base of the proximal anchor on Genome B (as alt)
    Output:
        variants:<list> - a list of variants, mismatches Xs, insertions I, and deletions D, 
        (1) linear coord on ref, 
        (2) offset on the intermediate string, 
        (3) linear coord on alt)
    '''
    variants = []
    cigar = processCigar(cigar)
    for i, s in enumerate(cigar): 
        if s == '=':
            alt_i += 1
            ref_i += 1

        if s == 'I':
            variants.append((numpy.nan, int_i+i, alt_i, s))
            alt_i += 1
        if s == 'D':
            variants.append((ref_i, int_i+i, numpy.nan, s))
            ref_i += 1
        if s == 'X':
            variants.append((ref_i, int_i+i, alt_i, s))
            alt_i += 1
            ref_i += 1
            
    return variants




# In[22]:


# run

# In[3]:

def load_alignment(args):
    '''Load pairwise alignment results between the two genome
    
    def load_alignment(args):
    Parameters:
        args: <input parser>
    Output:
        anchor_dict <dict>: 
        keys: anchor name
        values: cigar strings
    '''
    chromo = args.Chromosome[0]
    alignmentdirectory = str(args.Alignment[0])
    print(alignmentdirectory)
    fp = alignmentdirectory
    with open(fp % str(chromo), 'r') as fp:
        alignment = json.load(fp)
    anchor_dict = alignment[str(chromo)]
    
    return anchor_dict

def load_alignment_txt(args):
    '''Load pairwise alignment results between the two genome
    
    def load_alignment(args):
    Parameters:
        args: <input parser>
    Output:
        anchor_dict <dict>: 
        keys: anchor name
        values: cigar strings
    '''

    alignmentfile = str(args.Alignment[0])
    chromo = args.Chromosome[0]
    with open(alignmentfile, 'r') as f:
        data = f.read()
        data = data.split('\n')
        data.pop(0)
    alignment = {}
    for item in data:
        itemlist = item.split(', ')
        #print itemlist
        if len(itemlist) == 2:
            alignment[itemlist[0]] = itemlist[1]
    return alignment

def get_anchor_intermediate_pos(anchor_pos_19, anchor_pos_37, alignment):
    
    anchorinfo = []

    pos = 1
    for i, (node, pos_19) in enumerate(anchor_pos_19[:-1]):
        anchorinfo.append((node, pos))
        a = alignment.get(node, "?")
        if a == "?":
            d_019 = anchor_pos_19[i+1, 1] - anchor_pos_19[i, 1]
            d_037 = anchor_pos_37[i+1, 1] - anchor_pos_37[i, 1] # fix error June14， 2022
            d = max(d_019, d_037)
            pos += d
            continue
        cigar = processCigar("45=" + a)
        pos += len(cigar)

    node = anchor_pos_19[i+1, 0]
    assert node == "SINK", "incomplete"
    anchorinfo.append((node, pos))
    anchor_dict = dict(anchorinfo)
    return anchor_dict

def load_anchor(args):
    '''Load pairwise alignment results between the two genome
    
    def load_anchor(args):
    Parameters:
        args: <input parser>
    Output:
        anchor_pos_19 <array>: anchor name, linear coordinates on Genome A
        anchor_pos_37 <array>: anchor name, linear coordinates on Genome B
    '''
    file1 = args.Anchors[0]
    # file to copy upon
    file2 = args.Anchors[1]
    chromo = args.Chromosome[0]
    print(chromo, file1, file2)

    anchor_info_19 = pd.read_csv(file1, index_col = None)
    anchor_info_37 = pd.read_csv(file2, index_col = None)
    
    anchor_info_19 = anchor_info_19[anchor_info_19['Chromosome'] == chromo]
    anchor_info_37 = anchor_info_37[anchor_info_37['Chromosome'] == chromo]

    anchor_pos_19 = anchor_info_19[['Anchor', 'Pos']].values
    anchor_pos_37 = anchor_info_37[['Anchor', 'Pos']].values
    
    return (anchor_pos_19, anchor_pos_37)



def extracting_variants(alignment, anchor_pos_19, anchor_pos_37):
    '''Construct mapping table for the chromosome
    
    def construct_mapping_table(alignment <dict>, anchor_pos_19 <dict>, anchor_pos_37 <dict>):
    Parameters:
        alignment <dict>: keys - anchor name; values - cigar string
        anchor_pos_19 <dict>: keys - anchor name; values - linear coordinates on genome A
        anchor_pos_37 <dict>:keys - anchor name; values - linear coordinates on genome B
    Output:
        mapping_table <array>: mapping between two genome and the intermediate string
        (col_1) linear coord on ref, 
        (col_2) offset on the intermediate string, 
        (col_3) linear coord on alt        
    '''
    variants = []
    anchorlist = sorted(alignment.keys())
    anchor_intermediate_dict = get_anchor_intermediate_pos(anchor_pos_19, anchor_pos_37, alignment)
    anchor_pos_19 = dict(anchor_pos_19)
    anchor_pos_37 = dict(anchor_pos_37)
    
    for node in anchorlist:
        ref_coord = anchor_pos_19[node]
        alt_coord = anchor_pos_37[node]
        intermediate_coord = anchor_intermediate_dict[node]

        if alignment.get(node) == "?":
            continue

        cigar = "45=" + alignment[node]
        offsets = report_variants(cigar, ref_coord, alt_coord, intermediate_coord)
        variants += offsets


    variants = numpy.array(variants)
    #mapping_table = mapping_table.astype(int)
    return variants


def get_alleles(data, seq_19, seq_37):
    data = data.astype(str)
    result = []
    for row in data:
        ref, intermediate, alt, t = row
        #print ref, alt, type(alt), alt == 'nan'
        if ref == 'nan':
            ra = ''
        else:
            ra = seq_19[int(ref)]
        if alt == 'nan':
            aa = ''
        else:
            aa = seq_37[int(alt)]
        result.append([ref, intermediate, alt, t, ra, aa])
    return result


def main():
    
    # create parser object
    parser = argparse.ArgumentParser(description = "Extracting Variants")
     
    parser.add_argument("-a", "--Alignment", type = str, nargs = 1,
                        metavar = "path", default = None,
                        help = "alignment file name")
    
  
    parser.add_argument("-c", "--Chromosome", type = int, nargs = 1,
                        metavar = "chromosome", default = None,
                        help = "specify chromosome, 1-20")
     
    parser.add_argument("-p", "--Anchors", type = str, nargs = 2,
                        metavar = ('file1','file2'), help = "anchor_info file for genome A and genome B")
    
    parser.add_argument("-g", "--Genomes", type = str, nargs = 2,
                        metavar = ('genomefile1','genomefile2'), 
                        help = "genome fasta file for genome A and genome B")
    
    parser.add_argument("-o", "--Output", type = str, nargs = 1,
                    metavar = "output filename", help = "output filename")
    
    
    args = parser.parse_args()
    Chromosome = args.Chromosome[0]
    print(Chromosome)
        
    anchor_dict = load_alignment_txt(args)
    
    anchor_pos_19, anchor_pos_37 = load_anchor(args)
    
    CC019genome = args.Genomes[0]
    CC037genome = args.Genomes[1]
    headers_019, sequences_019 = loadFasta(CC019genome)
    headers_037, sequences_037 = loadFasta(CC037genome)


    seq_19 = "+" + sequences_019[int(Chromosome) - 1]
    seq_37 = "+" + sequences_037[int(Chromosome) - 1]

    
    variants = extracting_variants(anchor_dict, anchor_pos_19, anchor_pos_37)
    results = get_alleles(variants, seq_19, seq_37)
    
    outputfile = args.Output[0]
    results = numpy.array(results)
    numpy.save(outputfile, results)
#     with open(outputfile, 'w') as fp:       
#         for pos in variants:
#             fp.write(str(pos) + '\n')

if __name__ == "__main__":
    # calling the main function
    main()
