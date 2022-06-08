#!/usr/bin/env python
# coding: utf-8

# In[6]:


# import sys
# sys.path.append('../kmerapp/')
# sys.path.append('../FVB_B6_Graph/package/')
# sys.path.append('./kmerapp/package/')
# sys.path.append('./package/')
#import CCGG_extension_C as CCGG
import gzip
import Bio
from Bio import pairwise2
import json
import pandas as pd
import argparse
import Levenshtein
#import Levenshtein





def load_anchor_info(CC019):
    df_019 = pd.read_csv(CC019, index_col = None)

    dict_019 = {}
    chrlist = range(1,21)

    for chromo in chrlist:
        data = df_019[df_019['Chromosome'] == chromo]
        dict_019[chromo] = dict(data[['Anchor', 'Pos']].values)
    return dict_019





def makeCigar(seq, ref):
    if (len(seq) > 16384) or (len(ref) > 16384):
        rmid = len(ref)/2        
        smid = len(seq)/2
        prox = makeCigar(seq[:smid],ref[:rmid])
        dist = makeCigar(seq[smid:],ref[rmid:])
        return prox+dist
    ops = Levenshtein.editops(seq,ref)
    code = ['=' for i in range(len(seq))]
    offset = 0
    for op, si, di in ops:
        if (op == "replace"):
            code[si+offset] = 'X'
        elif (op == "insert"):
            code.insert(si+offset,'D')
            offset += 1
        elif (op == "delete"):
            code[si+offset] = 'I'# LM: fixed bug here 2019-04-15
    cigar = ''
    count = 1
    prev = code[0]
    for c in code[1:]:
        if (c == prev):
            count += 1
        else:
            cigar += "%d%c" % (count, prev)
            count = 1
            prev = c
    cigar += "%d%c" % (count, prev)
    return cigar

# mask Ns on the alternative path
def maskNs(seq):
    seq = seq.upper()
    newseq = ''
    for i in seq:
        if i == "N":
            newseq += 'n'
        else:
            newseq += i
    return newseq

def gen_cigar(ref, qry):
    if len(ref) != len(qry):
        raise Exception('unequal length')
    cigar = []
    for i in range(len(ref)):
        r, q = ref[i], qry[i];
        if r is '-' and q is '-':
            raise Exception('both gaps')
        op = '=' if r is q else 'I' if r is '-' else 'D' if q is '-' else 'X';
        if len(cigar) > 0 and cigar[-1][1] is op: # add to the last operation
            cigar[-1][0] += 1
        else: cigar.append([1, op]);              # a new operation
    return "".join(map(lambda x: str(x[0]) + x[1], cigar)); # turn to string

def GapOpenAligner(reference,sequence):
    alignments = pairwise2.align.globalms(sequence,reference,0, -6, -5, -3)
    #print(len(sequence), len(reference))
    alignment = alignments[0]
    seq, ref = alignment[0],alignment[1]
    cigar = gen_cigar(ref,seq)
    return cigar


# In[19]:


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

def get_ref_kmers(seq, k):
    kmerlist = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        kmerlist[kmer] = kmerlist.get(kmer, []) + [i]
    return kmerlist

def get_alt_kmer(alt_seq, k):
    alt_kmer = {}
    #alt_seq = seq_37
    for i in range(len(alt_seq) - k + 1):
        kmer = alt_seq[i:i+k]
        if kmer.count('N')>0:
            continue

        # uniqueness
        ref_indexlist = kmerlist.get(kmer, "")    
        if len(ref_indexlist) != 1:
            continue

        alt_kmer[kmer] = alt_kmer.get(kmer, []) + [i]

    alt_kmer = sorted(alt_kmer.items(), key=lambda item: item[1])    
    return alt_kmer

def get_indextable(alt_kmer, kmerlist):
    indextable = []
    for kmer, index in alt_kmer:
        assert len(index) == 1 and len(kmerlist[kmer]) == 1
        indextable.append([kmer, kmerlist[kmer][0], index[0]]) # kmer, reference_index, alternative_index
    indextable = numpy.array(indextable)
    # test monotonicity
    positionlist = indextable[:,1:].astype(int)
    assert numpy.all((positionlist[1:, :] - positionlist[:-1, :])>0)
    
    return positionlist

def find_blocks(Indexlist, k):
    count = 0
    Blocks = []
    for i,s in enumerate(Indexlist):
        if i == 0:
            if Indexlist[i+1] - Indexlist[i] <= k:
                count += 1
                start = i
        elif i > 0 and i < len(Indexlist)-1:
            if Indexlist[i] - Indexlist[i-1] > k and Indexlist[i+1] - Indexlist[i] <= k:
                count +=1
                start = i
            elif Indexlist[i] - Indexlist[i-1] <= k and Indexlist[i+1] - Indexlist[i] > k:
                end = i+1
                Blocks.append((start, end))
        else:
            if Indexlist[i] - Indexlist[i-1] <= k:
                end = i+1
                Blocks.append((start, end))
                
    return count, Blocks


def sort_by_kmerlength(poslist, k, index):
    indexlist = numpy.array(poslist[:,index]).astype(int)
    errorindex = []
    count, blocks = find_blocks(indexlist, k)
    for s,e in blocks:
        if e - s <3:
            for i in range(s+1,e):
                errorindex.append(i)
        else:
            for i in range(s+1,e-1):
                errorindex.append(i)
    mask = numpy.ones(len(indexlist), dtype=bool)
    mask[errorindex] = False
    result = poslist[mask,:]
    return result


# mask Ns on the alternative path
def maskNs(seq):
    seq = seq.upper()
    newseq = ''
    for i in seq:
        if i == "N":
            newseq += 'n'
        else:
            newseq += i
    return newseq

def processCigar(cigar):
    """Helper Function, may not be used directly, expand Cigar string

    Parameters:
        cigar: <str> - compressed cigar
    """
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
    """Helper Function, may not be used directly, compress Cigar string

    Parameters:
        cigar: <str> - expanded cigar
    """
    cigar = cigar +'$'
    out = ''
    N = 0
    start = 0
    for i in range(1,len(cigar)):
        if cigar[i-1] != cigar[i]:
            out += str(i-start) + cigar[i-1]
            start = i
    return out  


def sort_table(positionlist):
    a,v = numpy.where((positionlist[1:,:] - positionlist[:-1, :])<45)
    while len(a)>0 and len(v) > 0:
        positionlist = sort_by_kmerlength(positionlist, 45, 0)
        positionlist = sort_by_kmerlength(positionlist, 45, 1)
        a,v = numpy.where((positionlist[1:,:] - positionlist[:-1, :])<45)
    return positionlist

def get_alignment(p_dict, seq_19, seq_37):
    '''get alignment for gap over 2kb'''
    
    ref_poslist = sorted(p_dict.keys())

    ref_start = 0
    alt_start = 0

    cigar = ''
    for s in ref_poslist:
        alt = p_dict[s]

        if (s == ref_start) or (alt_start == alt):
            continue
        ref_end = s
        alt_end = alt

        assert ref_end > ref_start and alt_end > alt_start

        ref_seq = seq_19[ref_start:ref_end]
        alt_seq = maskNs(seq_37[alt_start:alt_end])


        #cigar = GapOpenAligner(ref_seq,alt_seq)
        c = makeCigar(alt_seq, ref_seq)    
        var = split(c)

        assert valid_size(var,alt_seq,ref_seq)
        assert valid_eqaul_and_mismatch(var, alt_seq, ref_seq)

        cigar += c

        ref_start = ref_end
        alt_start = alt_end


    ref_seq = seq_19[ref_start:]
    alt_seq = maskNs(seq_37[alt_start:])
    cigar += makeCigar(alt_seq, ref_seq)


    cigar = processCigar(cigar)
    cigar = combineCigar(cigar)

    var = split(cigar)
    #print var
    assert valid_size(var,seq_37,seq_19)
    assert valid_eqaul_and_mismatch(var, seq_37, seq_19)

    return cigar


# functions for validating cigars
def split(variant):
    splitVariants = []
    previous = 0
    for i in range(len(variant)):
        v = variant[i]
        if v.isdigit():
            continue
        else:
            splitVariants += [variant[previous:i]]
            splitVariants += [v]
            previous = i + 1
            
    numbers, types = [],[]
    for i in range(len(splitVariants)):
        v = splitVariants[i]
        if i %2 == 0:
            numbers += [v]
        else:
            types += [v]
            
    variant = []
    for i in range(len(numbers)):
        variant += [(numbers[i],types[i])]
    return variant

def search_valid_size(variants,seq):
    ref_length,alt_length = 0,0
    for v in variants:
        base, op = v
        base = int(base)
        if op == '=' or op == 'X':
            ref_length += base
            alt_length += base
        elif op == 'I':
            alt_length += base
        elif op == 'D':
            ref_length += base

    return len(seq) == alt_length 


# check the length of the cigar string is valid
def valid_size(variants,seq,reference):
    ref_length,alt_length = 0,0
    for v in variants:
        base, op = v
        base = int(base)
        if op == '=' or op == 'X':
            ref_length += base
            alt_length += base
        elif op == 'I':
            alt_length += base
        elif op == 'D':
            ref_length += base

    return len(seq) == alt_length and len(reference) == ref_length

def valid_eqaul_and_mismatch(variants,seq,reference):
    ref_pos,alt_pos = 0,0
    for v in variants:
        base,op = v
        base = int(base)
        if op == '=':
            if seq[alt_pos:alt_pos+base] != reference[ref_pos:ref_pos+base]:
                print seq[alt_pos:alt_pos+base], reference[ref_pos:ref_pos+base]
                return False
            ref_pos += base
            alt_pos += base
        elif op == 'X':
            if seq[alt_pos:alt_pos+base].count("N")>0 or reference[ref_pos:ref_pos+base].count("N")>0:
                pass
            elif seq[alt_pos:alt_pos+base] == reference[ref_pos:ref_pos+base]:
                print seq[alt_pos:alt_pos+base], reference[ref_pos:ref_pos+base]
                return False
            ref_pos += base
            alt_pos += base
        elif op == 'I':
            alt_pos += base
        elif op == 'D':
            ref_pos += base
    else:
        return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--Chromosome', help='chromosome used to do search',required=True)
    parser.add_argument("-a", "--Anchors", type = str, nargs = 2,
                        metavar = ('anchorfile1','anchorfile2'), help = "anchor_info file for genome A and genome B", required = True)
    parser.add_argument("-g", "--Genomes", type = str, nargs = 2,
                        metavar = ('genomefile1','genomefile2'), help = "genome fasta file for genome A and genome B", required = True)
    parser.add_argument("-m", "--MaxLength", type = int, nargs = 1,
                        metavar = "maximal length of the interval", help = "an integer for the maximal length between anchor pairs")
    parser.add_argument("-o", "--Output", type = str, nargs = 1,metavar = "output filename", help = "output filename")



    args = parser.parse_args()
    Chromosome = args.Chromosome
    CC019 = args.Anchors[0]
    CC037 = args.Anchors[1]
    print(Chromosome, CC019, CC037)

    dict_019 = load_anchor_info(CC019)
    dict_037 = load_anchor_info(CC037)


    # ref 019, alt 037
    CC019genome = args.Genomes[0]
    CC037genome = args.Genomes[1]
    headers_019, sequences_019 = loadFasta(CC019genome)
    headers_037, sequences_037 = loadFasta(CC037genome)


    # In[ ]:


    Cigars = {}
    maxlength = args.MaxLength[0]
    unmappable = 0

    chromo = int(Chromosome)

    keys = sorted(dict_019[chromo].keys())
    keys.remove("SOURCE")
    keys = ['SOURCE'] + keys
    genome_019 = "+" + sequences_019[chromo-1].upper()
    genome_037 = "+" + sequences_037[chromo-1].upper()

    Cigars[chromo] = {}

    for i in range(len(keys)-1):
        sanchor = keys[i]
        eanchor = keys[i+1]
        spos_019 = dict_019[chromo][sanchor] + 45
        epos_019 = dict_019[chromo][eanchor]
        seq_019 = genome_019[spos_019:epos_019]

        spos_037 = dict_037[chromo][sanchor] + 45
        epos_037 = dict_037[chromo][eanchor]
        seq_037 = genome_037[spos_037:epos_037]

        alt_seq = maskNs(seq_037) # mask Ns
        ref_seq = seq_019

        if len(ref_seq) > maxlength or len(alt_seq) > maxlength:
            # kmer profiling
            k = 45
            kmerlist = get_ref_kmers(ref_seq, k)
            alt_kmer = get_alt_kmer(alt_seq, k)
            try:
                positionlist = get_indextable(alt_kmer, kmerlist)
                passed += 1
            except:
                Cigars[chromo][sanchor] = '?'
                #print anchor, passed
                continue

            if len(positionlist) < 3:
                Cigars[chromo][sanchor] = '?'
                continue

            positionlist = sort_table(positionlist)
            p_dict = dict(positionlist)
            try:
                cigar = get_alignment(p_dict, seq_19, seq_37)
                Cigars[chromo][sanchor] = cigar
            except:
                Cigars[chromo][sanchor] = '?'
                print "wrong alignment", anchor
        else:
            if len(seq_019) == len(seq_037):
                if seq_019 == seq_037:
                    Cigars[chromo][sanchor] = str(len(seq_019)) + '='
                else:
                    cigar = makeCigar(alt_seq,ref_seq)
                    Cigars[chromo][sanchor] = cigar
            else:
                cigar = GapOpenAligner(ref_seq, alt_seq)
                Cigars[chromo][sanchor] = cigar

    outputfilename = args.Output[0]        
    with open(outputfilename + 'pairwisealignment_%s.json' % str(chromo), 'w') as outfile:
        json.dump(Cigars, outfile)




if __name__ == "__main__":
    # calling the main function
    main()


