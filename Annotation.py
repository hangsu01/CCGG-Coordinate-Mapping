import json
import pandas as pd
import numpy
import gzip
import argparse
import datetime

def load_alignment_txt(alignmentfile):
    '''Load pairwise alignment results between the two genome
    
    def load_alignment(args):
    Parameters:
        args: <input parser>
    Output:
        anchor_dict <dict>: 
        keys: anchor name
        values: cigar strings
    '''
    with open(alignmentfile, 'r') as f:
        data = f.read()
        data = data.split('\n')
        data.pop(0)
    alignment = {}
    for item in data:
        itemlist = item.split(',')
        #print itemlist
        if len(itemlist) == 2:
            alignment[itemlist[0]] = itemlist[1]
    return alignment


def processCigar(cigar):
    '''expand cigar string
    Input: 
    cigar, e.g. 2=
    Output:==
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
    '''compress cigar string'''
    cigar = cigar +'$'
    out = ''
    N = 0
    start = 0
    for i in range(1,len(cigar)):
        if cigar[i-1] != cigar[i]:
            out += str(i-start) + cigar[i-1]
            start = i
    return out 


def offsetexchange(cigar, ref_i, alt_i, int_i):
    '''Given a compressed cigar string between two sequences, return the base pair mapping between the two sequences.
    def offsetexchange(cigar <str>, ref_i <int>, alt_i <int>):
    Parameters:
        cigar: <str> - Samtool Compatible Cigar String. e.g.:2=
        ref_i: <int> - integer, linear coordinates of the first base of the proximal anchor on Genome A (as ref)
        alt_i: <int> - integer, linear coordinates of the first base of the proximal anchor on Genome B (as alt)
    Output:
        offsets:<list> - a list of identical coordinates, 
        (1) linear coord on ref, 
        (2) offset on the intermediate string, 
        (3) linear coord on alt)
    '''
    offsets = []
    cigar = processCigar(cigar)
    for i, s in enumerate(cigar): 
        if s == '=':
            offsets.append((ref_i, int_i+i, alt_i))
            alt_i += 1
            ref_i += 1

        if s == 'I':
            offsets.append((numpy.nan, int_i+i, alt_i))
            alt_i += 1
        if s == 'D':
            offsets.append((ref_i, int_i+i, numpy.nan))
            ref_i += 1
        if s == 'X':
            offsets.append((ref_i, int_i+i, alt_i))
            alt_i += 1
            ref_i += 1
            
    return offsets

def get_anchor_intermediate_pos(anchor_pos_19, anchor_pos_37, alignment):
    
    anchorinfo = []

    pos = 1
    for i, (node, pos_19) in enumerate(anchor_pos_19[:-1]):
        anchorinfo.append((node, pos))
        a = alignment.get(node, "?")
        if a == "?":
            d_019 = anchor_pos_19[i+1, 1] - anchor_pos_19[i, 1]
            d_037 = anchor_pos_37[i+1, 1] - anchor_pos_37[i, 1]
            d = max(d_019, d_037)
            pos += d
            continue
        cigar = "45=" + alignment[node]
        cigar = processCigar(cigar)
        pos += len(cigar)

    node = anchor_pos_19[i+1, 0]
    assert node == "SINK", "incomplete"
    anchorinfo.append((node, pos))
    anchor_dict = dict(anchorinfo)
    return anchor_dict

def construct_mapping_table(alignment, anchor_pos_ref, anchor_pos_19):
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
    mapping_table = []
    anchorlist = sorted(alignment.keys())
    anchor_intermediate_dict = get_anchor_intermediate_pos(anchor_pos_ref, anchor_pos_19, alignment)
    anchor_pos_ref = dict(anchor_pos_ref)
    anchor_pos_19 = dict(anchor_pos_19)
    
    for node in anchorlist:
        ref_coord = anchor_pos_ref[node]
        alt_coord = anchor_pos_19[node]
        intermediate_coord = anchor_intermediate_dict[node]

        if alignment.get(node) == "?":
            continue

        cigar = "45=" + alignment[node]
        offsets = offsetexchange(cigar, ref_coord, alt_coord, intermediate_coord)
        mapping_table += offsets


    mapping_table = numpy.array(mapping_table)
    mapping_table = mapping_table.astype(int)
    return mapping_table


def main():
    
    # create parser object
    parser = argparse.ArgumentParser(description = "Annotate Genome")
 
    # defining arguments for parser object
    parser.add_argument("-c", "--CCStrain", type = str, nargs = 1,
                        metavar = "ccstrain", default = None,
                        help = "three digit for a cc strain")
    parser.add_argument("-o", "--Output", type = str, nargs = 1,
                    metavar = "output filename", help = "output filename")
    parser.add_argument("-a", "--Anchor", type = str, nargs = 1,
                    metavar = "anchor filename", help = "anchor filename, csv")
    parser.add_argument("-r", "--RefAnnotation", type = str, nargs = 1,
                    metavar = "B38 annotation filename", help = "B38 annotation filename, gff3")
    parser.add_argument("-l", "--Alignment", type = str, nargs = 1,
                    metavar = "Directory of the alignment files", help = "Directory of the alignment files")

     
    
    args = parser.parse_args()
    cc = args.CCStrain[0]
    outputfile = args.Output[0]
    anchorfile =args.Anchor[0]
    refannotationfile = args.RefAnnotation[0]
    alignment_dir = args.Alignment[0]
     
    mapping = ['##gff-version 3\n',
     '#!genome-build CCGG.CC%s\n' % cc,
     '#!genome-version CCGG_Version1\n',
     '#!genome-date %s\n' % str(datetime.datetime.today())]

    with open(outputfile, 'w') as fp:
        for item in mapping:
            fp.write(item)
            
            
    
    chrlist = range(1,21)
    anchor_info = pd.read_csv(anchorfile, index_col = None)

    with gzip.open(refannotationfile, 'r') as fp:
        data = fp.readlines()

    for chromo in chrlist:

        anchor_info_37 = anchor_info[anchor_info['Chromosome'] == chromo]

        anchor_pos_ref = anchor_info_37[['Anchor', 'Pos38']].values
        anchor_pos_37 = anchor_info_37[['Anchor', 'Pos']].values

        alignment = load_alignment_txt(alignment_dir + '/pairwisealignment_%s.txt' % str(chromo))
        if "SOURCE" in alignment:
            del alignment["SOURCE"]
        mapping_table = construct_mapping_table(alignment, anchor_pos_ref, anchor_pos_37)  
        Mapping = dict(zip(mapping_table[:,0], mapping_table[:,2]))

        info = []
        total = 0
        avail = 0

        if chromo == 20:
            chromo = "X"

        for line in data:
            itemlist = line.decode().split('\t')

            if itemlist[0] == str(chromo) or itemlist[0] == chromo:
                total += 1

                start = int(itemlist[3])
                end = int(itemlist[4])

                cc_start = Mapping.get(start, "?")
                cc_end = Mapping.get(end, "?")

                if (cc_start != "?") and (cc_end != "?") and (int(cc_start)>0) and (int(cc_end) >0):
                    avail += 1
                    itemlist[3] = str(cc_start)
                    itemlist[4] = str(cc_end)
                    info.append("\t".join(itemlist))

        print(chromo, total, avail, avail/float(total), len(info))
        del mapping_table
        del Mapping

        with open(outputfile, 'a') as fp:
            for it in info:
                fp.write(it)        


if __name__ == "__main__":
    # calling the main function
    main()

