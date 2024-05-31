import argparse
import csv
from collections import defaultdict
import re
import glob



__author__ = 'Xiaoli Dong'
__version__ = '0.1.0'
__maintainer__ = 'Xiaoli Dong'
__email__ = 'xiaoli.dong@albertaprecisionlabs.ca'
__status__ = 'Dev'


# Utility function to create dictionary
def multi_dict(K, type):
    if K == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: multi_dict(K-1, type))
def readNames(nfile):
    name_dict = defaultdict(str)
    with open(nfile) as nf:
        for line in nf:
            match = re.findall("(^\S+?)\s+(\S+)", line)
            g = match[0]
            for name in g[1].split(','):
                name_dict[name] = g[0]
    return name_dict

def readFasta(fasta):
    seqs = defaultdict(str)
    
    with open (fasta, 'r+') as in_f:
        seqid = ""
        for line in in_f:   
            if line[0] == ">":
                seqid = line[1:].split()[0]
                seqs[seqid] = ""
                
            else:
                seqs[seqid] += line.strip()
    return seqs

def convert(seqdict, allseq_dict, name_dict, ref):
    #print(name_dict)
    # example key: 240112_S_I_008-S13-bwa_freebayes-segment_2-ref_accession_KX351456
    for seqid in seqdict:
        #[('240112_S_I_008-S13', 'bwa_freebayes', 'segment_2', 'KX351456')]
        if "ref_accession" not in seqid:
            ref[seqid] = seqdict[seqid]
            continue
        match = re.findall("(^\w+-\w+?)-(\w+?)-(\w+?)-ref_accession_(\w+)", seqid)
        g = match[0]
        #sampleid:ref:segid:tool
        sampleid = g[0]
        if sampleid in name_dict:
            sampleid = name_dict[sampleid]
        #print(sampleid)
        #sampleid:refid:segid:tool:
        allseq_dict[sampleid][g[3]][g[2]][g[1]] = seqdict[seqid]

    return allseq_dict, ref

# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("-i", "--input", help = "input directory containing fasta files", required=True)
parser.add_argument("-n", "--name-list", help = "name list in the format of name\tname1,name2", required=False)


# Read arguments from command line
args = parser.parse_args()

name_dict = defaultdict(str)

if args.name_list is not None:
    name_dict = readNames(args.name_list)
    #print(name_dict)

for filename in glob.glob(args.input + '/*.fasta'):
    allseqs_dict = multi_dict(4, str)
    ref = defaultdict(str)
    #print(filename)
    read_dict = readFasta(filename) 
    #print(read_dict)
    convert(read_dict, allseqs_dict, name_dict, ref)    
    #print(c_dict)  
    f = open(filename + ".diff.tsv", "w")
    
    #print(allseqs_dict.keys()) 
    #print(allseqs_dict)
    header = ["sample", "refid", "reflen", "segid", "pos", "ref", "tool1", "tool1_base", "tool1", "tool2_base"]
    rows = []
    #sampleid:ref:segid:tool
    for sample in allseqs_dict:
        for refid in allseqs_dict[sample]:
            refseq = ref[refid]
            reflen = len(refseq)
            for seqid in allseqs_dict[sample][refid]:
                tool_list = list(allseqs_dict[sample][refid][seqid].keys())
                #print(tool_list)
                for t1 in tool_list:
                    tool_list.remove(t1)
                    for t2 in tool_list: 
                        
                        seq1 = allseqs_dict[sample][refid][seqid][t1]
                        seq2 = allseqs_dict[sample][refid][seqid][t2]
                        
                        zip_object = zip(refseq, seq1, seq2)
                        for index, (r, b1, b2) in enumerate(zip_object):
                            if b1.lower() == 'n' or b2.lower() == 'n':
                                continue
                            elif b1 != b2:
                                rows.append([sample, refid, reflen, seqid, index+1, r, t1, b1,t2, b2])

    print("\t".join(header))
    f.write("\t".join(header) + "\n")
    for row in rows: 
        s = '\t'.join(str(x) for x in row)
        print(s)
        f.write(s + "\n")
    f.close()               