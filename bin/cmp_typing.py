import argparse
import csv
from collections import defaultdict
import re
import glob
import sys
import os


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

def read_summary(csvfile):

    """
    cid,gene_length,total_ATCG,total_nonATCG,total_N,pct_completeness,pct_Ns,numreads_mapped,covbases,coverage,meandepth,clade,clade_database,type,influenza_db_version,typing_sseqid,typing_pident,typing_mismatch,typing_gapopen,typing_eval
ue,typing_qcovs,ref_identity,ref_shared-hashes,ref_query-comment
240112_S_I_008-S10-bwa_bcftools-seg_1-ref_PP664072,2280,2280,0,0,100.0,0.0,23286,2280,100,1205.11,,,,influenzaDB-2024-08-26,,,,,,,0.988479,784/1000,Environment|1|PB2|H3N2|United_Kingdom|A/environment/Northern_Ireland/4/2022|na|na|na|2
280 
240112_S_I_008-S10-bwa_bcftools-seg_2-ref_PP302466,2277,2277,0,0,100.0,0.0,12587,2277,100,687.31,,,,influenzaDB-2024-08-26,,,,,,,0.990074,811/1000,Human|2|PB1|H3N2|USA|A/Georgia/47/2023|na|na|na|2277 
240112_S_I_008-S10-bwa_bcftools-seg_3-ref_MH346871,2154,2152,0,2,99.91,0.09,22468,2154,100,1230.54,,,,influenzaDB-2024-08-26,,,,,,,0.9857,739/1000,Human|3|PA|H3N2|Chile|A/Santiago/PUC-MVL_188/2017|na|na|na|2154 
    """

    summary_dict = {}

    if csvfile is None:
        return summary_dict

    if os.path.exists(csvfile) == False or os.path.isfile(csvfile) == False:
        print(f"\nERROR: typing file {csvfile} does not exist.\n", file=sys.stderr)
        exit(1)
    if os.path.getsize(csvfile) == 0:
        print(f"\nWarning: typing file {csvfile} is empty.\n", file=sys.stderr)
        # exit(1)
        return summary_dict

    print(csvfile)
    with open(csvfile, "r") as csvf:
            reader = csv.DictReader(csvf, delimiter=",")
            cols_to_keep = ["pct_completeness", "clade", "type"]
            
            for row in reader:
                #print(type(row))
                id = row.pop("cid")
                summary_dict[id] = {}
                #type example: CY163681~~7~~A, CY163682~~6~~N2
                for col in cols_to_keep:
                    summary_dict[id][col] = row[col]

    sorted_summary_dict = {i: summary_dict[i] for i in sorted(summary_dict.keys())}
    #print(sorted_summary_dict)
    return sorted_summary_dict

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

def convert(seqdict, allseq_dict, name_dict):
    #print(name_dict)
    # example key: 240112_S_I_008-S13-bwa_freebayes-segment_2-ref_accession_KX351456
   
    for seqid in seqdict:
        #[('240112_S_I_008-S13', 'bwa_freebayes', 'segment_2', 'KX351456')]
        print(seqid)
        match = re.findall("(^\w+-\w+?)-(\w+?)-(\w+?)-ref_accession_(\w+)", seqid)
        g = match[0]
        #sampleid:ref:segid:tool
        sampleid = g[0]
        if sampleid in name_dict:
            sampleid = name_dict[sampleid]
        #print(sampleid)
        #sampleid:refid:segid:tool:
        #print("ttttttttttttt=" + g[1])
        allseq_dict[sampleid][g[3]][g[2]][g[1]] = seqdict[seqid]

    return allseq_dict

# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("-i", "--input", help = "input directory contains cvs summary files", required=True)
parser.add_argument("-n", "--name-list", help = "name list in the format of name\tname1,name2", required=True)


# Read arguments from command line
args = parser.parse_args()

name_dict = defaultdict(str)

if args.name_list is not None:
    name_dict = readNames(args.name_list)
    #print(name_dict)
#ref = readFasta(args.ref_fasta)
allseqs_dict = multi_dict(2, str)
for filename in glob.glob(args.input + '/*.consensus_summary.csv'):
    #print(filename)
    summary_dict = read_summary(filename) 
    #print(read_dict)
    convert(summary_dict, allseqs_dict, name_dict)    
    #print(c_dict)  

"""     
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
            #print(len(tool_list))
            for t1 in tool_list:
                #print(t1)
                tool_list.remove(t1)
                for t2 in tool_list: 
                    
                    #print("xxxxxxxxxxxxxxx")
                    seq1 = allseqs_dict[sample][refid][seqid][t1]
                    seq2 = allseqs_dict[sample][refid][seqid][t2]
                    
                    zip_object = zip(refseq, seq1, seq2)
                    for index, (r, b1, b2) in enumerate(zip_object):
                        if b1 != b2:
                            rows.append([sample, refid, reflen, seqid, index+1, r, t1, b1,t2, b2])

print("\t".join(header))
for row in rows: 
    s = '\t'.join(str(x) for x in row)
    print(s)
                 """