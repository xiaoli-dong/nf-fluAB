import argparse
import csv
from collections import defaultdict
import re
import glob
import os
import subprocess

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

# pair the consensus with their own ref and write to files
def cat_seqs(allseq_dict, outdir):
    
    for sample in allseqs_dict:
        
        for refid in allseqs_dict[sample]:
            refseq = ref[refid]
            for seqid in allseqs_dict[sample][refid]:
                tool_list = list(allseqs_dict[sample][refid][seqid].keys())
                #print(tool_list)
                fname = sample + "-" + seqid + "." + "fasta"
                file_path = os.path.join(outdir, fname)
                f = open(file_path, "w")
                for t1 in tool_list:
                    tool_list.remove(t1)
                    for t2 in tool_list: 
                        seq1 = allseqs_dict[sample][refid][seqid][t1]
                        seq2 = allseqs_dict[sample][refid][seqid][t2]
                        # example key: 230721_S_I_314-S87-bwa_freebayes-seg_1-ref_OQ367282
                        #header = [sample, t1, seqid, refid]
                        f.write(">" + refid + "\n")
                        f.write(refseq + "\n")
                        f.write(">" + sample + "-" + t1 + "-" + seqid + "-ref_" + refid + "\n")
                        f.write(seq1 + "\n")
                        f.write(">" + sample + "-" + t2 + "-" + seqid + "-ref_" + refid + "\n")
                        f.write(seq2 + "\n")
                f.close()

def convert(seqdict, allseq_dict, name_dict):
    #print(name_dict)
    # example key: 230721_S_I_314-S87-bwa_freebayes-seg_1-ref_OQ367282
    for seqid in seqdict:
        #[('240112_S_I_008-S13', 'bwa_freebayes', 'segment_2', 'KX351456')]
       
        match = re.findall("(^\w+-\w+?)-(\w+?)-(\w+?)-ref_(\w+)", seqid)
        g = match[0]
        #sampleid:ref:segid:tool
        sampleid = g[0]
        if sampleid in name_dict:
            sampleid = name_dict[sampleid]
        #print(sampleid)
        #sampleid:refid:segid:tool:
        allseq_dict[sampleid][g[3]][g[2]][g[1]] = seqdict[seqid]

    return allseq_dict

# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("-i", "--indir", help = "input directory containing fasta files", required=True)
parser.add_argument("-o", "--outdir", help = "outdir directory containing mafft aligned fasta files", required=True)
parser.add_argument("-f", "--fasta", help = "reference fasta", required=True)
parser.add_argument("-n", "--name-list", help = "name list in the format of name\tname1,name2", required=False)


# Read arguments from command line
args = parser.parse_args()

name_dict = defaultdict(str)

if args.name_list is not None:
    name_dict = readNames(args.name_list)
    print(name_dict)

ref = readFasta(args.fasta)
allseqs_dict = multi_dict(4, str)
for filename in glob.glob(args.indir + '/*.fa'):
    #print(filename)
    read_dict = readFasta(filename) 
    #print(read_dict)
    convert(read_dict, allseqs_dict, name_dict)    
    #print(c_dict)  

#prepare output dir for split segment fasta file
seg_fasta = os.path.join(args.outdir, "fasta")
mafft_fasta = os.path.join(args.outdir, "mafft")

if not os.path.exists(args.outdir):
    os.mkdir(path = args.outdir)

if not os.path.exists(seg_fasta):
    os.mkdir(path = seg_fasta)
if not os.path.exists(mafft_fasta):
    os.mkdir(path = mafft_fasta)

cat_seqs(allseqs_dict, seg_fasta)   
#do msa over each segment
for file in os.listdir(seg_fasta):
    filename = os.fsdecode(file)
    if filename.endswith(".fasta"): 
        # mafft on the fasta file
        input_file = os.path.join(seg_fasta, filename)
        output_file = os.path.join(mafft_fasta, filename + ".mafft.fasta")
        cmd = "mafft --auto " + input_file + " > " + output_file
        print(cmd)
        returned_value = subprocess.call(cmd, shell=True)  # returns the exit code in unix
        #print('returned value:', returned_value)
    else:
        continue
