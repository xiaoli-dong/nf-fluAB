
# fluA, B
sequence data, xml, and csv file downloaded from: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=taxid:197911&VirusLineage_ss=taxid:197912&VirusLineage_ss=taxid:197913&VirusLineage_ss=taxid:1511083&LabHost_s=include on 2024-08-12

#meta data downloaded from: 2024-08-23
https://www.bv-brc.org/view/Taxonomy/11308#view_tab=genomes&filter=false
BVBRC_genome.csv

#get rid of too short, without required meta data fileds 
python ../bin/parseCSV.py --input BVBRC_genome.csv > BVBRC_genome.filter_and_reformat.csv

#only keep the sequences which contained in BVBRC_genome.filter_and_reformat.csv
python ../bin/build_influenza_db.py  --csv BVBRC_genome.filter_and_reformat.csv --fasta sequences.fasta > sequences.filter.fasta

mothur > trim.seqs(fasta=sequences.filter.fasta, maxambig=0, maxhomop=10)
mothur > summary.seqs(fasta=sequences.filter.trim.fasta)

Using 24 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       800     800     0       4       1
2.5%-tile:      1       838     838     0       5       25991
25%-tile:       1       1401    1401    0       6       259908
Median:         1       1683    1683    0       6       519816
75%-tile:       1       2215    2215    0       6       779723
97.5%-tile:     1       2341    2341    0       7       1013640
Maximum:        1       2600    2600    0       10      1039630
Mean:   1       1663    1663    0       5
# of Seqs:      1039630

#cluster away the near identical sequences to reduce the db size 119626
mmseqs easy-cluster ../sequences.filter.trim.fasta clusterRes tmp --min-seq-id 0.99 -c 0.8 --cov-mode 1

#filter out problematic sequences
fasta-trim-terminal-ambigs.pl mmseqs/clusterRes_rep_seq.fasta > mmseqs/clusterRes_rep_seq.vadr_trim.fasta
v-annotate.pl --split --cpu 8 -r --atgonly --xnocomp --nomisc --alt_fail extrant5,extrant3 --mkey flu --mdir /nfs/APL_Genomics/db/prod/vadr/vadr-models-flu-1.6.3-2 mmseqs/clusterRes_rep_seq.vadr_trim.fasta vadr
