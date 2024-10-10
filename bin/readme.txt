

python ../../../nf-fluAB/bin/split_and_msa.py  -i freebayes_vs_clair3 -o freebayes_vs_clair3_msa -f ref.fasta -n namelist.tsv
python ../../../../nf-fluAB/bin/cmp_alignedFasta.py -i mafft -n ../namelist.tsv 

2024-05-06 08:21:05.886106    INFO: running_client (basecall_manager)
    arguments: --port
               ipc:///tmp/.guppy/5555
               --server_file_load_timeout
               600--save_path
               
               /home/prl_taqmancalgary/Desktop/minknow_outdir/240501_S_N_075/16S_Noro_Reanalyzed/pod5_pass
               --config
               dna_r10.4.1_e8.2_400bps_5khz_sup.cfg
               --progress_stats_frequency
               2
               --input_path
               /home/prl_taqmancalgary/Desktop/minknow_outdir/240501_S_N_075/16S_Noro_Reanalyzed/pod5_pass
               --compress_fastq
               --recursive
               --barcode_kits
               SQK-NBD114-96
               --enable_trim_barcodes
               --require_barcodes_both_ends
    executable: /opt/ont/dorado/bin/ont_basecall_client


    redo the base call 
    ont_basecall_client --input_path /APL_Genomics/virus_influenza/fluA/2024/240112_S_N_008/FluA/20240112_1042_MN33158_FAX57730_663132be/pod5_pass --save_path  /APL_Genomics/virus_influenza/fluA/2024/240112_S_N_008/FluA/20240112_1042_MN33158_FAX57730_663132be/basecall_sup --config dna_r10.4.1_e8.2_400bps_5khz_sup.cfg --barcode_kits SQK-NBD114-24 --compress_fastq --recursive  --enable_trim_barcodes --require_barcodes_both_ends --port ipc:///tmp/.guppy/5555 --server_file_load_timeout 600


     ont_basecall_client --input_path /APL_Genomics/virus_influenza/fluA/2024/240112_S_N_008/FluA/20240112_1042_MN33158_FAX57730_663132be/pod5_pass --save_path  /APL_Genomics/virus_influenza/fluA/2024/240112_S_N_008/FluA/20240112_1042_MN33158_FAX57730_663132be/basecall_sup --config dna_r10.4.1_e8.2_400bps_5khz_sup.cfg --barcode_kits EXP-NBD196 --compress_fastq --recursive  --enable_trim_barcodes --require_barcodes_both_ends --port ipc:///tmp/.guppy/5555  --server_file_load_timeout 600 --use_tcp