#!/usr/bin/env python

import argparse
import vcfpy
import sys

def main():
    description = """
        program is used to set GT value (0/0, 1/1, 0/1, 1/2, 1/2/3) for the vcf recrod


        """
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=f"input vcf file\n",
    )
    parser.add_argument(
        "-t",
        "--tbi",
        required=True,
        help=f"input vcf tbi file\n",
    )

    parser.add_argument(
        "-o",
        "--output-vcf",
        required=True,
        help=f"the output file\n",
    )
    
    parser.add_argument(
        "-l",
        "--lower-allele-freq-limit",
        type=float,
        default=0.25,
        help=f"The lowest allele freq cutoff for a biallelic variant to be kept\n",
    )
    parser.add_argument(
        "-u",
        "--upper-allele-freq-limit",
        type=float,
        default=0.75,
        help=f"The lowest allele frequence limit for a variant to be included into the consensus\n",
    )
    
    parser.add_argument(
        "-d",
        "--min-depth",
        type=int,
        default=10,
        help=f"The min read depth requred to keep a site\n",
    )
    
    

    args = parser.parse_args()
    
    reader = vcfpy.Reader.from_path(args.input, tabix_path=args.tbi)
    writer = vcfpy.Writer.from_path(args.output_vcf, reader.header)
    #print(reader.header)
    for record in reader:
        #print(record.FORMAT, file=sys.stderr)
        #print(record.CHROM, file=sys.stderr)
        #print(record.POS, file=sys.stderr)

        #print(record.calls)
        for call in record.calls:
            #print(call.data.get('GT') or './.', file=sys.stderr)
            #print(call.data.get('AD') or '', file=sys.stderr)
            #print(call.data.get('DP') or '', file=sys.stderr)

             #ALT is a list
            #filter out the regions having no variants
            #print(record.ALT[0])
            
            
            if record.ALT[0].value == '*' or record.ALT[0].value == 'NON_REF':
                #print(len(record.REF), file=sys.stderr)
                break
            
            dp = call.data.get('DP')

            # create a bed file for masking the low depth region
            if dp < args.min_depth:
                break
           
            
            ad_freq_dict = {}
            ad_list = call.data.get('AD')
            
            reflen = len(record.REF[0])
            reffreq = ad_list[0]/dp
            
            if reffreq >= args.lower_allele_freq_limit:
                ad_freq_dict[0] = reffreq

            i = 1 #alt index started with 1 in AD tag array
            
            max_freq = reffreq
            max_index = 0

            while i < len(ad_list):
                altfreq = ad_list[i]/dp
                    # ignore low frequency variants
                if altfreq >= args.lower_allele_freq_limit:
                    ad_freq_dict[i] = altfreq
                    if altfreq > max_freq:
                        max_freq = altfreq
                        max_index = i
                i += 1
              
            keys = list(ad_freq_dict.keys())

            if max_freq >= float(args.upper_allele_freq_limit):
                call.data['GT'] = str(max_index) + "/" + str(max_index)
                writer.write_record(record) 
            elif len(keys) >= 1:
                #print(keys)
                #all.data['GT'] = "/".join(keys)
                if len(keys) == 1:
                    call.data['GT'] = str(keys[0]) + "/" + str(keys[0])
                    writer.write_record(record)
                else:
                    call.data['GT'] = "/".join([str(i) for i in keys])
                    writer.write_record(record) 
                         
if __name__ == "__main__":
    main()
