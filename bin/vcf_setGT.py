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
            print(record.ALT[0])
            
            
            if record.ALT[0].value == '*' or record.ALT[0].value == 'NON_REF':
                break

            dp = call.data.get('DP')

            # create a bed file for masking the low depth region
            if dp < args.min_depth:
                break
           
            ad_freq = []
            gt_value = []
            ad_list = call.data.get('AD')
            for var in ad_list:
                # calculate allele freq
                ad_freq.append(var/dp)
            max_freq = max(ad_freq)
            max_freq_index = ad_freq.index(max_freq)
            if max_freq >= float(args.upper_allele_freq_limit):
                call.data['GT'] = str(max_freq_index) + "/" + str(max_freq_index)
                writer.write_record(record) 
            else:
                for i in range(len(ad_freq)):
                    if ad_freq[i] >= args.lower_allele_freq_limit:
                        gt_value.append(i)
                if len(gt_value) == 1:
                    call.data['GT'] = str(gt_value[0]) + "/" + str(gt_value[0])
                    writer.write_record(record) 
                elif len(gt_value) >= 2:
                    call.data['GT'] = "/".join([str(i) for i in gt_value])
                    writer.write_record(record) 
                  
       
if __name__ == "__main__":
    main()
