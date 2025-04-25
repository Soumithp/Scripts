import sys
import argparse
import os
import gzip
import re

def main():
    args = parse_args()

    oheader = open(args.oheader, 'w')
    omat = open(args.omat, 'w')

    with gzip.open(args.input, mode = 'rt') as file:
        fields = {}
        count = 0
        for line in file:
            if line.startswith('##INFO='):
                field = {}
                info_header = line.replace('##INFO=<','').replace('>','')
                p = re.compile('ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="([^"]+)"')
                m = p.match(info_header)
                field = {
                    'ID': m.group(1),
                    'Number': m.group(2),
                    'Type': m.group(3),
                    'Description': '"' + m.group(4) + '"'
                }

                oheader.write(",".join(field.values()) + '\n')
                fields[field['ID']] = field
                 
            elif not line.startswith('#'):
                if count == 0:
                    omat.write('ID,REF,ALT,' + ",".join(fields.keys()) + '\n')
                count += 1
                
                # init info list for each record
                info = {key: "NA" for key in fields.keys()}

                line = line.strip()
                line = line.split('\t')
                rsid = line[2]
                ref_allele = '"{}"'.format(line[3]); alt_allele = '"{}"'.format(line[4])
                info_col = line[7].split(';')
                for i in info_col:
                    try:
                        a,b = i.split('=')
                        info[a] = b.strip('\n')
                    except:
                        if i != '':
                            info[i] = "yes"
                omat.write(rsid + ',' + ref_allele + ',' + alt_allele + ',' + ",".join(['"{}"'.format(value) for value in info.values()]) + '\n')

    oheader.close()
    omat.close()




def parse_args():
    parser = argparse.ArgumentParser(prog='multivcf2mat.py', description="Convert multi-sample vcf file into custom matrix of genotypes.")

    ## required ##
    essential_args = parser.add_argument_group("Required")
    essential_args.add_argument("--input", type=str, help="input dbSNP/COSMIC. Required.", required=True)
    essential_args.add_argument("--oheader", type=str, help="exlanation of headers. Required.", required=True)
    essential_args.add_argument("--omat", type=str, help="annotation in matrix CSV format. Required.", required=True)


    args = parser.parse_args()

    ## Required ##
    args.input = os.path.abspath(args.input)
    args.oheader = os.path.abspath(args.oheader)
    args.omat = os.path.abspath(args.omat)
    

    return args

if __name__ == "__main__":
    main()