import sys
import argparse
import os

def format_geno(x):
    return {
        '0/0': '00',
        '0/1': '01',
        '1/1': '11',
        '0/2': '02',
        '2/2': '22',
        '0/3': '03',
        '3/3': '33',
        '1/2': '12',
        '1/3': '13',
        '2/3': '23',
        './.': 'NA',
    }.get(x, x) 

def main():
    args = parse_args()

    out = open(args.out, 'w')

    with open(args.vcf, mode = 'r') as file:
        for line in file:
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                line = line.split('\t')
                nsamples = len(line) - 9
                outline = ['rsID']
                outline.extend([line[i] for i in range(9, len(line))])
                out.write(','.join(outline))
            else:
                line = line.split('\t')
                id = line[2]
                # skip missing rs ids
                if id == '.': continue

                # initialize a output line
                outline = [id]
                # initialize a list of genotypes
                genotypes = []
                for col in line[9:(9+nsamples)]:
                    format = col.split(':')
                    # skip malformatted FORMAT fields
                    try:
                        gt = str(format[0])
                        dp = int(format[2])
                    except:
                        continue
                    # label low depth sites as missing
                    if dp < args.DP: gt = './.'
                    # print(id, format_geno(gt), dp)
                    outline.append(format_geno(gt))
                    genotypes.append(gt)
                # skip variants with too many missing genotypes
                if genotypes.count('./.') < args.missing * nsamples:
                    out.write(','.join(outline) + '\n')
    out.close()




def parse_args():
    parser = argparse.ArgumentParser(prog='multivcf2mat.py', description="Convert multi-sample vcf file into custom matrix of genotypes.")

    ## required ##
    essential_args = parser.add_argument_group("Required")
    essential_args.add_argument("--vcf", type=str, help="Path to the multi-sample vcf generated with GATK GenotypeGVCFs. Required.", required=True)
    essential_args.add_argument("--out", type=str, help="Path to output matrix. Required.", required=True)

    ## optional ##
    optional_args = parser.add_argument_group('Optional optional')
    optional_args.add_argument("--DP", type=int, help="Read depth cutoff to determine a missing genotype. [default = 10]", default=10, required=False)
    optional_args.add_argument("--missing", type=float, help="The percentage of missing genotypes to ignore the variant. Default is remove the variants only when genotypes are missing acorss all samples. [default = 1]", default=1, required=False)
    # optional_args.add_argument("--memory", type=str, help="The number of memory in 'G' to use for individual cluster jobs. [default = '60G']", default="60G", required=False)

    args = parser.parse_args()

    ## Required ##
    args.vcf = os.path.abspath(args.vcf)
    args.out = os.path.abspath(args.out)

    ## Optional ##
    args.DP = int(args.DP)
    args.missing = float(args.missing)

    return args

if __name__ == "__main__":
    main()