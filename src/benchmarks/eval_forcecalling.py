import sys
import argparse
import logging
import time

def population_statistic(pop_merged_vcf, output_file):
    output = open(output_file, 'w')
    idx = 0
    with open(pop_merged_vcf, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            seq = line.strip().split('\t')
            svlen = abs(int(seq[7].split(';SVLEN=')[1].split(';')[0]))
            svtype = seq[7].split(';SVTYPE=')[1].split(';')[0]
            if (svtype != 'TRA' and svtype != 'BND') and svlen < 50:
                continue
            idx += 1
            af = float(seq[7].split(';AF=')[1].split(';')[0])
            hwe = float(seq[7].split(';HWE=')[1].split(';')[0])
            exchet = float(seq[7].split(';ExcHet=')[1])
            missing_cnt = 0
            for i in range(9, 109, 1):
                if seq[i][0] == '.':
                    missing_cnt += 1
                if seq[i][2] == '.':
                    missing_cnt += 1
            output.write('%d\t%f\t%f\t%f\t%f\n'%(idx, missing_cnt/200, af, hwe, exchet))

def compare_callsets(pop_vcf_file, base_vcf_file, output_file):
    def parse(file):
        svs_dict = dict()
        with open(file, 'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                seq = line.strip().split('\t')
                chrom = seq[0]
                pos = int(seq[1])
                svtype = seq[7].split('SVTYPE=')[1].split(';')[0]
                svlen = abs(int(seq[7].split(';SVLEN=')[1].split(';')[0]))
                if (svtype != 'TRA' and svtype != 'BND') and svlen < 50:
                    continue
                hwe = float(seq[7].split(';HWE=')[1].split(';')[0])
                exchet = float(seq[7].split(';ExcHet=')[1])
                af = float(seq[7].split(';AF=')[1].split(';')[0])
                missing_cnt = 0
                for gt in seq[9:]:
                    if gt[0] == '.':
                        missing_cnt += 1
                    if gt[2] == '.':
                        missing_cnt += 1
                if missing_cnt > 10 or hwe < 0.000001 or exchet < 0.000001:
                    continue
                if svtype == 'DEL' or svtype == 'INS':
                    svlen = abs(int(seq[7].split('SVLEN=')[1].split(';')[0]))
                    af = float(seq[7].split(';AF=')[1].split(';')[0])
                    if chrom not in svs_dict:
                        svs_dict[chrom] = list()
                    svs_dict[chrom].append([pos, svtype, svlen, af])
        return svs_dict
    def parse_base(file):
        svs_dict = dict()
        with open(file, 'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                seq = line.strip().split('\t')
                chrom = seq[0]
                pos = int(seq[1])
                svtype = seq[7].split('SVTYPE=')[1].split(';')[0]
                af = float(seq[7].split(';AF=')[1].split(';')[0])
                if svtype == 'DEL' or svtype == 'INS':
                    svlen = abs(int(seq[7].split('SVLEN=')[1].split(';')[0]))
                    if chrom not in svs_dict:
                        svs_dict[chrom] = list()
                    svs_dict[chrom].append([pos, svtype, svlen, af])
        return svs_dict
    base = parse_base(base_vcf_file)
    comp = parse(pop_vcf_file)
    output = open(output_file, 'w')
    pos_bias = 1000
    length_ratio = 0.7
    for chrom in base:
        if chrom in comp:
            for basesv in base[chrom]:
                for compsv in comp[chrom]:
                    if basesv[1] == compsv[1] and abs(basesv[0] - compsv[0]) <= pos_bias and min(basesv[2], compsv[2]) / max(basesv[2], compsv[2]) > length_ratio:
                        output.write('%s\t%f\t%f\t%f\n'%(basesv[1], basesv[3], compsv[3], basesv[3]-compsv[3]))
                        break
    output.close()

def pre_cmrg(input_vcf, output_vcf):
    # input_vcf -> HG002_GRCh37_CMRG_SV_v1.00.vcf
    output = open(output_vcf, 'w')
    with open(input_vcf, 'r') as f:
        for line in f:
            if line[0] == '#':
                if line[1] != '#':
                    output.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
                    output.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
                output.write(line)
            else:
                seq = line.strip().split('\t')
                ref = seq[3]
                alt = seq[4]
                for i in range(7):
                    output.write(seq[i] + '\t')
                if len(ref) > len(alt): # DEL
                    output.write('SVTYPE=DEL;SVLEN=%d'%(len(alt) - len(ref)))
                else: # INS
                    output.write('SVTYPE=INS;SVLEN=%d'%(len(alt) - len(ref)))
                output.write('\t%s\t%s\n'%(seq[8], seq[9]))
    output.close()

def main(argv):
    args = parseArgs(argv)
    setupLogging(False)
    starttime = time.time()
    if args.handle == 'POP':
        population_statistic(args.input, args.output)
    if args.handle == 'COMP':
        compare_callsets(args.input, args.base_vcf, args.output)
    if args.handel == 'CMRG':
        pre_cmrg(args.input, args.output)
    logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

USAGE="""\
    Processing and evaluation of force calling.
    Author: Shuqi Cao
    Email: sqcao@stu.hit.edu.cn
"""

def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="evaluation on population Statistics/callsets Compare worldwide cohort", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("handle", type=str, help="The aspect of evaluation, contains CMRG/POP/COMP.")
    parser.add_argument("--input", type=str, help="Input VCF, including vcf file to be preprocessed or merged population callsets from force calling methods.")
    parser.add_argument("--base_vcf", type=str, help="Worldwide population callsets.")
    parser.add_argument("--output", type=str, help="Output file.")
    args = parser.parse_args(argv)
    return args

def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
    logging.info("Running %s" % " ".join(sys.argv))

if __name__ == '__main__':
    main(sys.argv[1:])
