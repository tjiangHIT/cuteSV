import sys
import argparse
import logging
import time
import vcf

def phase_bnd(alt):
	if alt[0] in [']', '[']:
		chr2 = alt.split(':')[0][1:]
		pos2 = int(alt.split(':')[1][:-2])
	else:
		chr2 = alt.split(':')[0][2:]
		pos2 = int(alt.split(':')[1][:-1])

	return chr2, pos2

def main_ctrl(args):
	filein = vcf.Reader(open(args.vcf, 'r'))
	fileout = open(args.bedpe, 'w')

	fileout.write("#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstarnd1\tstrand2\tsvtype\tnumber_of_support_read\n")

	for record in filein:
		if record.INFO['SVTYPE'] in ["DEL", "INS", "INV", "DUP"]:
			fileout.write("{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{name}\t{score}\t{starnd1}\t{strand2}\t{svtype}\t{RE}\n".format(
				chrom1 = record.CHROM,
				start1 = record.POS+1,
				end1 = record.POS+1,
				chrom2 = record.CHROM,
				start2 = record.INFO["END"]+1,
				end2 = record.INFO["END"]+1,
				name = record.ID,
				score = record.QUAL,
				starnd1 = '+',
				strand2 = '-',
				svtype = record.INFO['SVTYPE'],
				RE = record.INFO["RE"]
				))

		else:
			chr2, pos2 = phase_bnd(str(record.ALT[0]))
			fileout.write("{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{name}\t{score}\t{starnd1}\t{strand2}\t{svtype}\t{RE}\n".format(
				chrom1 = record.CHROM,
				start1 = record.POS+1,
				end1 = record.POS+1,
				chrom2 = chr2,
				start2 = pos2+1,
				end2 = pos2+1,
				name = record.ID,
				score = record.QUAL,
				starnd1 = '+',
				strand2 = '-',
				svtype = record.INFO['SVTYPE'],
				RE = record.INFO["RE"]
				))
	fileout.close()

def main(argv):
	args = parseArgs(argv)
	setupLogging(False)
	starttime = time.time()
	main_ctrl(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

USAGE="""\
	Convert vcf file generated from cuteSV to bedpe file.
	Author: Tao Jiang
	Email: tjiang@hit.edu.cn
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="vcf2bedpe.py", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("vcf", type=str, help="A vcf file generated from cuteSV.")
	parser.add_argument("bedpe", type=str, help="The output bedpe file name.")
	args = parser.parse_args(argv)
	return args

def setupLogging(debug=False):
	logLevel = logging.DEBUG if debug else logging.INFO
	logFormat = "%(asctime)s [%(levelname)s] %(message)s"
	logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
	logging.info("Running %s" % " ".join(sys.argv))

if __name__ == '__main__':
	main(sys.argv[1:])

