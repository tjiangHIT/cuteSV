import sys
import argparse
import logging
import time
import vcf

def call_gt(tag):
	if sum(tag) == 2:
		return '1/1'
	else:
		if tag[0] == 1:
			return '1/0'
		if tag[1] == 1:
			return '0/1'

	return './.'

def main_ctrl(args):

	fileout = open(args.outvcf, 'w')

	header = open(args.invcf, 'r')
	for line in header:
		if line[0] == '#':
			fileout.write(line)
		else:
			continue
	header.close()

	filein = vcf.Reader(open(args.invcf, 'r'))
	for record in filein:
		if len(record.FILTER) == 0:
			filter_table = 'PASS'
		else:
			filter_table = record.FILTER[0]
		tag = [0, 0]
		for i in record.INFO['RNAMES']:
			if 'cutesvh1' in i:
				tag[0] = 1
			if 'cutesvh2' in i:
				tag[1] = 1

		try:
			fileout.write("{chr}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{gt}\n".format(
				chr = record.CHROM,
				pos = record.POS,
				id = record.ID,
				ref = record.REF,
				alt = record.ALT[0],
				qual = record.QUAL,
				filter = filter_table,
				info = "SVTYPE=%s;SVLEN=%d;END=%d;RE=%d;RNAMES=%s"%(record.INFO['SVTYPE'], 
					record.INFO['SVLEN'], 
					record.INFO['END'], 
					record.INFO['RE'], 
					','.join(record.INFO['RNAMES'])),
				format = 'GT',
				# gt = '1/1' if sum(tag) == 2 else '0/1'
				gt = call_gt(tag)
				))
		except:
			if 'TRA' in record.INFO['SVTYPE'] or 'BND' in record.INFO['SVTYPE']:
				fileout.write("{chr}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{gt}\n".format(
					chr = record.CHROM,
					pos = record.POS,
					id = record.ID,
					ref = record.REF,
					alt = record.ALT[0],
					qual = record.QUAL,
					filter = filter_table,
					info = "SVTYPE=%s;RE=%d;RNAMES=%s"%(record.INFO['SVTYPE'], 
						record.INFO['RE'], 
						','.join(record.INFO['RNAMES'])),
					format = 'GT',
					# gt = '1/1' if sum(tag) == 2 else '0/1'
					gt = call_gt(tag)
				))
			else:
				pass
	fileout.close()

def main(argv):
	args = parseArgs(argv)
	setupLogging(False)
	starttime = time.time()
	main_ctrl(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

USAGE="""\
	Convert the typical SV callsets generated from cuteSV to diploid based SV callsets.
	Author: Tao Jiang
	Email: tjiang@hit.edu.cn
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="diploid_calling.py", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("invcf", type=str, help="A vcf file generated from cuteSV.")
	parser.add_argument("outvcf", type=str, help="The output diploid based SV callsets.")
	args = parser.parse_args(argv)
	return args

def setupLogging(debug=False):
	logLevel = logging.DEBUG if debug else logging.INFO
	logFormat = "%(asctime)s [%(levelname)s] %(message)s"
	logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
	logging.info("Running %s" % " ".join(sys.argv))

if __name__ == '__main__':
	main(sys.argv[1:])

