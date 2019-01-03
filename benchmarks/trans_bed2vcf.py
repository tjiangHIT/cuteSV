import sys

def main(argv):
	file = open(argv[0], 'r')
	num = 0
	print("##fileformat=VCFv4.2")
	print("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">")
	print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">")
	print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">")
	print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	print("##INFO=<ID=RE,Number=1,Type=Integer,Description=\"read support\">")
	print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	print("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference reads\">")
	print("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant reads\">")
	print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002")
	for line in file:
		seq = line.strip("\n").split("\t")
		chr = seq[0]
		if int(seq[4]) < 10:
			continue
		svtype = seq[1]
		if svtype == "TRA":
			#pass
			continue
		elif svtype == "INV":
			pos = seq[2]
			svlen = str(int(seq[3])-int(seq[2]))
			num += 1
			re = seq[4]
		else:
			num += 1
			pos = seq[2]
			svlen = seq[3]
			re = seq[4]
		print("%s\t%s\t%d\tN\t<%s>\t.\tPASS\tPRECISE;SVTYPE=%s;SVLEN=%s;RE=%s\tGT:DR:DV\t./.:10:10"%(chr, pos, num, svtype, svtype, svlen, re))


if __name__ == '__main__':
	main(sys.argv[1:])
