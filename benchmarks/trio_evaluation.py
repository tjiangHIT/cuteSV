import sys
import argparse
import logging
import time

def load_parent_info(path):
	homo_P = dict()
	total_P = dict()
	linkage = 0
	file = open(path, "r")
	for line in file:
		seq = line.strip('\n').split('\t')
		if len(seq) == 7:
			chr = seq[0]
			pos = int(seq[1])
			svlen = int(seq[2])
			svtype = seq[5]
			GT = seq[6]

			if svtype not in total_P:
				total_P[svtype] = dict()
			if chr not in total_P[svtype]:
				total_P[svtype][chr] = dict()
			hash_1 = int(pos/10000)
			mod = pos % 10000
			hash_2 = int(mod / 50)
			if hash_1 not in total_P[svtype][chr]:
				total_P[svtype][chr][hash_1] = dict()
			if hash_2 not in total_P[svtype][chr][hash_1]:
				total_P[svtype][chr][hash_1][hash_2] = list()

			if GT == '1/1':
				linkage += 1
				total_P[svtype][chr][hash_1][hash_2].append([pos, svlen, 'Y', linkage])
				homo_P[linkage] = seq
			else:
				total_P[svtype][chr][hash_1][hash_2].append([pos, svlen, 'N', 0])

		else:
			# TRA
			chr_1 = seq[0]
			breakpoint_1 = int(seq[1])
			chr_2 = seq[2]
			breakpoint_2 = int(seq[3])
			svtype = seq[6]
			GT = seq[7]

			if svtype not in total_P:
				total_P[svtype] = dict()
			if chr_1 not in total_P[svtype]:
				total_P[svtype][chr_1] = dict()
			hash_1 = int(breakpoint_1/10000)
			mod = breakpoint_1 % 10000
			hash_2 = int(mod / 50)
			if hash_1 not in total_P[svtype][chr_1]:
				total_P[svtype][chr_1][hash_1] = dict()
			if hash_2 not in total_P[svtype][chr_1][hash_1]:
				total_P[svtype][chr_1][hash_1][hash_2] = list()
			# total_P[svtype][chr_1][hash_1][hash_2].append([breakpoint_1, chr_2, breakpoint_2])
			if GT == '1/1':
				linkage	+= 1
				total_P[svtype][chr_1][hash_1][hash_2].append([breakpoint_1, chr_2, breakpoint_2, 'Y', linkage])
				homo_P[linkage] = seq
			else:
				total_P[svtype][chr_1][hash_1][hash_2].append([breakpoint_1, chr_2, breakpoint_2, 'N', 0])

	file.close()
	return homo_P, total_P

def acquire_locus(down, up, keytype, chr, MainCandidate):
	if keytype not in MainCandidate:
		return []
	if chr not in MainCandidate[keytype]:
		return []
	if int(up/10000) == int(down/10000):
		key_1 = int(down/10000)
		if key_1 not in MainCandidate[keytype][chr]:
			return []
		for i in xrange(int((up%10000)/50)-int((down%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = int((down%10000)/50)+i
			if key_2 not in MainCandidate[keytype][chr][key_1]:
				continue
			for ele in MainCandidate[keytype][chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					return ele
	else:
		key_1 = int(down/10000)
		if key_1 in MainCandidate[keytype][chr]:
			for i in xrange(200-int((down%10000)/50)):
				# exist a bug ***********************************
				key_2 = int((down%10000)/50)+i
				if key_2 not in MainCandidate[keytype][chr][key_1]:
					continue
				for ele in MainCandidate[keytype][chr][key_1][key_2]:
					if ele[0] >= down and ele[0] <= up:
						return ele
		key_1 += 1
		if key_1 not in MainCandidate[keytype][chr]:
			return []
		for i in xrange(int((up%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = i
			if key_2 not in MainCandidate[keytype][chr][key_1]:
				continue
			for ele in MainCandidate[keytype][chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					return ele
	return []

def main_ctrl(args):
	logging.info("Loading Male parent callsets.")
	homo_MP, total_MP = load_parent_info(args.MP)
	# print len(homo_MP)

	logging.info("Loading Female parent callsets.")
	homo_FP, total_FP = load_parent_info(args.FP)
	# print len(homo_FP)
	# for svtype in total_FP:
	# 	# print svtype
	# 	for chr in total_FP[svtype]:
	# 		# print chr
	# 		for hash_1 in total_FP[svtype][chr]:
	# 			for hash_2 in total_FP[svtype][chr][hash_1]:
	# 				print total_FP[svtype][chr][hash_1][hash_2]
	# 				# break
	logging.info("Loading Offspring callsets.")
	recall = 0
	total_call = 0
	file = open(args.F1, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		total_call += 1
		if len(seq) == 7:
			flag = 0
			chr = seq[0]
			pos = int(seq[1])
			svlen = int(seq[2])
			svtype = seq[5]
			ans_MP = acquire_locus(pos-50, pos+50, svtype, chr, total_MP)
			if len(ans_MP) > 0:
				flag = 1
				if ans_MP[2] == 'Y':
					homo_MP[ans_MP[3]] = 1
			ans_FP = acquire_locus(pos-50, pos+50, svtype, chr, total_FP)
			if len(ans_FP) > 0:
				flag = 1
				if ans_FP[2] == 'Y':
					homo_FP[ans_FP[3]] = 1
			if flag == 1:
				recall += 1

		else:
			flag = 1
			chr_1 = seq[0]
			breakpoint_1 = int(seq[1])
			chr_2 = seq[2]
			breakpoint_2 = int(seq[3])
			svtype = seq[6]
			ans_MP = acquire_locus(pos-50, pos+50, svtype, chr, total_MP)
			if len(ans_MP) > 0:
				flag = 1
				# print ans_MP
				if ans_MP[3] == 'Y':
					homo_MP[ans_MP[4]] = 1
			ans_FP = acquire_locus(pos-50, pos+50, svtype, chr, total_FP)
			if len(ans_FP) > 0:
				flag = 1
				if ans_FP[3] == 'Y':
					homo_FP[ans_FP[4]] = 1
			if flag == 1:
				recall += 1

	file.close()
	logging.info("Accuracy rate %d/%d"%(recall, total_call))
	total_right = 0
	for key in homo_MP:
		if homo_MP[key] == 1:
			total_right += 1
		else:
			print homo_MP[key]
	for key in homo_FP:
		if homo_FP[key] == 1:
			total_right += 1
		else:
			print homo_FP[key]
	logging.info("Sensitivity rate %d/%d"%(total_right, len(homo_FP)+len(homo_MP)))

def main(argv):
	args = parseArgs(argv)
	setupLogging(False)
	# print args
	starttime = time.time()
	main_ctrl(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

USAGE="""\
	Evaluate SV callset generated by cuteSV
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="Trio_eval", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("MP", type=str, help="Male parent callsets")
	parser.add_argument('FP', type=str, help = "Female parent callsets")
	parser.add_argument('F1', type=str, help = "Offspring callsets")
	# parser.add_argument('-s', '--min_support', help = "Minimum number of reads that support a SV to be reported.[%(default)s]", default = 10, type = int)
	args = parser.parse_args(argv)
	return args

def setupLogging(debug=False):
	logLevel = logging.DEBUG if debug else logging.INFO
	logFormat = "%(asctime)s [%(levelname)s] %(message)s"
	logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
	logging.info("Running %s" % " ".join(sys.argv))

if __name__ == '__main__':
	main(sys.argv[1:])
