import sys
import argparse
import logging
import time

def pase_info(seq):
	info = {'SVLEN': 0, 'END': 0, "SVTYPE": '', "RE": 0, "CHR2": ''}
	for i in seq.split(';'):
		if i.split('=')[0] in ["SVLEN", "END", "RE"]:
			try:
				info[i.split('=')[0]] = abs(int(float(i.split('=')[1])))
			except:
				pass
		if i.split('=')[0] in ["CHR2"]:
			info[i.split('=')[0]] = i.split('=')[1]
		if i.split('=')[0] in ["SVTYPE"]:
			info[i.split('=')[0]] = i.split('=')[1][0:3]

	return info

def load_callset(path, filter=0):
	callset = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		pos = int(seq[1])
		info = pase_info(seq[7])

		if info['SVTYPE'] in ['DEL', 'INS', 'DUP', 'INV']:
			if info['SVTYPE'] not in callset:
				callset[info['SVTYPE']] = dict()
			if info['SVLEN'] == 0:
				info['SVLEN'] = info['END'] - pos + 1

			if chr not in callset[info['SVTYPE']]:
				callset[info['SVTYPE']][chr] = list()

			callset[info['SVTYPE']][chr].append([pos, info['END'], info['SVLEN'], [0,0,0]])

		if info['SVTYPE'] == "BND":
			if seq[4][0] == ']':
				form = ']]N'
				chr2 = seq[4].split(':')[0][1:]
				pos2 = int(seq[4].split(':')[1][:-2])
			elif seq[4][0] == '[':
				form = '[[N'
				chr2 = seq[4].split(':')[0][1:]
				pos2 = int(seq[4].split(':')[1][:-2])
			else:
				if seq[4][1] == ']':
					form = 'N]]'
					chr2 = seq[4].split(':')[0][2:]
					pos2 = int(seq[4].split(':')[1][:-1])
				else:
					form = 'N[['
					chr2 = seq[4].split(':')[0][2:]
					pos2 = int(seq[4].split(':')[1][:-1])
			if info['SVTYPE'] not in callset:
				callset[info['SVTYPE']] = dict()
			if info['END'] == 0:
				info['CHR2'] = chr2
				info['END'] = pos2

			if filter > 0:
				if int(seq[-1].split(":")[1].split(',')[1]) < filter:
					continue

			if chr not in callset[info['SVTYPE']]:
				callset[info['SVTYPE']][chr] = list()
			callset[info['SVTYPE']][chr].append([pos, info['CHR2'], info['END'], form, [0,0,0]])

	file.close()
	return callset

def eva_record(call_A, call_B, bias, offect, tag1, tag2):
	# call_A 0/1
	# call_B 1/1
	for svtype in call_A:
		if svtype not in call_B:
			continue

		for chr in call_A[svtype]:
			if chr not in call_B[svtype]:
				continue

			for i in call_A[svtype][chr]:
				for j in call_B[svtype][chr]:

					if svtype == 'INS':
						if abs(i[0]-j[0]) <= offect and float(min(i[2],j[2])/max(i[2],j[2])) >= bias:
							i[-1][tag1] = 1
							j[-1][tag2] = 1
					elif svtype == 'BND':
						if i[1] == j[1] and i[3] == j[3]:
							if abs(i[0]-j[0]) <= offect and abs(i[2]-j[2]) <= offect:
								i[-1][tag1] = 1
								j[-1][tag2] = 1
					else:
						if max(i[0]-offect, j[0]) <= min(i[1]+offect, j[1]) and float(min(i[2],j[2])/max(i[2],j[2])) >= bias:
							i[-1][tag1] = 1
							j[-1][tag2] = 1

def statistics(callset, a, b, c, d):
	for svtype in callset:
		record = 0
		record_000 = 0
		record_100 = 0
		record_010 = 0
		record_001 = 0
		record_110 = 0
		record_101 = 0
		record_011 = 0
		record_111 = 0

		for chr in callset[svtype]:
			for i in callset[svtype][chr]:
				record += 1
				if i[-1][0] == 0 and i[-1][1] == 0 and i[-1][2] == 0:
					record_000 += 1
				if i[-1][0] == 1 and i[-1][1] == 0 and i[-1][2] == 0:
					record_100 += 1
				if i[-1][0] == 0 and i[-1][1] == 1 and i[-1][2] == 0:
					record_010 += 1
				if i[-1][0] == 0 and i[-1][1] == 0 and i[-1][2] == 1:
					record_001 += 1
				if i[-1][0] == 1 and i[-1][1] == 1 and i[-1][2] == 0:
					record_110 += 1
				if i[-1][0] == 1 and i[-1][1] == 0 and i[-1][2] == 1:
					record_101 += 1
				if i[-1][0] == 0 and i[-1][1] == 1 and i[-1][2] == 1:
					record_011 += 1
				if i[-1][0] == 1 and i[-1][1] == 1 and i[-1][2] == 1:
					record_111 += 1
				
		logging.info("%s number of %s:\t%d"%(svtype, a, record))
		logging.info("Only %s:\t%d"%(a, record_000))
		logging.info("%s and %s:\t%d"%(a, b, record_100))
		logging.info("%s and %s:\t%d"%(a, c, record_010))
		logging.info("%s and %s:\t%d"%(a, d, record_001))
		logging.info("%s and %s and %s:\t%d"%(a, b, c, record_110))
		logging.info("%s and %s and %s:\t%d"%(a, b, d, record_101))
		logging.info("%s and %s and %s:\t%d"%(a, c, d, record_011))
		logging.info("%s and %s and %s and %s:\t%d"%(a, b, c, d, record_111))
		logging.info("-----")

def main_ctrl(args):	
	logging.info("Load SV callset of selected caller.")

	cuteSV_callset = load_callset(args.c1)
	# for key in cuteSV_callset:
	# 	print(key, len(cuteSV_callset[key]))
	Sniffles_callset = load_callset(args.c2)
	# for key in Sniffles_callset:
	# 	print(key, len(Sniffles_callset[key]))
	PBSV_callset = load_callset(args.c3, 3)
	# for key in PBSV_callset:
	# 	print(key, len(PBSV_callset[key]))
	SVIM_callset = load_callset(args.c4)
	# for key in SVIM_callset:
	# 	print(key, len(SVIM_callset[key]))

	logging.info("Comparing...")
	eva_record(cuteSV_callset, Sniffles_callset, args.bias, args.offect, 0, 0)
	eva_record(cuteSV_callset, PBSV_callset, args.bias, args.offect, 1, 0)
	eva_record(cuteSV_callset, SVIM_callset, args.bias, args.offect, 2, 0)
	eva_record(Sniffles_callset, PBSV_callset, args.bias, args.offect, 1, 1)
	eva_record(Sniffles_callset, SVIM_callset, args.bias, args.offect, 2, 1)
	eva_record(PBSV_callset, SVIM_callset, args.bias, args.offect, 2, 2)

	logging.info("Final results:")
	statistics(cuteSV_callset, "cuteSV", "Sniffles", "PBSV", "SVIM")
	logging.info("Final results:")
	statistics(Sniffles_callset, "Sniffles", "cuteSV", "PBSV", "SVIM")
	logging.info("Final results:")
	statistics(PBSV_callset, "PBSV", "cuteSV", "Sniffles", "SVIM")
	logging.info("Final results:")
	statistics(SVIM_callset, "SVIM", "cuteSV", "Sniffles", "PBSV")


def main(argv):
	args = parseArgs(argv)
	setupLogging(False)
	# print args
	starttime = time.time()
	main_ctrl(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

USAGE="""\
	Evaluate SV callset generated by cuteSV/Sniffles/PBSV/SVIM
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="Venn_eval", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("c1", type=str, help="cuteSV callset")
	parser.add_argument('c2', type=str, help = "Sniffles callset")
	parser.add_argument('c3', type=str, help = "PBSV callset")
	parser.add_argument('c4', type=str, help = "SVIM callset")
	parser.add_argument('-b', '--bias', help = "Bias of overlaping.[%(default)s]", default = 0.7, type = float)
	parser.add_argument('-o', '--offect', help = "Offect of translocation overlaping.[%(default)s]", default = 1000, type = int)
	args = parser.parse_args(argv)
	return args

def setupLogging(debug=False):
	logLevel = logging.DEBUG if debug else logging.INFO
	logFormat = "%(asctime)s [%(levelname)s] %(message)s"
	logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
	logging.info("Running %s" % " ".join(sys.argv))

if __name__ == '__main__':
	main(sys.argv[1:])