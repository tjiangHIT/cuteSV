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

def phase_GT(seq):
	i = seq.split(':')
	if i[0] in ['0/1', '1/0']:
		return 'het'
	elif i[0] == '1/1':
		return 'hom'
	else:
		return 'unknown'


def load_callset(path, svtype_list):
	callset = dict()
	abtype = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		pos = int(seq[1])
		info = pase_info(seq[7])

		if len(svtype_list) == 3 and info['SVTYPE'] == "DUP":
			info['SVTYPE'] = "INS"

		if info['SVTYPE'] in svtype_list:
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
					callset[info['SVTYPE']] = list()
				if info['END'] == 0:
					info['CHR2'] = chr2
					info['END'] = pos2
				try:
					if int(chr) <= int(info['CHR2']):
						if form == 'N[[':
							form = ']]N'
						if form == ']]N':
							form = 'N[['
						callset[info['SVTYPE']].append([chr, pos, info['CHR2'], info['END'], form, phase_GT(seq[9]), 0])
					else:
						callset[info['SVTYPE']].append([info['CHR2'], info['END'], chr, pos, form, phase_GT(seq[9]), 0])
				except:
					callset[info['SVTYPE']].append([chr, pos, info['CHR2'], info['END'], form, phase_GT(seq[9]), 0])

			else:
				if info['SVTYPE'] not in callset:
					callset[info['SVTYPE']] = list()
				if info['SVLEN'] == 0:
					info['SVLEN'] = info['END'] - pos + 1
				callset[info['SVTYPE']].append([chr, pos, info['END'], info['SVLEN'], phase_GT(seq[9]), 0])
		else:
			if info['SVTYPE'] not in abtype:
				abtype[info['SVTYPE']] = 0
			abtype[info['SVTYPE']] += 1

	file.close()
	return callset, abtype

def eval(call, ans, bias, offect, opt, genotype):
	for svtype in call:
		if svtype not in ans:
			# continue
			if svtype == 'INS':
				for i in call[svtype]:
					for key in ans:
						for j in ans[key]:
							if i[0] == j[0]:
								if abs(i[1] - j[1]) <= offect and float(min(i[3],j[3])/max(i[3],j[3])) >= bias:
									i[-1] = 1
									j[3+opt] = 1
									if i[4] == genotype[j[0]]:
										i[-1] = 2
										j[3+opt] = 2
		else:
			for i in call[svtype]:
				for j in ans[svtype]:
					if i[0] != j[0]:
						continue
					else:
						if svtype in ['INS']:
							if abs(i[1] - j[1]) <= offect and float(min(i[3],j[2])/max(i[3],j[2])) >= bias:
								j[2+opt] = 1
								i[-1] = 1
								if i[4] == genotype[j[0]]:
									j[2+opt] = 2
									i[-1] = 2
						elif svtype == 'BND':
							if i[2] != j[2]:
								continue
							else:
								# if i[4] == j[4]:
								if 1:
									if abs(i[1]-j[1]) <= offect and abs(i[3]-j[3]) <= offect:
										i[-1] = 1
										j[4+opt] = 1
										if i[5] == genotype[j[0]] or i[5] == genotype[j[2]]:
											i[-1] = 2
											j[4+opt] = 2
						else:
							if max(i[1]-offect, j[1]) <= min(i[2]+offect, j[2]) and float(min(i[3],j[3])/max(i[3],j[3])) >= bias:
								j[3+opt] = 1
								i[-1] = 1
								if i[4] == genotype[j[0]]:
									j[3+opt] = 2
									i[-1] = 2
						# genotype


def statistics(call, ans, opt, res):
	for svtype in call:
		tp = 0
		total = 0
		for ele in call[svtype]:
			total += 1
			if ele[-1] >= res:
				tp += 1
			# if ele[-1] < res and ele[-1] == 1:
			# 	print(ele)
		logging.info('TP-%d of %s:\t%d\t%d'%(res, svtype, tp, total))

	for svtype in ans:
		fn = 0
		total = 0
		for ele in ans[svtype]:
			total += 1
			if svtype == 'INS':
				if ele[2+opt] >= res:
					fn += 1
			elif svtype == 'BND':
				if ele[4+opt] >= res:
					fn += 1
			else:
				if ele[3+opt] >= res:
					fn += 1
		logging.info('TN-%d of %s:\t%d\t%d'%(res, svtype, fn, total))

typetrans = {'insertion':'INS', 
			'deletion':'DEL', 
			'inversion':'INV',
			'tandem duplication':'DUP',
			'reciprocal translocation':'BND'
			}

def load_ans(path):
	ansbed = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		svtype = typetrans[seq[3]]
		start = int(seq[1])
		end = int(seq[2])

		if svtype not in ansbed:
			ansbed[svtype] = list()

		if svtype in ['INS']:
			ansbed[svtype].append([chr, start, len(seq[4]), 0, 0, 0, 0])
		elif svtype in ['BND']:
			chr2 = seq[4].split(':')[1]
			start2 = int(seq[4].split(':')[2])
			strand1 = seq[4].split(':')[3]
			strand2 = seq[4].split(':')[4]

			if strand1[0] == 'f':
				if strand2[0] == 'f':
					ansbed[svtype].append([chr, start, chr2, start2, "N[[", 0, 0, 0, 0])
					# ansbed[svtype].append([chr, start, chr2, start2, "]]N", 0, 0, 0, 0])
					# ansbed[svtype].append([chr, end, chr2, start2+end-start, "]]N", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2+end-start, "N[[", 0, 0, 0, 0])
				else:
					ansbed[svtype].append([chr, start, chr2, start2, "N[[", 0, 0, 0, 0])
					ansbed[svtype].append([chr, start, chr2, start2+end-start, "[[N", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2, "N]]", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2+end-start, "]]N", 0, 0, 0, 0])
			else:
				if strand2[0] == 'f':
					ansbed[svtype].append([chr, start, chr2, start2+end-start, "N]]", 0, 0, 0, 0])
					ansbed[svtype].append([chr, start, chr2, start2, "]]N", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2, "[[N", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2+end-start, "N[[", 0, 0, 0, 0])
				else:
					ansbed[svtype].append([chr, start, chr2, start2+end-start, "N]]", 0, 0, 0, 0])
					ansbed[svtype].append([chr, end, chr2, start2, "N]]", 0, 0, 0, 0])
					# ansbed[svtype].append([chr, end, chr2, start2, "[[N", 0, 0, 0, 0])
					# ansbed[svtype].append([chr, start, chr2, start2+end-start, "[[N", 0, 0, 0, 0])
		else:
			ansbed[svtype].append([chr, start, end, end-start+1, 0, 0, 0, 0])

	file.close()
	return ansbed

def load_gt(path):
	GT = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		if chr not in GT:
			GT[chr] = ''
		if float(seq[-1]) > 80.0:
			GT[chr] = 'hom'
		elif 80.0 >= float(seq[-1]) > 20.0:
			GT[chr] = 'het'
		else:
			GT[chr] = 'None'
	return GT

def main_ctrl(args):

	# load ground truth set
	ans = load_ans(args.ans)
	genotype = load_gt(args.gt)

	if args.choice == "BND":
		cuteSV, abcuteSV = load_callset(args.cuteSV, ["BND"])
		logging.info("The number of calls within abnormal SV type in cuteSV:")
		for key in abcuteSV:
			logging.info("<cuteSV-%s>\t%d."%(key, abcuteSV[key]))
		logging.info("Evaluation on cuteSV callsets...")
		eval(cuteSV, ans, args.bias, args.offect, 1, genotype)
		statistics(cuteSV, ans, 1, 1)
		statistics(cuteSV, ans, 1, 2)

		sniffles, absniffles = load_callset(args.sniffles, ["BND"])
		logging.info("The number of calls within abnormal SV type in Sniffles:")
		for key in absniffles:
			logging.info("<Sniffles-%s>\t%d."%(key, absniffles[key]))
		logging.info("Evaluation on Sniffles callsets...")
		eval(sniffles, ans, args.bias, args.offect, 2, genotype)
		statistics(sniffles, ans, 2, 1)
		statistics(sniffles, ans, 2, 2)

		pbsv, abpbsv = load_callset(args.pbsv, ["BND"])
		logging.info("The number of calls within abnormal SV type in PBSV:")
		for key in abpbsv:
			logging.info("<PBSV-%s>\t%d."%(key, abpbsv[key]))
		logging.info("Evaluation on PBSV callsets...")
		eval(pbsv, ans, args.bias, args.offect, 3, genotype)
		statistics(pbsv, ans, 3, 1)
		statistics(pbsv, ans, 3, 2)

		svim, absvim = load_callset(args.svim, ["BND"])
		logging.info("The number of calls within abnormal SV type in SVIM:")
		for key in absvim:
			logging.info("<SVIM-%s>\t%d."%(key, absvim[key]))
		logging.info("Evaluation on SVIM callsets...")
		eval(svim, ans, args.bias, args.offect, 4, genotype)
		statistics(svim, ans, 4, 1)
		statistics(svim, ans, 4, 2)

	if args.choice == "DUP":
		cuteSV, abcuteSV = load_callset(args.cuteSV, ["INS", "DUP"])
		logging.info("The number of calls within abnormal SV type in cuteSV:")
		for key in abcuteSV:
			logging.info("<cuteSV-%s>\t%d."%(key, abcuteSV[key]))
		logging.info("Evaluation on cuteSV callsets...")
		eval(cuteSV, ans, args.bias, args.offect, 1, genotype)
		statistics(cuteSV, ans, 1, 1)
		statistics(cuteSV, ans, 1, 2)

		sniffles, absniffles = load_callset(args.sniffles, ["INS", "DUP"])
		logging.info("The number of calls within abnormal SV type in Sniffles:")
		for key in absniffles:
			logging.info("<Sniffles-%s>\t%d."%(key, absniffles[key]))
		logging.info("Evaluation on Sniffles callsets...")
		eval(sniffles, ans, args.bias, args.offect, 2, genotype)
		statistics(sniffles, ans, 2, 1)
		statistics(sniffles, ans, 2, 2)

		pbsv, abpbsv = load_callset(args.pbsv, ["INS", "DUP"])
		logging.info("The number of calls within abnormal SV type in PBSV:")
		for key in abpbsv:
			logging.info("<PBSV-%s>\t%d."%(key, abpbsv[key]))
		logging.info("Evaluation on PBSV callsets...")
		eval(pbsv, ans, args.bias, args.offect, 3, genotype)
		statistics(pbsv, ans, 3, 1)
		statistics(pbsv, ans, 3, 2)

		svim, absvim = load_callset(args.svim, ["INS", "DUP"])
		logging.info("The number of calls within abnormal SV type in SVIM:")
		for key in absvim:
			logging.info("<SVIM-%s>\t%d."%(key, absvim[key]))
		logging.info("Evaluation on SVIM callsets...")
		eval(svim, ans, args.bias, args.offect, 4, genotype)
		statistics(svim, ans, 4, 1)
		statistics(svim, ans, 4, 2)

	if args.choice == "IID":
		
		cuteSV, abcuteSV = load_callset(args.cuteSV, ["INS", "INV", "DEL"])
		logging.info("The number of calls within abnormal SV type in cuteSV:")
		for key in abcuteSV:
			logging.info("<cuteSV-%s>\t%d."%(key, abcuteSV[key]))
		logging.info("Evaluation on cuteSV callsets...")
		eval(cuteSV, ans, args.bias, args.offect, 1, genotype)
		statistics(cuteSV, ans, 1, 1)
		statistics(cuteSV, ans, 1, 2)

		sniffles, absniffles = load_callset(args.sniffles, ["INS", "INV", "DEL"])
		logging.info("The number of calls within abnormal SV type in Sniffles:")
		for key in absniffles:
			logging.info("<Sniffles-%s>\t%d."%(key, absniffles[key]))
		logging.info("Evaluation on Sniffles callsets...")
		eval(sniffles, ans, args.bias, args.offect, 2, genotype)
		statistics(sniffles, ans, 2, 1)
		statistics(sniffles, ans, 2, 2)

		pbsv, abpbsv = load_callset(args.pbsv, ["INS", "INV", "DEL"])
		logging.info("The number of calls within abnormal SV type in PBSV:")
		for key in abpbsv:
			logging.info("<PBSV-%s>\t%d."%(key, abpbsv[key]))
		logging.info("Evaluation on PBSV callsets...")
		eval(pbsv, ans, args.bias, args.offect, 3, genotype)
		statistics(pbsv, ans, 3, 1)
		statistics(pbsv, ans, 3, 2)

		svim, absvim = load_callset(args.svim, ["INS", "INV", "DEL"])
		logging.info("The number of calls within abnormal SV type in SVIM:")
		for key in absvim:
			logging.info("<SVIM-%s>\t%d."%(key, absvim[key]))
		logging.info("Evaluation on SVIM callsets...")
		eval(svim, ans, args.bias, args.offect, 4, genotype)
		statistics(svim, ans, 4, 1)
		statistics(svim, ans, 4, 2)

		# ***********optional output****************
		# '''
		# for key in ans:
		# 	for i in ans[key]:
		# 		if key == 'INS':
		# 			if i[3] == 0:
		# 				print('%s\t%d\t%s\t%d'%(i[0], i[1], key, i[2]))
		# 		elif key == 'DEL':
		# 			if i[4] == 0:
		# 				print('%s\t%d\t%s\t%d'%(i[0], i[1], key, i[3]))
		# 		else:
		# 			if i[4] == 0:
		# 				print(key, i)
		# '''
		# for key in cuteSV:
		# 	for i in cuteSV[key]:
		# 		if key == 'INS':
		# 			if i[4] == 0:
		# 				print('%s\t%d\t%s\t%d'%(i[0], i[1], key, i[3]))
		# 		else:
		# 			if i[4] == 0:
		# 				print('%s\t%d\t%s\t%d'%(i[0], i[1], key, i[3]))


def main(argv):
	args = parseArgs(argv)
	setupLogging(False)
	# print args
	starttime = time.time()
	main_ctrl(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

USAGE="""\
	Evaluate SV callset generated by simulations.
	Author: Tao Jiang
	Email: tjiang@hit.edu.cn
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="Trio_eval", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("choice", type=str, help="Chose specific SV type.[IID/DUP/BND]")
	parser.add_argument("ans", type=str, help="Ground truth of simulations.")
	parser.add_argument("gt", type=str, help="Genotype in each chromosome.")
	parser.add_argument('cuteSV', type=str, help = "cuteSV callsets")
	parser.add_argument('sniffles', type=str, help = "Sniffles callsets")
	parser.add_argument('pbsv', type=str, help = "PBSV callsets")
	parser.add_argument('svim', type=str, help = "SVIM callsets")
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