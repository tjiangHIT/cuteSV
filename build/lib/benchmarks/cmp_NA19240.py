import sys
import argparse
import logging
import time

callset = {1: "cuteSV", 2: "Sniflles", 3: "PBSV", 4: "SVIM"}

USAGE="""\
	Evaluate SV callset on NA19240 dataset
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="NA19240_eval", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("base", type=str, help="Base vcf file of NA19240.")
	parser.add_argument("cuteSV", type=str, help="CuteSV vcf file of NA19240.")
	parser.add_argument("sniffles", type=str, help="Sniffles vcf file of NA19240.")
	parser.add_argument("pbsv", type=str, help="PBSV vcf file of NA19240.")
	parser.add_argument("svim", type=str, help="SVIM vcf file of NA19240.")
	parser.add_argument('-b', '--bias', help = "Bias of overlaping.[%(default)s]", default = 0.7, type = float)
	parser.add_argument('-o', '--offect', help = "Offect of breakpoint distance.[%(default)s]", default = 1000, type = int)
	args = parser.parse_args(argv)
	return args

def pase_base_info(seq):
	info = {'SVLEN': 0, 'END': 0, "SVTYPE": '', "RE": 0}
	for i in seq.split(';'):
		if i.split('=')[0] in ["SVLEN", "END", "RE"]:
			try:
				info[i.split('=')[0]] = abs(int(i.split('=')[1]))
			except:
				pass
		if i.split('=')[0] == "SVTYPE":
			info[i.split('=')[0]] = i.split('=')[1][0:3]
	return info


def load_base(base_path):
	base_call = dict()
	file = open(base_path, 'r')
	for line in file:
		seq = line.strip('\n').split("\t")
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		pos = int(seq[1])
		ALT = seq[4][1:4]
		if ALT not in ["INS", "INV", "DEL", "DUP"]:
			continue
		if ALT == "DUP":
			ALT = "INS"
		info = pase_base_info(seq[7])
		if ALT not in base_call:
			base_call[ALT] = dict()

		if chr not in base_call[ALT]:
			base_call[ALT][chr] = list()

		if ALT == "INV":
			base_call[ALT][chr].append([pos, info["END"] - pos + 1, info["END"], 0])
		else:
			if info["SVLEN"] >= 50 and info["SVLEN"] <= 100000:
				base_call[ALT][chr].append([pos, info["SVLEN"], info["END"], 0])
	file.close()
	return base_call

def load_cuteSV(cuteSV_path):
	# inv_tag = 0
	last_inv = list()
	cuteSV_call = dict()
	file = open(cuteSV_path, 'r')
	for line in file:
		seq = line.strip('\n').split("\t")
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		pos = int(seq[1])
		ALT = seq[2][7:10]
		# if ALT == "DUP":
		# 	ALT = "INS"
		if ALT not in ["INS", "INV", "DEL", "DUP"]:
			continue

		info = pase_base_info(seq[7])
		if ALT not in cuteSV_call:
			cuteSV_call[ALT] = dict()

		if chr not in cuteSV_call[ALT]:
			cuteSV_call[ALT][chr] = list()

		if info["SVLEN"] >= 50 and info["SVLEN"] <= 100000:
			if ALT == "INV":
				last_inv.append([ALT, chr, pos, info["SVLEN"], info["END"], info["RE"]])
				# if inv_tag == 0
			else:
				cuteSV_call[ALT][chr].append([pos, info["SVLEN"], info["END"], 0])
				# inv_tag = 0
				if len(last_inv):
					last_inv = sorted(last_inv, key = lambda x:-x[3])
					cuteSV_call[last_inv[0][0]][last_inv[0][1]].append([last_inv[0][2], last_inv[0][3], last_inv[0][4], 0])
					last_inv = list()
	file.close()
	return cuteSV_call

def load_sniffles(sniffles_path):
	sniffles_call = dict()
	last_inv = list()
	file = open(sniffles_path, 'r')
	for line in file:
		seq = line.strip('\n').split("\t")
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		pos = int(seq[1])
		info = pase_base_info(seq[7])

		if info["SVTYPE"] not in ["INS", "INV", "DEL", "DUP"]:
			continue

		# if info["SVTYPE"] == "DUP":
		# 	info["SVTYPE"] = "INS"

		if info["SVTYPE"] not in sniffles_call:
			sniffles_call[info["SVTYPE"]] = dict()

		if chr not in sniffles_call[info["SVTYPE"]]:
			sniffles_call[info["SVTYPE"]][chr] = list()

		if info["SVLEN"] >= 50 and info["SVLEN"] <= 100000:
			if info["SVTYPE"] == "INV":
				last_inv.append([info["SVTYPE"], chr, pos, info["SVLEN"], info["END"], info["RE"]])
			else:
				sniffles_call[info["SVTYPE"]][chr].append([pos, info["SVLEN"], info["END"], 0])
				if len(last_inv):
					last_inv = sorted(last_inv, key = lambda x:-x[3])
					sniffles_call[last_inv[0][0]][last_inv[0][1]].append([last_inv[0][2], last_inv[0][3], last_inv[0][4], 0])
					last_inv = list()

	file.close()
	return sniffles_call

def load_pbsv(pbsv_path):
	pbsv_call = dict()
	file = open(pbsv_path, 'r')
	for line in file:
		seq = line.strip('\n').split("\t")
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		pos = int(seq[1])

		info = pase_base_info(seq[7])

		if info["SVTYPE"] not in ["INS", "INV", "DEL", "DUP"]:
			continue

		# if info["SVTYPE"] == "DUP":
		# 	info["SVTYPE"] = "INS"

		if info["SVTYPE"] not in pbsv_call:
			pbsv_call[info["SVTYPE"]] = dict()

		if chr not in pbsv_call[info["SVTYPE"]]:
			pbsv_call[info["SVTYPE"]][chr] = list()

		if info["SVTYPE"] == "INV":
			pbsv_call[info["SVTYPE"]][chr].append([pos, info["END"] - pos + 1, info["END"], 0])
		else:
			if info["SVLEN"] >= 50 and info["SVLEN"] <= 100000:
				pbsv_call[info["SVTYPE"]][chr].append([pos, info["SVLEN"], info["END"], 0])
	file.close()
	return pbsv_call

def load_svim(base_path):
	base_call = dict()
	file = open(base_path, 'r')
	for line in file:
		seq = line.strip('\n').split("\t")
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		pos = int(seq[1])
		ALT = seq[4][1:4]
		if ALT not in ["INS", "INV", "DEL", "DUP"]:
			continue
		# if ALT == "DUP":
		# 	ALT = "INS"
		info = pase_base_info(seq[7])
		if ALT not in base_call:
			base_call[ALT] = dict()

		if chr not in base_call[ALT]:
			base_call[ALT][chr] = list()

		if ALT == "INV":
			base_call[ALT][chr].append([pos, info["END"] - pos + 1, info["END"], 0])
		else:
			if info["SVLEN"] >= 50 and info["SVLEN"] <= 100000:
				base_call[ALT][chr].append([pos, info["SVLEN"], info["END"], 0])
	file.close()
	return base_call

def cmp_callsets(base, call, flag, Bias, Offect):
	for svtype in base:
		if svtype not in call:
			continue
		else:
			for chr in base[svtype]:
				if chr not in call[svtype]:
					continue
				else:
					for i in base[svtype][chr]:
						for j in call[svtype][chr]:
							if i[0] - Offect <= j[0] <= i[2] + Offect or i[0] - Offect <= j[2] <= i[2] + Offect or j[0] - Offect <= i[0] <= j[2] + Offect:
								if min(i[1], j[1])*1.0/max(i[1], j[1]) >= Bias:
									i[3] = flag
									j[3] = flag
								else:
									pass
	total_base = 0
	tp_base = 0
	# for svtype in ["INS"]:
	# for svtype in ["DUP"]:
	# for svtype in ["DEL"]:
	# for svtype in ["INS", "DEL"]:
	for svtype in ["INS", "DEL", "INV"]:
	# for svtype in ["INS", "DEL", "INV", "DUP"]:
	# for svtype in ["INV"]:
		for chr in base[svtype]:
			for i in base[svtype][chr]:
				total_base += 1
				if i[3] == flag:
					tp_base += 1
				# else:
				# 	print(flag, svtype, chr, i[0], i[1], i[2])
	# logging.info("Base count: %d"%(total_base))
	# logging.info("TP-base count: %d"%(tp_base))
	logging.info("====%s===="%(callset[flag]))
	total_call = 0
	tp_call = 0
	# for svtype in ["INS"]:
	# for svtype in ["DUP"]:
	# for svtype in ["DEL"]:
	# for svtype in ["INS", "DEL"]:
	for svtype in ["INS", "DEL", "INV"]:
	# for svtype in ["INS", "DEL", "INV", "DUP"]:
	# for svtype in ["INV"]:
		for chr in call[svtype]:
			for i in call[svtype][chr]:
				total_call += 1
				if i[3] == flag:
					tp_call += 1

	
	logging.info("Camp count: %d"%(total_call))
	logging.info("TP-call count: %d"%(tp_call))
	logging.info("Precision: %.2f"%(100.0*tp_call/total_call))
	logging.info("Recall: %.2f"%(100.0*tp_base/total_base))
	logging.info("F-measure: %.2f"%(200.0*tp_base*tp_call/(total_base*tp_call+tp_base*total_call)))



def main_ctrl(args):
	# pass
	base_call = load_base(args.base)
	cuteSV_call = load_cuteSV(args.cuteSV)
	sniffles_call = load_sniffles(args.sniffles)
	pbsv_call = load_pbsv(args.pbsv)
	svim_call = load_svim(args.svim)
	# for svtype in sniffles_call:
	# 	for chr in sniffles_call[svtype]:
	# 		for i in sniffles_call[svtype][chr]:
	# 			print(svtype, chr, i)

	cmp_callsets(base_call, cuteSV_call, 1, args.bias, args.offect)
	cmp_callsets(base_call, sniffles_call, 2, args.bias, args.offect)
	cmp_callsets(base_call, pbsv_call, 3, args.bias, args.offect)
	cmp_callsets(base_call, svim_call, 4, args.bias, args.offect)
	

def main(argv):
	args = parseArgs(argv)
	setupLogging(False)
	# print args
	starttime = time.time()
	main_ctrl(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

def setupLogging(debug=False):
	logLevel = logging.DEBUG if debug else logging.INFO
	logFormat = "%(asctime)s [%(levelname)s] %(message)s"
	logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
	logging.info("Running %s" % " ".join(sys.argv))

if __name__ == '__main__':
	main(sys.argv[1:])
