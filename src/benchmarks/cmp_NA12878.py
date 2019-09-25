import sys
import argparse
import logging
import time

def load_high_conf_bed(path):
	confbed = dict()
	file = open(path,'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		svtype = seq[1]
		istart = int(seq[2])
		iend = int(seq[3])

		if svtype not in confbed:
			confbed[svtype] = dict()
		if chr not in confbed[svtype]:
			confbed[svtype][chr] = list()
		confbed[svtype][chr].append([istart, iend, 0, 0, 0])
	file.close()
	return confbed

def pase_info(seq):
	info = {'SVLEN': 0, 'END': 0, "SVTYPE": '', "RE": 0, "CHR2": ''}
	for i in seq.split(';'):
		if i.split('=')[0] in ["SVLEN", "END", "RE"]:
			try:
				info[i.split('=')[0]] = abs(int(float(i.split('=')[1])))
			except:
				pass
		if i.split('=')[0] in ["SVTYPE", "CHR2"]:
			info[i.split('=')[0]] = i.split('=')[1]
	return info

def load_callset_cuteSV(path, filter):
	callset = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		pos = int(seq[1])
		ALT = seq[2][7:10]
		# if ALT == "DUP":
		# 	ALT = "INS"
		if ALT not in ["INS", "INV", "DEL", "DUP", "BND"]:
			continue

		info = pase_info(seq[7])

		if ALT not in callset:
			callset[ALT] = dict()
		if chr not in callset[ALT]:
			if ALT == 'BND':
				callset[ALT][chr] = dict()
			else:
				callset[ALT][chr] = list()

		if ALT == "BND":
			if info["CHR2"] not in callset[ALT][chr]:
				callset[ALT][chr][info["CHR2"]] = list()
			# if info["RE"] >= filter:
			callset[ALT][chr][info["CHR2"]].append([pos, info["END"], 0])
		else:
			# if info["SVLEN"] >= 50:
			if info["SVLEN"] >= 50 and info["SVLEN"] <= 100000:
				callset[ALT][chr].append([pos, info["SVLEN"], info["END"], 0])
	file.close()
	return callset

def eva_record(call_A, call_B, bias, offect):
	for svtype in call_A:
		if svtype not in call_B:
			continue
		else:
			if svtype == "BND":
				for chr1 in call_A[svtype]:
					if chr1 not in call_B[svtype]:
						continue
					else:
						for chr2 in call_A[svtype][chr1]:
							if chr2 not in call_B[svtype][chr1]:
								continue
							else:
								for i in call_A[svtype][chr1][chr2]:
									for j in call_B[svtype][chr1][chr2]:
										if abs(i[0] - j[0]) <= offect and abs(i[1] - j[1]) <= offect:
											i[2] += 1
											j[2] += 1
										else:
											pass
			else:
				for chr in call_A[svtype]:
					if chr not in call_B[svtype]:
						continue
					else:
						for i in call_A[svtype][chr]:
							for j in call_B[svtype][chr]:
								if i[0] - offect <= j[0] <= i[2] + offect or i[0] - offect <= j[2] <= i[2] + offect or j[0] - offect <= i[0] <= j[2] + offect:
									if min(i[1], j[1])*1.0/max(i[1], j[1]) >= bias:
										i[3] += 1
										j[3] += 1
									else:
										pass

def statistics_true_possitive(callset, SVTYPE="ALL"):
	record = 0
	non_record = 0
	if SVTYPE == "ALL":
		for svtype in callset:
			if svtype == "BND":
				for chr1 in callset[svtype]:
					for chr2 in callset[svtype][chr1]:
						for res in callset[svtype][chr1][chr2]:
							record += 1
							if res[2] == 0:
								non_record += 1
			else:
				for chr in callset[svtype]:
					for res in callset[svtype][chr]:
						record += 1
						if res[3] == 0:
							non_record += 1
	else:
		if SVTYPE == "BND":
			for chr1 in callset[SVTYPE]:
				for chr2 in callset[SVTYPE][chr1]:
					for res in callset[SVTYPE][chr1][chr2]:
						record += 1
						if res[2] == 0:
							non_record += 1
		else:
			for chr in callset[SVTYPE]:
				for res in callset[SVTYPE][chr]:
					record += 1
					if res[3] == 0:
						non_record += 1
	return record, non_record

def eva_hc(pac, ont, hc, bias, offect):
	for chr in hc['INS']:
		for ele in hc["INS"][chr]:
			if chr not in pac["INS"]:
				pass
			else:
				for i in pac["INS"][chr]:
					if ele[0] - offect <= i[0] <= ele[1] + offect:
						ele[2] = 1
					else:
						pass
			if chr in ont["INS"]:
				for j in ont["INS"][chr]:
					if ele[0] - offect <= j[0] <= ele[1] + offect:
						ele[3] = 1
					else:
						pass

	for chr in hc['DEL']:
		for ele in hc['DEL'][chr]:
			if chr in pac['DEL']:
				for i in pac['DEL'][chr]:
					if ele[0] - offect <= i[0] <= ele[1] + offect or ele[0] - offect <= i[2] <= ele[1] + offect or i[0] - offect <= ele[0] <= i[2] + offect:
						if min(ele[1]-ele[0], i[1])*1.0/max(ele[1]-ele[0], i[1]) >= bias:
							ele[2] = 1
			if chr in ont['DEL']:
				for i in ont['DEL'][chr]:
					if ele[0] - offect <= i[0] <= ele[1] + offect or ele[0] - offect <= i[2] <= ele[1] + offect or i[0] - offect <= ele[0] <= i[2] + offect:
						if min(ele[1]-ele[0], i[1])*1.0/max(ele[1]-ele[0], i[1]) >= bias:
							ele[3] = 1

	for SVTYPE in ["INS", 'DEL']:
		hc_pac = 0
		hc_ont = 0
		hc_pac_ont = 0
		hc_num = 0
		for chr in hc[SVTYPE]:
			for i in hc[SVTYPE][chr]:
				hc_num += 1
				if i[2] == 1:
					hc_pac += 1
				if i[3] == 1:
					hc_ont += 1
				if i[2]+i[3] == 2:
					i[4] = 1
					hc_pac_ont += 1
		logging.info("%s\thc:%d\tPacBio_hc:%d\tONT_hc:%d\tPacBio_ONT_hc:%d."%(SVTYPE, hc_num, hc_pac, hc_ont, hc_pac_ont))

def main_ctrl(args):	
	pac = load_callset_cuteSV(args.C1, 10)
	ont = load_callset_cuteSV(args.C2, 5)
	hc = load_high_conf_bed(args.C3)
	# for SVTYPE in hc:
	# 	for chr in hc[SVTYPE]:
	# 		for i in hc[SVTYPE][chr]:
	# 			print(SVTYPE, chr, i)
	logging.info("Calculate overlaps.")
	eva_record(pac, ont, args.bias, args.offect)
	eva_hc(pac, ont, hc, args.bias, args.offect)
	svtype = ["DEL", "DUP", "INS", "INV", "BND", "ALL"]
	for i in svtype:
		pac_r, pac_nol = statistics_true_possitive(pac, i)
		pac_ol = pac_r - pac_nol
		ont_r, ont_nol = statistics_true_possitive(ont, i)
		ont_ol = ont_r - ont_nol
		logging.info("%s\tPacBio:%d %d %.4f\tONT:%d %d %.4f."%(i, pac_r, pac_ol, 1.0*pac_ol/pac_r, ont_r, ont_ol, 1.0*ont_ol/ont_r))


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
	parser.add_argument("C1", type=str, help="PacBio callset")
	parser.add_argument('C2', type=str, help = "ONT callset")
	parser.add_argument('C3', type=str, help = "High confidence callset")
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