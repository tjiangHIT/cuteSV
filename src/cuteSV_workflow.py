
import pysam
import logging
import os

from multiprocessing import Pool
from cuteSV_utils import exe


def alignment_preprocessing(BAM_PATH, BATCHES):
	'''cuteSV preprocess alignments and establish task scheduling strategies.'''
	samfile = pysam.AlignmentFile(BAM_PATH)
	contig_num = len(samfile.get_index_statistics())
	logging.info("Chromosome total number is: %d"%(contig_num))

	Task_list = list()
	chr_name_list = list()
	for chr in samfile.get_index_statistics():
		chr_name_list.append(chr[0])
		local_ref_len = samfile.get_reference_length(chr[0])
		if local_ref_len < BATCHES:
			'''length of chromosome shorter than BATCHES'''
			Task_list.append([chr[0], 0, local_ref_len])
		else:
			'''length of chromosome longer than BATCHES'''
			local_pos = 0
			task_round = int(local_ref_len / BATCHES)
			for j in xrange(task_round):
				Task_list.append([chr[0], local_pos, local_pos + BATCHES])
				local_pos += BATCHES
			if local_pos < local_ref_len:
				Task_list.append([chr[0], local_pos, local_ref_len])
	samfile.close()
	return Task_list, chr_name_list

# def multi_run_wrapper(args):
# 	return single_pipe(*args)

def sigs_extract(args, Tasks):
	'''extraction of DEL/DUP/INS/INV/TRA signals'''
	temp_folder = args.temp_dir if args.temp_dir.endswith('/') else '%s/'%(args.temp_dir)
	os.mkdir("%ssignatures"%temp_folder)
	analysis_pools = Pool(processes = int(args.threads))
	for i in Tasks:
		para = [(args.input, args.min_length, args.min_mapq, args.max_split_parts, args.min_seq_size, temp_folder, i)]
		analysis_pools.map_async(single_pipe, *para)
	analysis_pools.close()
	analysis_pools.join()
	return temp_folder

def merge_sigs(temp_folder, args):
	'''signatures intergration'''
	logging.info("Merge SV signals by SV-subtypes.")
	analysis_pools = Pool(processes=int(args.threads))
	cmd_del = ("cat %ssignatures/*.bed | grep DEL | sort -u | sort -k 2,2 -k 3,3n > %sDEL.sigs"%(temp_folder, temp_folder))
	cmd_ins = ("cat %ssignatures/*.bed | grep INS | sort -u | sort -k 2,2 -k 3,3n > %sINS.sigs"%(temp_folder, temp_folder))
	cmd_inv = ("cat %ssignatures/*.bed | grep INV | sort -u | sort -k 2,2 -k 3,3n > %sINV.sigs"%(temp_folder, temp_folder))
	cmd_tra = ("cat %ssignatures/*.bed | grep TRA | sort -u | sort -k 2,2 -k 4,4 -k 3,3n > %sTRA.sigs"%(temp_folder, temp_folder))
	cmd_dup = ("cat %ssignatures/*.bed | grep DUP | sort -u | sort -k 1,1r -k 2,2 -k 3,4n > %sDUP.sigs"%(temp_folder, temp_folder))
	for i in [cmd_ins, cmd_del, cmd_dup, cmd_tra, cmd_inv]:
		analysis_pools.map_async(exe, (i,))
	analysis_pools.close()
	analysis_pools.join()

def clean_func(temp_folder):
	'''clean temporary folders and files'''
	logging.info("Cleaning temporary files.")
	cmd_remove_tempfile = ("rm -r %ssignatures %s*.sigs"%(temp_folder, temp_folder))
	exe(cmd_remove_tempfile)

def wrk_flow(args):
	'''cuteSV workflow'''

	# data preprocessing
	Tasks, chr_names = alignment_preprocessing(args.input, args.batches)

	# extraction of signals
	temp_folder = sigs_extract(args, Tasks)

	# merging signals by sv types
	merge_sigs(temp_folder, args)

	# cluster
	chr_names.sort()
	process_tra = list()
	result = list()
	for i in xrange(len(chr_names)-1):
		for j in chr_names[i+1:]:
			process_tra.append([chr_names[i], j])


	logging.info("Clustering structural variants.")
	analysis_pools = Pool(processes=int(args.threads))
	for chr in chr_names:
		for svtype in ["DEL", "INS", "INV"]:
			if svtype == "INV":
				# para = [("%s%s.sigs"%(temporary_dir, svtype), chr, svtype, args.min_support, 0.5, 20, args.min_length)]
				para = [("%s%s.sigs"%(temporary_dir, svtype), chr, svtype, args.min_support, 0.5, args.max_cluster_bias, args.min_length, args.max_distance)]
			else:
				# para = [("%s%s.sigs"%(temporary_dir, svtype), chr, svtype, args.min_support, 0.7, 10, args.min_length)]
				para = [("%s%s.sigs"%(temporary_dir, svtype), chr, svtype, args.min_support, 0.5, args.max_cluster_bias, args.min_length, args.max_distance)]
			# resolution_INDEL(sys.argv[1], "1", "INV", 10, 0.5, 20, 50)
			result.append(analysis_pools.map_async(run_indel_inv, para))
			# print chr, svtype
		# DUP
		# resolution_DUP(path, chr, read_count, max_cluster_bias, sv_size=50)
		# resolution_DUP(sys.argv[1], "1", 10, 50, 50)
		para = [("%s%s.sigs"%(temporary_dir, "DUP"), chr, args.min_support, args.max_cluster_bias, args.min_length, args.max_distance)]
		result.append(analysis_pools.map_async(run_dup, para))
	for i in process_tra:
		# TRA
		# resolution_TRA(path, chr_1, chr_2, read_count, overlap_size, max_cluster_bias)
		para = [("%s%s.sigs"%(temporary_dir, "TRA"), i[0], i[1], args.min_support, 0.6, args.max_cluster_bias)]
		result.append(analysis_pools.map_async(run_tra, para))
	analysis_pools.close()
	analysis_pools.join()

	semi_result = list()
	for res in result:
		semi_result += res.get()[0]

	logging.info("Writing into disk.")
	file = open(args.output, 'w')
	for i in semi_result:
		file.write(i)
	file.close()

	# clean files
	clean_func(temp_folder)
