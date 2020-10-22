from cuteSV.cuteSV_resolveINV import call_gt as call_gt_inv
from cuteSV.cuteSV_resolveTRA import call_gt as call_gt_tra
from cuteSV.cuteSV_resolveINDEL import call_gt as call_gt_indel
from cuteSV.cuteSV_resolveDUP import call_gt as call_gt_dup
from multiprocessing import Pool, Manager
import vcf
import time
import logging


def parse_sigs(var_type, work_dir):
    var_dict = dict()  # [chrom, pos, len, read_id] or [chrom1, pos1, chrom2, pos2, read_id]
    with open(work_dir + var_type + '.sigs', 'r') as f:
        for line in f:
           seq = line.strip('\n').split('\t')
           if seq[1] in var_dict:
                var_dict[seq[1]].append([seq[1], int(seq[2]), int(seq[3]), seq[4]])
           else:
                var_dict[seq[1]] = []
                var_dict[seq[1]].append([seq[1], int(seq[2]), int(seq[3]), seq[4]])
    return var_dict


# var_list read_id near pos in [start, end], the reads support the variant
def find_in_list(var_type, var_list, start, end, pos, sv_end):
    if len(var_list) == 0:
        return []
    left = 0
    right = len(var_list) - 1
    mid = 0
    while left < right:
        mid = (left + right) >> 1
        if var_list[mid][1] < start:
            left = mid + 1
        else:
            right = mid
    read_id_list = []
    for i in range(right, len(var_list)):
        if var_list[i][1] > end:
            break
        if abs(var_list[i][2] - sv_end) < 1000:
            read_id_list.append(var_list[i][3])
    return read_id_list
        

def call_gt_wrapper(call_gt_args, gt_list, idx, record, var_type):
    if var_type == 'INS' or var_type == 'DEL':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_indel(*call_gt_args)
        gt_list[idx] = [record.CHROM,
                        record.POS,
                        genotype,
                        record.INFO['SVTYPE'],
                        record.INFO['SVLEN'],
                        record.INFO['END'],
                        record.INFO['CIPOS'],
                        record.INFO['CIEND'],
                        gt_re,
                        ','.join(call_gt_args[3]),
                        record.ID,
                        record.REF, 
                        record.ALT[0],
                        record.QUAL,
                        record.INFO['STRANDS'],
                        record.FILTER
        ]
    if idx > 0 and idx % 4000 == 0:
        logging.info('Finished ' + str(idx - 3999) + ' to ' + str(idx))
        #print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Finished ' + str(idx))
    
    
def force_calling(bam_path, ivcf_path, output_path, sigs_dir, max_cluster_bias_dict, gt_round, threads):
    logging.info('Start force calling.')
    #print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Start force calling.')
    ins_dict = parse_sigs('INS', sigs_dir)
    del_dict = parse_sigs('DEL', sigs_dir)    
    max_cluster_bias_INS = max_cluster_bias_dict['INS']
    max_cluster_bias_DEL = max_cluster_bias_dict['DEL']
    vcf_reader = vcf.Reader(filename = ivcf_path)
    row_count = 0
    for record in vcf_reader:
        row_count += 1
    idx = -1
    gt_list = Manager().list([[] for x in range(row_count)])
    result = []
    process_pool = Pool(processes = threads)
    vcf_reader = vcf.Reader(filename = ivcf_path)
    for record in vcf_reader:
        idx += 1
        chrom = record.CHROM
        pos = record.POS
        var_type = record.INFO['SVTYPE']
        sv_len = record.INFO['SVLEN']
        if var_type == 'INS':
            search_start = max(int(pos) - max_cluster_bias_INS, 0)
            search_end = search_start + 2 * max_cluster_bias_INS
            read_id_list = find_in_list(var_type, ins_dict[chrom], search_start, search_end,
                    pos, sv_len)
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, chrom, read_id_list, max_cluster_bias_INS, gt_round], gt_list, idx, record, 'INS'))
        if var_type == 'DEL':
            search_start = max(int(pos) - max_cluster_bias_DEL, 0)
            search_end = search_start + 2 * max_cluster_bias_DEL
            read_id_list = find_in_list(var_type, del_dict[chrom], search_start, search_end,
                    pos, -sv_len) 
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, chrom, read_id_list, max_cluster_bias_DEL, gt_round], gt_list, idx, record, 'DEL'))
    process_pool.close()
    process_pool.join()
    logging.info('Finish calling')
    #print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ' Finish calling')
    #result = gt_list
    return gt_list


def run_fc(args):
    return force_calling(*args)