
from cuteSV.cuteSV_resolveINV import call_gt as call_gt_inv
from cuteSV.cuteSV_resolveTRA import call_gt as call_gt_tra
from cuteSV.cuteSV_resolveINDEL import call_gt as call_gt_indel
from cuteSV.cuteSV_resolveDUP import call_gt as call_gt_dup
from multiprocessing import Pool, Manager
import vcf
import math
import logging


def parse_sigs(var_type, work_dir):
    var_dict = dict()  # [chrom, start, end, read_id]
    with open(work_dir + var_type + '.sigs', 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            if seq[1] in var_dict:
                if var_type == 'INV':
                    var_dict[seq[1]].append([seq[1], int(seq[3]), int(seq[4]), seq[5]])
                else:
                    var_dict[seq[1]].append([seq[1], int(seq[2]), int(seq[3]), seq[4]])
            else:
                var_dict[seq[1]] = []
                if var_type == 'INV':
                    var_dict[seq[1]].append([seq[1], int(seq[3]), int(seq[4]), seq[5]])
                else:
                    var_dict[seq[1]].append([seq[1], int(seq[2]), int(seq[3]), seq[4]])
    return var_dict


def parse_trasigs(work_dir):
    var_dict = dict()  # var_dict[chrom1][tra_type][chrom2] = [[chrom2, start, end, read_id]]
    with open(work_dir + 'TRA.sigs', 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            chrom1 = seq[1]
            tra_type = seq[2]
            pos1 = int(seq[3])
            chrom2 = seq[4]
            pos2 = int(seq[5])
            read_id = seq[6]
            if chrom1 in var_dict:
                if tra_type in var_dict[chrom1]:
                    if chrom2 in var_dict[chrom1][tra_type]:
                        var_dict[chrom1][tra_type][chrom2].append([chrom2, pos1, pos2, read_id])
                    else:
                        var_dict[chrom1][tra_type][chrom2] = [[chrom2, pos1, pos2, read_id]]
                else:
                    var_dict[chrom1][tra_type] = dict()
                    var_dict[chrom1][tra_type][chrom2] = [[chrom2, pos1, pos2, read_id]]

            else:
                var_dict[chrom1] = dict()
                var_dict[chrom1][tra_type] = dict()
                var_dict[chrom1][tra_type][chrom2] = [[chrom2, pos1, pos2, read_id]]
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
        

def call_gt_wrapper(call_gt_args, gt_list, idx, row_count, record, var_type):
    if var_type == 'INS' or var_type == 'DEL':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_indel(*call_gt_args)
    if var_type == 'DUP':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_dup(*call_gt_args)
    if var_type == 'INV':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_inv(*call_gt_args)
    if var_type == 'TRA':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_tra(*call_gt_args)
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
    if idx > 0 and idx % 5000 == 0:
        logging.info(str(math.floor(idx / row_count * 100)) + '% SV calls of the given vcf has been processed.')


def force_calling(bam_path, ivcf_path, output_path, sigs_dir, max_cluster_bias_dict, gt_round, threads):
    logging.info('Check the parameter -Ivcf: OK.')
    logging.info('Enable to perform force calling.')
    #print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    sv_dict = dict()
    for sv_type in ["DEL", "INS", "INV", "DUP"]:
        sv_dict[sv_type] = parse_sigs(sv_type, sigs_dir)
    sv_dict['TRA'] = parse_trasigs(sigs_dir)
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
        sv_type = record.INFO['SVTYPE']
        sv_len = -record.INFO['SVLEN'] if sv_type == 'DEL' else record.INFO['SVLEN']
        if sv_type == 'TRA':
            tra_alt = str(record.ALT[0])
            if tra_alt[0] == 'N':
                if tra_alt[1] == '[':
                    tra_type = 'A'
                else:
                    tra_type = 'B'
            elif tra_alt[0] == '[':
                tra_type = 'C'
            else:
                tra_type = 'D'
            if tra_alt[0] == 'N':
                tra_alt = tra_alt[2:-1]
            else:
                tra_alt = tra_alt[1:-2]
            sv_chr2 = tra_alt.split(':')[0]
            sv_end = tra_alt.split(':')[1]
        search_start = max(int(pos) - max_cluster_bias_dict[sv_type], 0)
        search_end = search_start + 2 * max_cluster_bias_dict[sv_type]
        search_id_list = []
        if sv_type == 'TRA' and chrom in sv_dict['TRA'] and tra_type in sv_dict['TRA'][chrom] and sv_chr2 in sv_dict['TRA'][chrom][tra_type]:
            search_id_list = sv_dict['TRA'][chrom][tra_type][sv_chr2]
        elif sv_type != 'TRA' and chrom in sv_dict[sv_type]:
            search_id_list = sv_dict[sv_type][chrom]
        read_id_list = find_in_list(sv_type, search_id_list, search_start, search_end, pos, sv_len)

        if sv_type == 'INS':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, chrom, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, record, 'INS'))
        if sv_type == 'DEL':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, chrom, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, record, 'DEL'))
        if sv_type == 'INV':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, sv_end, chrom, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, record, 'INV'))
        if sv_type == 'DUP':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, sv_end, chrom, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, record, 'DUP'))
        if sv_type == 'TRA':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, sv_end, chrom, sv_chr2, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, record, 'TRA'))

    process_pool.close()
    process_pool.join()
    logging.info('Finished force calling.')
    return gt_list


def run_fc(args):
    return force_calling(*args)
