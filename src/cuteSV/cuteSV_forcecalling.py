
from cuteSV.cuteSV_resolveINV import call_gt as call_gt_inv
from cuteSV.cuteSV_resolveTRA import call_gt as call_gt_tra
from cuteSV.cuteSV_resolveINDEL import call_gt as call_gt_indel
from cuteSV.cuteSV_resolveDUP import call_gt as call_gt_dup
from multiprocessing import Pool, Manager
import vcf
from pysam import VariantFile
import math
import time
import numpy as np
import logging


class Para(object):
    def __init__(self, record):
        self.chrom = record.chrom
        self.pos = parse_to_int(record.pos)
        self.svlen = parse_to_int(record.info['SVLEN']) if 'SVLEN' in record.info else 0
        self.end = parse_to_int(record.stop)
        self.cipos = record.info['CIPOS'] if 'CIPOS' in record.info else ('.', '.')
        self.ciend = record.info['CIEND'] if 'CIEND' in record.info else ('.', '.')
        self.id = record.id
        self.ref = record.ref
        self.alts = record.alts[0]
        self.qual = record.qual


def parse_to_int(sth):
    if sth == None:
        return 0
    elif isinstance(sth, str):
        return int(sth)
    elif isinstance(sth, list):
        return parse_to_int(sth[0])
    elif isinstance(sth, tuple):
        return parse_to_int(sth[0])
    elif isinstance(sth, int):
        return sth
    else:
        return sth


def parse_record(record):
    sv_type = record.info['SVTYPE']
    if sv_type == 'BND':
        sv_type = 'TRA'
    chrom1 = record.chrom
    start = parse_to_int(record.pos)
    if (record.info['SVTYPE'] == 'INS' or record.info['SVTYPE'] == 'DEL') and 'SVLEN' in record.info:  ###
        end = abs(parse_to_int(record.info['SVLEN']))
    else:
        try:
            end = parse_to_int(record.stop)
        except:
            pass   
    if sv_type == 'TRA':
        tra_alt = str(record.alts[0])
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
        chrom2 = tra_alt.split(':')[0]
        end = int(tra_alt.split(':')[1])
    strand = '.'
    if sv_type != 'TRA':
        chrom2 = record.chrom
        if 'STRAND' in record.info:
            strand = record.info['STRAND']
        elif 'STRANDS' in record.info:
            strand = record.info['STRANDS']
    return sv_type, chrom1, chrom2, start, end, strand


def parse_sigs(var_type, work_dir):
    var_dict = dict()  #var_dict[chrom] = [chrom, start, end, read_id]
    with open(work_dir + var_type + '.sigs', 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            if seq[1] in var_dict:
                var_dict[seq[1]].append([seq[1], int(seq[2]), int(seq[3]), seq[4]])
            else:
                var_dict[seq[1]] = []
                var_dict[seq[1]].append([seq[1], int(seq[2]), int(seq[3]), seq[4]])
    return var_dict


def parse_invsigs(work_dir):
    var_dict = dict() # var_dict[chrom][strand] = [chrom, start, end, read_id]
    with open(work_dir + 'INV.sigs', 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            chrom = seq[1]
            strand = seq[2]
            if chrom in var_dict:
                if strand in var_dict[chrom]:
                    var_dict[chrom][strand].append([chrom, int(seq[3]), int(seq[4]), seq[5]])
                else:
                    var_dict[chrom][strand] = [[chrom, int(seq[3]), int(seq[4]), seq[5]]]
            else:
                var_dict[chrom] = dict()
                var_dict[chrom][strand] = [[chrom, int(seq[3]), int(seq[4]), seq[5]]]
    return var_dict


def parse_trasigs(work_dir):
    var_dict = dict()  #var_dict[chrom1][chrom2] = [[chrom2, start, end, read_id]]
    with open(work_dir + 'TRA.sigs', 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            chrom1 = seq[1]
            tra_type = seq[2]
            pos1 = int(seq[3])
            chrom2 = seq[4]
            pos2 = int(seq[5])
            read_id = seq[6]
            '''
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
            '''
            if chrom1 in var_dict:
                if chrom2 in var_dict[chrom1]:
                    var_dict[chrom1][chrom2].append([chrom2, pos1, pos2, read_id])
                else:
                    var_dict[chrom1][chrom2] = [[chrom2, pos1, pos2, read_id]]

            else:
                var_dict[chrom1] = dict()
                var_dict[chrom1][chrom2] = [[chrom2, pos1, pos2, read_id]]
    for chr1 in var_dict:
        for chr2 in var_dict[chr1]:
            var_dict[chr1][chr2].sort(key=lambda x:x[1])
    return var_dict


def check_same_variant(sv_type, end1, end2):
    if sv_type == 'INS' or sv_type == 'DEL':
        return 0.7 < min(end1, end2) / max(end1, end2) <= 1
    return abs(end1 - end2) < 1000


# var_list read_id which is similar to (pos, sv_end), the reads support the variant
def find_in_list(var_type, var_list, bias, pos, sv_end):
    '''
    if interval_start == 44058795:
        print('start find in list')
        print('interval_start=%d, interval_end=%d, pos=%d, sv_end=%d'%(interval_start, interval_end, pos, sv_end))
    '''
    if len(var_list) == 0:
        return [], 0
    left = 0
    right = len(var_list) - 1
    mid = 0
    while left < right:
        mid = (left + right) >> 1
        if var_list[mid][1] < pos:
            left = mid + 1
        else:
            right = mid
    read_id_list = set()
    search_start = 0
    search_end = 0
    if right > 0 and pos - var_list[right - 1][1] <= 1000:
        for i in range(right - 1, -1, -1):
            if check_same_variant(var_type, var_list[i][2], sv_end):
                read_id_list.add(var_list[i][3])
                search_start = var_list[i][1]
            if i > 0 and var_list[i][1] - var_list[i - 1][1] > bias:
                break
    if var_list[right][1] - pos <= 1000:
        for i in range(right, len(var_list)):
            if check_same_variant(var_type, var_list[i][2], sv_end):  # if abs(var_list[i][2] - sv_end) < 1000:
                read_id_list.add(var_list[i][3])
                search_end = var_list[i][1]
            if i < len(var_list) - 1 and var_list[i + 1][1] - var_list[i][1] > bias:
                break
    search_threshold = min(abs(pos - search_start), abs(pos - search_end))
    return list(read_id_list), search_threshold


def find_in_indel_list(var_type, var_list, bias, pos, sv_end, threshold_gloab):
    if len(var_list) == 0:
        return [], 0
    left = 0
    right = len(var_list) - 1
    mid = 0
    while left < right:
        mid = (left + right) >> 1
        if var_list[mid][1] < pos:
            left = mid + 1
        else:
            right = mid
    candidates = []
    if right > 0 and pos - var_list[right - 1][1] <= 1000:
        for i in range(right - 1, -1, -1):
            candidates.append(var_list[i])
            if i > 0 and var_list[i][1] - var_list[i - 1][1] > bias and pos - var_list[i - 1][1] > 1000:
                break
    if var_list[right][1] - pos <= 1000:
        for i in range(right, len(var_list)):
            candidates.append(var_list[i])
            if i < len(var_list) - 1 and var_list[i + 1][1] - var_list[i][1] > bias and var_list[i + 1][1] - pos > 1000:
                break
    if len(candidates) == 0:
        return [], 0
    read_tag = dict()
    for element in candidates:
        if element[3] not in read_tag:
            read_tag[element[3]] = element
        else:
            if element[2] > read_tag[element[3]][2]:
                read_tag[element[3]] = element
    read_tag2SortedList = sorted(list(read_tag.values()), key = lambda x:x[2])
    global_len = [i[2] for i in read_tag2SortedList]
    DISCRETE_THRESHOLD_LEN_CLUSTER_TEMP = threshold_gloab * np.mean(global_len)

    last_len = read_tag2SortedList[0][2]
    allele_collect = list()
    allele_collect.append([[read_tag2SortedList[0][1]],  # start
                            [read_tag2SortedList[0][2]],  # len
                            [],  # support
                            [read_tag2SortedList[0][3]]])  # read_id
    for i in read_tag2SortedList[1:]:
        if i[2] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_TEMP:
            allele_collect[-1][2].append(len(allele_collect[-1][0]))
            allele_collect.append([[],[],[],[]])

        allele_collect[-1][0].append(i[1])
        allele_collect[-1][1].append(i[2])
        allele_collect[-1][3].append(i[3])
        last_len = i[2]
    allele_collect[-1][2].append(len(allele_collect[-1][0]))
    allele_sort = sorted(allele_collect, key = lambda x:x[2])

    allele_idx = 0
    nearest_gap = 0x3f3f3f3f
    for i in range(len(allele_sort)):
        allele = allele_sort[i]
        signalLen = np.mean(allele[1])
        if abs(signalLen - sv_end) < nearest_gap:
            allele_idx = i
            nearest_gap = abs(signalLen - sv_end)
    read_id_set = set(allele_sort[allele_idx][3])
    search_start = min(allele_sort[allele_idx][0])
    search_end = max(allele_sort[allele_idx][0])
    search_threshold = min(abs(pos - search_start), abs(pos - search_end))
    return list(read_id_set), search_threshold
        

def call_gt_wrapper(call_gt_args, gt_list, idx, row_count, para, strands, var_type):
    rname_list = []
    if var_type == 'INS' or var_type == 'DEL':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_indel(*call_gt_args)
        rname_list = call_gt_args[3]
    if var_type == 'DUP':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_dup(*call_gt_args)
        rname_list = call_gt_args[4]
    if var_type == 'INV':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_inv(*call_gt_args)
        rname_list = call_gt_args[4]
    if var_type == 'TRA':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_tra(*call_gt_args)
        rname_list = call_gt_args[5]
    rname = ','.join(rname_list)
    if rname == '':
        rname = 'NULL'
    gt_list[idx] = [para.chrom,
                    para.pos,
                    genotype,
                    var_type,
                    para.svlen,
                    para.end,
                    para.cipos,
                    para.ciend,
                    [gt_re, DR, GL, GQ, QUAL],
                    rname,
                    para.id,
                    para.ref, 
                    para.alts,
                    para.qual,
                    strands
    ]
    if idx > 0 and idx % 5000 == 0:
        logging.info(str(math.floor(idx / row_count * 100)) + '% SV calls of the given vcf has been processed.')


def force_calling(bam_path, ivcf_path, output_path, sigs_dir, max_cluster_bias_dict, threshold_gloab_dict, gt_round, threads):
    logging.info('Check the parameter -Ivcf: OK.')
    logging.info('Enable to perform force calling.')
    #print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    sv_dict = dict()
    #'''
    for sv_type in ["DEL", "INS", "DUP"]:
        sv_dict[sv_type] = parse_sigs(sv_type, sigs_dir)
    sv_dict['INV'] = parse_invsigs(sigs_dir)
    sv_dict['TRA'] = parse_trasigs(sigs_dir)
    #'''
    vcf_reader = VariantFile(ivcf_path, 'r')
    row_count = 0
    for record in vcf_reader.fetch():
        row_count += 1
    idx = -1
    gt_list = Manager().list([[] for x in range(row_count)])
    result = []
    process_pool = Pool(processes = threads)
    vcf_reader = VariantFile(ivcf_path, 'r')
    for record in vcf_reader.fetch():
        idx += 1
        sv_type, chrom, sv_chr2, pos, sv_end, sv_strand = parse_record(record)
        search_id_list = []
        if sv_type == 'TRA' and 'TRA' in sv_dict and chrom in sv_dict['TRA'] and sv_chr2 in sv_dict['TRA'][chrom]:
            search_id_list = sv_dict['TRA'][chrom][sv_chr2]
        elif sv_type == 'INV' and 'INV' in sv_dict and chrom in sv_dict['INV']:
            if sv_strand in sv_dict['INV'][chrom]:
                search_id_list = sv_dict['INV'][chrom][sv_strand]
            else:
                for strand_iter in sv_dict['INV'][chrom]:
                    sv_strand = strand_iter
                    search_id_list = sv_dict['INV'][chrom][strand_iter]
                    break
        elif sv_type != 'TRA' and sv_type != 'INV' and sv_type in sv_dict and chrom in sv_dict[sv_type]:
            search_id_list = sv_dict[sv_type][chrom]
        max_cluster_bias = 0
        if sv_type == 'INS' or sv_type == 'DEL':
            read_id_list, max_cluster_bias = find_in_indel_list(sv_type, search_id_list, max_cluster_bias_dict[sv_type], pos, sv_end, threshold_gloab_dict[sv_type])
        else:
            read_id_list, max_cluster_bias = find_in_list(sv_type, search_id_list, max_cluster_bias_dict[sv_type], pos, sv_end)
        if sv_type == 'INV' and len(read_id_list) == 0:
            for strand_iter in sv_dict['INV'][chrom]:
                if strand_iter != sv_strand:
                    search_id_list = sv_dict['INV'][chrom][strand_iter]
                    read_id_list, max_cluster_bias = find_in_list(sv_type, search_id_list, max_cluster_bias_dict[sv_type], pos, sv_end)
                    if len(read_id_list) != 0:
                        sv_strand = strand_iter
                        break
        #print(read_id_list)
        if sv_type == 'INS':
            max_cluster_bias = max(1000, max_cluster_bias)
        else:
            max_cluster_bias = max(max_cluster_bias_dict[sv_type], max_cluster_bias)

        para = Para(record)
        '''
        if sv_type == 'INS':
            call_gt_wrapper([bam_path, pos, chrom, read_id_list, max_cluster_bias, gt_round], gt_list, idx, row_count, para, sv_strand, 'INS')
        if sv_type == 'DEL':
            call_gt_wrapper([bam_path, pos, chrom, read_id_list, max_cluster_bias, gt_round], gt_list, idx, row_count, para, sv_strand, 'DEL')
        if sv_type == 'INV':
            call_gt_wrapper([bam_path, pos, sv_end, chrom, read_id_list, max_cluster_bias, gt_round], gt_list, idx, row_count, para, sv_strand, 'INV')
        if sv_type == 'DUP':
            call_gt_wrapper([bam_path, pos, sv_end, chrom, read_id_list, max_cluster_bias, gt_round], gt_list, idx, row_count, para, sv_strand, 'DUP')
        if sv_type == 'TRA':
            call_gt_wrapper([bam_path, pos, sv_end, chrom, sv_chr2, read_id_list, max_cluster_bias, gt_round], gt_list, idx, row_count, para, sv_strand, 'TRA')
        '''
        #'''
        if sv_type == 'INS':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, chrom, read_id_list, max_cluster_bias, gt_round], gt_list, idx, row_count, para, sv_strand, 'INS'))
        if sv_type == 'DEL':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, chrom, read_id_list, max_cluster_bias, gt_round], gt_list, idx, row_count, para, sv_strand, 'DEL'))
        if sv_type == 'INV':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, sv_end, chrom, read_id_list, max_cluster_bias, gt_round], gt_list, idx, row_count, para, sv_strand, 'INV'))
        if sv_type == 'DUP':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, sv_end, chrom, read_id_list, max_cluster_bias, gt_round], gt_list, idx, row_count, para, sv_strand, 'DUP'))
        if sv_type == 'TRA':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, sv_end, chrom, sv_chr2, read_id_list, max_cluster_bias, gt_round], gt_list, idx, row_count, para, sv_strand, 'TRA'))
        #'''
    process_pool.close()
    process_pool.join()
    logging.info('Finished force calling.')
    return gt_list


def run_fc(args):
    return force_calling(*args)

