
from cuteSV.cuteSV_resolveINV import call_gt as call_gt_inv
from cuteSV.cuteSV_resolveTRA import call_gt as call_gt_tra
from cuteSV.cuteSV_resolveINDEL import call_gt as call_gt_indel
from cuteSV.cuteSV_resolveDUP import call_gt as call_gt_dup
from multiprocessing import Pool, Manager
import vcf
from pysam import VariantFile
import math
import time
import logging


class Para(object):
    def __init__(self, record):
        self.chrom = record.chrom
        self.pos = parse_to_int(record.pos)
        self.svlen = parse_to_int(record.info['SVLEN']) if 'SVLEN' in record.info else 0
        self.end = parse_to_int(record.stop)
        self.cipos = record.info['CIPOS'] if 'CIPOS' in record.info else (0, 0)
        self.ciend = record.info['CIEND'] if 'CIEND' in record.info else (0, 0)
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
    else:
        return sth


def parse_record(record):
    sv_type = record.info['SVTYPE']
    chrom1 = record.chrom
    start = parse_to_int(record.pos)
    if record.info['SVTYPE'] == 'INS' and 'SVLEN' in record.info:
        end = parse_to_int(record.info['SVLEN'])
    else:
        try:
            end = parse_to_int(record.info['END'])
        except:
            try:
                end = parse_to_int(record.stop)
            except:
                pass   
    if record.info['SVTYPE'] == 'BND' or record.info['SVTYPE'] == 'TRA':
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
    if record.info['SVTYPE'] != 'TRA' and record.info['SVTYPE'] != 'BND':
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
                if var_type == 'DEL':
                    var_dict[seq[1]].append([seq[1], int(seq[2]), int(seq[2]) + int(seq[3]), seq[4]])
                else:
                    var_dict[seq[1]].append([seq[1], int(seq[2]), int(seq[3]), seq[4]])
            else:
                var_dict[seq[1]] = []
                if var_type == 'DEL':
                    var_dict[seq[1]].append([seq[1], int(seq[2]), int(seq[2]) + int(seq[3]), seq[4]])
                else:
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
    #var_dict[chrom1][tra_type][chrom2] = [[chrom2, start, end, read_id]]
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
    if sv_type == 'INS':
        return 0.7 < min(end1, end2) / max(end1, end2) <= 1
    return abs(end1 - end2) < 1000


# var_list read_id which is similar to (pos, sv_end) in [interval_start, interval_end], the reads support the variant
def find_in_list(var_type, var_list, bias, pos, sv_end):
    '''
    if interval_start == 44058795:
        print('start find in list')
        print('interval_start=%d, interval_end=%d, pos=%d, sv_end=%d'%(interval_start, interval_end, pos, sv_end))
    '''
    if len(var_list) == 0:
        return []
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
    if right > 0 and pos - var_list[right - 1][1] <= bias:
        for i in range(right - 1, 0, -1):
            if check_same_variant(var_type, var_list[i][2], sv_end):
                read_id_list.add(var_list[i][3])
            if i > 0 and var_list[i][1] - var_list[i - 1][1] > bias:
                break
    if var_list[right][1] - pos <= bias:
        for i in range(right, len(var_list)):
            if check_same_variant(var_type, var_list[i][2], sv_end):  # if abs(var_list[i][2] - sv_end) < 1000:
                read_id_list.add(var_list[i][3])
            if var_list[i][1] - var_list[i - 1][1] > bias:
                break
    return list(read_id_list)
        

def call_gt_wrapper(call_gt_args, gt_list, idx, row_count, para, strands, var_type):
    if var_type == 'INS' or var_type == 'DEL':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_indel(*call_gt_args)
    if var_type == 'DUP':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_dup(*call_gt_args)
    if var_type == 'INV':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_inv(*call_gt_args)
    if var_type == 'TRA':
        gt_re, DR, genotype, GL, GQ, QUAL = call_gt_tra(*call_gt_args)
    gt_list[idx] = [para.chrom,
                    para.pos,
                    genotype,
                    var_type,
                    para.svlen,
                    para.end,
                    para.cipos,
                    para.ciend,
                    [gt_re, DR, GL, GQ, QUAL],
                    ','.join(call_gt_args[3]),
                    para.id,
                    para.ref, 
                    para.alts,
                    para.qual,
                    strands
    ]
    if idx > 0 and idx % 5000 == 0:
        logging.info(str(math.floor(idx / row_count * 100)) + '% SV calls of the given vcf has been processed.')


def force_calling(bam_path, ivcf_path, output_path, sigs_dir, max_cluster_bias_dict, gt_round, threads):
    logging.info('Check the parameter -Ivcf: OK.')
    logging.info('Enable to perform force calling.')
    #print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    sv_dict = dict()
    for sv_type in ["DEL", "INS", "DUP"]:
        sv_dict[sv_type] = parse_sigs(sv_type, sigs_dir)
    sv_dict['INV'] = parse_invsigs(sigs_dir)
    sv_dict['TRA'] = parse_trasigs(sigs_dir)
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
        read_id_list = find_in_list(sv_type, search_id_list, max_cluster_bias_dict[sv_type], pos, sv_end)
        if sv_type == 'INV' and len(read_id_list) == 0:
            for strand_iter in sv_dict['INV'][chrom]:
                if strand_iter != sv_strand:
                    search_id_list = sv_dict['INV'][chrom][strand_iter]
                    read_id_list = find_in_list(sv_type, search_id_list, max_cluster_bias_dict[sv_type], pos, sv_end)
                    if len(read_id_list) != 0:
                        break
        #print(read_id_list)
        '''
        if sv_type == 'INS':
            call_gt_wrapper([bam_path, pos, chrom, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, record, sv_strand, 'INS')
        if sv_type == 'DEL':
            call_gt_wrapper([bam_path, pos, chrom, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, record, sv_strand, 'DEL')
        if sv_type == 'INV':
            call_gt_wrapper([bam_path, pos, sv_end, chrom, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, record, sv_strand, 'INV')
        if sv_type == 'DUP':
            call_gt_wrapper([bam_path, pos, sv_end, chrom, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, record, sv_strand, 'DUP')
        if sv_type == 'TRA':
            call_gt_wrapper([bam_path, pos, sv_end, chrom, sv_chr2, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, record, sv_strand, 'TRA')
        '''
        #'''
        para = Para(record)
        if sv_type == 'INS':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, chrom, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, para, sv_strand, 'INS'))
        if sv_type == 'DEL':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, chrom, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, para, sv_strand, 'DEL'))
        if sv_type == 'INV':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, sv_end, chrom, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, para, sv_strand, 'INV'))
        if sv_type == 'DUP':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, sv_end, chrom, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, para, sv_strand, 'DUP'))
        if sv_type == 'TRA':
            process_pool.apply_async(call_gt_wrapper, 
                args=([bam_path, pos, sv_end, chrom, sv_chr2, read_id_list, max_cluster_bias_dict[sv_type], gt_round], gt_list, idx, row_count, para, sv_strand, 'TRA'))
        #'''
    process_pool.close()
    process_pool.join()
    logging.info('Finished force calling.')
    return gt_list


def run_fc(args):
    return force_calling(*args)

