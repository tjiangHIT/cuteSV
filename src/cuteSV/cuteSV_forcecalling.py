from cuteSV.cuteSV_genotype import cal_CIPOS, overlap_cover, assign_gt_fc
from multiprocessing import Pool
from pysam import VariantFile
import math
import time
import logging
import numpy as np
import pickle
from sklearn.cluster import KMeans

def parse_svtype(sv_type):
    if 'DEL' in sv_type:
        return 'DEL'
    if 'INS' in sv_type:
        return 'INS'
    if 'INV' in sv_type:
        return 'INV'
    if 'DUP' in sv_type:
        return 'DUP'
    if 'TRA' in sv_type:
        return 'TRA'
    if 'BND' in sv_type:
        return 'BND'
    return 'NA'

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

# INDEL end->len, others:end
def parse_record(record):
    sv_type = parse_svtype(record.info['SVTYPE'])
    chrom1 = record.chrom
    start = parse_to_int(record.pos)
    chrom2 = record.chrom
    end = None
    ref = record.ref
    alts = record.alts[0]
    # read sv length from 1.SVLEN= 2.alt-ref 3.record as 0
    if 'SVLEN' in record.info:
        svlen = abs(parse_to_int(record.info['SVLEN']))
    elif alts[0] != '<' and (sv_type != 'TRA' and sv_type != 'BND'):
        svlen = abs(len(alts) - len(ref))
    else:
        svlen = 0
    # read end from 1.BND-alt 2.END= 3.start+svlen
    if sv_type == 'TRA' or sv_type == 'BND':
        try:
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
            if ':' in tra_alt:
                chrom2 = tra_alt.split(':')[0]
                end = int(tra_alt.split(':')[1])
        except:
            pass
    if end == None:
        if 'END' in record.info:
            end = parse_to_int(record.info['END'])
        else:
            end = start + svlen
    # chrom2 for BND
    if 'CHR2' in record.info:
        chrom2 = record.info['CHR2']
    # strand
    strand = '.'
    if 'STRAND' in record.info:
        strand = record.info['STRAND']  
    elif 'STRANDS' in record.info:
        strand = record.info['STRANDS']
    if isinstance(strand, tuple) or isinstance(strand, list):
        strand = strand[0]
    svid = record.id
    # update ref and alt with SEQ
    if 'SEQ' in record.info:
        if record.info['SVTYPE'] == 'INS' and record.alts[0] == '<INS>':
            alts = record.info['SEQ']
        if record.info['SVTYPE'] == 'DEL' and record.alts[0] == '<DEL>':
            ref = record.info['SEQ']
    return sv_type, chrom1, chrom2, start, end, svlen, strand, svid, ref, alts

def parse_sigs_chrom(var_type, work_dir, chrom_list, index):
    var_dict={}
    if False and var_type != 'TRA':
        with open(work_dir + var_type + '.pickle', 'rb') as f:
            for chrom in chrom_list:
                if chrom not in index[var_type].keys():
                    continue
                f.seek(index[var_type][chrom])
                var_dict[chrom]=pickle.load(f)
    else:
        with open(work_dir + var_type + '.pickle', 'rb') as f:
            for chrom in chrom_list:
                if chrom not in index[var_type].keys():
                    continue
                f.seek(index[var_type][chrom])
                sigs=pickle.load(f)
                #file_sig[]=mem_sig[DEL: -2,-1,0,1,2, INS: -2,-1,0,1,2,3, DUP: -2,-1,0,1,2, INV: -2,-1,0,1,2,3, TRA: -2,-1,0,1,2,3,4, reads: -1,0,1,2,3]
                #from file: DEL,DUP:  1,2,3,4, INS: 1,2,3,4,5, INV:1,3,4,5, TRA[1][4]:4,3,5,6
                #from mem: DEL,DUP: -1,0,1,2, INS: -1,0,1,2,3, INV: -1,1,2,3, TRA[-1][2]:2,1,3,4
                if var_type == 'DEL' or var_type == 'DUP':
                    for seq in sigs:
                        if chrom not in var_dict:
                            var_dict[chrom] = []
                        var_dict[chrom].append([seq[-1], int(seq[0]), int(seq[1]), seq[2]])
                elif var_type == 'INS':
                    for seq in sigs:
                        if chrom not in var_dict:
                            var_dict[chrom] = []
                        if len(seq) < 6:
                            cigar = '<INS>'
                        else:
                            cigar = seq[3]
                        cigar = '<INS>'
                        var_dict[chrom].append([seq[-1], int(seq[0]), int(seq[1]), seq[2], cigar])
                elif var_type == 'INV':
                    for seq in sigs:
                        if chrom not in var_dict:
                            var_dict[chrom] = []
                        var_dict[chrom].append([seq[-1], int(seq[1]), int(seq[2]), seq[3]])
                else:
                    for seq in sigs:
                        chrom1 = seq[-1]
                        tra_type = seq[0]
                        pos1 = int(seq[1])
                        chrom2 = seq[2]
                        pos2 = int(seq[3])
                        read_id = seq[4]
                        if chrom1 not in var_dict:
                            var_dict[chrom1] = dict()
                        if chrom2 not in var_dict[chrom1]:
                            var_dict[chrom1][chrom2] = []
                        var_dict[chrom1][chrom2].append([chrom2, pos1, pos2, read_id])
                    for chr1 in var_dict:
                        for chr2 in var_dict[chr1]:
                            var_dict[chr1][chr2].sort(key=lambda x:x[1])
            return var_dict

def check_same_variant(sv_type, end1, end2, bias):
    if sv_type == 'INS' or sv_type == 'DEL':
        return 0.7 < min(end1, end2) / max(end1, end2) <= 1
    return abs(end1 - end2) < bias

# var_list read_id which is similar to (pos, sv_end), the reads support the variant
def find_in_list(var_type, var_list, bias, pos, sv_end):
    if len(var_list) == 0:
        return [], pos, pos
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
    search_start = -1
    search_end = -1
    if right > 0 and pos - var_list[right - 1][1] <= bias: ##
        for i in range(right - 1, -1, -1):
            if check_same_variant(var_type, var_list[i][2], sv_end, bias):
                read_id_list.add(var_list[i][3])
                search_start = var_list[i][1]
            if i > 0 and (var_list[i][1] - var_list[i - 1][1] > bias or pos - var_list[i - 1][1] > bias):
                break
    if var_list[right][1] - pos <= bias:
        for i in range(right, len(var_list)):
            if check_same_variant(var_type, var_list[i][2], sv_end, bias):  # if abs(var_list[i][2] - sv_end) < 1000:
                read_id_list.add(var_list[i][3])
                search_end = var_list[i][1]
            if i < len(var_list) - 1 and (var_list[i + 1][1] - var_list[i][1] > bias or var_list[i + 1][1] - pos > bias):
                break
    if search_start == -1:
        search_start = pos
    if search_end == -1:
        search_end = pos
    search_threshold = max(abs(pos - search_start), abs(pos - search_end))
    if search_start > search_end:
        search_start, search_end = search_end, search_start
    if search_start == search_end:
        search_end += 1
    return list(read_id_list), search_start, search_end

def find_in_indel_list(var_type, var_list, bias_origin, pos, sv_end, threshold_gloab, multi_allele):
    # bias = min(bias_origin * 10, 2000)
    bias = bias_origin
    debug_pos = -1
    if len(var_list) == 0:
        return [], pos, pos, '.,.', '.,.'
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
    if right > 0 and pos - var_list[right - 1][1] <= bias:
        for i in range(right - 1, -1, -1):
            candidates.append(var_list[i]) # [chrom, start, len, read_id] or [chrom, start, len, read_id, seq]
            if i > 0 and (var_list[i][1] - var_list[i - 1][1] > bias or pos - var_list[i - 1][1] > 2 * bias):
                break
    if var_list[right][1] - pos <= bias:
        for i in range(right, len(var_list)):
            candidates.append(var_list[i])
            if i < len(var_list) - 1 and (var_list[i + 1][1] - var_list[i][1] > bias or var_list[i + 1][1] - pos > 2 * bias):
                break
    if len(candidates) == 0:
        return [], pos, pos, '.,.', '.,.'
    read_tag = dict()
    for element in candidates:
        if element[3] not in read_tag:
            read_tag[element[3]] = []
        read_tag[element[3]].append(element)

    # merge sigs on the same read
    read_tag2SortedList = []

    for read_id in read_tag: # 40
        for i in range(len(read_tag[read_id])):
            read_tag2SortedList.append(read_tag[read_id][i])
            if i + 1 < len(read_tag[read_id]):
                j = i + 1
                if var_type == 'DEL':
                    read_tag2SortedList.append(
                        [read_tag[read_id][i][0], int((read_tag[read_id][i][1] + read_tag[read_id][j][1]) / 2),
                         read_tag[read_id][i][2] + read_tag[read_id][j][2], read_tag[read_id][i][3]])
                else:
                    read_tag2SortedList.append(
                        [read_tag[read_id][i][0], int((read_tag[read_id][i][1] + read_tag[read_id][j][1]) / 2),
                         read_tag[read_id][i][2] + read_tag[read_id][j][2], read_tag[read_id][i][3],
                         read_tag[read_id][i][4]])
                if j + 1 < len(read_tag[read_id]):
                    k = j + 1
                    if var_type == 'DEL':
                        read_tag2SortedList.append([read_tag[read_id][i][0], int((read_tag[read_id][i][1] +
                                                                                read_tag[read_id][j][1] +
                                                                                read_tag[read_id][k][1]) / 3),
                                                    read_tag[read_id][i][2] + read_tag[read_id][j][2] +
                                                    read_tag[read_id][k][2], read_tag[read_id][i][3]])
                    else:
                        read_tag2SortedList.append([read_tag[read_id][i][0], int((read_tag[read_id][i][1] +
                                                                                read_tag[read_id][j][1] +
                                                                                read_tag[read_id][k][1]) / 3),
                                                    read_tag[read_id][i][2] + read_tag[read_id][j][2] +
                                                    read_tag[read_id][k][2], read_tag[read_id][i][3],
                                                    read_tag[read_id][i][4]])

    read_tag2SortedList = sorted(read_tag2SortedList, key = lambda x:x[2])
    if pos == debug_pos:
        print(read_tag2SortedList)
    global_len = [i[2] for i in read_tag2SortedList]
    DISCRETE_THRESHOLD_LEN_CLUSTER_TEMP = threshold_gloab * np.mean(global_len)
    if var_type == 'DEL':
        last_len = read_tag2SortedList[0][2]
        cur_bias = last_len * threshold_gloab
        allele_collect = list()
        allele_collect.append([[read_tag2SortedList[0][1]],  # start
                                [read_tag2SortedList[0][2]],  # len
                                [],  # support
                                [read_tag2SortedList[0][3]]])  # read_id
        for i in read_tag2SortedList[1:]:
            if i[2] - last_len > cur_bias:
                allele_collect[-1][2].append(len(allele_collect[-1][0]))
                allele_collect.append([[],[],[],[]])
            allele_collect[-1][0].append(i[1])
            allele_collect[-1][1].append(i[2])
            allele_collect[-1][3].append(i[3])
            last_len = (last_len * (len(allele_collect[-1][0]) - 1) + i[2]) / len(allele_collect[-1][0])
            cur_bias = last_len * threshold_gloab
        if pos == debug_pos:
            print(allele_collect)
            print('above is all allele collect')
        allele_collect[-1][2].append(len(allele_collect[-1][0]))
        allele_idx = -1
        nearest_gap = 0x3f3f3f3f
        for i in range(len(allele_collect)):
            allele = allele_collect[i]
            signalLen = np.mean(allele[1])
            if min(signalLen, sv_end) / max(signalLen, sv_end) > 0.7:
                if abs(signalLen - sv_end) < nearest_gap:
                    allele_idx = i
                    nearest_gap = abs(signalLen - sv_end)
        if allele_idx == -1:
            lower_length = sv_end * 0.7
            upper_length = sv_end / 0.7
            final_alleles = [[],[],[],[]]
            for i in range(len(allele_collect[allele_idx][0])):
                if lower_length <= allele_collect[allele_idx][1][i] <= upper_length:
                    final_alleles[0].append(allele_collect[allele_idx][0][i])
                    final_alleles[1].append(allele_collect[allele_idx][1][i])
                    final_alleles[3].append(allele_collect[allele_idx][3][i])
        else:
            final_alleles = allele_collect[allele_idx]

        if pos == debug_pos:
            print(allele_idx)
            print(final_alleles[:2])
        if multi_allele:
            bimodal_data = np.array(final_alleles[1])
            if len(bimodal_data) > 1 and final_alleles[1][0] != final_alleles[1][-1]:
                model = KMeans(n_clusters=2, init=np.array([int(len(bimodal_data)/4), int(len(bimodal_data)/4*3)]).reshape(-1, 1), n_init=1)
                # try:
                model.fit(bimodal_data.reshape(-1, 1))
                # except:
                #     print('reshape error')
                labels = model.labels_
                if pos == debug_pos:
                    print(labels)
                cate = 0
                for i in range(0, len(labels) - 1):
                    if labels[i] != labels[i + 1]:
                        cate = i + 1
                        break
                delta0 = math.ceil(cate / 8) if cate >= 3 else 0
                delta1 = math.ceil((len(labels) - cate + 1) / 8) if len(labels) - cate >= 3 else 0
                min_alleles = [final_alleles[1][delta0], final_alleles[1][cate + delta1]]
                max_alleles = [final_alleles[1][cate - delta0 - 1], final_alleles[1][len(labels) - delta1 - 1]]
                final_alleles_filter = [[],[],[],[]]
                if abs(max_alleles[0]-max_alleles[1]) >= max(3*max(max_alleles[0]-min_alleles[0], max_alleles[1]-min_alleles[1]), 6):
                    allele0 = np.mean(final_alleles[1][delta0:(cate-delta0)])
                    allele1 = np.mean(final_alleles[1][cate+delta1:]) if delta1 == 0 else np.mean(final_alleles[1][cate+delta1:-delta1])
                    if min(allele0, sv_end) / max(allele0, sv_end) >= min(allele1, sv_end) / max(allele1, sv_end): # choose front
                        if min(min_alleles[0], sv_end) / max(min_alleles[0], sv_end) > 0.9 and min(max_alleles[0], sv_end) / max(max_alleles[0], sv_end) > 0.9:
                            if cate >= max(3, len(labels) / 5):
                                for i in range(cate):
                                    # if min(final_alleles[1][i], sv_end) / max(final_alleles[1][i], sv_end) > 0.9:
                                    for j in [0, 1, 3]:
                                        final_alleles_filter[j].append(final_alleles[j][i])
                    elif min(min_alleles[1], sv_end) / max(min_alleles[1], sv_end) > 0.9 and min(max_alleles[1], sv_end) / max(max_alleles[1], sv_end) > 0.9:
                        if len(labels) - cate >= max(3, len(labels) / 5):
                            for i in range(cate, len(labels)):
                                # if min(final_alleles[1][i], sv_end) / max(final_alleles[1][i], sv_end) > 0.9:
                                for j in [0, 1, 3]:
                                    final_alleles_filter[j].append(final_alleles[j][i])
                if len(final_alleles_filter[0]) > 0:
                    final_alleles = final_alleles_filter
                if pos == debug_pos:
                    print('final alleles:')
                    print(final_alleles[:2])
        if len(final_alleles[3]) > 0:
            read_id_set = set(final_alleles[3])
            CIPOS = cal_CIPOS(np.std(final_alleles[0]), len(final_alleles[0]))
            CILEN = cal_CIPOS(np.std(final_alleles[1]), len(final_alleles[1]))
            seq = '<DEL>'
            # search_start = np.median(final_alleles[0])
            # search_end = np.median(final_alleles[0])
            search_start = min(final_alleles[0])
            search_end = max(final_alleles[0])
            search_threshold = min(abs(pos - search_start), abs(pos - search_end))
            if pos == debug_pos:
                print('search_start:%f'%(search_start))
                print(final_alleles[0])
        else:
            read_id_set = set()
            CIPOS = "-0,0"
            CILEN = "-0,0"
            seq = '<DEL>'
            search_start = pos
            search_end = pos
            search_threshold = 0
    if var_type == 'INS':
        last_len = read_tag2SortedList[0][2]
        cur_bias = last_len * threshold_gloab
        allele_collect = list()
        allele_collect.append([[read_tag2SortedList[0][1]],  # start
                                [read_tag2SortedList[0][2]],  # len
                                [],  # support
                                [read_tag2SortedList[0][3]],  # read_id
                                [read_tag2SortedList[0][4]]])  # ins_seq
        for i in read_tag2SortedList[1:]:
            if i[2] - last_len > cur_bias:
                allele_collect[-1][2].append(len(allele_collect[-1][0]))
                allele_collect.append([[],[],[],[],[]])
            allele_collect[-1][0].append(i[1])
            allele_collect[-1][1].append(i[2])
            allele_collect[-1][3].append(i[3])
            allele_collect[-1][4].append(i[4])
            last_len = (last_len * (len(allele_collect[-1][0]) - 1) + i[2]) / len(allele_collect[-1][0])
            cur_bias = last_len * threshold_gloab
        if pos == debug_pos:
            print(allele_collect)
            print('above is all allele collect')
        allele_collect[-1][2].append(len(allele_collect[-1][0]))
        allele_idx = -1
        nearest_gap = 0x3f3f3f3f
        for i in range(len(allele_collect)):
            allele = allele_collect[i]
            signalLen = np.mean(allele[1])
            if min(signalLen, sv_end) / max(signalLen, sv_end) > 0.7:
                if abs(signalLen - sv_end) < nearest_gap:
                    allele_idx = i
                    nearest_gap = abs(signalLen - sv_end)
        if allele_idx == -1:
            lower_length = sv_end * 0.7
            upper_length = sv_end / 0.7
            final_alleles = [[],[],[],[]]
            for i in range(len(allele_collect[allele_idx][0])):
                if lower_length <= allele_collect[allele_idx][1][i] <= upper_length:
                    final_alleles[0].append(allele_collect[allele_idx][0][i])
                    final_alleles[1].append(allele_collect[allele_idx][1][i])
                    final_alleles[3].append(allele_collect[allele_idx][3][i])
                    # final_alleles[4].append(allele_collect[allele_idx][4][i])
        else:
            final_alleles = allele_collect[allele_idx]

        if pos == debug_pos:
            print(allele_idx)
            print(final_alleles[:2])
        if multi_allele:
            bimodal_data = np.array(final_alleles[1])
            if len(bimodal_data) > 1 and final_alleles[1][0] != final_alleles[1][-1]:
                model = KMeans(n_clusters=2, init=np.array([int(len(bimodal_data)/4), int(len(bimodal_data)/4*3)]).reshape(-1, 1), n_init=1)
                # try:
                model.fit(bimodal_data.reshape(-1, 1))
                # except:
                #     print('reshape error')
                labels = model.labels_
                if pos == debug_pos:
                    print(labels)
                cate = 0
                for i in range(0, len(labels) - 1):
                    if labels[i] != labels[i + 1]:
                        cate = i + 1
                        break
                delta0 = math.ceil(cate / 8) if cate >= 5 else 0
                delta1 = math.ceil((len(labels) - cate) / 8) if len(labels) - cate >= 5 else 0
                min_alleles = [final_alleles[1][delta0], final_alleles[1][cate + delta1]]
                max_alleles = [final_alleles[1][cate - delta0 - 1], final_alleles[1][len(labels) - delta1 - 1]]
                final_alleles_filter = [[],[],[],[]]
                if abs(max_alleles[0]-max_alleles[1]) >= max(3*max(max_alleles[0]-min_alleles[0], max_alleles[1]-min_alleles[1]), 6): # or > 3
                    allele0 = np.mean(final_alleles[1][delta0:(cate-delta0)])
                    allele1 = np.mean(final_alleles[1][cate+delta1:]) if delta1 == 0 else np.mean(final_alleles[1][cate+delta1:-delta1])
                    if min(allele0, sv_end) / max(allele0, sv_end) >= min(allele1, sv_end) / max(allele1, sv_end): # choose front
                        if min(min_alleles[0], sv_end) / max(min_alleles[0], sv_end) > 0.9 and min(max_alleles[0], sv_end) / max(max_alleles[0], sv_end) > 0.9:
                            if cate >= max(3, len(labels) / 5):
                                for i in range(cate):
                                    # if min(final_alleles[1][i], sv_end) / max(final_alleles[1][i], sv_end) > 0.9:
                                    for j in [0, 1, 3]:
                                        final_alleles_filter[j].append(final_alleles[j][i])
                    elif min(min_alleles[1], sv_end) / max(min_alleles[1], sv_end) > 0.9 and min(max_alleles[1], sv_end) / max(max_alleles[1], sv_end) > 0.9:
                        if len(labels) - cate >= max(3, len(labels) / 5):
                            for i in range(cate, len(labels)):
                                # if min(final_alleles[1][i], sv_end) / max(final_alleles[1][i], sv_end) > 0.9:
                                for j in [0, 1, 3]:
                                    final_alleles_filter[j].append(final_alleles[j][i])
                if len(final_alleles_filter[0]) > 0:
                    final_alleles = final_alleles_filter
                if pos == debug_pos:
                    print('final alleles:')
                    print(final_alleles[:2])
        if len(final_alleles[3]) > 0:
            read_id_set = set(final_alleles[3])
            CIPOS = cal_CIPOS(np.std(final_alleles[0]), len(final_alleles[0]))
            CILEN = cal_CIPOS(np.std(final_alleles[1]), len(final_alleles[1]))
            seq = '<INS>'
            # search_start = np.median(final_alleles[0])
            # search_end = np.median(final_alleles[0])
            search_start = min(final_alleles[0])
            search_end = max(final_alleles[0])
            search_threshold = min(abs(pos - search_start), abs(pos - search_end))
        else:
            read_id_set = set()
            CIPOS = "-0,0"
            CILEN = "-0,0"
            seq = '<INS>'
            search_start = pos
            search_end = pos
            search_threshold = 0
    # return list(read_id_set), search_threshold, CIPOS, CILEN
    return list(read_id_set), search_start, search_end, CIPOS, CILEN

def generate_dispatch(reads_count, chrom_list):
    dispatch = [[]]
    cur_count = 0
    for item in reads_count:
        if cur_count >= 10000:
            dispatch.append([])
            cur_count = 0
        cur_count = item[1]
        dispatch[-1].append(item[0])
    # chrom without reads
    for chrom in chrom_list:
        flag = 0
        for x in reads_count:
            if chrom == x[0]:
                flag = 1
        if flag == 0:
            dispatch[0].append(chrom)
    return dispatch


def force_calling_chrom(ivcf_path, temporary_dir, max_cluster_bias_dict, threshold_gloab_dict, gt_round, read_range, threads, sigs_index):
    logging.info('Check the parameter -Ivcf: OK.')
    logging.info('Enable to perform force calling.')
    if sigs_index==None:
        with open("%s/sigindex.pickle"%temporary_dir,"rb") as f:
            sigs_index=pickle.load(f)
    # parse svs tobe genotyped
    vcf_reader = VariantFile(ivcf_path, 'r')
    svs_tobe_genotyped = dict()
    svs_pre = dict()
    svs_multi = dict()
    for record in vcf_reader.fetch():
        sv_type, chrom, sv_chr2, pos, sv_end, sv_len, sv_strand, svid, ref, alts = parse_record(record)
        if sv_type not in ["DEL", "INS", "DUP", "INV", "TRA", "BND"]:
            continue
        if chrom not in svs_tobe_genotyped:
            svs_tobe_genotyped[chrom] = list()
            svs_pre[chrom] = dict()
        svs_tobe_genotyped[chrom].append([sv_type, sv_chr2, pos, sv_end, sv_len, svid, ref, alts, sv_strand, chrom])
        if pos not in svs_pre[chrom]:
            svs_pre[chrom][pos] = 0
        svs_pre[chrom][pos] += 1
    for c in svs_pre:
        for s in svs_pre[c]:
            if svs_pre[c][s] == 2:
                if c not in svs_multi:
                    svs_multi[c] = set()
                svs_multi[c].add(s)
    start_time = time.time()
    # parse reads in alignment
    reads_count = sigs_index["reads_count"]
    reads_count = sorted(reads_count.items(), key=lambda x:x[1])
    dispatch = generate_dispatch(reads_count, svs_tobe_genotyped.keys())
    logging.info('finish dispatch in {}.'.format(time.time() - start_time))

    # force calling
    pool_result = list()
    result= {}
    process_pool = Pool(processes = threads)
    # dispatch = [['MT']]
    for chroms in dispatch:
        genotype_sv_list = dict()
        for chrom in chroms:
            if chrom in svs_tobe_genotyped:
                genotype_sv_list[chrom] = svs_tobe_genotyped[chrom]
        if len(genotype_sv_list) == 0:
            continue
        fx_para = [(chroms, genotype_sv_list, temporary_dir, max_cluster_bias_dict, threshold_gloab_dict, gt_round, sigs_index, read_range, svs_multi)]
        pool_result.append(process_pool.map_async(solve_fc_wrapper, fx_para))
    process_pool.close()
    process_pool.join()

    for x in pool_result:
        result.update(x.get()[0])
    return result

def solve_fc_wrapper(args):
    return solve_fc(*args)
def solve_fc(chrom_list, svs_dict, temporary_dir, max_cluster_bias_dict, threshold_gloab_dict, gt_round, sigs_index, read_range, svs_multi):
    reads_info = dict() # [10000, 10468, 0, 'm54238_180901_011437/52298335/ccs']
    readsfile = open("%sreads.pickle"%(temporary_dir), 'rb')
    for chrom in chrom_list:
        try:
            readsfile.seek(sigs_index["reads"][chrom])
            reads_info[chrom]=pickle.load(readsfile)
        except:
            reads_info[chrom] = []
    readsfile.close()
    sv_dict = dict()
    for sv_type in ["DEL", "DUP", "INS", "INV", "TRA"]:
        sv_dict[sv_type] = parse_sigs_chrom(sv_type, temporary_dir, chrom_list, sigs_index)
    
    gt_list = {}
    for chrom in svs_dict:
        gt_list[chrom]=[]
        start_time = time.time()
        read_id_dict = dict()
        svtype_id_dict = dict()
        ci_dict = dict()
        search_list = list()
        for i in range(len(svs_dict[chrom])):
            record = svs_dict[chrom][i]
            # [sv_type, sv_chr2, pos, sv_end, sv_len, svid, ref, alts, sv_strand, chrom]
            sv_type = record[0]
            sv_chr2 = record[1]
            sv_start = record[2]
            sv_end = record[3]
            sv_len = record[4]
            chrom = record[9]
            search_id_list = list()
            # rewrite!
            if (sv_type == 'TRA' or sv_type == 'BND') and 'TRA' in sv_dict and chrom in sv_dict['TRA'] and sv_chr2 in sv_dict['TRA'][chrom]:
                search_id_list = sv_dict['TRA'][chrom][sv_chr2]
            elif sv_type != 'TRA' and sv_type in sv_dict and chrom in sv_dict[sv_type]:
                search_id_list = sv_dict[sv_type][chrom]
            max_cluster_bias = 0
            if sv_type == 'INS' or sv_type == 'DEL':
                sigs_bias = max_cluster_bias_dict[sv_type]
                if chrom in svs_multi and sv_start in svs_multi[chrom]:
                    multi_allele = True
                else:
                    multi_allele = False
                read_id_list, search_start, search_end, CIPOS, CILEN = find_in_indel_list(sv_type, search_id_list, sigs_bias, sv_start, sv_len, threshold_gloab_dict[sv_type], multi_allele)
            else:
                sv_temp = sv_type
                sigs_bias = max_cluster_bias_dict[sv_type if sv_type != 'BND' else 'TRA']
                if sv_len / 2 > sigs_bias:
                    sigs_bias = sv_len / 2
                # read_id_list, max_cluster_bias = find_in_list(sv_type, search_id_list, sigs_bias, sv_start, sv_end)
                read_id_list, search_start, search_end = find_in_list(sv_type, search_id_list, sigs_bias, sv_start, sv_end)
                CIPOS = '.'
                CILEN = '.'
            max_cluster_bias = max(abs(sv_start - search_start), abs(sv_start - search_end))
            max_cluster_bias = max(read_range, max_cluster_bias)
            if sv_type == 'INS' or sv_type == 'TRA' or sv_type == 'BND':
                search_list.append((max(sv_start - max_cluster_bias, 0), sv_start + max_cluster_bias))
            elif sv_type == 'DEL':
                if read_range < 500:
                    search_list.append((max(sv_start - max_cluster_bias, 0), sv_start + max_cluster_bias))
                else:
                    search_list.append((max(sv_start + abs(sv_len) / 5, 0), sv_start + abs(sv_len) - abs(sv_len) / 5))
            elif sv_type == 'INV':
                search_list.append((search_start, search_end+1))
            elif sv_type == 'DUP':
                search_list.append((sv_start, sv_end))

            read_id_dict[i] = read_id_list
            svtype_id_dict[i] = sv_type
            ci_dict[i] = (CIPOS, CILEN)
        if chrom in reads_info:
            iteration_dict, primary_num_dict, cover_dict, overlap_dict = overlap_cover(search_list, reads_info[chrom]) # both key(sv idx), value(set(read id))
        else:
            iteration_dict = dict()
            primary_num_dict = dict()
            cover_dict = dict()
            overlap_dict = dict()
            for i in read_id_dict:
                iteration_dict[i] = 0
                primary_num_dict[i] = 0
                cover_dict[i] = set()
                overlap_dict[i] = set()
        assert len(iteration_dict) == len(read_id_dict), "overlap length error"
        assert len(cover_dict) == len(svs_dict[chrom]), "cover length error"
        assert len(overlap_dict) == len(svs_dict[chrom]), "overlap length error"
        assign_list = assign_gt_fc(iteration_dict, primary_num_dict, cover_dict, overlap_dict, read_id_dict, svtype_id_dict)
        for i in range(len(svs_dict[chrom])):
            assert len(assign_list[i]) == 6, "assign genotype error"
            record = svs_dict[chrom][i]
            rname = ','.join(read_id_dict[i])
            if rname == '':
                rname = 'NULL'
            if record[7] == '<TRA>' or record[7] == '<BND>':
                seq = str(record[1]) + ':' + str(record[3])
            else:
                seq = '<' + record[0] + '>'
            # [sv_type, sv_chr2, pos, sv_end, sv_len, svid, ref, alts, sv_strand, chrom]
            gt_list[record[9]].append([record[9], record[2], assign_list[i][2], record[0], record[3],
                        ci_dict[i][0], ci_dict[i][1], assign_list[i], rname, record[5],
                        record[6], record[7],
                        record[8], seq, record[4]])
        logging.info("Finished calling %s."%(chrom))
    return gt_list
