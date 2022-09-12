import numpy as np
import logging
from cuteSV.cuteSV_genotype import overlap_cover, assign_gt

'''
*******************************************
                TO DO LIST
*******************************************
    1. Identify DP with samfile pointer;
    2. Add CIPOS, CILEN and/or CIEND;
    3. Determine (IM)PRECISE type.
    4. Filter DUP to improve INS FN rate.
*******************************************
'''

def resolution_DUP(path, chr, read_count, max_cluster_bias, sv_size, 
    bam_path, action, MaxSize, gt_round):
    semi_dup_cluster = list()
    semi_dup_cluster.append([0, 0, ''])
    candidate_single_SV = list()

    file = open("%s%s.sigs"%(path, "DUP"), 'r')
    for line in file:
        seq = line.strip('\n').split('\t')
        if seq[1] != chr:
            continue

        pos_1 = int(seq[2])
        pos_2 = int(seq[3])
        read_id = seq[4]
        
        if pos_1 - semi_dup_cluster[-1][0] > max_cluster_bias or pos_2 - semi_dup_cluster[-1][1] > max_cluster_bias:
            if len(semi_dup_cluster) >= read_count:
                if semi_dup_cluster[-1][0] == semi_dup_cluster[-1][1] == 0:
                    pass
                else:
                    generate_dup_cluster(semi_dup_cluster, 
                                        chr, 
                                        read_count, 
                                        max_cluster_bias, 
                                        sv_size, 
                                        candidate_single_SV,
                                        action,
                                        MaxSize,
                                        gt_round)
            semi_dup_cluster = []
            semi_dup_cluster.append([pos_1, pos_2, read_id])
        else:
            if semi_dup_cluster[-1][0] == semi_dup_cluster[-1][1] == 0:
                semi_dup_cluster = []
                semi_dup_cluster.append([pos_1, pos_2, read_id])
            else:
                semi_dup_cluster.append([pos_1, pos_2, read_id])

    if len(semi_dup_cluster) >= read_count:
        if semi_dup_cluster[-1][0] == semi_dup_cluster[-1][1] == 0:
            pass
        else:
            generate_dup_cluster(semi_dup_cluster, 
                                chr, 
                                read_count, 
                                max_cluster_bias, 
                                sv_size, 
                                candidate_single_SV,
                                action,
                                MaxSize,
                                gt_round)
    file.close()
    if action:
        candidate_single_SV_gt = call_gt(path, chr, candidate_single_SV, max_cluster_bias)
        logging.info("Finished %s:%s."%(chr, "DUP"))
        return candidate_single_SV_gt
    else:
        logging.info("Finished %s:%s."%(chr, "DUP"))
        return candidate_single_SV

def generate_dup_cluster(semi_dup_cluster, chr, read_count, max_cluster_bias, 
    sv_size, candidate_single_SV, action, MaxSize, gt_round):
    # calculate support reads
    support_read = list(set([i[2] for i in semi_dup_cluster]))
    if len(support_read) < read_count:
        return

    low_b = int(len(semi_dup_cluster)*0.4)
    up_b = int(len(semi_dup_cluster)*0.6)

    if low_b == up_b:
        breakpoint_1 = semi_dup_cluster[low_b][0]
        breakpoint_2 = semi_dup_cluster[low_b][1]
    else:
        breakpoint_1 = [i[0] for i in semi_dup_cluster[low_b:up_b]]
        breakpoint_2 = [i[1] for i in semi_dup_cluster[low_b:up_b]]
        breakpoint_1 = int(sum(breakpoint_1)/len(semi_dup_cluster[low_b:up_b]))
        breakpoint_2 = int(sum(breakpoint_2)/len(semi_dup_cluster[low_b:up_b]))


    if sv_size <= breakpoint_2 - breakpoint_1 <= MaxSize or (sv_size <= breakpoint_2 - breakpoint_1 and MaxSize == -1):
        if action:
            candidate_single_SV.append([chr,
                                        'DUP', 
                                        breakpoint_1, 
                                        breakpoint_2,
                                        support_read])
            # print("DUP", chr, int(breakpoint_1), int(breakpoint_2), DR, DV, QUAL, "%.4f"%cost_time)
        else:
            candidate_single_SV.append([chr,
                                        'DUP', 
                                        str(breakpoint_1), 
                                        str(breakpoint_2 - breakpoint_1), 
                                        str(len(support_read)),
                                        '.',
                                        './.',
                                        '.,.,.',
                                        '.',
                                        '.',
                                        str(','.join(support_read))])


def run_dup(args):
    return resolution_DUP(*args)

def call_gt(temporary_dir, chr, candidate_single_SV, max_cluster_bias):
    reads_list = list() # [(10000, 10468, 0, 'm54238_180901_011437/52298335/ccs'), ...]
    readsfile = open("%sreads.sigs"%(temporary_dir), 'r')
    for line in readsfile:
        seq = line.strip().split('\t')
        if seq[0] != chr:
            continue
        reads_list.append([int(seq[1]), int(seq[2]), int(seq[3]), seq[4]])
    readsfile.close()
    svs_list = list()
    for item in candidate_single_SV:
        new_cluster_bias = min(max_cluster_bias, item[3] - item[2])
        svs_list.append((max(item[2] - new_cluster_bias/2, 0), item[2] + new_cluster_bias/2))
    for item in candidate_single_SV:
        new_cluster_bias = min(max_cluster_bias, item[3] - item[2])
        svs_list.append((max(item[3] - new_cluster_bias/2, 0), item[3] + new_cluster_bias/2))
    iteration_dict, primary_num_dict, cover_dict = overlap_cover(svs_list, reads_list) # both key(sv idx), value(set(read id))
    assert len(cover_dict) == 2 * len(candidate_single_SV), "overlap length error"
    candidate_single_SV_length = len(candidate_single_SV)
    for idx in range(candidate_single_SV_length):
        for item in cover_dict[idx + candidate_single_SV_length]:
            cover_dict[idx].add(item)
    for idx in range(candidate_single_SV_length, candidate_single_SV_length * 2, 1):
        cover_dict.pop(idx)
    assert len(cover_dict) == len(candidate_single_SV), "overlap length error"

    read_id_dict = dict()
    for i in range(len(candidate_single_SV)):
        read_id_dict[i] = candidate_single_SV[i][4]
    assign_list = assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict)
    # [[DV, DR, GT, GL, GQ, QUAL] ...]
    assert len(candidate_single_SV) == len(assign_list), "assign error"
    candidate_single_SV_gt = list()
    for i in range(len(candidate_single_SV)):
        candidate_single_SV_gt.append([candidate_single_SV[i][0], 
                                    candidate_single_SV[i][1], 
                                    str(candidate_single_SV[i][2]), 
                                    str(candidate_single_SV[i][3] - candidate_single_SV[i][2]), 
                                    str(len(candidate_single_SV[i][4])),
                                    str(assign_list[i][1]),
                                    str(assign_list[i][2]),
                                    str(assign_list[i][3]),
                                    str(assign_list[i][4]),
                                    str(assign_list[i][5]),
                                    ','.join(candidate_single_SV[i][4])])
    return candidate_single_SV_gt	