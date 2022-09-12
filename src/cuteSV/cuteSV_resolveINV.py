import numpy as np
import logging
from cuteSV.cuteSV_genotype import overlap_cover, assign_gt

def resolution_INV(path, chr, svtype, read_count, max_cluster_bias, sv_size, 
    bam_path, action, MaxSize, gt_round):
    '''
    cluster INV
    ************************************************************************
    path:	INV.sigs
    chr:	chromosome id
    svtype:	<INV>
    
    SEQTYPE		read_count 	max_cluster_bias 	sv_size 
    ------------------------------------------------------------------------
    CCS			5			10 bp (<500 bp)		20 bp 	
    CLR			5			20 bp (<500 bp)		50 bp 	
    ------------------------------------------------------------------------
    
    Input file format
    ------------------------------------------------------------------------
    column	#1	#2	#3	#4	#5
            INV	CHR	BP1	BP2	ID	
    #1	inversion type
    #2	chromosome number
    #3	breakpoint_1 in each read
    #4	breakpoint_2 in each read
    #5	read ID
    ************************************************************************
    '''

    # Initialization of some temporary variables
    semi_inv_cluster = list()
    semi_inv_cluster.append([0,0,'',''])
    candidate_single_SV = list()

    # Load inputs & cluster breakpoint from each signature read 
    file = open("%s%s.sigs"%(path, "INV"), 'r')
    for line in file:
        seq = line.strip('\n').split('\t')
        if seq[1] != chr:
            continue

        strand = seq[2]
        breakpoint_1_in_read = int(seq[3])
        breakpoint_2_in_read = int(seq[4])
        read_id = seq[5]

        # print("new")
        # print(seq[1], seq[2], seq[3], seq[4], seq[5])
        # print(semi_inv_cluster)

        if breakpoint_1_in_read - semi_inv_cluster[-1][0] > max_cluster_bias or breakpoint_2_in_read - semi_inv_cluster[-1][1] > max_cluster_bias or strand != semi_inv_cluster[-1][-1]:
            if len(semi_inv_cluster) >= read_count:
                if semi_inv_cluster[-1][0] == semi_inv_cluster[-1][1] == 0:
                    pass
                else:
                    generate_semi_inv_cluster(semi_inv_cluster, 
                                            chr, 
                                            svtype, 
                                            read_count, 
                                            sv_size, 
                                            candidate_single_SV, 
                                            max_cluster_bias,
                                            action,
                                            MaxSize,
                                            gt_round)
            semi_inv_cluster = []
            semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id, strand])
        else:
            if semi_inv_cluster[-1][0] == semi_inv_cluster[-1][1] == 0:
                semi_inv_cluster = []
                semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id, strand])
            else:
                semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id, strand])

    if len(semi_inv_cluster) >= read_count:
        if semi_inv_cluster[-1][0] == semi_inv_cluster[-1][1] == 0:
            pass
        else:
            generate_semi_inv_cluster(semi_inv_cluster, 
                                    chr, svtype, 
                                    read_count, 
                                    sv_size, 
                                    candidate_single_SV, 
                                    max_cluster_bias,
                                    action,
                                    MaxSize,
                                    gt_round)
    file.close()
    if action:
        candidate_single_SV_gt = call_gt(path, chr, candidate_single_SV, max_cluster_bias)
        logging.info("Finished %s:%s."%(chr, "INV"))
        return candidate_single_SV_gt
    else:
        logging.info("Finished %s:%s."%(chr, "INV"))
        return candidate_single_SV

def generate_semi_inv_cluster(semi_inv_cluster, chr, svtype, read_count, sv_size, 
    candidate_single_SV, max_cluster_bias, action, MaxSize, gt_round):

    strand = semi_inv_cluster[0][-1]

    read_id = [i[2] for i in semi_inv_cluster]
    support_read = len(list(set(read_id)))
    if support_read < read_count:
        return

    inv_cluster_b2 = sorted(semi_inv_cluster, key = lambda x:x[1])

    # breakpoint_1 = np.mean(breakpoint_1_candidate)
    last_bp = inv_cluster_b2[0][1]
    temp_count = 1
    # max_count = 0
    temp_sum_b1 = inv_cluster_b2[0][0]
    temp_sum_b2 = last_bp

    # max_sum = 0
    temp_id = dict()
    temp_id[inv_cluster_b2[0][2]] = 0

    for i in inv_cluster_b2[1:]:
        if i[1] - last_bp > max_cluster_bias:
            if temp_count >= read_count:
                max_count_id = len(temp_id)

                breakpoint_1 = round(temp_sum_b1 / temp_count)
                breakpoint_2 = round(temp_sum_b2 / temp_count)
                inv_len = breakpoint_2 - breakpoint_1
                if inv_len >= sv_size and max_count_id >= read_count:
                    # candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\n'%(chr, svtype, breakpoint_1, breakpoint_2, max_count_id))
                    if inv_len <= MaxSize or MaxSize == -1:
                        if action:
                            candidate_single_SV.append([chr, 
                                                        svtype, 
                                                        breakpoint_1, 
                                                        inv_len, 
                                                        max_count_id,
                                                        strand,
                                                        list(temp_id.keys()),
                                                        breakpoint_2])
                        else:
                            candidate_single_SV.append([chr, 
                                                        svtype, 
                                                        str(int(breakpoint_1)), 
                                                        str(int(inv_len)), 
                                                        str(max_count_id),
                                                        '.',
                                                        './.',
                                                        strand,
                                                        '.,.,.',
                                                        ".",
                                                        ".",
                                                        str(','.join(list(temp_id.keys())))])
                        # print(chr, svtype, str(int(breakpoint_1)), str(int(inv_len)), str(max_count_id), str(DR), str(GT), strand)

            temp_id = dict()
            temp_count = 1
            temp_sum_b1 = i[0]
            temp_sum_b2 = i[1]
            temp_id[i[2]] = 0
        else:
            if i[2] not in temp_id:
                temp_id[i[2]] = 0
            else:
                temp_id[i[2]] += 1
            temp_count += 1
            temp_sum_b1 += i[0]
            temp_sum_b2 += i[1]
        last_bp = i[1]
    if temp_count >= read_count:
        max_count_id = len(temp_id)
        breakpoint_1 = round(temp_sum_b1 / temp_count)
        breakpoint_2 = round(temp_sum_b2 / temp_count)
        inv_len = breakpoint_2 - breakpoint_1
        if inv_len >= sv_size and max_count_id >= read_count:
            # candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\n'%(chr, svtype, breakpoint_1, breakpoint_2, max_count_id))
            if inv_len <= MaxSize or MaxSize == -1:
                if action:
                    candidate_single_SV.append([chr, 
                                                svtype, 
                                                breakpoint_1, 
                                                inv_len, 
                                                max_count_id,
                                                strand,
                                                list(temp_id.keys()),
                                                breakpoint_2])
                else:
                    candidate_single_SV.append([chr, 
                                                svtype, 
                                                str(int(breakpoint_1)), 
                                                str(int(inv_len)), 
                                                str(max_count_id),
                                                '.',
                                                './.',
                                                strand,
                                                '.,.,.',
                                                ".",
                                                ".",
                                                str(','.join(list(temp_id.keys())))])
                # print(chr, svtype, str(int(breakpoint_1)), str(int(inv_len)), str(max_count_id), str(DR), str(GT), strand)

def run_inv(args):
    return resolution_INV(*args)

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
        svs_list.append((max(item[2] - max_cluster_bias/2, 0), item[2] + max_cluster_bias/2))
    for item in candidate_single_SV:
        svs_list.append((max(item[7] - max_cluster_bias/2, 0), item[7] + max_cluster_bias/2))
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
        read_id_dict[i] = candidate_single_SV[i][6]
    assign_list = assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict)
    # [[DV, DR, GT, GL, GQ, QUAL] ...]
    assert len(candidate_single_SV) == len(assign_list), "assign error"
    candidate_single_SV_gt = list()
    for i in range(len(candidate_single_SV)):
        candidate_single_SV_gt.append([candidate_single_SV[i][0], 
                                    candidate_single_SV[i][1], 
                                    str(int(candidate_single_SV[i][2])), 
                                    str(int(candidate_single_SV[i][3])), 
                                    str(candidate_single_SV[i][4]), 
                                    str(assign_list[i][1]),
                                    str(assign_list[i][2]),
                                    candidate_single_SV[i][5],
                                    str(assign_list[i][3]),
                                    str(assign_list[i][4]),
                                    str(assign_list[i][5]),
                                    ','.join(candidate_single_SV[i][6])])
    return candidate_single_SV_gt