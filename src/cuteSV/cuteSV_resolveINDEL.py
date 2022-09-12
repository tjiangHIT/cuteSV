import numpy as np
from cuteSV.cuteSV_genotype import cal_CIPOS, overlap_cover, assign_gt
import logging

'''
*******************************************
                TO DO LIST
*******************************************
    1. Identify DP with samfile pointer;
    2. Add CIPOS, CILEN and/or CIEND;
    3. Determine (IM)PRECISE type.
*******************************************

'''

def resolution_DEL(path, chr, svtype, read_count, threshold_gloab, max_cluster_bias,
                 minimum_support_reads, bam_path, action, gt_round, remain_reads_ratio):

    '''
    cluster DEL
    ********************************************************************************************
    path:	DEL.sigs
    chr:	chromosome id
    svtype:	<DEL>
    
    SEQTYPE		read_count 	max_cluster_bias 	sv_size		threshold_gloab 	threshold_local 
    --------------------------------------------------------------------------------------------
    CCS			3			200 bp (<500 bp)	30 bp 		0.4					0.5
    CLR			5/10		200 bp (<500 bp)	50 bp 		0.3					0.7
    --------------------------------------------------------------------------------------------
    
    Input file format
    --------------------------------------------------------------------------------------------
    column	#1	#2	#3	#4	#5
            DEL	CHR	BP	LEN	ID	
    #1	deletion type
    #2	chromosome number
    #3	breakpoint in each read
    #4	DEL_len in each read
    #5	read ID
    ********************************************************************************************
    '''
    if remain_reads_ratio > 1:
        remain_reads_ratio = 1
    semi_del_cluster = list()
    semi_del_cluster.append([0,0,''])
    candidate_single_SV = list()

    file = open("%s%s.sigs"%(path, "DEL"), 'r')
    for line in file:
        seq = line.strip('\n').split('\t')
        if seq[1] != chr:
            continue

        pos = int(seq[2])
        indel_len = int(seq[3])
        read_id = seq[4]
        
        if pos - semi_del_cluster[-1][0] > max_cluster_bias:
            if len(semi_del_cluster) >= read_count:
                if semi_del_cluster[-1][0] == semi_del_cluster[-1][1] == 0:
                    pass
                else:
                    generate_del_cluster(semi_del_cluster, 
                                        chr, 
                                        svtype, 
                                        read_count, 
                                        threshold_gloab, 
                                        # threshold_local, 
                                        minimum_support_reads, 
                                        candidate_single_SV,
                                        action,
                                        gt_round,
                                        remain_reads_ratio)
            semi_del_cluster = []
            semi_del_cluster.append([pos, indel_len, read_id])
        else:
            if semi_del_cluster[-1][0] == semi_del_cluster[-1][1] == 0:
                semi_del_cluster = []
                semi_del_cluster.append([pos, indel_len, read_id])
            else:
                semi_del_cluster.append([pos, indel_len, read_id])

    if len(semi_del_cluster) >= read_count:
        if semi_del_cluster[-1][0] == semi_del_cluster[-1][1] == 0:
            pass
        else:
            generate_del_cluster(semi_del_cluster, 
                                chr, 
                                svtype, 
                                read_count, 
                                threshold_gloab, 
                                # threshold_local, 
                                minimum_support_reads, 
                                candidate_single_SV,
                                action,
                                gt_round,
                                remain_reads_ratio)
    file.close()
    if action:
        candidate_single_SV_gt = call_gt(path, chr, candidate_single_SV, max_cluster_bias, 'DEL')
        logging.info("Finished %s:%s."%(chr, "DEL"))
        return candidate_single_SV_gt
    else:
        logging.info("Finished %s:%s."%(chr, "DEL"))
        return candidate_single_SV

def generate_del_cluster(semi_del_cluster, chr, svtype, read_count, 
    threshold_gloab, minimum_support_reads, candidate_single_SV, 
    action, gt_round, remain_reads_ratio):

    '''
    generate deletion
    *************************************************************
    threshold_gloab 	threshold_local 	minimum_support_reads
    -------------------------------------------------------------
        0.3					0.7 					5		CLR
        0.4					0.5 				  <=5		CCS
    *************************************************************
    '''

    # Remove duplicates
    read_tag = dict()
    for element in semi_del_cluster:
        if element[2] not in read_tag:
            read_tag[element[2]] = element
        else:
            if element[1] > read_tag[element[2]][1]:
                read_tag[element[2]] = element

    if len(read_tag) < read_count:
        return

    read_tag2SortedList = sorted(list(read_tag.values()), key = lambda x:x[1])
    global_len = [i[1] for i in read_tag2SortedList]
    DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP = threshold_gloab * np.mean(global_len)

    last_len = read_tag2SortedList[0][1]

    allele_collect = list()
    '''
    *************************************************************
        #1 				#2			#3			#4
    -------------------------------------------------------------
        del-breakpoint	del-len		#support 	read-id
    *************************************************************
    '''
    allele_collect.append([[read_tag2SortedList[0][0]],[read_tag2SortedList[0][1]],[],
        [read_tag2SortedList[0][2]]])

    for i in read_tag2SortedList[1:]:
        if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP:
            allele_collect[-1][2].append(len(allele_collect[-1][0]))
            allele_collect.append([[],[],[],[]])

        allele_collect[-1][0].append(i[0])
        allele_collect[-1][1].append(i[1])
        allele_collect[-1][3].append(i[2])
        last_len = i[1]
    allele_collect[-1][2].append(len(allele_collect[-1][0]))
    allele_sort = sorted(allele_collect, key = lambda x:x[2])

    for allele in allele_sort:
        if allele[2][0] >= minimum_support_reads:
            allele_list = list()
            var_list = list()
            remain_allele_num = max(int(remain_reads_ratio * allele[2][0]), 1)
            pos_mean = np.mean(allele[0])
            for i in range(len(allele[0])):
                var_list.append((abs(allele[0][i] - pos_mean), i))
            var_list.sort(key=lambda x:x[0])
            for i in range(remain_allele_num):
                allele_list.append(allele[0][var_list[i][1]])
            breakpointStart = np.mean(allele_list)
            search_threshold = allele_list[0]

            allele_list = list()
            var_list = list()
            len_mean = np.mean(allele[1])
            for i in range(len(allele[1])):
                var_list.append((abs(allele[1][i] - len_mean), i))
            var_list.sort(key=lambda x:x[0])
            for i in range(remain_allele_num):
                allele_list.append(allele[1][var_list[i][1]])
            signalLen = np.mean(allele_list)

            # breakpointStart = np.mean(allele[0])
            # search_threshold = np.min(allele[0])
            CIPOS = cal_CIPOS(np.std(allele[0]), len(allele[0]))
            # signalLen = np.mean(allele[1])
            signalLen_STD = np.std(allele[1])
            CILEN = cal_CIPOS(np.std(allele[1]), len(allele[1]))

            if action:
                candidate_single_SV.append([chr, 
                                            svtype, 
                                            int(breakpointStart), 
                                            int(-signalLen), 
                                            allele[2][0], 
                                            str(CIPOS),
                                            str(CILEN),
                                            int(search_threshold),
                                            allele[3]])
            else:
                candidate_single_SV.append([chr, 
                                            svtype, 
                                            str(int(breakpointStart)), 
                                            str(int(-signalLen)), 
                                            str(allele[2][0]), 
                                            str(CIPOS),
                                            str(CILEN),
                                            '.',
                                            './.',
                                            '.,.,.',
                                            '.',
                                            '.',
                                            str(','.join(allele[3]))])
    

def resolution_INS(path, chr, svtype, read_count, threshold_gloab, 
    max_cluster_bias, minimum_support_reads, bam_path, action, gt_round, remain_reads_ratio):
    
    '''
    cluster INS
    ********************************************************************************************
    path:	INS.sigs
    chr:	chromosome id
    svtype:	<INS>
    
    SEQTYPE		read_count 	max_cluster_bias 	sv_size		threshold_gloab 	threshold_local 
    --------------------------------------------------------------------------------------------
    CCS			3			200 bp (<500 bp)	30 bp 		0.65				0.7
    CLR			5/10		100 bp (<500 bp)	50 bp 		0.2					0.6
    --------------------------------------------------------------------------------------------
    
    Input file format
    --------------------------------------------------------------------------------------------
    column	#1	#2	#3	#4	#5 #6
            INS	CHR	BP	LEN	ID SEQ
    #1	insertion type
    #2	chromosome number
    #3	breakpoint in each read
    #4	INS_len in each read
    #5	read ID
    #6  INS sequence
    ********************************************************************************************
    '''
    if remain_reads_ratio > 1:
        remain_reads_ratio = 1
    semi_ins_cluster = list()
    semi_ins_cluster.append([0,0,'',''])
    candidate_single_SV = list()

    file = open("%s%s.sigs"%(path, "INS"), 'r')
    for line in file:
        seq = line.strip('\n').split('\t')
        if seq[1] != chr:
            continue

        pos = int(seq[2])
        indel_len = int(seq[3])
        read_id = seq[4]
        try:
            ins_seq = seq[5]
        except:
            ins_seq = ''
        
        if pos - semi_ins_cluster[-1][0] > max_cluster_bias:
            if len(semi_ins_cluster) >= read_count:
                if semi_ins_cluster[-1][0] == semi_ins_cluster[-1][1] == 0:
                    pass
                else:
                    generate_ins_cluster(semi_ins_cluster, 
                                        chr, 
                                        svtype, 
                                        read_count, 
                                        threshold_gloab, 
                                        # threshold_local, 
                                        minimum_support_reads, 
                                        candidate_single_SV,
                                        action,
                                        gt_round,
                                        remain_reads_ratio)
            semi_ins_cluster = []
            semi_ins_cluster.append([pos, indel_len, read_id, ins_seq])
        else:
            if semi_ins_cluster[-1][0] == semi_ins_cluster[-1][1] == 0:
                semi_ins_cluster = []
                semi_ins_cluster.append([pos, indel_len, read_id, ins_seq])
            else:
                semi_ins_cluster.append([pos, indel_len, read_id, ins_seq])

    if len(semi_ins_cluster) >= read_count:
        if semi_ins_cluster[-1][0] == semi_ins_cluster[-1][1] == 0:
            pass
        else:
            generate_ins_cluster(semi_ins_cluster, 
                                chr, 
                                svtype, 
                                read_count, 
                                threshold_gloab, 
                                # threshold_local, 
                                minimum_support_reads, 
                                candidate_single_SV,
                                action,
                                gt_round,
                                remain_reads_ratio)
    file.close()
    if action:
        candidate_single_SV_gt = call_gt(path, chr, candidate_single_SV, 1000, 'INS') # max_cluster_bias
        logging.info("Finished %s:%s."%(chr, "INS"))
        return candidate_single_SV_gt
    else:
        logging.info("Finished %s:%s."%(chr, "INS"))
        return candidate_single_SV

def generate_ins_cluster(semi_ins_cluster, chr, svtype, read_count, 
    threshold_gloab, minimum_support_reads, candidate_single_SV, 
    action, gt_round, remain_reads_ratio):
        
    '''
    generate insertion
    *************************************************************
    threshold_gloab 	threshold_local 	minimum_support_reads
    -------------------------------------------------------------
        0.2					0.6 					5		CLR
        0.65				0.7 				  <=5		CCS
    *************************************************************
    '''
    # Remove duplicates
    read_tag = dict()
    for element in semi_ins_cluster:
        if element[2] not in read_tag:
            read_tag[element[2]] = element
        else:
            if element[1] > read_tag[element[2]][1]:
                read_tag[element[2]] = element

    if len(read_tag) < read_count:
        return

    read_tag2SortedList = sorted(list(read_tag.values()), key = lambda x:x[1])
    # start&end breakpoint
    global_len = [i[1] for i in read_tag2SortedList]
    DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * np.mean(global_len)
    last_len = read_tag2SortedList[0][1]

    allele_collect = list()
    allele_collect.append([[read_tag2SortedList[0][0]],
                            [read_tag2SortedList[0][1]],
                            [], 
                            [read_tag2SortedList[0][2]],
                            [read_tag2SortedList[0][3]]])

    for i in read_tag2SortedList[1:]:
        if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP:
            allele_collect[-1][2].append(len(allele_collect[-1][0]))
            allele_collect.append([[],[],[],[],[]])

        allele_collect[-1][0].append(i[0])
        allele_collect[-1][1].append(i[1])
        allele_collect[-1][3].append(i[2])
        allele_collect[-1][4].append(i[3])
        last_len = i[1]
    allele_collect[-1][2].append(len(allele_collect[-1][0]))
    allele_sort = sorted(allele_collect, key = lambda x:x[2])

    for allele in allele_sort:
        if allele[2][0] >= minimum_support_reads:
            allele_list = list()
            var_list = list()
            remain_allele_num = max(int(remain_reads_ratio * allele[2][0]), 1)
            pos_mean = np.mean(allele[0])
            for i in range(len(allele[0])):
                var_list.append((abs(allele[0][i] - pos_mean), i))
            var_list.sort(key=lambda x:x[0])
            for i in range(remain_allele_num):
                allele_list.append(allele[0][var_list[i][1]])
            breakpointStart = np.mean(allele_list)

            allele_list = list()
            var_list = list()
            len_mean = np.mean(allele[1])
            for i in range(len(allele[1])):
                var_list.append((abs(allele[1][i] - len_mean), i))
            var_list.sort(key=lambda x:x[0])
            for i in range(remain_allele_num):
                allele_list.append(allele[1][var_list[i][1]])
            signalLen = np.mean(allele_list)

            # breakpointStart = np.mean(allele[0])
            CIPOS = cal_CIPOS(np.std(allele[0]), len(allele[0]))
            # signalLen = np.mean(allele[1])
            signalLen_STD = np.std(allele[1])
            CILEN = cal_CIPOS(np.std(allele[1]), len(allele[1]))
            ideal_ins_seq = '<INS>'
            for pos,i in zip(allele[0],allele[4]):
                if len(i) >= int(signalLen):
                    breakpointStart = pos
                    ideal_ins_seq = i[0:int(signalLen)]
                    break
            if ideal_ins_seq == '<INS>':
                continue

            if action:
                candidate_single_SV.append([chr, 
                                            svtype, 
                                            int(breakpointStart), 
                                            int(signalLen), 
                                            allele[2][0], 
                                            str(CIPOS),
                                            str(CILEN),
                                            int(breakpointStart), 
                                            allele[3],
                                            ideal_ins_seq])
            else:
                candidate_single_SV.append([chr, 
                                            svtype, 
                                            str(int(breakpointStart)), 
                                            str(int(signalLen)), 
                                            str(allele[2][0]), 
                                            str(CIPOS),
                                            str(CILEN),
                                            '.',
                                            './.',
                                            '.,.,.',
                                            ".",
                                            ".",
                                            str(','.join(allele[3])),
                                            ideal_ins_seq])


def run_del(args):
    return resolution_DEL(*args)

def run_ins(args):
    return resolution_INS(*args)

def call_gt(temporary_dir, chr, candidate_single_SV, max_cluster_bias, svtype):
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
        svs_list.append((max(item[7] - max_cluster_bias, 0), item[7] + max_cluster_bias))
    iteration_dict, primary_num_dict, cover_dict = overlap_cover(svs_list, reads_list) # both key(sv idx), value(set(read id))
    assert len(cover_dict) == len(candidate_single_SV), "overlap length error"

    read_id_dict = dict()
    for i in range(len(candidate_single_SV)):
        read_id_dict[i] = candidate_single_SV[i][8]
    assign_list = assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict)
    # [[DV, DR, GT, GL, GQ, QUAL] ...]
    assert len(candidate_single_SV) == len(assign_list), "assign error"
    candidate_single_SV_gt = list()
    for i in range(len(candidate_single_SV)):
        assert len(assign_list[i]) == 6, "assign genotype error"
        candidate_single_SV_gt.append([candidate_single_SV[i][0], 
                                    candidate_single_SV[i][1], 
                                    str(candidate_single_SV[i][2]), 
                                    str(candidate_single_SV[i][3]), 
                                    str(candidate_single_SV[i][4]), 
                                    candidate_single_SV[i][5],
                                    candidate_single_SV[i][6],
                                    str(assign_list[i][1]),
                                    str(assign_list[i][2]),
                                    str(assign_list[i][3]),
                                    str(assign_list[i][4]),
                                    str(assign_list[i][5]),
                                    ','.join(candidate_single_SV[i][8])])
        if svtype == 'INS':
            candidate_single_SV_gt[i].append(candidate_single_SV[i][9])
    return candidate_single_SV_gt