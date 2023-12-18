import joblib
import os, sys
import pickle
import file_handle

SG_FA = sys.argv[1] # sgRNA文件，输入的系数1。
NCRNA_FA = 'Annotation/ncrna.fa'  # 非编码基因文件，输入系数2。
PROMOTER_FILE = 'Annotation/ncRNA_withPromoter.fa'  # 非编码带有启动子序列的序列文件，输入系数3。
LOCATION_DICT = 'dict/location_dict.pkl'  # 位置区域文件，输入系数4.
CHR_BED_DICT = 'dict/chr_dict.pkl'  # 染色体文件，输入系数5
genome_fa = sys.argv[2] # 用于脱靶的基因组序列文件
distance = int(sys.argv[3])
needoff = int(sys.argv[4])
device = int(sys.argv[5])
CRISPRko_SVM = joblib.load('SVM_model/CRISPRko.pkl')# open the file in the write mode

if __name__ == '__main__':
    needoff = needoff
    fa_dict = file_handle.read_fa(NCRNA_FA)
    oldsg_dict = file_handle.read_sg(SG_FA)
    promoterfa_dict = file_handle.read_promoterfa(PROMOTER_FILE, distance)
    chr_read = open(CHR_BED_DICT, 'rb')
    location_read = open(LOCATION_DICT, 'rb')
    chr_dict = pickle.load(chr_read)
    location_dict = pickle.load(location_read)

    #脱靶输入文件
    sgRNAfile = open('sgRNA_input/sgRNA_input.txt', 'w')
    sgRNAfile.write(genome_fa + '\n')
    # sgRNAfile.write(sys.path[0]+'/fasta/human_hg38.fa'+'\n')
    # second line
    sgRNAfile.write('NNNNNNNNNNNNNNNNNNNNNRG' + '\n')

    # 输出文件
    f = open('CRISPRlnc_CRISPRko_result.tsv', 'w')  # create the csv writer
    head = 'transcript\tpair id\tpair score\tsgRNA1 sequence\tsgRNA1 PAM\tsgRNA1 GC content\tsgRNA1 position\tsgRNA 1 strand\tsgRNA1 Cutting effectiveness\tsgRNA1 region\tsgRNA1 location\tsgRNA1 off score\tsgRNA2 sequence\tsgRNA2 PAM\tsgRNA2 GC content\tsgRNA2 position\tsgRNA2 strand\tsgRNA2 Cutting effectiveness\tsgRNA2 region\tsgRNA2 location\tsgRNA2 off score\n'
    f.write(head)
    result = {}
    input = {}
    pair_dict = {}
    # sg_list
    for gene in oldsg_dict.keys():
        sequence_num = oldsg_dict[gene]['sequence_num']
        sg_wholeseq = oldsg_dict[gene]['sequence'].replace('\n', '')
        input[sequence_num] = []
#         result[sequence_num] = {}
        # 先来匹配整段基因
        whole_dict = fa_dict
        transcript_list = []
        for transcript in whole_dict:
            matchseq = whole_dict[transcript]['sequence'].replace('\n', '').find(sg_wholeseq)
            if matchseq != -1:
                transcript_list.append(transcript)

        # 如果输入的序列匹配不到基因
        if transcript_list == []:
            sg_dict = file_handle.Designsg(sg_wholeseq)
            #将得到的sgRNA写入文件，计算脱靶
            for gene in sg_dict:
                for key1 in sg_dict[gene]:
                    sgRNA = str(key1[0][0:20]) + 'NNN'
                    sgRNAfile.write(sgRNA + ' ' + '2' + '\n')
            sgRNAfile.close()


            # pachong
            if needoff == 1:
                # deal with off-target data
                if device == 0:
                    cmd = " ./cas-offinder sgRNA_input/sgRNA_input.txt C off-target/off-target_output.txt"
                elif device == 1:
                    cmd = " ./cas-offinder sgRNA_input/sgRNA_input.txt G off-target/off-target_output.txt"
                os.system(cmd)
                print(cmd)
                seqs = []
                for key1 in sg_dict:
                    # sg_list.append(key1)
                    seqs.append((key1[0:20], key1[20:23]))
                # 计算出脱靶情况和脱靶得分
                try:
                    fh = open("off-target/offtarget_output.txt", "r")
                    offsituation = file_handle.offline_off('off-target/offtarget_output.txt')
                except IOError:
                    print("Error: No off-target result files were found, please check if Cas-offinder is working properly.")
                else:
                    print("The off-target result was calculated successfully.")
                    fh.close()
#                 offsituation = file_handle.offline_off('offtarget_result/output.txt')
                offscore_list = file_handle.geoff_score(offsituation)

            # 计算pair
            pair_dict[sequence_num] = {}
            pair_count = 1
            new_sg_dict = list(sg_dict.values())
            for sg in new_sg_dict:
                sg = sg[0]
                first_start = int(sg[1][0])
                for sg2 in new_sg_dict:
                    sg2 = sg2[0]
                    start = int(sg2[1][0])
                    if start > first_start + 50:
                        pair_dict[sequence_num][pair_count] = (first_start, start)
                        pair_count += 1

            for key1 in sg_dict:
                CRISPRko_Predcit = str(CRISPRko_SVM.predict(file_handle.CRISPRkofeature(key1[0:20]))[0])
#                 CRISPRko_Predcit = str(1)
                if ('TTTT' in key1[0:20]) or (
                        (key1[0:20].count('G') + key1[0:20].count('C')) / len(key1[0:20]) < 0.20) or (
                        (key1[0:20].count('G') + key1[0:20].count('C')) / len(key1[0:20]) > 0.80):
                    CRISPRko_Predcit = str(0)
                if CRISPRko_Predcit == '0':
                    CRISPRko_Predcit = 'inefficient'
                    predict_score = 0
                if CRISPRko_Predcit == '1':
                    CRISPRko_Predcit = 'efficient'
                    predict_score = 1
                GCcontent = (key1[0:20].count('G') + key1[0:20].count('C')) / len(key1[0:20])
                for sg in sg_dict[key1]:
                    start = sg[1][0]
                    end = sg[1][1]
                    sg_strand = sg[1][2]
                    strand = sg_strand
                    sg_sequence = sg[0][0:20]
                    pam = sg[0][20:23]
                    chr_position = 'sequence:' + str(start) + '-' + str(end)
                    region = 'NA'
                    location = 'NA'
                    transcript = 'NA'
                    if needoff == 0:  # 不需要脱靶的情况
                        offscore = 'NA'
                        sgRNA_off = []
                        all_score = predict_score * 0.875
                        sgRNA = (
                            sg_sequence, pam, GCcontent, chr_position, strand, CRISPRko_Predcit, region, location,
                            offscore,
                            sgRNA_off, all_score)
                        result[sequence_num][start] = sgRNA

                    if needoff == 1:  # 需要脱靶的情况
                        offscore = offscore_list[sg[0]]
                        sgRNA_off = offsituation[sg[0]]
                        all_score = (predict_score - offscore) * 0.875
                        sgRNA = (
                            sg_sequence, pam, GCcontent, chr_position, strand, CRISPRko_Predcit, region, location,
                            offscore,
                            sgRNA_off, all_score)
                        result[sequence_num][start] = sgRNA

        else:
            for transcript in transcript_list:
                # 定义返回的第二个字典，匹配上的输入序列的所有转录本信息
                transcript_info = (transcript, promoterfa_dict[transcript]['gene'], promoterfa_dict[transcript]['type'],
                                   promoterfa_dict[transcript]['symbol'], promoterfa_dict[transcript]['chr_name'],
                                   promoterfa_dict[transcript]['start'], promoterfa_dict[transcript]['end'],
                                   promoterfa_dict[transcript]['strand'])
                input[sequence_num].append(transcript_info)

                sg_list = []
                result[transcript] = {}
                chr_name = promoterfa_dict[transcript]['chr_name']
                sg_dict = file_handle.Designsg(promoterfa_dict[transcript]['sequence'])  # 从启动子区域开始设计sgRNA
                # 将得到的sgRNA写入文件，计算脱靶
                for gene in sg_dict:
                    for key1 in sg_dict[gene]:
                        sgRNA = str(key1[0][0:20]) + 'NNN'
                        sgRNAfile.write(sgRNA + ' ' + '2' + '\n')
                sgRNAfile.close()


                gene_start = promoterfa_dict[transcript]['start']
                gene_end = promoterfa_dict[transcript]['end']
                type = promoterfa_dict[transcript]['type']
                gene_name = promoterfa_dict[transcript]['gene']
                symbol = promoterfa_dict[transcript]['symbol']
                transcript_infomation = (gene_name, chr_name, gene_start, gene_end, symbol, type)
                result[transcript]['info'] = transcript_infomation
                chr_sequence = chr_dict[chr_name].upper()

                # pachong
                if needoff == 1:
                    # deal with off-target data
                    if device == 0:
                        cmd = " ./cas-offinder sgRNA_input/sgRNA_input.txt C off-target/offtarget_output.txt"
                    elif device == 1:
                        cmd = " ./cas-offinder sgRNA_input/sgRNA_input.txt G off-target/offtarget_output.txt"
                    os.system(cmd)
                    seqs = []
                    for key1 in sg_dict:
                        sg_list.append(key1)
                        seqs.append((key1[0:20], key1[20:23]))
                    # 计算出脱靶情况和脱靶得分
                    try:
                        fh = open("off-target/offtarget_output.txt", "r")
                        offsituation = file_handle.offline_off('off-target/offtarget_output.txt')
                    except IOError:
                        print("Error: No off-target result files were found, please check if Cas-offinder is working properly.")
                    else:
                        print("The off-target result was calculated successfully.")
                        fh.close()
    #                 offsituation = file_handle.offline_off('offtarget_result/output.txt')
                    offscore_list = file_handle.geoff_score(offsituation)

                # 计算pair
                pair_dict[transcript] = {}
                pair_count = 1
                new_sg_dict = list(sg_dict.values())
                for sg in new_sg_dict:
                    sg = sg[0]
                    first_start = int(sg[1][0])
                    for sg2 in new_sg_dict:
                        sg2 = sg2[0]
                        start = int(sg2[1][0])
                        if start > first_start + 200:
                            pair_dict[transcript][pair_count] = (first_start, start)
                            pair_count += 1
            

                for key1 in sg_dict:
                    CRISPRko_Predcit = str(CRISPRko_SVM.predict(file_handle.CRISPRkofeature(key1[0:20]))[0])
                    if ('TTTT' in key1[0:20]) or (
                            (key1[0:20].count('G') + key1[0:20].count('C')) / len(key1[0:20]) < 0.20) or (
                            (key1[0:20].count('G') + key1[0:20].count('C')) / len(key1[0:20]) > 0.80):
                        CRISPRko_Predcit = str(0)
                    if CRISPRko_Predcit == '0':
                        CRISPRko_Predcit = 'inefficient'
                        predict_score = 0
                    if CRISPRko_Predcit == '1':
                        CRISPRko_Predcit = 'efficient'
                        predict_score = 0.5
                    GCcontent = (key1[0:20].count('G') + key1[0:20].count('C')) / len(key1[0:20])
                    for sg in sg_dict[key1]:
                        start = sg[1][0]
                        end = sg[1][1]
                        sg_strand = sg[1][2]
                        strand = sg_strand
                        if (promoterfa_dict[transcript]['strand'] == '-') and (sg_strand == '+'):
                            strand = '-'  # 最终相对于染色体的正负链
                        if (promoterfa_dict[transcript]['strand'] == '-') and (sg_strand == '-'):
                            strand = '+'
                        sg_sequence = sg[0][0:20]
                        pam = sg[0][20:23]
                        if strand == '+':
                            position = chr_sequence.find(sg[0])
                            chr_position = chr_name + ':' + str(position + 1) + '-' + str(position + 23)
                            location = file_handle.getlocation(chr_name, strand, position + 1, position + 23,
                                                               location_dict)
                        if strand == '-':
                            position = chr_sequence.find(file_handle.reverse_complement(sg[0]))
                            chr_position = chr_name + ':' + str(position + 1) + '-' + str(position + 23)
                            location = file_handle.getlocation(chr_name, strand, position + 1, position + 23,
                                                               location_dict)
                        if (promoterfa_dict[transcript]['strand'] == '+') and (position < int(gene_start)):
                            region = 'Promoter'
                            region_score = 0.1
                        elif (promoterfa_dict[transcript]['strand'] == '+') and (position >= int(gene_start)):
                            region = 'gene body'
                            region_score = 0.5
                        elif (promoterfa_dict[transcript]['strand'] == '-') and (position > int(gene_end)):
                            region = 'Promoter'
                            region_score = 0.1
                        elif (promoterfa_dict[transcript]['strand'] == '-') and (position <= int(gene_end)):
                            region = 'gene body'
                            region_score = 0.5

                        if int(gene_start) - distance <= position <= int(
                                gene_end) + distance and needoff == 0:  # 不需要脱靶的情况
                            offscore = 'NA'
                            sgRNA_off = []
                            all_score = (region_score + predict_score) * 0.875
                            sgRNA = (
                            sg_sequence, pam, GCcontent, chr_position, strand, CRISPRko_Predcit, region, location,
                            offscore,
                            sgRNA_off, all_score)
                            result[transcript][start] = sgRNA
                            

                        if int(gene_start) - distance <= position <= int(
                                gene_end) + distance and needoff == 1:  # 需要脱靶的情况
                            offscore = offscore_list[sg[0][0:20]]
                            sgRNA_off = offsituation[sg[0][0:20]]
                            all_score = region_score + predict_score - offscore
                            sgRNA = (
                                sg_sequence, pam, GCcontent, chr_position, strand, CRISPRko_Predcit, region, location,
                                offscore,
                                sgRNA_off, all_score)
                            result[transcript][start] = sgRNA
    
    res_list_format = []
    for transcript in result:
        for pair_id in pair_dict[transcript]:
            first_sgRNA = pair_dict[transcript][pair_id][0]
            second_sgRNA = pair_dict[transcript][pair_id][1]
            if first_sgRNA in result[transcript] and second_sgRNA in result[transcript]:
                pair_distance = str(second_sgRNA - first_sgRNA) + 'bp'
                # print(pair_distance)
                first_sgRNA_info = result[transcript][first_sgRNA]
                second_sgRNA_info = result[transcript][second_sgRNA]
                sgRNA1_sequence = first_sgRNA_info[0]
                sgRNA1_pam = first_sgRNA_info[1]
                sgRNA1_GC = first_sgRNA_info[2]
                sgRNA1_position = first_sgRNA_info[3]
                sgRNA1_strand = first_sgRNA_info[4]
                sgRNA1_efficience = first_sgRNA_info[5]
                sgRNA1_region = first_sgRNA_info[6]
                sgRNA1_location = first_sgRNA_info[7]
                sgRNA1_offscore = first_sgRNA_info[8]
                sgRNA1_offsituation = first_sgRNA_info[9]
                sgRNA1_allscore = first_sgRNA_info[10]
                sgRNA2_sequence = second_sgRNA_info[0]
                sgRNA2_pam = second_sgRNA_info[1]
                sgRNA2_GC = second_sgRNA_info[2]
                sgRNA2_position = second_sgRNA_info[3]
                sgRNA2_strand = second_sgRNA_info[4]
                sgRNA2_efficience = second_sgRNA_info[5]
                sgRNA2_region = second_sgRNA_info[6]
                sgRNA2_location = second_sgRNA_info[7]
                sgRNA2_offscore = second_sgRNA_info[8]
                sgRNA2_offsituation = second_sgRNA_info[9]
                sgRNA2_allscore = second_sgRNA_info[10]
                pair_score = (sgRNA1_allscore + sgRNA2_allscore) * 0.5
                pair = {
                    'transcript_id': transcript,
                    'pair_id': pair_id,
                    'pair_score': pair_score,
                    'sgRNA1_sequence': sgRNA1_sequence,
                    'sgRNA1_pam': sgRNA1_pam,
                    'sgRNA1_GC': sgRNA1_GC,
                    'sgRNA1_position': sgRNA1_position,
                    'sgRNA1_strand': sgRNA1_strand,
                    'sgRNA1_efficience': sgRNA1_efficience,
                    'sgRNA1_region': sgRNA1_region,
                    'sgRNA1_location': sgRNA1_location,
                    'sgRNA1_offscore': sgRNA1_offscore,
                    'sgRNA1_offsituation': sgRNA1_offsituation,
                    'sgRNA2_sequence': sgRNA2_sequence,
                    'sgRNA2_pam': sgRNA2_pam,
                    'sgRNA2_GC': sgRNA2_GC,
                    'sgRNA2_position': sgRNA2_position,
                    'sgRNA2_strand': sgRNA2_strand,
                    'sgRNA2_efficience': sgRNA2_efficience,
                    'sgRNA2_region': sgRNA2_region,
                    'sgRNA2_location': sgRNA2_location,
                    'sgRNA2_offscore': sgRNA2_offscore,
                    'sgRNA2_offsituation': sgRNA2_offsituation,
                }
                pair_detail = f'{transcript}\t{pair_id}\t{pair_score}\t{sgRNA1_sequence}\t{sgRNA1_pam}\t{sgRNA1_GC}\t{sgRNA1_position}\t{sgRNA1_strand}\t{sgRNA1_efficience}\t{sgRNA1_region}\t{sgRNA1_location}\t{sgRNA1_offscore}\t{sgRNA2_sequence}\t{sgRNA2_pam}\t{sgRNA2_GC}\t{sgRNA2_position}\t{sgRNA2_strand}\t{sgRNA2_efficience}\t{sgRNA2_region}\t{sgRNA2_location}\t{sgRNA2_offscore}\n'
                f.write(pair_detail)
                res_list_format.append(pair)

    res_list_format = sorted(res_list_format, key=lambda k: k['pair_score'], reverse=True)