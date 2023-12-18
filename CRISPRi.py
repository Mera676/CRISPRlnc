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
CRISPRi_SVM = joblib.load('SVM_model/CRISPRi.pkl')# open the file in the write mode

if __name__ == '__main__':
    needoff = needoff
    fa_dict = file_handle.read_fa(NCRNA_FA)
    oldsg_dict = file_handle.read_sg(SG_FA)
    promoterfa_dict = file_handle.read_promoterfa(PROMOTER_FILE, distance)
    chr_read = open(CHR_BED_DICT, 'rb')
    location_read = open(LOCATION_DICT, 'rb')
    chr_dict = pickle.load(chr_read)
    location_dict = pickle.load(location_read)
    # 脱靶输入文件
    sgRNAfile = open('sgRNA_input/sgRNA_input.txt', 'w')
    sgRNAfile.write(sys.path[0] + '/' + genome_fa + '\n')
    # sgRNAfile.write(sys.path[0]+'/fasta/human_hg38.fa'+'\n')
    # second line
    sgRNAfile.write('NNNNNNNNNNNNNNNNNNNNNRG' + '\n')

    # 输出文件
    f = open('CRISPRlnc_CRISPRi_result.tsv', 'w')  # create the csv writer
    head = 'transcript\tsequence\tPAM\tGC content\tposition\tstrand\tCutting effectiveness\tregion\tlocation\toff score\tComprehensive score\n'
    f.write(head)
    result = {}
    input = {}
    res_list_format = []
    # sg_list
    for gene in oldsg_dict.keys():
        sequence_num = oldsg_dict[gene]['sequence_num']
        sg_wholeseq = oldsg_dict[gene]['sequence'].replace('\n', '')
        input[sequence_num] = []
        # 先来匹配整段基因
        whole_dict = fa_dict
        transcript_list = []
        for transcript in whole_dict:
            matchseq = whole_dict[transcript]['sequence'].replace('\n', '').find(sg_wholeseq)
            if matchseq != -1:
                transcript_list.append(transcript)
        # print(transcript_list)

        # 如果输入的序列匹配不到基因
        if transcript_list == []:
            sg_dict = file_handle.Designsg(sg_wholeseq)
            # 将得到的sgRNA写入文件，计算脱靶
            for gene in sg_dict:
                for key1 in sg_dict[gene]:
                    sgRNA = str(key1[0][0:20]) + 'NNN'
                    sgRNAfile.write(sgRNA + ' ' + '2' + '\n')
            sgRNAfile.close()
            # pachong
            if needoff == 1:
                # deal with off-target data
                if device == 0:
                    cmd = " ./cas-offinder sgRNA_input/sgRNA_input.txt C off-target/output.txt"
                elif device == 1:
                    cmd = " ./cas-offinder sgRNA_input/sgRNA_input.txt G off-target/output.txt"
                os.system(cmd)
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

            for key1 in sg_dict:
                CRISPRi_Predcit = str(CRISPRi_SVM.predict(file_handle.CRISPRifeature(key1[0:20]))[0])
                if ('TTTT' in key1[0:20]) or (
                        (key1[0:20].count('G') + key1[0:20].count('C')) / len(key1[0:20]) < 0.20) or (
                        (key1[0:20].count('G') + key1[0:20].count('C')) / len(key1[0:20]) > 0.80):
                    CRISPRi_Predcit = str(0)
                if CRISPRi_Predcit == '0':
                    CRISPRi_Predcit = 'inefficient'
                    predict_score = 0
                if CRISPRi_Predcit == '1':
                    CRISPRi_Predcit = 'efficient'
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
                        sgRNA = {
                            'transcript_id': sequence_num,
                            'sg_sequence': sg_sequence,
                            'pam': pam,
                            'GCcontent': GCcontent,
                            'chr_position': chr_position,
                            'strand': strand,
                            'CRISPRi_Predcit': CRISPRi_Predcit,
                            'region': region,
                            'location': location,
                            'offscore': offscore,
                            'sgRNA_off': sgRNA_off,
                            'all_score': all_score
                        }
                        sg_detail = f'{sequence_num}\t{sg_sequence}\t{pam}\t{GCcontent}\t{chr_position}\t{strand}\t{CRISPRi_Predcit}\t{region}\t{location}\t{offscore}\t{all_score}\n'
                        f.write(sg_detail)
                        res_list_format.append(sgRNA)

                    if needoff == 1:  # 需要脱靶的情况
                        offscore = offscore_list[sg[0]]
                        sgRNA_off = offsituation[sg[0]]
                        all_score = (predict_score - offscore) * 0.875
                        sgRNA = {
                            'transcript_id': sequence_num,
                            'sg_sequence': sg_sequence,
                            'pam': pam,
                            'GCcontent': GCcontent,
                            'chr_position': chr_position,
                            'strand': strand,
                            'CRISPRi_Predcit': CRISPRi_Predcit,
                            'region': region,
                            'location': location,
                            'offscore': offscore,
                            'sgRNA_off': sgRNA_off,
                            'all_score': all_score
                        }
                        sg_detail = f'{sequence_num}\t{sg_sequence}\t{pam}\t{GCcontent}\t{chr_position}\t{strand}\t{CRISPRi_Predcit}\t{region}\t{location}\t{offscore}\t{all_score}\n'
                        f.write(sg_detail)
                        res_list_format.append(sgRNA)

        # 如果可以匹配到gene
        else:
            for transcript in transcript_list:
                # 定义返回的短第二个字典，匹配上的输入序列的所有转录本信息
                transcript_info = (transcript, promoterfa_dict[transcript]['gene'], promoterfa_dict[transcript]['type'],
                                   promoterfa_dict[transcript]['symbol'], promoterfa_dict[transcript]['chr_name'],
                                   promoterfa_dict[transcript]['start'], promoterfa_dict[transcript]['end'],
                                   promoterfa_dict[transcript]['strand'])
                input[sequence_num].append(transcript_info)

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
                result[transcript]['sgRNAs'] = []
                chr_sequence = chr_dict[chr_name].upper()

                # pachong
                if needoff == 1:
                    # deal with off-target data
                    if device == 0:
                        cmd = " ./cas-offinder sgRNA_input/sgRNA_input.txt C off-target/output.txt"
                    elif device == 1:
                        cmd = " ./cas-offinder sgRNA_input/sgRNA_input.txt G off-target/output.txt"
                    os.system(cmd)
                    seqs = []
                    for key1 in sg_dict:
                        seqs.append((key1[0:20], key1[20:23]))
                    # 计算出脱靶情况和脱靶得分
                    try:
                        fh = open("off-target/offtarget_output.txt", "r")
                        offsituation = file_handle.offline_off('off-target/offtarget_output.txt')
                    except IOError:
                        print("Error: No off-target result file found, please check if Cas-offinder is working correctly!")
                    else:
                        print("The off-target result was calculated successfully.")
                        fh.close()
    #                 offsituation = file_handle.offline_off('offtarget_result/output.txt')
                    offscore_list = file_handle.geoff_score(offsituation)

                for key1 in sg_dict:
                    CRISPRi_Predcit = str(CRISPRi_SVM.predict(file_handle.CRISPRifeature(key1[0:20]))[0])
                    if ('TTTT' in key1[0:20]) or (
                            (key1[0:20].count('G') + key1[0:20].count('C')) / len(key1[0:20]) < 0.20) or (
                            (key1[0:20].count('G') + key1[0:20].count('C')) / len(key1[0:20]) > 0.80):
                        CRISPRi_Predcit = str(0)
                    if CRISPRi_Predcit == '0':
                        CRISPRi_Predcit = 'inefficient'
                        predict_score = 0
                    if CRISPRi_Predcit == '1':
                        CRISPRi_Predcit = 'efficient'
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
                            region_score = 0.5
                        elif (promoterfa_dict[transcript]['strand'] == '+') and (position >= int(gene_start)):
                            region = 'gene body'
                            region_score = 0.1
                        elif (promoterfa_dict[transcript]['strand'] == '-') and (position > int(gene_end)):
                            region = 'Promoter'
                            region_score = 0.5
                        elif (promoterfa_dict[transcript]['strand'] == '-') and (position <= int(gene_end)):
                            region = 'gene body'
                            region_score = 0.1

                        if int(gene_start) - distance <= position <= int(
                                gene_end) + distance and needoff == 0:  # 不需要脱靶的情况
                            offscore = 'NA'
                            sgRNA_off = []
                            all_score = (predict_score + region_score) * 0.875
                            sg_detail = f'{transcript}\t{sg_sequence}\t{pam}\t{GCcontent}\t{chr_position}\t{strand}\t{CRISPRi_Predcit}\t{region}\t{location}\t{offscore}\t{all_score}\n'
                            f.write(sg_detail)


                        if int(gene_start) - distance <= position <= int(
                                gene_end) + distance and needoff == 1:  # 需要脱靶的情况
                            offscore = offscore_list[sg[0]]
                            sgRNA_off = offsituation[sg[0]]
                            all_score = (predict_score + region_score - offscore) * 0.875
                            sg_detail = f'{transcript}\t{sg_sequence}\t{pam}\t{GCcontent}\t{chr_position}\t{strand}\t{CRISPRi_Predcit}\t{region}\t{location}\t{offscore}\t{all_score}\n'
                            f.write(sg_detail)
