# coding:utf-8
import re
import joblib
import os, sys
import pickle
import file_handle
import csv
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.svm import SVC
from numpy import loadtxt
from numpy import sort
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import SelectFromModel
import warnings
warnings.filterwarnings("ignore")
import RNA

#Reading file
#fa_dict,anti_fa_dict=file_handle.read_fa('Annotation/Homo_sapiens.GRCh38.ncrna.fa')    
#bed_dict=file_handle.read_Promoter('Annotation/GRCh38transcript.promoter.bed')
#sg_dict=file_handle.read_sg('sg.fa')


ncrna_fa = sys.argv[1]
bed_fa = sys.argv[2]
sg_fa = sys.argv[3]
genome_fa=sys.argv[4]

fa_dict,anti_fa_dict=file_handle.read_fa(ncrna_fa)    
bed_dict=file_handle.read_Promoter(bed_fa)
sg_dict=file_handle.read_sg(sg_fa)






CRISPRi_SVM= joblib.load('SVM_model/CRISPRi.pkl')
CRISPRko_SVM = joblib.load('SVM_model/CRISPRko.pkl')# open the file in the write mode
f = open('result.csv', 'w')# create the csv writer
sgRNAfile=open('sgRNA_input/sgRNA_input.txt','w')
writer = csv.writer(f)# write a row to the csv file
#Define the output file header
head='sgRNA_id'+','+'source'+','+'sgRNA_sequence'+','+'PAM'+','+'position'+','+'strand'+','+'CRISPRko predict'+','+'CRISPRi predict'+','+'chrom'+','+'start'+','+'end'+','+'gene'+','+'transcript'+','+'biotype'+','+'gene symbol'+','+'Region'+','+'Off-targets counts for0-1-2-3-4 mismatches'+'\n'
f.write(head)
#生成sg_list
sg_list=[]
sg_count=1
#first line
sgRNAfile.write(sys.path[0]+'/'+genome_fa+'\n')
#sgRNAfile.write(sys.path[0]+'/fasta/human_hg38.fa'+'\n')
#second line
sgRNAfile.write('NNNNNNNNNNNNNNNNNNNNNRG'+'\n')
for gene in sg_dict:
    for key1 in sg_dict[gene]:
        sgRNA=key1[0:20]+'NNN'
        sgRNAfile.write(sgRNA+' '+'3'+'\n')
sgRNAfile.close()

#deal with off-target data
cmd = " ./cas_offinder sgRNA_input/sgRNA_input.txt G  offtarget_result/output.txt"
os.system(cmd)

mismatch0_dict,mismatch1_dict,mismatch2_dict,mismatch3_dict=file_handle.mismatch('offtarget_result/output.txt')

for gene in sg_dict:
    sg_wholeseq = sg_dict[gene]['sequence']
    whole_dict = fa_dict
    whole_dict.update(anti_fa_dict)
    current_transcript = ''
    transcript_list=[]
    for transcript in whole_dict:
         matchseq = whole_dict[transcript]['sequence'].replace('\n', '').find(sg_wholeseq)
         if matchseq!=-1:
             transcript_list.append(transcript)
    del sg_dict[gene]['sequence']
    for key1 in sg_dict[gene]:
        muti_target=0
        sg_list.append(key1)
        if key1[0:20] in mismatch0_dict:
            off0_count=len(mismatch0_dict[key1[0:20]])
        elif key1[0:20] not in mismatch0_dict:
            off0_count=0
        if key1[0:20] in mismatch1_dict:
            off1_count=len(mismatch1_dict[key1[0:20]])
        elif key1[0:20] not in mismatch1_dict:
            off1_count=0
        if key1[0:20] in mismatch2_dict:
            off2_count=len(mismatch2_dict[key1[0:20]])
        elif key1[0:20] not in mismatch2_dict:
            off2_count=0
        if key1[0:20] in mismatch3_dict:
            off3_count=len(mismatch3_dict[key1[0:20]])
        elif key1[0:20] not in mismatch3_dict:
            off3_count=0
        off_situation=str(off0_count)+'-'+str(off1_count)+'-'+str(off2_count)+'-'+str(off3_count)
        CRISPRi_Predcit = str(CRISPRi_SVM.predict(file_handle.CRISPRifeature(key1[0:20]))[0])
        CRISPRko_Predcit = str(CRISPRko_SVM.predict(file_handle.CRISPRkofeature(key1[0:20]))[0])
        if ('TTTT' in key1[0:20]) or ((key1[0:20].count('G') + key1[0:20].count('C')) / len(key1[0:20])<0.20) or( (key1[0:20].count('G') + key1[0:20].count('C')) / len(key1[0:20])>0.80):
            CRISPRi_Predcit = str(0)
            CRISPRko_Predcit = str(0)
        if len(sg_dict[gene][key1])>1:
            muti_target=1
        for sg in sg_dict[gene][key1]:
            results=[]
            start=sg[1][0]
            end=sg[1][1]
            strand=sg[1][2]
            str_position=str(start)+'-'+str(end)
            sg_sequence=sg[0][0:20]
            pam=sg[0][20:23]
            for transcript in transcript_list:
                matchseq = whole_dict[transcript]['sequence'].replace('\n', '').find(sg_wholeseq)
                flag = 0
                if matchseq != -1:
                    position = int(whole_dict[transcript]['start']) + matchseq + start + 1
                    if transcript not in bed_dict:
                        region = 'Undocumented Promoter'
                    elif int(bed_dict[transcript]['start']) < position < int(bed_dict[transcript]['end']):
                        region = 'Promoter(<1kb)'
                    elif int(bed_dict[transcript]['start']) + 1000 < position < int(
                            bed_dict[transcript]['end']) + 500:
                        region = 'Promoter(1-2kb)'
                    else:
                        region = 'gene body'
                    sg_detail = str(
                        sg_count) + ',' + gene.split(' ')[0] + ',' + sg_sequence + ',' + pam + ',' + str_position + ',' + strand + ',' + CRISPRko_Predcit + ',' + CRISPRi_Predcit + ',' + 'Chr' + \
                                whole_dict[transcript][
                                    'chrom'] + ',' + str(position) + ',' + str(position + 19) + ',' + \
                                whole_dict[transcript]['gene'] + ',' + transcript + ',' + \
                                whole_dict[transcript]['biotype'] + ',' + whole_dict[transcript][
                                    'gene_symbol'] + ',' + region + '\n'
                    results.append(matchseq)
                    return_value += sg_detail
                    f.write(sg_detail)
                    flag += 1
            if flag == 0 and strand == '-':
                for transcript in anti_fa_dict:
                    sequence=anti_fa_dict[transcript]['sequence']
                    result=sequence.find(sg_sequence)
                    if result!=-1:
                        position=int(anti_fa_dict[transcript]['start']) + result+1
                        if transcript not in bed_dict:
                            region='Undocumented Promoter'
                        elif int(bed_dict[transcript]['start'])<position<int(bed_dict[transcript]['end']):
                            region='Promoter(<1kb)'
                        elif int(bed_dict[transcript]['start'])+1000<position<int(bed_dict[transcript]['end'])+500:
                            region = 'Promoter(1-2kb)'
                        else:
                            region = 'gene body'
                        sg_detail = str(sg_count)+','+gene.split(' ')[0] + ',' + sg_sequence + ',' +pam+','+str_position+','+ strand +','+CRISPRko_Predcit+','+CRISPRi_Predcit+ ',' +'Chr'+ anti_fa_dict[transcript][
                            'chrom'] + ',' + str(position)+',' + str(position+19) +','+anti_fa_dict[transcript]['gene'] + ','+transcript+','+ \
                                    anti_fa_dict[transcript]['biotype'] + ',' + anti_fa_dict[transcript][
                                        'gene_symbol'] + ',' + region+','+off_situation+'\n'
                        results.append(result)
                        f.write(sg_detail)

            if flag == 0 and strand == '+':
                for transcript in fa_dict:
                    sequence=fa_dict[transcript]['sequence']
                    result=sequence.find(sg_sequence)
                    if result!=-1:
                        position=int(fa_dict[transcript]['start']) + result+1
                        #print('Chrom'+fa_dict[transcript]['chrom']+':'+str(position))
                        if transcript not in bed_dict:
                            region='Undocumented Promoter'
                        elif int(bed_dict[transcript]['start'])<position<int(bed_dict[transcript]['end']):
                            region='Promoter(<1kb)'
                        elif int(bed_dict[transcript]['start'])+1000<position<int(bed_dict[transcript]['end'])+500:
                            region = 'Promoter(1-2kb)'
                        else:
                            region = 'gene body'
                        sg_detail = str(sg_count) +','+ gene.split(' ')[0] + ',' + sg_sequence + ',' +pam+','+str_position+','+ strand + ',' +CRISPRko_Predcit+','+CRISPRi_Predcit+ ','+ 'Chr'+str(fa_dict[transcript]['chrom']) + ',' + str(position)+','+str(position+19) +','+ fa_dict[transcript]['gene'] + ','+transcript+',' + \
                                    fa_dict[transcript]['biotype'] + ',' + fa_dict[transcript][
                                        'gene_symbol'] + ',' + region+','+off_situation+'\n'
                        f.write(sg_detail)
                        results.append(result)
            if results==[]:
                sg_detail=str(sg_count)+','+gene+','+sg_sequence+','+pam+','+str_position+','+strand+','+CRISPRko_Predcit+','+CRISPRi_Predcit+ ','+'Not Found'+','+'Not Found'+','+'Not Found'+','+'Not Found'+','+'Not Found'+','+'Not Found'+','+'Not Found'+','+'Not Found'+','+off_situation+'\n'
                f.write(sg_detail)
            sg_count+=1