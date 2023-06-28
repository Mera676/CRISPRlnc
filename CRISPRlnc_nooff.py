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

fa_dict,anti_fa_dict=file_handle.read_fa(ncrna_fa)    
bed_dict=file_handle.read_Promoter(bed_fa)
sg_dict=file_handle.read_sg(sg_fa)


print(ncrna_fa+'\n')
print(bed_fa+'\n')
print(sg_fa)




CRISPRi_SVM= joblib.load('SVM_model/CRISPRi.pkl')
CRISPRko_SVM = joblib.load('SVM_model/CRISPRko.pkl')# open the file in the write mode
f = open('CRISPRlnc_nooff_result.csv', 'w')# create the csv writer
writer = csv.writer(f)# write a row to the csv file
#Define the output file header
head='sgRNA_id'+','+'source'+','+'sgRNA_sequence'+','+'PAM'+','+'position'+','+'strand'+','+'CRISPRko predict'+','+'CRISPRi predict'+','+'chrom'+','+'start'+','+'end'+','+'gene'+','+'transcript'+','+'biotype'+','+'gene symbol'+','+'Region'+'\n'
f.write(head)
#sg_list
sg_list=[]
sg_count=1


for gene in sg_dict:
    sg_wholeseq = sg_dict[gene]['sequence'].replace('\n', '')
    find_flag = 0#表示该基因是否能在fa文件中查找到序列
    #先来匹配整段基因
    whole_dict = fa_dict
    whole_dict.update(anti_fa_dict)
    current_transcript = ''
    transcript_list=[]
    for transcript in whole_dict:
         matchseq = whole_dict[transcript]['sequence'].replace('\n', '').find(sg_wholeseq)
         if matchseq!=-1:
             transcript_list.append(transcript)
    print(transcript_list)
    del sg_dict[gene]['sequence']
    for key1 in sg_dict[gene]:
        muti_target=0
        sg_list.append(key1)
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
            #先来匹配整段基因
            whole_dict = fa_dict
            whole_dict.update(anti_fa_dict)
            flag = 0
            for transcript in transcript_list:
                matchseq = whole_dict[transcript]['sequence'].replace('\n', '').find(sg_wholeseq)#改成只执行一次，这段代码的含义是查找一整段的基因，如果找到了，整段都不用重新再找
                if matchseq != -1:
                    find_flag = 1
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
                    #return_value += sg_detail
                    f.write(sg_detail)
                    flag += 1
            
            if flag== 0 and strand == '-':
                print('hello')
                for transcript in anti_fa_dict:
                    sequence=anti_fa_dict[transcript]['sequence']
                    result=sequence.find(sg_sequence)
                    if result!=-1:
                        position=int(anti_fa_dict[transcript]['start']) + result+1
                        if transcript not in bed_dict:
                            region='Not Found'
                        elif int(bed_dict[transcript]['start'])<position<int(bed_dict[transcript]['end']):
                            region='Promoter(<1kb)'
                        elif int(bed_dict[transcript]['start'])+1000<position<int(bed_dict[transcript]['end'])+500:
                            region = 'Promoter(1-2kb)'
                        else:
                            region = 'gene body'
                        sg_detail = str(sg_count)+','+gene.split(' ')[0] + ',' + sg_sequence + ',' +pam+','+str_position+','+ strand +','+CRISPRko_Predcit+','+CRISPRi_Predcit+ ',' +'Chr'+ anti_fa_dict[transcript][
                            'chrom'] + ',' + str(position)+',' + str(position+19) +','+anti_fa_dict[transcript]['gene'] + ','+transcript+','+ \
                                    anti_fa_dict[transcript]['biotype'] + ',' + anti_fa_dict[transcript][
                                        'gene_symbol'] + ',' + region+'\n'
                        results.append(result)
                        f.write(sg_detail)

            if flag == 0 and strand == '+':
                print('hi')
                for transcript in fa_dict:
                    sequence=fa_dict[transcript]['sequence']
                    result=sequence.find(sg_sequence)
                    if result!=-1:
                        position=int(fa_dict[transcript]['start']) + result+1
                        #print('Chrom'+fa_dict[transcript]['chrom']+':'+str(position))
                        if transcript not in bed_dict:
                            region='Not Found'
                        elif int(bed_dict[transcript]['start'])<position<int(bed_dict[transcript]['end']):
                            region='Promoter(<1kb)'
                        elif int(bed_dict[transcript]['start'])+1000<position<int(bed_dict[transcript]['end'])+500:
                            region = 'Promoter(1-2kb)'
                        else:
                            region = 'gene body'
                        sg_detail = str(sg_count) +','+ gene.split(' ')[0] + ',' + sg_sequence + ',' +pam+','+str_position+','+ strand + ',' +CRISPRko_Predcit+','+CRISPRi_Predcit+ ','+ 'Chr'+str(fa_dict[transcript]['chrom']) + ',' + str(position)+','+str(position+19) +','+ fa_dict[transcript]['gene'] + ','+transcript+',' + \
                                    fa_dict[transcript]['biotype'] + ',' + fa_dict[transcript][
                                        'gene_symbol'] + ',' + region+'\n'
                        f.write(sg_detail)
                        results.append(result)
            if results==[]:
                sg_detail=str(sg_count)+','+gene+','+sg_sequence+','+pam+','+str_position+','+strand+','+CRISPRko_Predcit+','+CRISPRi_Predcit+ ','+'Not Found'+','+'Not Found'+','+'Not Found'+','+'Not Found'+','+'Not Found'+','+'Not Found'+','+'Not Found'+','+'Not Found'+'\n'
                f.write(sg_detail)
            sg_count+=1