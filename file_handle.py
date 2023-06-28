# coding:utf-8
import os, sys
import pickle
import regex
import numpy as np
import RNA
import re

#deal with fa file
def read_fa(filename):
    fa_dict={}
    anti_fa_dict={}
    transcript_list=[]
    with open (filename,'r') as fa_file:
        content=fa_file.read()
    fa_list = content.split('>')
    for i in range(1, len(fa_list)):
        #print(fa_list[i])
        info=fa_list[i]
        info_list = info.split('\n')[0].split(' ')
        transcript=info_list[0].split('.')[0]
        chrom=info_list[2].split(':')[2]
        start=info_list[2].split(':')[3]
        end=info_list[2].split(':')[4]
        strand=info_list[2].split(':')[5]
        gene=info_list[3].split(':')[1].split('.')[0]
        biotype=info_list[4].split(':')[1]
        if len(info_list) < 7:
            gene_symbol = 'No symbol'
        else:
            if (info_list[6].find('description') != -1):
                gene_symbol = 'No symbol'
            else:
                gene_symbol = info_list[6].split(':')[1]
        sequence=info.split('\n',1)[1].replace('\n','')
        transcript_list=[chrom,start,end,strand,gene,biotype,gene_symbol,sequence]
        if str(strand) == '-1':
            anti_fa_dict[transcript]={}
            anti_fa_dict[transcript]['chrom']=chrom
            anti_fa_dict[transcript]['start']=start
            anti_fa_dict[transcript]['end']=end
            anti_fa_dict[transcript]['strand']=strand
            anti_fa_dict[transcript]['gene']=gene
            anti_fa_dict[transcript]['biotype']=biotype
            anti_fa_dict[transcript]['gene_symbol']=gene_symbol
            anti_fa_dict[transcript]['sequence']=sequence
        elif str(strand) == '1':
            fa_dict[transcript]={}
            fa_dict[transcript]['chrom']=chrom
            fa_dict[transcript]['start']=start
            fa_dict[transcript]['end']=end
            fa_dict[transcript]['strand']=strand
            fa_dict[transcript]['gene']=gene
            fa_dict[transcript]['biotype']=biotype
            fa_dict[transcript]['gene_symbol']=gene_symbol
            fa_dict[transcript]['sequence']=sequence
    return fa_dict,anti_fa_dict
#deal with promoter bed_file
def read_Promoter(filename):
    promoter_dict={}
    with open (filename,'r') as bed_file:
        content=bed_file.read()
    promoter_list=content.split('\n')
    for promoter in promoter_list:
        info_list=promoter.split('\t')
        transcript=info_list[3]
        chrom=info_list[0]
        start=info_list[1]
        end=info_list[2]
        gene=info_list[4]
        strand=info_list[5]
        promoter_dict[transcript]={}
        promoter_dict[transcript]['chrom']=chrom
        promoter_dict[transcript]['gene']=gene
        promoter_dict[transcript]['start']=start
        promoter_dict[transcript]['end']=end
        promoter_dict[transcript]['strand']=strand
    return promoter_dict
def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)
#deal with sg_txt      
def read_sg(file_name):
    length=21
    with open (file_name,'r') as sg_file:
        content=sg_file.read()  
    sg_dict={}
    gene_list = content.split('>')
    del gene_list[0]
    for gene in gene_list:
        gene_name=gene.split('\n',1)[0]
        gene_sequence=gene.split('\n',1)[1].replace('\n','')
        if gene_name not in sg_dict:
            sg_dict[gene_name] = {}
            sg_dict[gene_name]['sequence'] = gene_sequence
            gene_sequence2=reverse_complement(gene_sequence)
            match_objs = regex.finditer('[ATCG]{{{0}}}[AG]G'.format(length), gene_sequence)
            match2_objs = regex.finditer('[ATCG]{{{0}}}[AG]G'.format(length), gene_sequence2)
            for match_obj in match_objs:
                match_start, match_end = match_obj.span()
                sg_info=(match_obj.group(), (match_start, match_end, '+'))
                if match_obj.group() not in sg_dict[gene_name]:
                    sg_dict[gene_name][match_obj.group()]=[sg_info]
                else:
                    sg_dict[gene_name][match_obj.group()].append(sg_info)
            for match_obj2 in match2_objs:
                match_start, match_end = match_obj2.span()
                sg_info=(match_obj2.group(), (len(gene_sequence)-match_end, len(gene_sequence)-match_start, '-'))
                if match_obj2.group() not in sg_dict[gene_name]:
                    sg_dict[gene_name][match_obj2.group()]=[sg_info]
                else:
                    sg_dict[gene_name][match_obj2.group()].append(sg_info)
    return sg_dict
def reverse_complement(sequence):
    trantab = str.maketrans('ACGTacgt', 'TGCAtgca')
    sequence = sequence.translate(trantab)
    return sequence[::-1]
def position(base):
    if base=='A':
        return 1
    elif base=='C':
        return 2
    elif base=='G':
        return 3
    elif base=='T':
        return 4
def GC(sg):
    seq=sg
    return (seq.count('G')+seq.count('C'))/len(seq)
def GC0_10(sg):
    seq=sg[0:10]
    return (seq.count('G')+seq.count('C'))/len(seq)
def GC10_20(sg):
    seq = sg[10:20]
    return (seq.count('G') + seq.count('C')) / len(seq)
def GC0_5(sg):
    seq=sg[0:5]
    return (seq.count('G')+seq.count('C'))/len(seq)
def GC5_10(sg):
    seq=sg[5:10]
    return (seq.count('G')+seq.count('C'))/len(seq)
def GC10_15(sg):
    seq=sg[10:15]
    return (seq.count('G')+seq.count('C'))/len(seq)
def GC15_20(sg):
    seq=sg[15:20]
    return (seq.count('G')+seq.count('C'))/len(seq)
def gestructure(sg):
    seq=sg.replace('T','U')
    fc = RNA.fold_compound(seq)
    (mfe_struct, mfe) = fc.mfe()
    return mfe_struct
def geenergy(sg):
    seq=sg.replace('T','U')
    fc = RNA.fold_compound(seq)
    (mfe_struct, mfe) = fc.mfe()
    return mfe
def Stem(sg):
    seq=sg.replace('T','U')
    fc = RNA.fold_compound(seq)
    (mfe_struct, mfe) = fc.mfe()
    stru = mfe_struct
    stem= stru.count('(')*2
    return stem
def Hairpin(sg):
    seq=sg.replace('T','U')
    fc = RNA.fold_compound(seq)
    (mfe_struct, mfe) = fc.mfe()
    stru = mfe_struct
    Hairpin = 0
    pattern = re.compile(r'\(\.+\)')
    results = pattern.findall(stru)
    for result in results:
        Hairpin += result.count('.')
    return Hairpin
def Free(sg):
    seq=sg.replace('T','U')
    fc = RNA.fold_compound(seq)
    (mfe_struct, mfe) = fc.mfe()
    stru = mfe_struct
    cycle = 0
    pattern = re.compile(r'\(+\.+\)+')
    results = pattern.findall(stru)
    for result in results:
        cycle += len(result)
    free = 20 - cycle
    return free

def CRISPRifeature(sg):
    sg=sg
    for i in range(0,len(sg)):
        value=position(sg[i])
        exec(f'feature_{i+1} = value')
    gc = GC(sg)
    gc0_10 = GC0_10(sg)
    gc10_20 = GC10_20(sg)
    gc0_5 = GC0_5(sg)
    gc5_10 = GC5_10(sg)
    gc10_15 = GC10_15(sg)
    gc15_20 = GC15_20(sg)
    energy = geenergy(sg)
    stem = Stem(sg)
    hairpin = Hairpin(sg)
    free = Free(sg)
    feature_list = [position(sg[9]), position(sg[19]), position(sg[14]), hairpin, position(sg[5]), gc10_15,
                    position(sg[12]), position(sg[11]), gc10_20, position(sg[3]), stem, position(sg[15]),
                    position(sg[0]), energy, position(sg[18]), gc, position(sg[8]), position(sg[4])]
    return (np.array(feature_list)).reshape(1, 18)

def CRISPRkofeature(sg):
    sg=sg
    gc = GC(sg)
    gc0_10 = GC0_10(sg)
    gc10_20 = GC10_20(sg)
    gc0_5 = GC0_5(sg)
    gc5_10 = GC5_10(sg)
    gc10_15 = GC10_15(sg)
    gc15_20 = GC15_20(sg)
    energy = geenergy(sg)
    stem = Stem(sg)
    hairpin = Hairpin(sg)
    free = Free(sg)
    feature_list = [energy, position(sg[3]), position(sg[19]), position(sg[4]), position(sg[6]), gc10_20, gc0_10,
                    position(sg[13]), position(sg[11]), hairpin, position(sg[1]), gc10_15, gc0_5, position(sg[7]), stem,
                    position(sg[16])]
    return (np.array(feature_list)).reshape(1, 16)
    
def mismatch(filename):
    mismatch0_dict={}
    mismatch1_dict={}
    mismatch2_dict={}
    mismatch3_dict={}
    with open (filename,'r')as offtarget_file:
        content=offtarget_file.read()
    off_list=content.split('\n')
    del off_list[-1]
    for off in off_list:
        sequence=off[0:20]
        if off[-1]=='0':
            mis=off.split('\t')[-3]
            mismatch0_dict.setdefault(sequence,[]).append(mis)
            #mismatch0_dict[sequence].append(mis)
        elif off[-1]=='1':
            mis=off.split('\t')[-3]
            mismatch1_dict.setdefault(sequence,[]).append(mis)
        elif off[-1]=='2':
            mis=off.split('\t')[-3]
            mismatch2_dict.setdefault(sequence,[]).append(mis)
        elif off[-1]=='3':
            mis=off.split('\t')[-3]
            mismatch3_dict.setdefault(sequence,[]).append(mis)
    return mismatch0_dict,mismatch1_dict,mismatch2_dict,mismatch3_dict