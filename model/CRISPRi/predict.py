import warnings
warnings.filterwarnings("ignore")
import joblib
import pandas as pd
import numpy as np
from numpy import loadtxt
from xgboost import XGBClassifier
from xgboost import plot_importance
import RNA
import os, sys
import re



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
    seq=sg[10:20]
    return (seq.count('G')+seq.count('C'))/len(seq)
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
    feature_list=[position(sg[9]),position(sg[19]),position(sg[14]),hairpin,position(sg[5]),gc10_15,position(sg[12]),position(sg[11]),gc10_20,position(sg[3]),stem,position(sg[15]),position(sg[0]),energy,position(sg[18]),gc,position(sg[8]),position(sg[4])]
    return (np.array(feature_list)).reshape(1,18)
def wholefeature(sg):
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
    feature_list=[position(sg[0]),position(sg[1]),position(sg[2]),position(sg[3]),position(sg[4]),position(sg[5]),position(sg[6]),position(sg[7]),position(sg[8]),position(sg[9]),position(sg[10]),position(sg[11]),position(sg[12]),position(sg[13]),position(sg[14]),position(sg[15]),position(sg[16]),position(sg[17]),position(sg[18]),position(sg[19]),gc,gc0_10,gc10_20,gc0_5,gc5_10,gc10_15,gc15_20,energy,stem,hairpin,free]
    return (np.array(feature_list)).reshape(1,31)
    


with open('./CRISPRi_result.tsv','w') as outputfile:
    outputfile.write('sgRNA'+'\t'+'predict'+'\n')
sg_txt = sys.argv[1]
CRISPRko_SVM= joblib.load('./save_model/CRISPRi.pkl')
with open(sg_txt,'r')as input_file:
    sgRNAs = input_file.read().split('\n')
    for sgRNA in sgRNAs:
        sgRNA_emb = CRISPRifeature(sgRNA[0:20])
        CRISPRi_Predcit = str(CRISPRko_SVM.predict(sgRNA_emb)[0])
        with open('./CRISPRi_results.tsv','a') as outputfile:
            outputfile.write(sgRNA+'\t'+CRISPRi_Predcit+'\n')
        print('CRISPRlnc software prediction is completed, please view the results on the CRISPRi_results.tsv file.')
        
        
