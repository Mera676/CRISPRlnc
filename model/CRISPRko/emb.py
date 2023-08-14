from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from imblearn.over_sampling import SMOTE, ADASYN
from sklearn.svm import SVC
import xgboost as xgb
from numpy import loadtxt
from numpy import sort
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score,precision_score,recall_score,f1_score,cohen_kappa_score
from sklearn.feature_selection import SelectFromModel
import warnings
warnings.filterwarnings("ignore")
import joblib
import pickle
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
def CRISPRkofeature(sg):
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
    feature_list=[energy , position(sg[3]), position(sg[19]), position(sg[4]) ,position(sg[6]) , gc10_20 , gc0_10, position(sg[13]), position(sg[11]), hairpin, position(sg[1]), gc10_15, gc0_5, position(sg[7]), stem, position(sg[16])]
    return (np.array(feature_list)).reshape(1,16)



train_txt = sys.argv[1]
x_emb = []
y_emb = []
train_data = dict()
with open(train_txt,'r') as inputfile:
    contents = inputfile.read().split('\n')
for content in contents:
    x_emb.append(sefeature(content.split(',')[0]).ravel().tolist())
    y_emb.append([int(content.split(',')[1])])
x_emb = np.array(x_emb)
y_emb = np.array(y_emb)
oversample = SMOTE()
x_emb, y_emb = oversample.fit_resample(x_emb, y_emb)
train_data['X'] =x_emb
train_data['Y'] =y_emb
f_save = open('data/CRISPRko_train_data.pkl', 'wb')
pickle.dump(train_data, f_save)
f_save.close()
        
        
        
