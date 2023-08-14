from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
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
import pandas as pd
import numpy as np
# plot feature importance using built-in function
from numpy import loadtxt
from xgboost import XGBClassifier
from xgboost import plot_importance
from matplotlib import pyplot
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score
import RNA
import os, sys
import pickle

def featuresort(X_train, X_test, y_train, y_test):
    #X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.3, random_state=18)
    # fit model on all training data
    model = XGBClassifier()
    feature = {}
    model.fit(X_train, y_train)
    # make predictions for test data and evaluate
    y_pred = model.predict(X_test)
    predictions = [round(value) for value in y_pred]
    accuracy = accuracy_score(y_test, predictions)
    print("Accuracy: %.2f%%" % (accuracy * 100.0))
    # Fit model using each importance as a threshold
    thresholds = sort(model.feature_importances_)
    print(model.feature_importances_)
    count = 0
    for fea in model.feature_importances_:
        feature['f'+str(count)] = fea
        count += 1
    for thresh in thresholds:
        # select features using threshold
        selection = SelectFromModel(model, threshold=thresh, prefit=True)
        select_X_train = selection.transform(X_train)
        # train model
        selection_model = XGBClassifier()
        selection_model.fit(select_X_train, y_train)
        # eval model
        select_X_test = selection.transform(X_test)
        y_pred = selection_model.predict(select_X_test)
        predictions = [round(value) for value in y_pred]
        accuracy = accuracy_score(y_test, predictions)
        print("Thresh=%.3f, n=%d, Accuracy: %.2f%%" % (thresh, select_X_train.shape[1], accuracy*100.0))
    fig, ax = plt.subplots(figsize=(16,8))
    pyplot.bar(range(len(model.feature_importances_)), model.feature_importances_)
    pyplot.show()
  
    return feature

def featureselect(feature,n):
    feature1 = sorted(feature.items(), key=lambda x: x[1], reverse=True)
    return feature1[:n]

def featurename(fea_dict):
    feature_list = []
    for item in fea_dict:
        print(item[0])
        name = feature_importance[feature_importance['num'] == item[0]]['features'].tolist()[0]
        feature_list.append(name)
    return feature_list 
    
def select_c_function(i,X,y):
    svm_model = SVC(kernel='rbf',C=i,gamma = 0.125)
    precision_score = cross_val_score(svm_model, X,y, scoring='precision', cv=10)
    recall_score = cross_val_score(svm_model, X,y, scoring='recall', cv=10)
    accuracy_score = cross_val_score(svm_model, X,y, scoring='accuracy', cv=10)
    f1_score = cross_val_score(svm_model,X,y, scoring='f1', cv=10)
    return accuracy_score.mean()
    
def select_gamma_function(i,X,y,c):
    svm_model = SVC(kernel='rbf',C=c,gamma = i)
    precision_score = cross_val_score(svm_model, X,y, scoring='precision', cv=10)
    recall_score = cross_val_score(svm_model, X,y, scoring='recall', cv=10)
    accuracy_score = cross_val_score(svm_model, X,y, scoring='accuracy', cv=10)
    f1_score = cross_val_score(svm_model,X,y, scoring='f1', cv=10)
    return accuracy_score.mean()
    
def bestC(X,y):
    max_score = select_c_function(1,X,y)
    max_count =1 
    c_range =range(2,10000)   
    for i in c_range:
        #print(i)
        avg_score = select_c_function(i,X,y)
        if avg_score > max_score:
            max_score = avg_score
            max_count = i
    return max_count,max_score


def bestGamma(X,y,c):
    max_score = select_gamma_function(1,X,y,c)
    max_count =1 
    gamma_2d_range = [3.01454,0.0001,0.00125,0.001, 0.01,0.005,0.0125,0.125,0.1,0.5,0.15,0.25,0.35,0.5,0.45,0.95,0.42]   
    for i in gamma_2d_range:
        avg_score = select_gamma_function(i,X,y,c)
        if avg_score > max_score:
            max_score = avg_score
            max_count = i
    return max_count,max_score
    
train_data = sys.argv[1]
f_read = open(train_data, 'rb')
data_dict = pickle.load(f_read)
f_read.close()
X = data_dict['X']
y = data_dict['Y']
print('Searching for optimal parameters.....')
bestc,maxc = bestC(X,y)
bestgamma,maxgamma =  bestGamma(X,y,bestc)
svm_model = SVC(kernel='rbf',C=bestc,gamma =bestgamma,decision_function_shape='ovr')
svm_model.fit(X,y)
accuracy_score = cross_val_score(svm_model,X,y, scoring='accuracy', cv=10).mean()
precision_score = cross_val_score(svm_model,X,y, scoring='precision', cv=10).mean()
recall_score = cross_val_score(svm_model,X,y, scoring='recall', cv=10).mean()
f1_score = cross_val_score(svm_model,X,y, scoring='f1', cv=10).mean()
print(f"The average accuracy of ten-fold cross-validation:{accuracy_score}")
joblib.dump(svm_model,'save_model/CRISPRi.pkl')
print('The trained model is stored in <./save_model/CRISPRi.pkl>')