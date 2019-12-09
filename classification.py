# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# Support Vector Machine (SVM)

# Importing the libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Importing the dataset
dataset = pd.read_csv('feature.csv')
X = dataset.iloc[:, :12].values
y = dataset.iloc[:, 12].values

from sklearn.model_selection import LeaveOneOut
loo=LeaveOneOut()
loo.get_n_splits(X)
k=0
for train_index,test_index in loo.split(X):
    X_train,X_test=X[train_index],X[test_index]
    y_train,y_test=y[train_index],y[test_index]
    from sklearn.preprocessing import StandardScaler
    sc = StandardScaler()
    X_train = sc.fit_transform(X_train)
    X_test = sc.transform(X_test)
    '''from sklearn.ensemble import RandomForestClassifier
    classifier = RandomForestClassifier(n_estimators = 1000, criterion = 'entropy')'''
    from sklearn.svm import SVC
    classifier = SVC(kernel = 'rbf')
    
    
    classifier.fit(X_train, y_train)
    y_pred = classifier.predict(X_test)
    print("Test:",y_test)
    print("Pred:",y_pred)
  
    if y_test==y_pred:
        
        k+=1
    
   
    

