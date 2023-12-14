#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on March 30 2023

#@author: Wenjie Ma

#check package version
#pip3 list

python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import numpy

from sklearn.ensemble import RandomForestClassifier
from sklearn import ensemble
from sklearn import svm, datasets
from sklearn.metrics import auc
from sklearn.metrics import plot_roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import ShuffleSplit


from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_auc_score, roc_curve
import pylab as plot
plt.style.use('seaborn-white')
import seaborn as sns
sns.set()



from sklearn.metrics import RocCurveDisplay



#adapted from https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html



##############################################################################3\
#set working directory
os.chdir('/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/')

#############################################
# Load in data
df = pd.read_csv("./mbx/data/rf/all_python.csv")
#print(df.shape)
#list(df.columns.values)

# Extract the labels  - do this before setting up as array
# y = np.array(input_data.pop('caco'))
y=np.where(df['caco']=='Divert', 1, 0)

X1=np.array(df.iloc[:,15:204]) # taxa 
X2=np.array(df.iloc[:,[9,13,14]]) # age,fiber,BMI
X3=np.array(df.iloc[:,204:756]) # metablites
cols = ([9] + list(range(13,756)))
X4=np.array(df.iloc[:,cols]) # age,fiber,BMI.taxa,metabolites
#set up as array


# Full dataset: https://www.kaggle.com/cdc/behavioral-risk-factor-surveillance-system

# #############################################################################
###########################Microbiome##########################################
###############################################################################
# Classification and ROC analysis
# RANDOM FOREST

# Run classifier with cross-validation and plot ROC curves
cv = ShuffleSplit(n_splits = 100, test_size = 0.25, random_state = 0)
#classifier = svm.SVC(kernel=‘linear’, probability=True,
#                     random_state=random_state)
classifier = ensemble.RandomForestClassifier(class_weight="balanced")

tprs = []
aucs = []
imp =  []

mean_fpr = np.linspace(0, 1, 100)

plt.figure(figsize=(6.5, 6.5)).clf()
plt.rcParams["axes.edgecolor"] = "black"
sns.set_style('white')

fig, ax = plt.subplots(figsize=(6, 6))
for i, (train, test) in enumerate(cv.split(X4, y)):
    classifier.fit(X4[train], y[train])
    viz = plot_roc_curve(
        classifier,
        X4[test],
        y[test],
        name='ROC fold {}'.format(i),
        alpha=0.3,
        lw=1,
        ax=ax,
    )
    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    interp_tpr[0] = 0.0
    interp_importance = classifier.feature_importances_
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)
    imp.append(interp_importance)
    
ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)

mean_importance = np.mean(imp, axis = 0)

ax.plot(
    mean_fpr,
    mean_tpr,
    color="b",
    label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
    lw=2,
    alpha=0.8,
)

std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
ax.fill_between(
    mean_fpr,
    tprs_lower,
    tprs_upper,
    color="grey",
    alpha=0.2,
    label=r"$\pm$ 1 std. dev.",
)
handles, labels = ax.get_legend_handles_labels()

ax.set(
    xlim=[-0.05, 1.05],
    ylim=[-0.05, 1.05],
    xlabel="False Positive Rate",
    ylabel="True Positive Rate",
    title="ROC of diverticulitis vs. controls",
)
ax.axis("square")
ax.legend(loc="lower right")
ax.legend(handles[-3:], labels[-3:])
#plt.show()


plt.savefig(r'/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/mbx/results/rf/roc_all_no_label.png',dpi=300)



