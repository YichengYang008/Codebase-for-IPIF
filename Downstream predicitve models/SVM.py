# coding=utf-8
import numpy as np
import pandas as pd
from sklearn import svm
import math
from sklearn.metrics import f1_score

#train
train_UP_FHDI=pd.read_csv("data/UP_FHDI.csv")
train_NAIVE=pd.read_csv("data/NAIVE.csv")
train_DELETE=pd.read_csv("data/DELETE.csv")
train_GAIN=pd.read_csv("data/GAIN.csv")
train_HI_VAE=pd.read_csv("data/HI_VAE.csv")

Y=train_UP_FHDI.iloc[:,-1]
Y_DELETE=train_DELETE.iloc[:,-1]


train_UP_FHDI=train_UP_FHDI.iloc[:,:-1]
train_NAIVE=train_NAIVE.iloc[:,:-1]
train_DELETE=train_DELETE.iloc[:,:-1]
train_GAIN=train_GAIN.iloc[:,:-1]
train_HI_VAE=train_HI_VAE.iloc[:,:-1]



#original_train
original_UP_FHDI=pd.read_csv("original/original_UP_FHDI.csv",header=None)
original_NAIVE=pd.read_csv("original/original_NAIVE.csv",header=None)
original_DELETE=pd.read_csv("original/original_DELETE.csv",header=None)
original_GAIN=pd.read_csv("original/original_GAIN.csv",header=None)
original_HI_VAE=pd.read_csv("original/original_HI_VAE.csv",header=None)



Y_original=original_UP_FHDI.iloc[:,-1]
Y_original_DELETE=original_DELETE.iloc[:,-1]

original_UP_FHDI=original_UP_FHDI.iloc[:,:-1]
original_NAIVE=original_NAIVE.iloc[:,:-1]
original_DELETE=original_DELETE.iloc[:,:-1]
original_GAIN=original_GAIN.iloc[:,:-1]
original_HI_VAE=original_HI_VAE.iloc[:,:-1]





#fix
model_UP_FHDI = svm.SVC()
model_UP_FHDI.fit(train_UP_FHDI, Y)
model_original_UP_FHDI = svm.SVC()
model_original_UP_FHDI.fit(original_UP_FHDI, Y_original)

model_NAIVE = svm.SVC()
model_NAIVE.fit(train_NAIVE, Y)
model_original_NAIVE = svm.SVC()
model_original_NAIVE.fit(original_NAIVE, Y_original)

model_DELETE = svm.SVC()
model_DELETE.fit(train_DELETE, Y_DELETE)
model_original_DELETE = svm.SVC()
model_original_DELETE.fit(original_DELETE, Y_original_DELETE)

model_GAIN = svm.SVC()
model_GAIN.fit(train_GAIN, Y)
model_original_GAIN = svm.SVC()
model_original_GAIN.fit(original_GAIN, Y_original)

model_HI_VAE = svm.SVC()
model_HI_VAE.fit(train_HI_VAE, Y)
model_original_HI_VAE = svm.SVC()
model_original_HI_VAE.fit(original_HI_VAE, Y_original)


#test
test_UP_FHDI=pd.read_csv("test/test_UP_FHDI.csv",header=None)
test_NAIVE=pd.read_csv("test/test_NAIVE.csv",header=None)
test_DELETE=pd.read_csv("test/test_DELETE.csv",header=None)
test_GAIN=pd.read_csv("test/test_GAIN.csv",header=None)
test_HI_VAE=pd.read_csv("test/test_HI_VAE.csv",header=None)

Y_true=test_UP_FHDI.iloc[:,-1]


test_UP_FHDI=test_UP_FHDI.iloc[:,:-1]
test_NAIVE=test_NAIVE.iloc[:,:-1]
test_DELETE=test_DELETE.iloc[:,:-1]
test_GAIN=test_GAIN.iloc[:,:-1]
test_HI_VAE=test_HI_VAE.iloc[:,:-1]


#pred
Y_UP_FHDI = model_UP_FHDI.predict(test_UP_FHDI)
Y_UP_FHDI_original= model_original_UP_FHDI.predict(test_UP_FHDI)

Y_NAIVE = model_NAIVE.predict(test_NAIVE)
Y_NAIVE_original = model_original_NAIVE.predict(test_NAIVE)

Y_DELETE = model_DELETE.predict(test_DELETE)
Y_DELETE_original = model_original_DELETE.predict(test_DELETE)

Y_GAIN = model_GAIN.predict(test_GAIN)
Y_GAIN_original = model_original_GAIN.predict(test_GAIN)

Y_HI_VAE = model_HI_VAE.predict(test_HI_VAE)
Y_HI_VAE_original = model_original_HI_VAE.predict(test_HI_VAE)


#F1
print("SVM")
f1_UP_FHDI = f1_score(Y_true, Y_UP_FHDI)
f1_UP_FHDI_original = f1_score(Y_true, Y_UP_FHDI_original)

f1_NAIVE = f1_score(Y_true, Y_NAIVE)
f1_NAIVE_original = f1_score(Y_true, Y_NAIVE_original)

f1_DELETE = f1_score(Y_true, Y_DELETE)
f1_DELETE_original = f1_score(Y_true, Y_DELETE_original)

f1_GAIN = f1_score(Y_true, Y_GAIN)
f1_GAIN_original = f1_score(Y_true, Y_GAIN_original)

f1_HI_VAE = f1_score(Y_true, Y_HI_VAE)
f1_HI_VAE_original = f1_score(Y_true, Y_HI_VAE_original)

print("F1-UP_FHDI:%f" %f1_UP_FHDI)
print("F1-UP_FHDI_original:%f" %f1_UP_FHDI_original)

print("F1-NAIVE:%f" %f1_NAIVE)
print("F1-NAIVE_origianl:%f" %f1_NAIVE_original)

print("F1-DELETE:%f" %f1_DELETE)
print("F1-DELETE_original:%f" %f1_DELETE_original)

print("F1-GAIN:%f" %f1_GAIN)
print("F1-GAIN_original:%f" %f1_GAIN_original)

print("F1-HI_VAE:%f" %f1_HI_VAE)
print("F1-HI_VAE_origianl:%f" %f1_HI_VAE_original)



nf1_UP_FHDI=f1_UP_FHDI/f1_UP_FHDI_original
nfi_naive=f1_NAIVE/f1_NAIVE_original
nf1_DELETE=f1_DELETE/f1_DELETE_original
nf1_GAIN=f1_GAIN/f1_GAIN_original
nf1_HI_VAE=f1_HI_VAE/f1_HI_VAE_original


print("NF1-UP_FHDI:%f" %nf1_UP_FHDI)
print("NF1-NAIVE:%f" %nfi_naive)
print("NF1-DELETE:%f" %nf1_DELETE)
print("NF1-GAIN:%f" %nf1_GAIN)
print("NF1-HI_VAE:%f" %nf1_HI_VAE)