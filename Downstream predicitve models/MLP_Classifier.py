# coding=utf-8
import pandas as pd
from numpy import array
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import f1_score

#train
train_UP_FHDI=pd.read_csv("data/final_daty_UP_FHDI_binary.csv",header=None)
#train_GAIN=pd.read_csv("data/imputed.txt",sep=" ",header=None)
train_GAIN=pd.read_csv("data/GAIN.csv",header=None)
train_NAIVE=pd.read_csv("data/NAIVE.csv",header=None)
train_DELETE=pd.read_csv("data/DELETE.csv",header=None)
train_HI_VAE=pd.read_csv("data/model_HIVAE_inputDropout_swarm_Missing20_1_z2_y5_s10_batch100_data_reconstruction.csv",header=None)
train_original=pd.read_csv("data/0.7data.csv",header=None)


Y=train_UP_FHDI.iloc[:,-1]
Y_delete=train_DELETE.iloc[:,-1]


train_UP_FHDI=train_UP_FHDI.iloc[:,:-1]
train_NAIVE=train_NAIVE.iloc[:,:-1]
train_DELETE=train_DELETE.iloc[:,:-1]
train_GAIN=train_GAIN.iloc[:,:-1]
train_HI_VAE=train_HI_VAE.iloc[:,:-1]
train_original=train_original.iloc[:,:-1]

n_features=train_UP_FHDI.shape[1]


#test
test=pd.read_csv("data/0.3data.csv",header=None)
Y_true=test.iloc[:,-1]
test=test.iloc[:,:-1]

#fix
model_UP_FHDI=MLPClassifier(hidden_layer_sizes=(50, 50, 50,25),random_state=7)
model_UP_FHDI.fit(train_UP_FHDI, Y)

model_NAIVE=MLPClassifier(hidden_layer_sizes=(50, 50, 50,25),random_state=7)
model_NAIVE.fit(train_NAIVE, Y)

model_DELETE=MLPClassifier(hidden_layer_sizes=(50, 50, 50,25),random_state=7)
model_DELETE.fit(train_DELETE, Y_delete)

model_GAIN=MLPClassifier(hidden_layer_sizes=(50, 50, 50,25),random_state=7)
model_GAIN.fit(train_GAIN, Y)

model_HI_VAE=MLPClassifier(hidden_layer_sizes=(50, 50, 50,25),random_state=7)
model_HI_VAE.fit(train_HI_VAE, Y)

model_original=MLPClassifier(hidden_layer_sizes=(50, 50, 50,25),random_state=7)
model_original.fit(train_original, Y)


#pred
Y_UP_FHDI = model_UP_FHDI.predict(test)
Y_GAIN = model_GAIN.predict(test)
Y_NAIVE=model_NAIVE.predict(test)
Y_DELETE=model_DELETE.predict(test)
Y_HI_VAE=model_HI_VAE.predict(test)
Y_original = model_original.predict(test)


#F1
print("MLP")
f1_UP_FHDI = f1_score(Y_true, Y_UP_FHDI)
f1_NAIVE = f1_score(Y_true, Y_NAIVE)
f1_DELETE = f1_score(Y_true, Y_DELETE)
f1_GAIN = f1_score(Y_true, Y_GAIN)
f1_HI_VAE = f1_score(Y_true, Y_HI_VAE)
f1_original = f1_score(Y_true, Y_original)
print("F1-UP_FHDI:%f" %f1_UP_FHDI)
print("F1-NAIVE:%f" %f1_NAIVE)
print("F1-DELETE:%f" %f1_DELETE)
print("F1-GAIN:%f" %f1_GAIN)
print("F1-HI_VAE:%f" %f1_HI_VAE)
print("F1-original:%f" %f1_original)


nf1_UP_FHDI=f1_UP_FHDI/f1_original
nf1_NAIVE=f1_NAIVE/f1_original
nf1_DELETE=f1_DELETE/f1_original
nf1_GAIN=f1_GAIN/f1_original
nf1_HI_VAE=f1_HI_VAE/f1_original
print("NF1-UP_FHDI:%f" %nf1_UP_FHDI)
print("NF1-NAIVE:%f" %nf1_NAIVE)
print("NF1-DELETE:%f" %nf1_DELETE)
print("NF1-GAIN:%f" %nf1_GAIN)
print("NF1-HI_VAE:%f" %nf1_HI_VAE)
