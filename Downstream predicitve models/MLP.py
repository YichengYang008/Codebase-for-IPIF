# coding=utf-8
import pandas as pd
import numpy as np
from sklearn.neural_network import MLPRegressor
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import math

#train
train_UP_FHDI=pd.read_csv("data/final_daty_UP_FHDI_binary.csv",header=None)
train_GAIN=pd.read_csv("data/imputed.txt",sep=" ",header=None)
#train_GAIN=pd.read_csv("data/GAIN.csv",header=None)
train_NAIVE=pd.read_csv("data/NAIVE.csv",header=None)
train_DELETE=pd.read_csv("data/DELETE.csv",header=None)
train_HI_VAE=pd.read_csv("data/model_HIVAE_inputDropout_fma_Missing20_1_z2_y5_s10_batch100_data_reconstruction.csv",header=None)
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
model_UP_FHDI=MLPRegressor(hidden_layer_sizes=(50, 50, 50,25),random_state=7)
model_UP_FHDI.fit(train_UP_FHDI, Y)

model_NAIVE=MLPRegressor(hidden_layer_sizes=(50, 50, 50,25),random_state=7)
model_NAIVE.fit(train_NAIVE, Y)

model_DELETE=MLPRegressor(hidden_layer_sizes=(50, 50, 50,25),random_state=7)
model_DELETE.fit(train_DELETE, Y_delete)

model_GAIN=MLPRegressor(hidden_layer_sizes=(50, 50, 50,25),random_state=7)
model_GAIN.fit(train_GAIN, Y)

model_HI_VAE=MLPRegressor(hidden_layer_sizes=(50, 50, 50,25),random_state=7)
model_HI_VAE.fit(train_HI_VAE, Y)

model_original=MLPRegressor(hidden_layer_sizes=(50, 50, 50,25),random_state=7)
model_original.fit(train_original, Y)

#pred
Y_UP_FHDI = model_UP_FHDI.predict(test)
Y_GAIN = model_GAIN.predict(test)
Y_NAIVE=model_NAIVE.predict(test)
Y_DELETE=model_DELETE.predict(test)
Y_HI_VAE=model_HI_VAE.predict(test)
Y_original = model_original.predict(test)

#RMSE
sum1 = 0.0
for i in range(len(Y_true)):
    a = (Y_DELETE[i] - Y_true[i]) / len(Y_true)
    b = a * a
    sum1 += b

RMSE_delete = math.sqrt(sum1)

sum2 = 0.0
for i in range(len(Y_true)):
    c = (Y_NAIVE[i] - Y_true[i]) / len(Y_true)
    d = c * c
    sum2 += d

RMSE_naive = math.sqrt(sum2)


sum3 = 0.0
for i in range(len(Y_true)):
    e = (Y_UP_FHDI[i] - Y_true[i]) / len(Y_true)
    f = e * e
    sum3 += f

RMSE_UP_FHDI = math.sqrt(sum3)

sum4 = 0.0
for i in range(len(Y_true)):
    g = (Y_GAIN[i] - Y_true[i]) / len(Y_true)
    h = g * g
    sum4 += h

RMSE_GAIN = math.sqrt(sum4)

sum5 = 0.0
for i in range(len(Y_true)):
    j = (Y_HI_VAE[i] - Y_true[i]) / len(Y_true)
    k = j * j
    sum5 += k

RMSE_HI_VAE = math.sqrt(sum5)

sum6 = 0.0
for i in range(len(Y_true)):
    l = (Y_original[i] - Y_true[i]) / len(Y_true)
    m = l * l
    sum6 += m

RMSE_orignal = math.sqrt(sum6)

print("MLP")
print(RMSE_UP_FHDI)
print(RMSE_naive)
print(RMSE_delete)
print(RMSE_GAIN)
print(RMSE_HI_VAE)
print(RMSE_orignal)

#nRMASE
nRMSE_delete = RMSE_delete/RMSE_orignal
nRMSE_naive = RMSE_naive/RMSE_orignal
nRMSE_UP_FHDI = RMSE_UP_FHDI/RMSE_orignal
nRMSE_GAIN = RMSE_GAIN/RMSE_orignal
nRMSE_HI_VAE = RMSE_HI_VAE/RMSE_orignal

print(nRMSE_UP_FHDI)
print(nRMSE_naive)
print(nRMSE_delete)
print(nRMSE_GAIN)
print(nRMSE_HI_VAE)