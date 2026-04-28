# coding=utf-8
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import f1_score
import math

# ======================
# train
# ======================
train_UP_FHDI = pd.read_csv("data/final_daty_UP_FHDI_binary.csv", header=None)
train_GAIN = pd.read_csv("data/imputed.txt", sep=" ", header=None)
train_NAIVE = pd.read_csv("data/NAIVE.csv", header=None)
train_DELETE = pd.read_csv("data/DELETE.csv", header=None)
train_HI_VAE = pd.read_csv(
    "data/model_HIVAE_inputDropout_fma_Missing20_1_z2_y5_s10_batch100_data_reconstruction.csv",
    header=None
)
train_original = pd.read_csv("data/0.7data.csv", header=None)

# ======================
# 每个数据集用自己的标签
# ======================
X_UP_FHDI = train_UP_FHDI.iloc[:, :-1]
Y_UP_FHDI = train_UP_FHDI.iloc[:, -1]

X_GAIN = train_GAIN.iloc[:, :-1]
Y_GAIN_train = train_GAIN.iloc[:, -1]

X_NAIVE = train_NAIVE.iloc[:, :-1]
Y_NAIVE_train = train_NAIVE.iloc[:, -1]

X_DELETE = train_DELETE.iloc[:, :-1]
Y_DELETE_train = train_DELETE.iloc[:, -1]

X_HI_VAE = train_HI_VAE.iloc[:, :-1]
Y_HI_VAE_train = train_HI_VAE.iloc[:, -1]

X_original = train_original.iloc[:, :-1]
Y_original_train = train_original.iloc[:, -1]

# ======================
# test
# ======================
test = pd.read_csv("data/0.3data.csv", header=None)
Y_true = test.iloc[:, -1].values
test = test.iloc[:, :-1]

# ======================
# 训练模型
# ======================
model_UP_FHDI = RandomForestRegressor()
model_UP_FHDI.fit(X_UP_FHDI, Y_UP_FHDI)

model_NAIVE = RandomForestRegressor()
model_NAIVE.fit(X_NAIVE, Y_NAIVE_train)

model_DELETE = RandomForestRegressor()
model_DELETE.fit(X_DELETE, Y_DELETE_train)

model_GAIN = RandomForestRegressor()
model_GAIN.fit(X_GAIN, Y_GAIN_train)

model_HI_VAE = RandomForestRegressor()
model_HI_VAE.fit(X_HI_VAE, Y_HI_VAE_train)

model_original = RandomForestRegressor()
model_original.fit(X_original, Y_original_train)



#pred
Y_UP_FHDI = model_UP_FHDI.predict(test)
Y_GAIN = model_GAIN.predict(test)
Y_NAIVE=model_NAIVE.predict(test)
Y_DELETE=model_DELETE.predict(test)
Y_HI_VAE=model_HI_VAE.predict(test)
Y_original = model_original.predict(test)

# RMSE
sum1 = 0.0
for i in range(len(Y_true)):
    a = (Y_DELETE[i] - Y_true[i])
    sum1 += a * a
RMSE_delete = math.sqrt(sum1 / len(Y_true))


sum2 = 0.0
for i in range(len(Y_true)):
    c = (Y_NAIVE[i] - Y_true[i])
    sum2 += c * c
RMSE_naive = math.sqrt(sum2 / len(Y_true))


sum3 = 0.0
for i in range(len(Y_true)):
    e = (Y_UP_FHDI[i] - Y_true[i])
    sum3 += e * e
RMSE_UP_FHDI = math.sqrt(sum3 / len(Y_true))


sum4 = 0.0
for i in range(len(Y_true)):
    g = (Y_GAIN[i] - Y_true[i])
    sum4 += g * g
RMSE_GAIN = math.sqrt(sum4 / len(Y_true))


sum5 = 0.0
for i in range(len(Y_true)):
    j = (Y_HI_VAE[i] - Y_true[i])
    sum5 += j * j
RMSE_HI_VAE = math.sqrt(sum5 / len(Y_true))


sum6 = 0.0
for i in range(len(Y_true)):
    l = (Y_original[i] - Y_true[i])
    sum6 += l * l
RMSE_orignal = math.sqrt(sum6 / len(Y_true))

print("Randomforest")
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