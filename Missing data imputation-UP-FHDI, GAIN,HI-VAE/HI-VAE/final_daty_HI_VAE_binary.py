# coding=utf-8
import numpy as np
import pandas as pd

data=pd.read_csv("Results/model_HIVAE_inputDropout_swarm_Missing20_1_z2_y5_s10_batch100/model_HIVAE_inputDropout_swarm_Missing20_1_z2_y5_s10_batch100_data_reconstruction.csv",header=None)
data=data.values
n,p=data.shape
print(data.shape)
daty_total_row = np.zeros((n * p),dtype=np.float64)

q = 0

# 复制daty到daty_total_column
for j in range(p):
    for i in range(n):
        daty_total_row[q] = data[i, j]
        q += 1

 # 写入二进制文件
with open("Results/final_daty_HI_VAE_binary.bin", "wb") as write_file:
    write_file.write(daty_total_row.tobytes())