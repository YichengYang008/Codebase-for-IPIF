# coding=utf-8
import pandas as pd
import numpy as np


data=pd.read_csv("Dataset/0.7data.csv",header=None)

data=data.values
n,p=data.shape
print(data.shape)
daty_total_column = np.zeros((n * p),dtype=np.float64)
q = 0

# 复制daty到daty_total_column
for j in range(p):
    for i in range(n):
        daty_total_column[q] = data[i, j]
        q += 1

 # 写入二进制文件
with open("data/full_daty_column_binary.bin", "wb") as write_file:
    write_file.write(daty_total_column.tobytes())