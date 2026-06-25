# 参数设置
n <- 6000  # 总行数
rho <- 0.5  # 相关系数

# 生成32个误差项（16列数据，每组4列，共4组，需要16个误差项，但为了严谨匹配分布，我们生成对应的误差）
# 这里沿用原逻辑：每4列一组，需要4个误差项。16列需要4组。
set.seed(1345)

# 第一组误差项 (用于y1-y4)
e1 = rnorm(n, 0, 1)
e2 = rnorm(n, 0, 1)
e3 = rgamma(n, 1, 1)
e4 = rnorm(n, 0, sd = sqrt(1.5))

# 第二组误差项 (用于y5-y8)
e5 = rnorm(n, 0, 1)
e6 = rnorm(n, 0, 1)
e7 = rgamma(n, 1, 1)
e8 = rnorm(n, 0, sd = sqrt(1.5))

# 第三组误差项 (用于y9-y12) - 新增
e9 = rnorm(n, 0, 1)
e10 = rnorm(n, 0, 1)
e11 = rgamma(n, 1, 1)
e12 = rnorm(n, 0, sd = sqrt(1.5))

# 第四组误差项 (用于y13-y16) - 新增
e13 = rnorm(n, 0, 1)
e14 = rnorm(n, 0, 1)
e15 = rgamma(n, 1, 1)
e16 = rnorm(n, 0, sd = sqrt(1.5))

# 生成完整的16列数据（无缺失）
# 第1-4列（第一组）
y1 = 1 + e1
y2 = y1 + 2.0 + rho*e1 + sqrt(0.75)*e2
y3 = y2 + e3
y4 = -1.0 + y1 + 0.25*y3 + e4

# 第5-8列（第二组）
y5 = 1 + e5
y6 = y5 + 2.0 + rho*e5 + sqrt(0.75)*e6
y7 = y6 + e7
y8 = -1.0 + y5 + 0.25*y7 + e8

# 第9-12列（第三组）- 新增
y9 = 1 + e9
y10 = y9 + 2.0 + rho*e9 + sqrt(0.75)*e10
y11 = y10 + e11
y12 = -1.0 + y9 + 0.25*y11 + e12

# 第13-16列（第四组）- 新增
y13 = 1 + e13
y14 = y13 + 2.0 + rho*e13 + sqrt(0.75)*e14
y15 = y14 + e15
y16 = -1.0 + y13 + 0.25*y15 + e16

# 完整的无缺失数据 (16列)
daty_full_all = cbind(y1, y2, y3, y4, y5, y6, y7, y8,
                      y9, y10, y11, y12, y13, y14, y15, y16)

# 将数据分成70%训练集和30%测试集
set.seed(123)
train_indices <- sample(1:n, size = round(n * 0.7))
test_indices <- setdiff(1:n, train_indices)

train_data <- daty_full_all[train_indices, ]
test_data <- daty_full_all[test_indices, ]

n_train <- nrow(train_data)
n_test <- nrow(test_data)
n_cols <- 16 # 列数

cat("训练集大小:", n_train, "行\n")
cat("测试集大小:", n_test, "行\n")

# 将测试集写入test.csv
write.csv(test_data, file = "test.csv", row.names = FALSE)
cat("测试集已保存为 test.csv\n")

# 处理训练集：在前50%的行中制造50%的缺失
n_missing_train <- round(n_train * 0.5)  # 训练集前40%的行
train_missing_rows <- 1:n_missing_train

# 创建缺失指示矩阵（训练集）
# 为了简化16列的处理，使用矩阵方式初始化
datr_train <- matrix(1, nrow = n_train, ncol = n_cols)
colnames(datr_train) <- paste0("r", 1:n_cols)

set.seed(456)
for (row in train_missing_rows) {
  # 为每一行随机选择哪些列缺失，每行大约40%的列缺失 (16 * 0.4 = 6.4 -> 约6列)
  missing_cols <- sample(1:n_cols, size = round(n_cols * 0.5))
  
  # 设置缺失指示符为0
  datr_train[row, missing_cols] <- 0
}

# 根据缺失指示矩阵设置训练数据缺失值
train_data_with_missing <- train_data

# 使用矩阵索引快速赋值NA
# 逻辑：datr_train == 0 的位置，在 train_data_with_missing 中设为 NA
train_data_with_missing[datr_train == 0] <- NA

# ========== 修改部分：在训练集的任何位置随机注入0.5%的异常值 ==========

# 计算训练集总单元格数 (现在是 16列)
total_cells <- n_train * n_cols

# 计算需要注入的异常值数量（0.5%的总单元格数）
n_outliers <- round(total_cells * 0.005)

cat("训练集总单元格数:", total_cells, "\n")
cat("需要注入的异常值数量:", n_outliers, "\n")

# 随机选择异常值的位置（行和列）
set.seed(789)
outlier_rows <- sample(1:n_train, n_outliers, replace = TRUE)
outlier_cols <- sample(1:n_cols, n_outliers, replace = TRUE) # 范围 1-16

# 创建包含异常值的训练数据（有缺失）
train_data_outlier <- train_data_with_missing

# 注入异常值
for(i in 1:n_outliers) {
  row_idx <- outlier_rows[i]
  col_idx <- outlier_cols[i]
  # 注入异常值：为了区分列，依旧使用 20 + col_idx (范围 21-36)
  train_data_outlier[row_idx, col_idx] <- 20 + col_idx 
}

# ========== 创建daty_full_outlier（无缺失，有异常值） ==========

# 创建完整的训练集数据（无缺失）并注入异常值
train_full_data_with_outliers <- train_data  # 先复制训练集的完整数据

# 在相同的位置注入异常值
for(i in 1:n_outliers) {
  row_idx <- outlier_rows[i]
  col_idx <- outlier_cols[i]
  train_full_data_with_outliers[row_idx, col_idx] <- 20 + col_idx  # 注入异常值
}

# 将训练集（包含异常值）和测试集合并为完整的daty_full_outlier
daty_full_outlier <- rbind(train_full_data_with_outliers, test_data)

# ========== 结束修改部分 ==========

# 创建训练集的完整数据（无缺失，无异常值，用于RMSE计算）
train_full_data <- train_data  # 训练集的完整版本

# 写入所有文件
# 确保输出目录存在，如果不存在可能会报错，这里假设目录已存在或直接写入当前目录
# 如果没有 "Outlier" 文件夹，请先创建： dir.create("Outlier", showWarnings = FALSE)
if(!dir.exists("Outlier-Chinese")) dir.create("Outlier-Chinese")

write.csv(as.data.frame(datr_train), file = "Outlier-Chinese/datr_train.csv", row.names = FALSE)
write.csv(train_data_outlier, file = "Outlier-Chinese/daty_train_outlier.csv", row.names = FALSE)
write.csv(train_full_data, file = "Outlier-Chinese/daty_train_full.csv", row.names = FALSE)
write.csv(daty_full_outlier, file = "Outlier-Chinese/daty_full_outlier.csv", row.names = FALSE)

# 输出统计信息
cat("\n========== 数据生成统计 ==========\n")
cat("总样本量:", n, "行\n")
cat("数据维度:", n_cols, "列\n")
cat("训练集大小:", n_train, "行\n")
cat("测试集大小:", n_test, "行\n")
cat("训练集前", n_missing_train, "行制造了缺失值\n")
cat("训练集缺失值总数:", sum(is.na(train_data_with_missing)), "\n")
cat("训练集总单元格数:", total_cells, "\n")
cat("训练集注入的异常值数量:", n_outliers, "（占总单元格的0.5%）\n")
cat("完整数据集中异常值数量:", n_outliers, "\n\n")

cat("========== 文件输出 ==========\n")
cat("1. 测试集: test.csv (", n_test, "行，16列)\n")
cat("2. 训练集缺失指示矩阵: Outlier-Chinese/datr_train.csv (", n_train, "行，16列)\n")
cat("3. 训练集含缺失和异常值数据: Outlier-Chinese/daty_train_outlier.csv\n")
cat("4. 训练集完整数据（无缺失）: Outlier-Chinese/daty_train_full.csv\n")
cat("5. 完整数据集（无缺失，有异常值）: Outlier-Chinese/daty_full_outlier.csv\n\n")

# 检查异常值情况
cat("========== 异常值检查 ==========\n")
# 找出训练数据中所有大于20的值（异常值）
outlier_cells_train <- which(train_data_outlier > 20, arr.ind = TRUE)
cat("训练集中异常值数量（值>20）:", nrow(outlier_cells_train), "\n")

if (nrow(outlier_cells_train) > 0) {
  cat("训练集中异常值分布（值 vs 频数）：\n")
  outlier_table_train <- table(train_data_outlier[outlier_cells_train])
  print(outlier_table_train)
}

# 数据摘要
cat("\n========== 数据摘要 (前8列展示) ==========\n")
cat("测试集摘要:\n")
print(summary(test_data[,1:8]))

cat("\n训练集（含缺失和异常值）摘要 (前8列展示):\n")
print(summary(train_data_outlier[,1:8]))