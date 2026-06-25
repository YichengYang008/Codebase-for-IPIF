# 参数设置
n <- 6000  # 总行数 (修改为6000)
rho <- 0.5  # 相关系数
n_cols <- 16 # 总列数

# 生成误差项（16列数据，每组4列，共4组）
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

# 第三组误差项 (用于y9-y12) - [新增]
e9 = rnorm(n, 0, 1)
e10 = rnorm(n, 0, 1)
e11 = rgamma(n, 1, 1)
e12 = rnorm(n, 0, sd = sqrt(1.5))

# 第四组误差项 (用于y13-y16) - [新增]
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

# 第9-12列（第三组）- [新增]
y9 = 1 + e9
y10 = y9 + 2.0 + rho*e9 + sqrt(0.75)*e10
y11 = y10 + e11
y12 = -1.0 + y9 + 0.25*y11 + e12

# 第13-16列（第四组）- [新增]
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

cat("训练集大小:", n_train, "行\n")
cat("测试集大小:", n_test, "行\n")

# 创建输出目录（可选，防止当前目录混乱）
if(!dir.exists("Noise_Data")) dir.create("Noise_Data")

# 将测试集写入文件
write.csv(test_data, file = "Noise_Data/0.3data.csv", row.names = FALSE)
cat("测试集已保存为 Noise_Data/0.3data.csv\n")

# ========== 第一步：先生成完整的加噪训练集 ==========

# 定义高斯噪声的标准差
noise_sigmas <- c(0.1, 0.5, 1.0)

# 为每个噪声标准差生成完整的加噪训练集
for (sigma in noise_sigmas) {
  cat(paste0("\n========== 生成完整的加噪训练集 σ = ", sigma, " ==========\n"))
  
  # 创建训练数据的完整副本
  train_data_with_noise_full <- train_data
  
  # 生成与训练集形状相同的高斯噪声 (n_train * 16)
  noise <- matrix(rnorm(n_train * n_cols, mean = 0, sd = sigma), 
                  nrow = n_train, ncol = n_cols)
  
  # 将噪声添加到所有值上
  train_data_with_noise_full <- train_data_with_noise_full + noise
  
  # 创建文件名
  sigma_name <- gsub("\\.", "_", as.character(sigma)) 
  full_filename <- paste0("Noise_Data/daty_", sigma_name, "_full.csv")
  
  # 保存完整的加噪训练集
  write.csv(train_data_with_noise_full, file = full_filename, row.names = FALSE)
  
  cat(paste0("已保存完整加噪训练集: ", full_filename, 
             " (", n_train, "行，", n_cols, "列，无缺失，σ = ", sigma, ")\n"))
}

# ========== 第二步：创建缺失指示矩阵 ==========

# 处理训练集：在前40%的行中制造40%的缺失
n_missing_train <- round(n_train * 0.5)  # 训练集前40%的行
train_missing_rows <- 1:n_missing_train

# 创建缺失指示矩阵（初始化为全1）
# 使用矩阵操作比单独定义r1-r16更简洁高效
datr_train <- matrix(1, nrow = n_train, ncol = n_cols)
colnames(datr_train) <- paste0("r", 1:n_cols)

set.seed(456)
for (row in train_missing_rows) {
  # 为每一行随机选择哪些列缺失，每行大约40%的列缺失 (16 * 0.4 = 6.4 -> 约6列)
  missing_cols <- sample(1:n_cols, size = round(n_cols * 0.5))
  
  # 设置缺失指示符为0
  datr_train[row, missing_cols] <- 0
}

# 写入datr文件
write.csv(datr_train, file = "Noise_Data/datr.csv", row.names = FALSE)
cat("\n缺失指示矩阵已保存为 Noise_Data/datr.csv\n")

# ========== 第三步：为每个噪声水平创建有缺失的版本 ==========

# 保存原始训练数据（无噪声，无缺失）
write.csv(train_data, file = "Noise_Data/daty_train_original.csv", row.names = FALSE)
cat("原始训练数据（无缺失，无噪声）已保存为 Noise_Data/daty_train_original.csv\n")

for (sigma in noise_sigmas) {
  cat(paste0("\n========== 为噪声水平 σ = ", sigma, " 创建有缺失的版本 ==========\n"))
  
  # 读取对应的完整加噪训练集
  sigma_name <- gsub("\\.", "_", as.character(sigma))
  full_filename <- paste0("Noise_Data/daty_", sigma_name, "_full.csv")
  train_data_noise_full <- read.csv(full_filename)
  
  # 创建有缺失的版本
  train_data_with_missing_and_noise <- train_data_noise_full
  
  # 根据缺失指示矩阵设置缺失值
  # 使用矩阵索引快速赋值：datr_train == 0 的位置设为 NA
  train_data_with_missing_and_noise[datr_train == 0] <- NA
  
  # 创建文件名
  filename <- paste0("Noise_Data/daty_", sigma_name, ".csv")
  
  # 保存有缺失的加噪训练集
  write.csv(train_data_with_missing_and_noise, file = filename, row.names = FALSE)
  
  cat(paste0("已保存有缺失的加噪训练集: ", filename, 
             " (", n_train, "行，", n_cols, "列，有缺失，σ = ", sigma, ")\n"))
}

# ========== 输出统计信息 ==========

cat("\n========== 数据生成统计 ==========\n")
cat("总样本量:", n, "行\n")
cat("数据维度:", n_cols, "列\n")
cat("训练集大小:", n_train, "行\n")
cat("测试集大小:", n_test, "行\n")
cat("训练集前", n_missing_train, "行制造了缺失值\n")
cat("训练集缺失值总数:", sum(datr_train == 0), "个\n")
cat("训练集总单元格数:", n_train * n_cols, "\n\n")

cat("========== 文件输出 (保存在 Noise_Data/ 文件夹中) ==========\n")
cat("1. 测试集: 0.3data.csv\n")
cat("2. 缺失指示矩阵: datr.csv\n")
cat("3. 原始训练数据: daty_train_original.csv\n")
cat("4. 完整加噪训练集 (daty_0_1_full.csv 等)\n")
cat("5. 有缺失的加噪训练集 (daty_0_1.csv 等)\n\n")

# ========== 数据摘要 (展示前8列) ==========

cat("========== 数据摘要 (前8列) ==========\n")
cat("测试集摘要:\n")
print(summary(test_data[, 1:8]))

# 对于每个噪声水平，显示简要摘要
for (sigma in noise_sigmas) {
  sigma_name <- gsub("\\.", "_", as.character(sigma))
  filename <- paste0("Noise_Data/daty_", sigma_name, ".csv")
  train_data_noise_missing <- read.csv(filename)
  
  cat(paste0("\n有缺失的加噪训练集（σ = ", sigma, "）摘要 (前8列):\n"))
  print(summary(train_data_noise_missing[, 1:8]))
}

# ========== 缺失模式验证 ==========

cat("\n========== 缺失模式验证 ==========\n")
cat("验证所有噪声水平是否有相同的缺失模式:\n")

# 读取所有有缺失的加噪训练集
daty_0_1 <- read.csv("Noise_Data/daty_0_1.csv")
daty_0_5 <- read.csv("Noise_Data/daty_0_5.csv")
daty_1_0 <- read.csv("Noise_Data/daty_1_0.csv") # 注意这里使用了统一的命名格式

# 检查缺失位置是否相同
missing_positions_0_1 <- which(is.na(daty_0_1), arr.ind = TRUE)
missing_positions_0_5 <- which(is.na(daty_0_5), arr.ind = TRUE)
missing_positions_1_0 <- which(is.na(daty_1_0), arr.ind = TRUE)

cat(paste0("daty_0_1 缺失数: ", nrow(missing_positions_0_1), "\n"))
cat(paste0("daty_0_5 缺失数: ", nrow(missing_positions_0_5), "\n"))
cat(paste0("daty_1_0 缺失数: ", nrow(missing_positions_1_0), "\n"))

if (identical(missing_positions_0_1, missing_positions_0_5) && 
    identical(missing_positions_0_1, missing_positions_1_0)) {
  cat("✓ 所有噪声水平的缺失模式完全相同\n")
} else {
  cat("✗ 缺失模式不完全相同，请检查代码\n")
}