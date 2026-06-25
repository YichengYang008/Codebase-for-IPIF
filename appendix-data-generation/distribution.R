# 参数设置
n <- 5000  # 总行数
rho <- 0.5  # 相关系数

# 定义三个级别的分布参数（越来越不均匀）
levels <- list(
  level1 = list(
    name = "relatively_even",
    e_shape = c(2, 2, 2, 2, 2, 2, 2, 2),  # 形状参数较大，分布相对均匀
    e_scale = c(1, 1, 1, 1, 1, 1, 1, 1)   # 尺度参数
  ),
  level2 = list(
    name = "moderately_skewed", 
    e_shape = c(1, 1, 0.8, 1, 1, 0.8, 0.8, 1),  # 形状参数减小，开始偏态
    e_scale = c(1, 1, 1.2, 1, 1, 1.2, 1.2, 1)   # 尺度参数适当调整
  ),
  level3 = list(
    name = "highly_skewed",
    e_shape = c(0.5, 0.5, 0.3, 0.5, 0.5, 0.3, 0.3, 0.5),  # 形状参数很小，高度偏态
    e_scale = c(1.5, 1.5, 2, 1.5, 1.5, 2, 2, 1.5)          # 尺度参数增大以保持方差
  )
)

# 为每个级别生成数据
for (level_idx in 1:length(levels)) {
  level <- levels[[level_idx]]
  level_name <- level$name
  
  cat("\n========== 生成数据级别 ", level_idx, ": ", level_name, " ==========\n", sep = "")
  
  # 设置随机种子确保可重复性（每个级别不同）
  set.seed(123 + level_idx * 100)
  
  # 生成Gamma分布的误差项（全部使用Gamma分布）
  e1 = rgamma(n, shape = level$e_shape[1], scale = level$e_scale[1])
  e2 = rgamma(n, shape = level$e_shape[2], scale = level$e_scale[2])
  e3 = rgamma(n, shape = level$e_shape[3], scale = level$e_scale[3])
  e4 = rgamma(n, shape = level$e_shape[4], scale = level$e_scale[4])
  e5 = rgamma(n, shape = level$e_shape[5], scale = level$e_scale[5])
  e6 = rgamma(n, shape = level$e_shape[6], scale = level$e_scale[6])
  e7 = rgamma(n, shape = level$e_shape[7], scale = level$e_scale[7])
  e8 = rgamma(n, shape = level$e_shape[8], scale = level$e_scale[8])
  
  # 标准化：使均值为0，方差为1（可选，如果需要）
  # e1 = (e1 - mean(e1)) / sd(e1)
  # 类似地处理e2-e8...
  
  # 生成完整的8列数据
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
  
  # 完整的无缺失数据
  daty_full_all = cbind(y1, y2, y3, y4, y5, y6, y7, y8)
  
  # 将数据分成70%训练集和30%测试集
  set.seed(456)  # 固定分割种子，确保所有级别的分割相同
  train_indices <- sample(1:n, size = round(n * 0.7))
  test_indices <- setdiff(1:n, train_indices)
  
  train_data <- daty_full_all[train_indices, ]
  test_data <- daty_full_all[test_indices, ]
  
  n_train <- nrow(train_data)
  n_test <- nrow(test_data)
  
  cat("训练集大小:", n_train, "行\n")
  cat("测试集大小:", n_test, "行\n")
  
  # 将测试集写入test.csv（每个级别不同的文件）
  test_filename <- paste0("test_", level_name, ".csv")
  write.csv(test_data, file = test_filename, row.names = FALSE)
  cat("测试集已保存为", test_filename, "\n")
  
  # 处理训练集：在前40%的行中制造40%的缺失
  n_missing_train <- round(n_train * 0.4)  # 训练集前40%的行
  train_missing_rows <- 1:n_missing_train
  
  # 创建缺失指示矩阵（训练集）
  r1_train = rep(1, n_train)
  r2_train = rep(1, n_train)
  r3_train = rep(1, n_train)
  r4_train = rep(1, n_train)
  r5_train = rep(1, n_train)
  r6_train = rep(1, n_train)
  r7_train = rep(1, n_train)
  r8_train = rep(1, n_train)
  
  set.seed(789)  # 固定缺失模式种子，确保所有级别的缺失模式相同
  for (row in train_missing_rows) {
    # 为每一行随机选择哪些列缺失，每行大约40%的列缺失
    missing_cols <- sample(1:8, size = round(8 * 0.4))
    
    # 设置缺失指示符
    if (1 %in% missing_cols) r1_train[row] <- 0
    if (2 %in% missing_cols) r2_train[row] <- 0
    if (3 %in% missing_cols) r3_train[row] <- 0
    if (4 %in% missing_cols) r4_train[row] <- 0
    if (5 %in% missing_cols) r5_train[row] <- 0
    if (6 %in% missing_cols) r6_train[row] <- 0
    if (7 %in% missing_cols) r7_train[row] <- 0
    if (8 %in% missing_cols) r8_train[row] <- 0
  }
  
  # 训练集的缺失指示矩阵
  datr_train = cbind(r1_train, r2_train, r3_train, r4_train, 
                     r5_train, r6_train, r7_train, r8_train)
  
  # 根据缺失指示矩阵设置训练数据缺失值
  train_data_with_missing <- train_data
  
  train_data_with_missing[r1_train == 0, 1] = NA
  train_data_with_missing[r2_train == 0, 2] = NA
  train_data_with_missing[r3_train == 0, 3] = NA
  train_data_with_missing[r4_train == 0, 4] = NA
  train_data_with_missing[r5_train == 0, 5] = NA
  train_data_with_missing[r6_train == 0, 6] = NA
  train_data_with_missing[r7_train == 0, 7] = NA
  train_data_with_missing[r8_train == 0, 8] = NA
  
  # 保存文件（每个级别不同的文件）
  datr_filename <- paste0("datr_", level_name, ".csv")
  train_full_filename <- paste0("daty_train_full_", level_name, ".csv")
  train_missing_filename <- paste0("daty_train_", level_name, ".csv")
  
  write.csv(datr_train, file = datr_filename, row.names = FALSE)
  write.csv(train_data, file = train_full_filename, row.names = FALSE)
  write.csv(train_data_with_missing, file = train_missing_filename, row.names = FALSE)
  
  cat("已保存文件:\n")
  cat("  - ", datr_filename, " (缺失指示矩阵)\n", sep = "")
  cat("  - ", train_full_filename, " (训练集完整数据)\n", sep = "")
  cat("  - ", train_missing_filename, " (训练集有缺失数据)\n", sep = "")
  
  # 计算并显示分布的偏度
  cat("\n分布偏度统计:\n")
  for (i in 1:8) {
    col_data <- train_data[, i]
    skewness <- sum((col_data - mean(col_data))^3) / (length(col_data) * sd(col_data)^3)
    cat(paste0("  y", i, " 偏度: ", round(skewness, 4), "\n"))
  }
  
  # 显示分布摘要
  cat("\n训练集完整数据摘要:\n")
  print(summary(train_data))
  
  cat("\n训练集有缺失数据摘要:\n")
  print(summary(train_data_with_missing))
}

# 汇总信息
cat("\n========== 数据生成汇总 ==========\n")
cat("已生成三个不同偏态级别的数据集:\n")
for (level_idx in 1:length(levels)) {
  level <- levels[[level_idx]]
  level_name <- level$name
  
  cat("\n级别 ", level_idx, ": ", level_name, "\n", sep = "")
  cat("  Gamma分布形状参数: ", paste(round(level$e_shape, 2), collapse = ", "), "\n", sep = "")
  cat("  Gamma分布尺度参数: ", paste(round(level$e_scale, 2), collapse = ", "), "\n", sep = "")
  
  # 计算Gamma分布的偏度：偏度 = 2/√shape
  theoretical_skewness <- 2 / sqrt(level$e_shape)
  cat("  理论偏度范围: ", round(min(theoretical_skewness), 2), " 到 ", 
      round(max(theoretical_skewness), 2), "\n", sep = "")
  
  cat("  生成的文件:\n")
  cat("    - test_", level_name, ".csv\n", sep = "")
  cat("    - datr_", level_name, ".csv\n", sep = "")
  cat("    - daty_train_full_", level_name, ".csv\n", sep = "")
  cat("    - daty_train_", level_name, ".csv\n", sep = "")
}

# Gamma分布参数说明
cat("\n========== Gamma分布参数说明 ==========\n")
cat("Gamma分布的概率密度函数：f(x) = x^(k-1) * e^(-x/θ) / (θ^k * Γ(k))\n")
cat("其中：k = 形状参数，θ = 尺度参数\n")
cat("Gamma分布的特性:\n")
cat("1. 均值 = k * θ\n")
cat("2. 方差 = k * θ^2\n")
cat("3. 偏度 = 2 / √k\n")
cat("\n形状参数(k)对分布的影响:\n")
cat("- k > 1: 分布相对对称，接近正态分布\n")
cat("- k = 1: 指数分布，右偏\n")
cat("- k < 1: 高度右偏，大部分质量集中在0附近，有长尾\n")
cat("\n三个级别的设计:\n")
cat("1. 级别1 (relatively_even): k=2，分布相对均匀\n")
cat("2. 级别2 (moderately_skewed): k=0.8-1，开始出现明显偏态\n")
cat("3. 级别3 (highly_skewed): k=0.3-0.5，高度偏态，长尾分布\n")