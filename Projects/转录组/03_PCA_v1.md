## PCA 逻辑

1. 读三个时间点归一化矩阵
2. 取共同基因
3. 合并成大矩阵
4. **转置成“样本 × 基因”**
5. 做 PCA
6. 画图

# 1-联合PCA

```R
# ==========================================================
# 联合PCA分析：1h / 4h / 24h 所有样本一起做 PCA
# 适配当前项目结构：
# F:/thermal_transcriptome/results/multi_timepoint_analysis
#
# 输入文件：
#   output_1h/03_1h_normalized_expression.csv
#   output_4h/03_4h_normalized_expression.csv
#   output_24h/03_24h_normalized_expression.csv
#
# 输出文件：
#   05_all_timepoints_PCA_plot.png
#   05_all_timepoints_PCA_scores.csv
# ==========================================================


# ==========================================================
# 第 1 步：加载所需 R 包
# 说明：
# - ggplot2：画 PCA 图
# - stringr：提取样本名中的时间和重复信息
# ==========================================================

if (!require("BiocManager")) install.packages("BiocManager")

required_pkgs <- c("ggplot2", "stringr")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, update = FALSE)
    library(pkg, character.only = TRUE)
  }
}


# ==========================================================
# 第 2 步：设置项目目录
# ==========================================================

project_dir <- "F:/thermal_transcriptome"
result_dir  <- file.path(project_dir, "results", "multi_timepoint_analysis")


# ==========================================================
# 第 3 步：定义时间点参数
# 说明：
# - 和差异分析、热图脚本保持一致
# ==========================================================

timepoint_params <- list(
  "1h" = list(
    count_cols = c("C-1-1", "C-1-2", "C-1-3", "P-1-1", "P-1-2", "P-1-3"),
    output_subdir = "output_1h"
  ),
  "4h" = list(
    count_cols = c("C-4-1", "C-4-2", "C-4-3", "P-4-1", "P-4-2", "P-4-3"),
    output_subdir = "output_4h"
  ),
  "24h" = list(
    count_cols = c("C-24-1", "C-24-2", "C-24-3", "P-24-1", "P-24-2", "P-24-3"),
    output_subdir = "output_24h"
  )
)


# ==========================================================
# 第 4 步：读取三个时间点的归一化表达矩阵
# 说明：
# - 这里读取的是你已经统一命名后的文件：
#   03_1h_normalized_expression.csv
#   03_4h_normalized_expression.csv
#   03_24h_normalized_expression.csv
# - 只保留 Geneid 和当前时间点样本列
# - 最后取三个时间点共有基因
# ==========================================================

norm_data_list <- list()
all_genes <- NULL

for (timepoint in names(timepoint_params)) {
  
  cat("\n------------------------------------------\n")
  cat("正在读取 ", timepoint, " 的归一化表达矩阵...\n", sep = "")
  
  params <- timepoint_params[[timepoint]]
  current_output_dir <- file.path(result_dir, params$output_subdir)
  
  norm_file <- file.path(
    current_output_dir,
    paste0("03_", timepoint, "_normalized_expression.csv")
  )
  
  if (!file.exists(norm_file)) {
    warning(paste0("⚠️ 未找到文件：", norm_file, "，跳过 ", timepoint))
    next
  }
  
  norm_data <- read.csv(
    norm_file,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM",
    na.strings = ""
  )
  
  if (!"Geneid" %in% colnames(norm_data)) {
    warning(paste0("⚠️ 文件缺少 Geneid 列：", norm_file, "，跳过 ", timepoint))
    next
  }
  
  missing_cols <- setdiff(params$count_cols, colnames(norm_data))
  if (length(missing_cols) > 0) {
    warning(
      paste0(
        "⚠️ ", timepoint, " 缺少样本列：",
        paste(missing_cols, collapse = ", "),
        "，跳过该时间点"
      )
    )
    next
  }
  
  norm_data <- norm_data[, c("Geneid", params$count_cols), drop = FALSE]
  norm_data_list[[timepoint]] <- norm_data
  
  # 逐步取共有基因
  if (is.null(all_genes)) {
    all_genes <- norm_data$Geneid
  } else {
    all_genes <- intersect(all_genes, norm_data$Geneid)
  }
  
  cat("✅ 成功读取：", norm_file, "\n")
  cat("✅ 当前文件基因数：", nrow(norm_data), "\n")
}


# ==========================================================
# 第 5 步：检查数据是否可用
# ==========================================================

if (length(norm_data_list) == 0) {
  stop("❌ 没有读取到任何可用的归一化表达矩阵，无法做 PCA！")
}

if (is.null(all_genes) || length(all_genes) == 0) {
  stop("❌ 三个时间点没有共同基因，无法做 PCA！")
}

cat("\n✅ 三个时间点共同基因数：", length(all_genes), "\n")


# ==========================================================
# 第 6 步：构建联合表达矩阵
# 说明：
# - 行：共同基因
# - 列：所有样本
# - 顺序：1h -> 4h -> 24h
# ==========================================================

merged_expr_list <- list()

for (timepoint in c("1h", "4h", "24h")) {
  
  if (!timepoint %in% names(norm_data_list)) next
  
  norm_data <- norm_data_list[[timepoint]]
  
  # 按共有基因顺序重新排列
  norm_data_sub <- norm_data[match(all_genes, norm_data$Geneid), ]
  
  # 只取样本列
  expr_part <- norm_data_sub[, timepoint_params[[timepoint]]$count_cols, drop = FALSE]
  
  merged_expr_list[[timepoint]] <- expr_part
}

expr_matrix <- do.call(cbind, merged_expr_list)
rownames(expr_matrix) <- all_genes

# 转成纯数值矩阵
expr_matrix <- apply(expr_matrix, 2, as.numeric)
rownames(expr_matrix) <- all_genes

cat("✅ 联合表达矩阵维度：", nrow(expr_matrix), " 个基因 × ", ncol(expr_matrix), " 个样本\n")


# ==========================================================
# 第 7 步：处理缺失值
# 说明：
# - 删除整行全是 NA 的基因
# - 其他 NA 用该基因在其他样本中的均值补齐
# ==========================================================

expr_matrix <- expr_matrix[!apply(is.na(expr_matrix), 1, all), , drop = FALSE]

for (i in 1:nrow(expr_matrix)) {
  row_na <- is.na(expr_matrix[i, ])
  if (any(row_na)) {
    expr_matrix[i, row_na] <- mean(expr_matrix[i, !row_na], na.rm = TRUE)
  }
}

cat("✅ 缺失值处理完成，当前矩阵维度：", nrow(expr_matrix), " × ", ncol(expr_matrix), "\n")


# ==========================================================
# 第 8 步：准备 PCA 输入矩阵
# 说明：
# - 原始表达矩阵：行=基因，列=样本
# - 但 prcomp 需要：行=样本，列=变量
# - 所以必须转置
# ==========================================================

expr_matrix_t <- t(expr_matrix)

# 执行 PCA
# center = TRUE：先减均值
# scale. = TRUE：再按标准差缩放
# 这是 PCA 常见做法，可以避免某些基因方差过大主导结果
pca_result <- prcomp(expr_matrix_t, center = TRUE, scale. = TRUE)


# ==========================================================
# 第 9 步：提取 PCA 结果
# 说明：
# - pca_result$x：每个样本在主成分空间中的坐标
# - PC1 / PC2 最常用于画图
# ==========================================================

pca_scores <- as.data.frame(pca_result$x[, 1:2])

# 计算解释率
pca_variance <- pca_result$sdev^2
pca_var_percent <- round(pca_variance / sum(pca_variance) * 100, 1)


# ==========================================================
# 第 10 步：构建样本注释表
# 说明：
# - Group：Control / Treatment
# - Time：1h / 4h / 24h
# - Replicate：1 / 2 / 3
# ==========================================================

sample_names <- rownames(pca_scores)

sample_anno <- data.frame(
  Sample = sample_names,
  Group = ifelse(grepl("^C-", sample_names), "Control", "Treatment"),
  Time = ifelse(
    grepl("-1-", sample_names), "1h",
    ifelse(grepl("-4-", sample_names), "4h", "24h")
  ),
  Replicate = stringr::str_extract(sample_names, "\\d+$"),
  stringsAsFactors = FALSE
)

sample_anno$Group <- factor(sample_anno$Group, levels = c("Control", "Treatment"))
sample_anno$Time <- factor(sample_anno$Time, levels = c("1h", "4h", "24h"))
sample_anno$Replicate <- factor(sample_anno$Replicate, levels = c("1", "2", "3"))

# 合并 PCA 得分和样本注释
pca_df <- cbind(pca_scores, sample_anno)


# ==========================================================
# 第 11 步：绘制 PCA 图
# 说明：
# - 颜色区分 Group
# - 点形区分 Time
# - 这样图上可以同时看“处理效应”和“时间效应”
# ==========================================================

color_palette <- c(
  "Control" = "#2166AC",
  "Treatment" = "#B2182B"
)

shape_palette <- c(
  "1h" = 16,
  "4h" = 17,
  "24h" = 18
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Group, shape = Time), size = 4, alpha = 0.85) +
  xlab(paste0("PC1 (", pca_var_percent[1], "%)")) +
  ylab(paste0("PC2 (", pca_var_percent[2], "%)")) +
  ggtitle("PCA of All Samples (1h / 4h / 24h)") +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5)
  )


# ==========================================================
# 第 12 步：保存 PCA 图和得分表
# 说明：
# - 放在总结果目录里
# - 文件名加编号，方便排序
# ==========================================================

pca_plot_file <- file.path(result_dir, "05_all_timepoints_PCA_plot.png")
pca_score_file <- file.path(result_dir, "05_all_timepoints_PCA_scores.csv")

ggsave(
  filename = pca_plot_file,
  plot = pca_plot,
  width = 10,
  height = 8,
  dpi = 300,
  device = "png"
)

write.csv(pca_df, pca_score_file, row.names = FALSE, na = "")


# ==========================================================
# 第 13 步：输出完成信息
# ==========================================================

cat("\n🎉 联合 PCA 绘制完成！\n")
cat("✅ PCA 图：", pca_plot_file, "\n")
cat("✅ PCA 得分表：", pca_score_file, "\n")
cat("✅ PC1解释率：", pca_var_percent[1], "%\n")
cat("✅ PC2解释率：", pca_var_percent[2], "%\n")
cat("✅ 参与 PCA 的样本数：", nrow(pca_df), "\n")
```


# 2-分开PCA

```R
# ==========================================================
# 各时间点独立 PCA
#
# 输入文件：
#   output_1h/03_1h_normalized_expression.csv
#   output_4h/03_4h_normalized_expression.csv
#   output_24h/03_24h_normalized_expression.csv
#
# 输出文件：
#   output_1h/06_1h_PCA_plot.png
#   output_1h/06_1h_PCA_scores.csv
#   output_4h/06_4h_PCA_plot.png
#   output_4h/06_4h_PCA_scores.csv
#   output_24h/06_24h_PCA_plot.png
#   output_24h/06_24h_PCA_scores.csv
# ==========================================================


if (!require("BiocManager")) install.packages("BiocManager")

required_pkgs <- c("ggplot2", "stringr")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, update = FALSE)
    library(pkg, character.only = TRUE)
  }
}

project_dir <- "F:/thermal_transcriptome"
result_dir  <- file.path(project_dir, "results", "multi_timepoint_analysis")

timepoint_params <- list(
  "1h" = list(
    count_cols = c("C-1-1", "C-1-2", "C-1-3", "P-1-1", "P-1-2", "P-1-3"),
    output_subdir = "output_1h"
  ),
  "4h" = list(
    count_cols = c("C-4-1", "C-4-2", "C-4-3", "P-4-1", "P-4-2", "P-4-3"),
    output_subdir = "output_4h"
  ),
  "24h" = list(
    count_cols = c("C-24-1", "C-24-2", "C-24-3", "P-24-1", "P-24-2", "P-24-3"),
    output_subdir = "output_24h"
  )
)

color_palette <- c(
  "Control" = "#2166AC",
  "Treatment" = "#B2182B"
)

shape_palette <- c(
  "1" = 16,
  "2" = 17,
  "3" = 18
)

cat("\n=== 开始绘制各时间点独立 PCA 图 ===\n")

for (timepoint in names(timepoint_params)) {
  
  cat("\n------------------------------------------\n")
  cat("正在处理 ", timepoint, " ...\n", sep = "")
  
  params <- timepoint_params[[timepoint]]
  current_output_dir <- file.path(result_dir, params$output_subdir)
  
  norm_file <- file.path(
    current_output_dir,
    paste0("03_", timepoint, "_normalized_expression.csv")
  )
  
  if (!file.exists(norm_file)) {
    warning(paste0("⚠️ 未找到文件：", norm_file, "，跳过 ", timepoint))
    next
  }
  
  norm_data <- read.csv(
    norm_file,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM",
    na.strings = ""
  )
  
  if (!"Geneid" %in% colnames(norm_data)) {
    warning(paste0("⚠️ 文件缺少 Geneid 列：", norm_file, "，跳过 ", timepoint))
    next
  }
  
  missing_cols <- setdiff(params$count_cols, colnames(norm_data))
  if (length(missing_cols) > 0) {
    warning(
      paste0(
        "⚠️ ", timepoint, " 缺少样本列：",
        paste(missing_cols, collapse = ", "),
        "，跳过该时间点"
      )
    )
    next
  }
  
  # 只保留当前时间点的 Geneid 和样本列
  norm_data <- norm_data[, c("Geneid", params$count_cols), drop = FALSE]
  
  # 转成表达矩阵
  expr_matrix <- as.matrix(norm_data[, -1])
  rownames(expr_matrix) <- norm_data$Geneid
  expr_matrix <- apply(expr_matrix, 2, as.numeric)
  rownames(expr_matrix) <- norm_data$Geneid
  
  # 删除全 NA 基因
  expr_matrix <- expr_matrix[!apply(is.na(expr_matrix), 1, all), , drop = FALSE]
  
  # 按行均值填补 NA
  for (i in 1:nrow(expr_matrix)) {
    row_na <- is.na(expr_matrix[i, ])
    if (any(row_na)) {
      expr_matrix[i, row_na] <- mean(expr_matrix[i, !row_na], na.rm = TRUE)
    }
  }
  
  # PCA 输入需要转置：行=样本，列=基因
  expr_matrix_t <- t(expr_matrix)
  
  # 做 PCA
  pca_result <- prcomp(expr_matrix_t, center = TRUE, scale. = TRUE)
  
  pca_scores <- as.data.frame(pca_result$x[, 1:2])
  
  pca_variance <- pca_result$sdev^2
  pca_var_percent <- round(pca_variance / sum(pca_variance) * 100, 1)
  
  # 样本注释
  sample_names <- rownames(pca_scores)
  
  sample_anno <- data.frame(
    Sample = sample_names,
    Group = ifelse(grepl("^C-", sample_names), "Control", "Treatment"),
    Replicate = stringr::str_extract(sample_names, "\\d+$"),
    stringsAsFactors = FALSE
  )
  
  sample_anno$Group <- factor(sample_anno$Group, levels = c("Control", "Treatment"))
  sample_anno$Replicate <- factor(sample_anno$Replicate, levels = c("1", "2", "3"))
  
  pca_df <- cbind(pca_scores, sample_anno)
  
  # 画图
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Group, shape = Replicate), size = 4, alpha = 0.85) +
    xlab(paste0("PC1 (", pca_var_percent[1], "%)")) +
    ylab(paste0("PC2 (", pca_var_percent[2], "%)")) +
    ggtitle(paste0("PCA: ", timepoint, " (Control vs Treatment)")) +
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = shape_palette) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5)
    )
  
  # 保存
  pca_plot_file <- file.path(current_output_dir, paste0("06_", timepoint, "_PCA_plot.png"))
  pca_score_file <- file.path(current_output_dir, paste0("06_", timepoint, "_PCA_scores.csv"))
  
  ggsave(
    filename = pca_plot_file,
    plot = pca_plot,
    width = 9,
    height = 7,
    dpi = 300,
    device = "png"
  )
  
  write.csv(pca_df, pca_score_file, row.names = FALSE, na = "")
  
  cat("✅ PCA 图：", pca_plot_file, "\n")
  cat("✅ PCA 得分表：", pca_score_file, "\n")
  cat("✅ PC1解释率：", pca_var_percent[1], "% | PC2解释率：", pca_var_percent[2], "%\n")
}

cat("\n🎉 所有时间点独立 PCA 完成！\n")
```