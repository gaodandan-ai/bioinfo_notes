联合版
```R
# ==========================================================
# 多时间点转录组 PCA 分析（教学完整版）
#
# 适配项目：
# F:/thermal_transcriptome
#
# 输入：
#   results/multi_timepoint_analysis/output_xx/03_*_normalized_expression.csv
#
# 输出：
#   05_all_timepoints_PCA_plot.png
#   05_all_timepoints_PCA_scores.csv
#   06_xh_PCA_plot.png
#   06_xh_PCA_scores.csv
#
# ==========================================================


# ==========================================================
# 第 1 步：加载必要包
# ==========================================================

if (!require("BiocManager")) install.packages("BiocManager")

required_pkgs <- c("ggplot2", "stringr")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


# ==========================================================
# 第 2 步：设置路径（和你前面完全统一）
# ==========================================================

project_dir <- "F:/thermal_transcriptome"
result_dir  <- file.path(project_dir, "results", "multi_timepoint_analysis")


# ==========================================================
# 第 3 步：定义时间点参数
# 说明：
# - count_cols：该时间点的样本列名
# - output_subdir：该时间点结果目录
# ==========================================================

timepoint_params <- list(
  "1h" = list(
    count_cols = c("C-1-1","C-1-2","C-1-3","P-1-1","P-1-2","P-1-3"),
    output_subdir = "output_1h"
  ),
  "4h" = list(
    count_cols = c("C-4-1","C-4-2","C-4-3","P-4-1","P-4-2","P-4-3"),
    output_subdir = "output_4h"
  ),
  "24h" = list(
    count_cols = c("C-24-1","C-24-2","C-24-3","P-24-1","P-24-2","P-24-3"),
    output_subdir = "output_24h"
  )
)


# ==========================================================
# 第 4 步：读取归一化表达矩阵（核心步骤）
#
# 思想：
# - 每个时间点一张表达矩阵
# - 最后要合并 → 必须保证基因一致
# - 所以要取“共同基因交集”
# ==========================================================

norm_data_list <- list()
all_genes <- NULL

for (timepoint in names(timepoint_params)) {
  
  cat("\n=== 读取 ", timepoint, " 数据 ===\n", sep = "")
  
  params <- timepoint_params[[timepoint]]
  current_dir <- file.path(result_dir, params$output_subdir)
  
  file_path <- file.path(
    current_dir,
    paste0("03_", timepoint, "_normalized_expression.csv")
  )
  
  if (!file.exists(file_path)) {
    warning("文件不存在：", file_path)
    next
  }
  
  # 读取数据
  df <- read.csv(file_path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  
  # 只保留 Geneid + 样本列
  df <- df[, c("Geneid", params$count_cols), drop = FALSE]
  
  norm_data_list[[timepoint]] <- df
  
  # 逐步取交集
  if (is.null(all_genes)) {
    all_genes <- df$Geneid
  } else {
    all_genes <- intersect(all_genes, df$Geneid)
  }
}

cat("\n共同基因数：", length(all_genes), "\n")


# ==========================================================
# 第 5 步：合并表达矩阵
#
# 为什么不用 merge？
# - merge 容易打乱顺序
# - 这里直接 match + cbind 更安全
# ==========================================================

merged_list <- list()

for (tp in c("1h","4h","24h")) {
  
  if (!tp %in% names(norm_data_list)) next
  
  df <- norm_data_list[[tp]]
  
  df <- df[match(all_genes, df$Geneid), ]
  
  merged_list[[tp]] <- df[, timepoint_params[[tp]]$count_cols]
}

expr_matrix <- do.call(cbind, merged_list)
rownames(expr_matrix) <- all_genes

# 转数值
expr_matrix <- apply(expr_matrix, 2, as.numeric)
rownames(expr_matrix) <- all_genes


# ==========================================================
# 第 6 步：缺失值处理
#
# 原则：
# - 删除全 NA 行
# - 剩下的 NA 用“该基因均值”补
# ==========================================================

expr_matrix <- expr_matrix[!apply(is.na(expr_matrix), 1, all), ]

for (i in 1:nrow(expr_matrix)) {
  idx <- is.na(expr_matrix[i, ])
  if (any(idx)) {
    expr_matrix[i, idx] <- mean(expr_matrix[i, !idx])
  }
}


# ==========================================================
# 第 7 步：准备 PCA 输入
#
# ❗关键点：
# PCA 要求：
#   行 = 样本
#   列 = 变量（基因）
# ==========================================================

expr_matrix_t <- t(expr_matrix)


# ==========================================================
# 第 8 步：执行 PCA
#
# center = TRUE：
#   每个基因减去均值
#
# scale. = TRUE：
#   每个基因除以标准差
#
# 👉 避免高表达基因主导结果
# ==========================================================

pca_res <- prcomp(expr_matrix_t, center = TRUE, scale. = TRUE)


# ==========================================================
# 第 9 步：提取 PCA 结果
# ==========================================================

pca_scores <- as.data.frame(pca_res$x[, 1:2])

# 解释率
var <- pca_res$sdev^2
var_percent <- round(var / sum(var) * 100, 1)


# ==========================================================
# 第 10 步：构建样本注释
#
# 从样本名解析：
# C-1-1 → Control / 1h / replicate1
# ==========================================================

samples <- rownames(pca_scores)

anno <- data.frame(
  Sample = samples,
  Group = ifelse(grepl("^C-", samples), "Control", "Treatment"),
  Time = ifelse(grepl("-1-", samples), "1h",
                ifelse(grepl("-4-", samples), "4h", "24h")),
  Rep = stringr::str_extract(samples, "\\d+$")
)

anno$Group <- factor(anno$Group, levels = c("Control","Treatment"))
anno$Time  <- factor(anno$Time, levels = c("1h","4h","24h"))

pca_df <- cbind(pca_scores, anno)


# ==========================================================
# 第 11 步：绘图（联合 PCA）
# ==========================================================

p <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(aes(color = Group, shape = Time), size = 4) +
  xlab(paste0("PC1 (", var_percent[1], "%)")) +
  ylab(paste0("PC2 (", var_percent[2], "%)")) +
  ggtitle("PCA (All Timepoints)") +
  scale_color_manual(values = c("Control"="#2166AC","Treatment"="#B2182B")) +
  theme_bw()


# ==========================================================
# 第 12 步：保存结果
# ==========================================================

ggsave(
  file.path(result_dir, "05_all_timepoints_PCA_plot.png"),
  p, width = 10, height = 8, dpi = 300
)

write.csv(
  pca_df,
  file.path(result_dir, "05_all_timepoints_PCA_scores.csv"),
  row.names = FALSE
)


cat("\n🎉 PCA 完成！\n")
```