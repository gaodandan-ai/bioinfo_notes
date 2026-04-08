```R
# ==========================================================
# 批量绘制 1h / 4h / 24h 的联合表达热图（详细注释学习版）
#
# 适配场景：
# - 使用你前面差异分析输出的归一化表达矩阵
# - 文件名采用序号前缀风格，例如：
#   03_1h_normalized_expression.csv
#   03_4h_normalized_expression.csv
#   03_24h_normalized_expression.csv
#
# 热图内容：
# - 把 1h / 4h / 24h 三个时间点的归一化表达矩阵合并
# - 取三个时间点共同存在的基因
# - 对每个基因做行标准化（Z-score）
# - 绘制联合热图
# ==========================================================


# ==========================================================
# 第 1 步：加载所需 R 包
# 说明：
# - 如果某个包不存在，则自动安装
# - ComplexHeatmap 是画高级热图最常用的包之一
# - circlize 用来定义连续型颜色映射
# ==========================================================

if (!require("BiocManager")) install.packages("BiocManager")

required_pkgs <- c(
  "ComplexHeatmap",
  "circlize",
  "RColorBrewer",
  "stringr"
)

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("ComplexHeatmap", "circlize")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg, update = FALSE)
    }
    library(pkg, character.only = TRUE)
  }
}


# ==========================================================
# 第 2 步：设置项目目录
# 说明：
# - 这里尽量和你前面的差异分析脚本保持一致
# - 这样整个项目路径统一，不容易乱
# ==========================================================

project_dir <- "F:/thermal_transcriptome"

# 差异分析结果总目录
result_dir <- file.path(project_dir, "results", "multi_timepoint_analysis")


heatmap_file <- file.path(
  heatmap_dir,
  "04_all_timepoints_expression_heatmap.png"
)

# ==========================================================
# 第 3 步：定义时间点参数
# 说明：
# - 和你前面差异分析脚本尽量保持一致
# - 这样列名、样本顺序、输出子目录都统一
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
# 第 4 步：读取每个时间点的归一化表达矩阵
# 说明：
# - 这里读取的是你前面分析输出的：
#   03_1h_normalized_expression.csv
#   03_4h_normalized_expression.csv
#   03_24h_normalized_expression.csv
# - 每个文件至少应包含：
#   Geneid + 当前时间点对应的6个样本列
# ==========================================================

norm_data_list <- list()   # 用来保存每个时间点的数据
all_genes <- NULL          # 用来逐步取三个时间点的共同基因

for (timepoint in names(timepoint_params)) {
  
  cat("\n------------------------------------------\n")
  cat("正在读取 ", timepoint, " 的归一化表达矩阵...\n", sep = "")
  
  params <- timepoint_params[[timepoint]]
  current_output_dir <- file.path(result_dir, params$output_subdir)
  
  # 这里使用你前面统一过的“03_”前缀文件名
  norm_file <- file.path(
    current_output_dir,
    paste0("03_", timepoint, "_normalized_expression.csv")
  )
  
  # 如果文件不存在，就跳过
  if (!file.exists(norm_file)) {
    warning(paste0("⚠️ 未找到文件：", norm_file, "，跳过 ", timepoint))
    next
  }
  
  # 读取归一化表达矩阵
  norm_data <- read.csv(
    norm_file,
    check.names = FALSE,       # 保留原始列名，不自动改名
    fileEncoding = "UTF-8-BOM",
    na.strings = ""
  )
  
  # 检查 Geneid 列是否存在
  if (!"Geneid" %in% colnames(norm_data)) {
    warning(paste0("⚠️ 文件缺少 Geneid 列：", norm_file, "，跳过 ", timepoint))
    next
  }
  
  # 检查样本列是否存在
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
  
  # 只保留 Geneid + 当前时间点的样本列
  norm_data <- norm_data[, c("Geneid", params$count_cols), drop = FALSE]
  
  # 保存到列表中
  norm_data_list[[timepoint]] <- norm_data
  
  # 逐步求三个时间点共有的基因
  if (is.null(all_genes)) {
    all_genes <- norm_data$Geneid
  } else {
    all_genes <- intersect(all_genes, norm_data$Geneid)
  }
  
  cat("✅ 成功读取：", norm_file, "\n")
  cat("✅ 当前文件基因数：", nrow(norm_data), "\n")
}


# ==========================================================
# 第 5 步：检查是否有可用数据
# ==========================================================

if (length(norm_data_list) == 0) {
  stop("❌ 没有读取到任何可用的归一化表达矩阵，无法绘图！")
}

if (is.null(all_genes) || length(all_genes) == 0) {
  stop("❌ 三个时间点没有共同基因，无法绘图！")
}

cat("\n✅ 三个时间点共同基因数：", length(all_genes), "\n")


# ==========================================================
# 第 6 步：构建合并表达矩阵
# 说明：
# - 行：共同基因
# - 列：三个时间点所有样本
# - 样本顺序按 1h -> 4h -> 24h 排列
#
# 为什么不用 merge 反复合并？
# - 因为这里每个时间点样本列都不同
# - 直接按共同基因顺序 cbind 更清晰，也更不容易出错
# ==========================================================

merged_expr_list <- list()

for (timepoint in c("1h", "4h", "24h")) {
  
  if (!timepoint %in% names(norm_data_list)) next
  
  norm_data <- norm_data_list[[timepoint]]
  
  # 按共同基因的顺序重新排列
  norm_data_sub <- norm_data[match(all_genes, norm_data$Geneid), ]
  
  # 只取表达矩阵部分，不要 Geneid
  expr_part <- norm_data_sub[, timepoint_params[[timepoint]]$count_cols, drop = FALSE]
  
  merged_expr_list[[timepoint]] <- expr_part
}

# 横向拼接所有时间点的样本列
expr_matrix <- do.call(cbind, merged_expr_list)

# 设置行名为基因ID
rownames(expr_matrix) <- all_genes

# 强制转成数值型矩阵（防止读入时混成字符）
expr_matrix <- apply(expr_matrix, 2, as.numeric)
rownames(expr_matrix) <- all_genes

cat("✅ 合并表达矩阵维度：", nrow(expr_matrix), " 个基因 × ", ncol(expr_matrix), " 个样本\n")


# ==========================================================
# 第 7 步：处理缺失值
# 说明：
# - 热图函数通常不喜欢 NA 太多
# - 这里采用“每个基因自己的均值”来填补该基因的 NA
# - 如果某个基因整行都是 NA，就删掉
# ==========================================================

# 删除整行全是 NA 的基因
expr_matrix <- expr_matrix[!apply(is.na(expr_matrix), 1, all), , drop = FALSE]

# 用每一行的均值填补 NA
for (i in 1:nrow(expr_matrix)) {
  row_na <- is.na(expr_matrix[i, ])
  if (any(row_na)) {
    expr_matrix[i, row_na] <- mean(expr_matrix[i, !row_na], na.rm = TRUE)
  }
}

cat("✅ 缺失值处理完成，当前矩阵维度：", nrow(expr_matrix), " × ", ncol(expr_matrix), "\n")


# ==========================================================
# 第 8 步：按行标准化（Z-score）
# 说明：
# - 热图里常见的做法是“每个基因单独标准化”
# - 这样可以突出该基因在不同样本之间的相对高低变化
#
# 原理：
#   Z = (x - mean) / sd
#
# 注意：
# - 画热图时，标准化后的值更适合观察模式
# - 但它不代表原始表达量高低
# ==========================================================

expr_matrix_scaled <- t(scale(t(expr_matrix)))

# 去掉标准化后出现 NA 的行
# （通常是因为某些基因在所有样本中表达完全一样，标准差为0）
expr_matrix_scaled <- expr_matrix_scaled[complete.cases(expr_matrix_scaled), , drop = FALSE]

cat("✅ 行标准化完成，最终热图矩阵维度：", nrow(expr_matrix_scaled), " × ", ncol(expr_matrix_scaled), "\n")


# ==========================================================
# 第 9 步：构建样本注释信息
# 说明：
# - 给每个样本标注：
#   1) Group：Control / Treatment
#   2) Time：1h / 4h / 24h
#
# 例如：
#   C-1-1 -> Control, 1h
#   P-24-3 -> Treatment, 24h
# ==========================================================

sample_names <- colnames(expr_matrix_scaled)

sample_anno <- data.frame(
  Sample = sample_names,
  Group = ifelse(grepl("^C-", sample_names), "Control", "Treatment"),
  Time = ifelse(
    grepl("-1-", sample_names), "1h",
    ifelse(grepl("-4-", sample_names), "4h", "24h")
  ),
  row.names = sample_names,
  stringsAsFactors = FALSE
)

# 设置因子顺序，保证图中显示顺序符合实验设计
sample_anno$Group <- factor(sample_anno$Group, levels = c("Control", "Treatment"))
sample_anno$Time  <- factor(sample_anno$Time,  levels = c("1h", "4h", "24h"))


# ==========================================================
# 第 10 步：设置注释颜色
# 说明：
# - Group 用蓝/红区分
# - Time 用 Set2 调色板区分
# ==========================================================

anno_colors <- list(
  Group = c(
    "Control" = "#2166AC",
    "Treatment" = "#B2182B"
  ),
  Time = c(
    "1h" = "#66C2A5",
    "4h" = "#FC8D62",
    "24h" = "#8DA0CB"
  )
)

# 构建顶部注释条
top_anno <- HeatmapAnnotation(
  df = sample_anno[, c("Group", "Time")],
  col = anno_colors,
  show_annotation_name = TRUE,
  annotation_name_side = "left"
)


# ==========================================================
# 第 11 步：设置热图配色
# 说明：
# - 常用蓝-白-红表示低表达 / 中间 / 高表达
# - 这里用固定断点比用分位数更直观
# - 对 Z-score 热图，常见范围是 -2 到 2
# ==========================================================

heatmap_col <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#053061", "#FFFFFF", "#67001F")
)


# ==========================================================
# 第 12 步：输出热图文件
# 说明：
# - 文件名前也统一加编号
# - 例如：
#   04_all_timepoints_expression_heatmap.png
# ==========================================================

heatmap_file <- file.path(
  result_dir,
  "04_all_timepoints_expression_heatmap.png"
)

# 根据基因数动态调整高度，避免基因太多时图太挤
plot_width <- max(10, ncol(expr_matrix_scaled) * 0.6)
plot_height <- max(8, min(20, nrow(expr_matrix_scaled) * 0.02))

png(
  filename = heatmap_file,
  width = plot_width,
  height = plot_height,
  units = "in",
  res = 300,
  type = "cairo"
)


# ==========================================================
# 第 13 步：绘制热图
# 说明：
# - cluster_rows = TRUE：基因聚类
# - cluster_columns = FALSE：样本不聚类，保持实验顺序
# - show_row_names = FALSE：基因太多时一般不显示行名
# ==========================================================

ht <- Heatmap(
  expr_matrix_scaled,
  name = "Z-score",
  col = heatmap_col,
  
  # 行列名称显示设置
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  
  # 聚类设置
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  # 分列显示：按时间点切分，更容易看
  column_split = sample_anno$Time,
  
  # 顶部注释
  top_annotation = top_anno,
  
  # 标题
  column_title = "Combined Heatmap of 1h / 4h / 24h Samples",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # 图例设置
  heatmap_legend_param = list(
    title = "Expression\n(Z-score)",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()


# ==========================================================
# 第 14 步：输出完成信息
# ==========================================================

cat("\n🎉 热图绘制完成！\n")
cat("✅ 输出文件：", heatmap_file, "\n")
cat("✅ 热图矩阵维度：", nrow(expr_matrix_scaled), " 个基因 × ", ncol(expr_matrix_scaled), " 个样本\n")
```