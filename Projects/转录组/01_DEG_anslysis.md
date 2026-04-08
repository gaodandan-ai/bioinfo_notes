```R
# ==========================================================
# 批量分析 1h / 4h / 24h 转录组差异表达
# 适用环境：Windows 本地 R 环境
# 分析方法：
#   edgeR + limma-voom
#
# 输出内容：
#   1) 每个时间点的全基因差异结果
#   2) 每个时间点的显著差异基因
#   3) 每个时间点的归一化表达矩阵
#   4) 三个时间点的DEG汇总表
# ==========================================================


# ==========================================================
# 第 1 步：加载所需 R 包
# 说明：
# - 如果某个包不存在，则自动安装
# - limma 和 edgeR 属于 Bioconductor 包，所以用 BiocManager 安装
# ==========================================================

if (!require("BiocManager")) install.packages("BiocManager")

required_pkgs <- c("limma", "edgeR", "jsonlite", "dplyr")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("limma", "edgeR")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg, update = FALSE)
    }
    library(pkg, character.only = TRUE)
  }
}


# ==========================================================
# 第 2 步：设置项目根目录和子目录
# 说明：
# - 以后如果项目整体移动位置，只需要修改 project_dir 这一行
# ==========================================================

project_dir <- "F:/thermal_transcriptome"

# 原始数据目录
data_dir <- file.path(project_dir, "data", "raw")

# 脚本目录（只是为了项目结构清晰，这里不直接使用）
script_dir <- file.path(project_dir, "scripts")

# 结果总目录
result_dir <- file.path(project_dir, "results", "multi_timepoint_analysis")

# 日志目录（当前脚本未实际写日志，但提前预留）
log_dir <- file.path(project_dir, "logs")

# 输入文件路径
counts_file <- file.path(data_dir, "PCA-all-wqz-ori.csv")


# ==========================================================
# 第 3 步：自动创建需要的目录
# 说明：
# - 如果目录不存在，就自动创建
# - recursive = TRUE 表示如果上级目录不存在，也一起创建
# ==========================================================

dir_list <- c(data_dir, script_dir, result_dir, log_dir)

for (d in dir_list) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    cat("✅ 已创建目录：", d, "\n")
  }
}


# ==========================================================
# 第 4 步：定义三个时间点的分析参数
# 说明：
# - 每个时间点都需要告诉程序：
#   1) 使用哪些样本列
#   2) 对照组名字
#   3) 处理组名字
#   4) 结果输出到哪个子目录
#   5) 对比名称是什么
#
# 这样做的好处：
# - 以后新增 48h 或 72h，只需要在这里加一段参数
# - 主分析循环不用重复写三遍
# ==========================================================

timepoint_params <- list(
  "1h" = list(
    count_cols = c("C-1-1", "C-1-2", "C-1-3", "P-1-1", "P-1-2", "P-1-3"),
    ctrl_group = "Control-1h",
    treat_group = "Treatment-1h",
    output_subdir = "output_1h",
    contrast_name = "Treat1h_vs_Control1h"
  ),
  
  "4h" = list(
    count_cols = c("C-4-1", "C-4-2", "C-4-3", "P-4-1", "P-4-2", "P-4-3"),
    ctrl_group = "Control-4h",
    treat_group = "Treatment-4h",
    output_subdir = "output_4h",
    contrast_name = "Treat4h_vs_Control4h"
  ),
  
  "24h" = list(
    count_cols = c("C-24-1", "C-24-2", "C-24-3", "P-24-1", "P-24-2", "P-24-3"),
    ctrl_group = "Control-24h",
    treat_group = "Treatment-24h",
    output_subdir = "output_24h",
    contrast_name = "Treat24h_vs_Control24h"
  )
)


# ==========================================================
# 第 5 步：检查输入文件是否存在
# 说明：
# - 如果原始数据文件不存在，就直接停止程序
# - 这样能避免后面报一堆不容易理解的错误
# ==========================================================

if (!file.exists(counts_file)) {
  stop(paste0("❌ 未找到输入文件：", counts_file))
}


# ==========================================================
# 第 6 步：读取输入文件
# 说明：
# - 仅读取一次，后面三个时间点都复用这份数据
# - check.names = FALSE 是为了保留原始列名
# - colClasses = "character" 先全部按字符读入，避免某些列自动类型转换出问题
# - UTF-8-BOM 是为了更稳妥处理 Windows 下带 BOM 的 CSV 编码
# ==========================================================

x <- read.delim(
  file = counts_file,
  sep = ",",
  check.names = FALSE,
  colClasses = "character",
  na.strings = "",
  skip = 0,
  fileEncoding = "UTF-8-BOM"
)

# 单独重新读取第一行作为列名，尽量完整保留原始格式
colnames(x) <- scan(
  file = counts_file,
  what = "",
  sep = ",",
  nlines = 1,
  strip.white = FALSE,
  quote = "\"",
  fileEncoding = "UTF-8-BOM"
)


# ==========================================================
# 第 7 步：初始化最终汇总表
# 说明：
# - 用于记录每个时间点：
#   总基因数、过滤后基因数、显著DEGs数、上调数、下调数
# - 最后会输出一个总表
# ==========================================================

deg_summary <- data.frame(
  Timepoint = character(),
  Total_Genes = integer(),
  Filtered_Genes = integer(),
  Significant_DEGs = integer(),
  Upregulated_DEGs = integer(),
  Downregulated_DEGs = integer(),
  stringsAsFactors = FALSE
)


# ==========================================================
# 第 8 步：循环分析每一个时间点
# 说明：
# - 这里是整个脚本的核心部分
# - 程序会依次分析 1h、4h、24h
# ==========================================================

for (timepoint in names(timepoint_params)) {
  
  cat("\n==================================================\n")
  cat("=== 开始分析 ", timepoint, " 转录组数据 ===\n", sep = "")
  cat("==================================================\n")
  
  # 取出当前时间点对应的参数
  params <- timepoint_params[[timepoint]]
  
  
  # ------------------------------------------------------
  # 8.1 检查当前时间点需要的样本列是否都存在
  # 说明：
  # - 比如 1h 分析必须要有 C-1-1 ~ P-1-3 这些列
  # - 如果缺列，就跳过这个时间点，避免程序报错中断
  # ------------------------------------------------------
  
  missing_cols <- setdiff(params$count_cols, colnames(x))
  
  if (length(missing_cols) > 0) {
    warning(
      paste0(
        "⚠️ ", timepoint, " 缺少样本列：",
        paste(missing_cols, collapse = ", "),
        "，跳过该时间点！"
      )
    )
    next
  }
  

  # ------------------------------------------------------
  # 8.2 创建当前时间点的输出目录
  # 例如：
  # - results/multi_timepoint_analysis/output_1h
  # - results/multi_timepoint_analysis/output_4h
  # ------------------------------------------------------
  
  current_output_dir <- file.path(result_dir, params$output_subdir)
  
  if (!dir.exists(current_output_dir)) {
    dir.create(current_output_dir, recursive = TRUE, showWarnings = FALSE)
    cat("✅ 创建输出子目录：", current_output_dir, "\n")
  }
  
  
  # ------------------------------------------------------
  # 8.3 提取当前时间点的 count 矩阵
  # 说明：
  # - 这里只取当前时间点需要的 6 个样本列
  # - 例如 1h 就取 C-1-* 和 P-1-* 这6列
  # ------------------------------------------------------
  
  counts <- x[, params$count_cols]
  
  
  # ------------------------------------------------------
  # 8.4 检查是否有非数值数据
  # 说明：
  # - 转录组 count 数据理论上应该是数值
  # - 如果有字符、空值、异常内容，转换数值时会变成 NA
  # - 这里先检查一下，方便发现数据问题
  # ------------------------------------------------------
  
  non_numeric_rows <- apply(counts, 1, function(row) {
    any(is.na(suppressWarnings(as.numeric(row))))
  })
  
  if (any(non_numeric_rows)) {
    warning(paste0("⚠️ 过滤 ", sum(non_numeric_rows), " 行非数值数据"))
  }
  
  # 把每一列转成数值型
  counts[,] <- apply(counts, 2, function(v) as.numeric(v))
  
  # 记录原始基因数
  original_gene_num <- nrow(counts)
  
  
  # ------------------------------------------------------
  # 8.5 进行基因过滤
  # 说明：
  # - keepCpm：要求 CPM > 1 的样本数至少有 3 个
  # - keepMin：要求该基因最大表达值 >= 0
  # - keep：两个条件同时满足才保留
  #
  # 为什么要过滤？
  # - 低表达基因噪音大
  # - 过滤后分析结果更稳
  # ------------------------------------------------------
  
  keepCpm <- rowSums(cpm(counts) > 1.0) >= 3
  keepMin <- apply(counts, 1, max) >= 0
  keep <- keepMin & keepCpm
  
  counts_filtered <- counts[keep, ]
  x_filtered <- x[keep, ]
  
  filtered_gene_num <- nrow(counts_filtered)
  
  cat("✅ 基因过滤：原始 ", original_gene_num, " 个 → 保留 ", filtered_gene_num, " 个\n", sep = "")
  
  
  # ------------------------------------------------------
  # 8.6 构建设计矩阵（design matrix）
  # 说明：
  # - 前3列是对照组
  # - 后3列是处理组
  #
  # 构造后类似：
  #                Control   Treatment
  # C-1-1             1         0
  # C-1-2             1         0
  # C-1-3             1         0
  # P-1-1             0         1
  # P-1-2             0         1
  # P-1-3             0         1
  # ------------------------------------------------------
  
  design <- matrix(
    data = c(rep(c(1, 0), 3), rep(c(0, 1), 3)),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(
      rownames = params$count_cols,
      colnames = c(params$ctrl_group, params$treat_group)
    )
  )
  
  
  # ------------------------------------------------------
  # 8.7 构建对比矩阵（contrast matrix）
  # 说明：
  # - 我们关心的是：处理组 - 对照组
  # - 所以系数写成 (-1, 1)
  # ------------------------------------------------------
  
  cont.matrix <- matrix(
    data = c(-1, 1),
    ncol = 1,
    dimnames = list(
      rownames = colnames(design),
      colnames = params$contrast_name
    )
  )
  
  
  # ------------------------------------------------------
  # 8.8 进行 limma-voom 差异分析
  # 步骤说明：
  # 1) calcNormFactors：计算归一化因子
  # 2) voom：把 count 数据转换为适合线性模型分析的形式
  # 3) lmFit：拟合线性模型
  # 4) contrasts.fit：应用组间对比
  # 5) eBayes：经验贝叶斯修正，提高统计稳定性
  # ------------------------------------------------------
  
  nf <- calcNormFactors(counts_filtered)
  
  y <- voom(
    counts_filtered,
    design,
    plot = FALSE,
    lib.size = colSums(counts_filtered) * nf
  )
  
  fit <- lmFit(y, design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  
  # ------------------------------------------------------
  # 8.9 提取差异分析结果
  # 说明：
  # - topTable 得到每个基因的统计结果
  # - 这里设置 n = Inf，表示保留全部基因
  # - sort.by = "none" 表示不重新排序，尽量保持原始顺序
  # ------------------------------------------------------
  
  out <- topTable(fit2, n = Inf, sort.by = "none")
  
  
  # ------------------------------------------------------
  # 8.10 整理输出表
  # 说明：
  # - 合并 log2FC、P值、FDR、平均表达量、Geneid 和原始 counts
  # - 第一列 log2FC 改成带时间点名字的列，便于后续识别
  # ------------------------------------------------------
  
  out2 <- cbind(
    log2FC = fit2$coef,
    P.Value = out$P.Value,
    adj.P.Val = out$adj.P.Val,
    AveExpr = out$AveExpr,
    x_filtered[, c("Geneid", params$count_cols)]
  )
  
  # 把第一列重命名得更明确
  colnames(out2)[1] <- paste0("log2FC_", params$contrast_name)
  
  # 同时保留一个统一列名 log2FC，方便后续筛选和统计
  # 这是为了避免后面写 out2$log2FC 时找不到列
  out2$log2FC <- out2[[1]]
  
  
  # ------------------------------------------------------
  # 8.11 保存全基因差异分析结果
  # ------------------------------------------------------
  
  result_file <- file.path(
    current_output_dir,
    paste0(timepoint, "_differential_result.csv")
  )
  
  write.csv(out2, file = result_file, row.names = FALSE, na = "")
  cat("✅ 保存差异结果：", result_file, "\n")
  
  
  # ------------------------------------------------------
  # 8.12 提取显著差异基因
  # 筛选标准：
  # - FDR < 0.05
  # - |log2FC| > 1
  # ------------------------------------------------------
  
  sig_degs <- out2[out2$adj.P.Val < 0.05 & abs(out2$log2FC) > 1, ]
  
  sig_deg_file <- file.path(
    current_output_dir,
    paste0(timepoint, "_significant_DEGs.csv")
  )
  
  write.csv(sig_degs, file = sig_deg_file, row.names = FALSE, na = "")
  cat("✅ 保存显著DEGs：", sig_deg_file, "\n")
  
  
  # ------------------------------------------------------
  # 8.13 保存归一化表达矩阵
  # 说明：
  # - y$E 是 voom 之后的表达矩阵
  # - 后续做 PCA、聚类、热图会比较常用
  # ------------------------------------------------------
  
  normalized_expr <- cbind(Geneid = x_filtered$Geneid, y$E)
  
  normalized_file <- file.path(
    current_output_dir,
    paste0(timepoint, "_normalized_expression.csv")
  )
  
  write.csv(normalized_expr, file = normalized_file, row.names = FALSE, na = "")
  cat("✅ 保存归一化矩阵：", normalized_file, "\n")
  
  
  # ------------------------------------------------------
  # 8.14 统计上调和下调基因数量
  # 说明：
  # - 上调：log2FC > 1
  # - 下调：log2FC < -1
  # ------------------------------------------------------
  
  up_degs <- sig_degs[sig_degs$log2FC > 1, ]
  down_degs <- sig_degs[sig_degs$log2FC < -1, ]
  
  
  # ------------------------------------------------------
  # 8.15 把当前时间点统计信息写入汇总表
  # ------------------------------------------------------
  
  deg_summary <- rbind(
    deg_summary,
    data.frame(
      Timepoint = timepoint,
      Total_Genes = original_gene_num,
      Filtered_Genes = filtered_gene_num,
      Significant_DEGs = nrow(sig_degs),
      Upregulated_DEGs = nrow(up_degs),
      Downregulated_DEGs = nrow(down_degs),
      stringsAsFactors = FALSE
    )
  )
  
  cat(
    "=== ",
    timepoint,
    " 分析完成！共得到 ",
    nrow(sig_degs),
    " 个显著差异基因（上调 ",
    nrow(up_degs),
    " 个，下调 ",
    nrow(down_degs),
    " 个） ===\n",
    sep = ""
  )
}


# ==========================================================
# 第 9 步：保存所有时间点的汇总结果
# ==========================================================

summary_file <- file.path(result_dir, "01_DEG_summary_all_timepoints.csv")
write.csv(deg_summary, file = summary_file, row.names = FALSE, na = "")

cat("\n✅ 保存所有时间点DEG汇总表：", summary_file, "\n")

cat("\n📊 三个时间点差异基因统计汇总：\n")
print(deg_summary, row.names = FALSE)

cat("\n🎉 所有时间点（1h / 4h / 24h）转录组差异分析全部完成！\n")
cat("📁 结果目录：", result_dir, "\n")
```