```R
# ==========================================================
# 批量进行 1h / 4h / 24h 的 KEGG 富集分析
#
# 适配项目：
# F:/thermal_transcriptome
#
# 输入文件：
#   output_1h/02_1h_significant_DEGs.csv
#   output_4h/02_4h_significant_DEGs.csv
#   output_24h/02_24h_significant_DEGs.csv
#
# 输出文件（以 1h 为例）：
#   output_1h/09_1h_KEGG_enrichment_results.csv
#   output_1h/09_1h_KEGG_dotplot.png
#   output_1h/10_1h_KEGG_barplot.png
#
# ==========================================================


# ==========================================================
# 第 1 步：加载需要的 R 包
# 说明：
# - clusterProfiler：做 KEGG 富集分析
# - dplyr：整理数据
# - ggplot2：clusterProfiler 的可视化依赖它
# ==========================================================

if (!require("BiocManager")) install.packages("BiocManager")

required_pkgs <- c("clusterProfiler", "dplyr", "ggplot2")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}


# ==========================================================
# 第 2 步：设置项目路径
# ==========================================================

project_dir <- "F:/thermal_transcriptome"
result_dir  <- file.path(project_dir, "results", "multi_timepoint_analysis")


# ==========================================================
# 第 3 步：定义时间点参数
# 说明：
# - sig_deg_file：显著差异基因文件
# - contrast_name：仅用于打印信息，便于查看
# ==========================================================

timepoint_params <- list(
  "1h" = list(
    output_subdir = "output_1h",
    sig_deg_file = "02_1h_significant_DEGs.csv",
    contrast_name = "Treat1h_vs_Control1h"
  ),
  "4h" = list(
    output_subdir = "output_4h",
    sig_deg_file = "02_4h_significant_DEGs.csv",
    contrast_name = "Treat4h_vs_Control4h"
  ),
  "24h" = list(
    output_subdir = "output_24h",
    sig_deg_file = "02_24h_significant_DEGs.csv",
    contrast_name = "Treat24h_vs_Control24h"
  )
)


# ==========================================================
# 第 4 步：设置 KEGG 物种参数
#
# 重要说明：
# - Corynebacterium glutamicum ATCC 13032 在 KEGG 中可见 cgl 
#
# keyType = "kegg" 的前提：
# - 输入基因 ID 本身就是 KEGG 可识别的基因 ID
# - 对于 C. glutamicum，常见就是 cglxxxx 这种
# ==========================================================

kegg_organism <- "cgl"
kegg_keytype  <- "kegg"


# ==========================================================
# 第 5 步：循环处理每个时间点
# ==========================================================

for (timepoint in names(timepoint_params)) {
  
  cat("\n==================================================\n")
  cat("开始进行 ", timepoint, " 的 KEGG 富集分析\n", sep = "")
  cat("==================================================\n")
  
  params <- timepoint_params[[timepoint]]
  current_output_dir <- file.path(result_dir, params$output_subdir)
  sig_deg_path <- file.path(current_output_dir, params$sig_deg_file)
  
  # --------------------------------------------------------
  # 5.1 检查输入文件是否存在
  # --------------------------------------------------------
  if (!file.exists(sig_deg_path)) {
    warning(paste0("⚠️ 未找到文件：", sig_deg_path, "，跳过 ", timepoint))
    next
  }
  
  cat("✅ 读取显著差异基因文件：", sig_deg_path, "\n")
  
  # --------------------------------------------------------
  # 5.2 读取显著差异基因数据
  # --------------------------------------------------------
  deg_df <- read.csv(
    sig_deg_path,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM",
    na.strings = ""
  )
  
  # --------------------------------------------------------
  # 5.3 检查关键列
  # 说明：
  # - Geneid：做 KEGG 富集一定要有基因 ID
  # --------------------------------------------------------
  if (!"Geneid" %in% colnames(deg_df)) {
    warning(paste0("⚠️ ", timepoint, " 文件缺少 Geneid 列，跳过"))
    next
  }
  
  # 去掉缺失值并去重
  gene_list <- deg_df %>%
    filter(!is.na(Geneid), Geneid != "") %>%
    pull(Geneid) %>%
    unique()
  
  cat("✅ ", timepoint, " 显著差异基因数（去重后）：", length(gene_list), "\n", sep = "")
  
  # 基因太少时富集通常没有统计意义
  if (length(gene_list) < 5) {
    warning(paste0("⚠️ ", timepoint, " 基因数少于 5，通常不足以进行稳定的 KEGG 富集，跳过"))
    next
  }
  
  
  # ========================================================
  # 第 6 步：进行 KEGG 富集分析
  #
  # 说明：
  # - enrichKEGG() 适合做 ORA（过表达分析）
  # - pvalueCutoff / qvalueCutoff：显著性阈值
  # - minGSSize：通路最少基因数
  # - maxGSSize：通路最多基因数
  # ========================================================
  
  cat("🔧 开始 KEGG 富集分析...\n")
  
  kegg_res <- tryCatch(
    {
      enrichKEGG(
        gene          = gene_list,
        organism      = kegg_organism,
        keyType       = kegg_keytype,
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        pAdjustMethod = "BH",
        minGSSize     = 5,
        maxGSSize     = 500
      )
    },
    error = function(e) {
      warning(paste0("⚠️ ", timepoint, " KEGG 富集失败：", e$message))
      return(NULL)
    }
  )
  
  # 如果没有结果，跳过绘图
  if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0) {
    warning(paste0("⚠️ ", timepoint, " 未获得显著 KEGG 富集结果"))
    next
  }
  
  kegg_df <- as.data.frame(kegg_res)
  
  cat("✅ ", timepoint, " 获得显著 KEGG 通路数：", nrow(kegg_df), "\n", sep = "")
  
  
  # ========================================================
  # 第 7 步：保存 KEGG 富集结果表
  # ========================================================
  
  kegg_csv_file <- file.path(
    current_output_dir,
    paste0("09_", timepoint, "_KEGG_enrichment_results.csv")
  )
  
  write.csv(kegg_df, kegg_csv_file, row.names = FALSE, na = "")
  cat("✅ KEGG 富集结果已保存：", kegg_csv_file, "\n")
  
  
  # ========================================================
  # 第 8 步：绘制 KEGG dotplot
  # 说明：
  # - dotplot 常用于展示前若干条显著通路
  # - showCategory = 15：显示前 15 条
  # ========================================================
  
  dotplot_file <- file.path(
    current_output_dir,
    paste0("09_", timepoint, "_KEGG_dotplot.png")
  )
  
  p_dot <- dotplot(
    kegg_res,
    showCategory = min(15, nrow(kegg_df)),
    title = paste0(timepoint, " KEGG Enrichment")
  )
  
  ggsave(
    filename = dotplot_file,
    plot = p_dot,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  cat("✅ KEGG dotplot 已保存：", dotplot_file, "\n")
  
  
  # ========================================================
  # 第 9 步：绘制 KEGG barplot
  # 说明：
  # - barplot 更直观地展示富集条目和显著性
  # ========================================================
  
  barplot_file <- file.path(
    current_output_dir,
    paste0("10_", timepoint, "_KEGG_barplot.png")
  )
  
  p_bar <- barplot(
    kegg_res,
    showCategory = min(15, nrow(kegg_df)),
    title = paste0(timepoint, " KEGG Enrichment")
  )
  
  ggsave(
    filename = barplot_file,
    plot = p_bar,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  cat("✅ KEGG barplot 已保存：", barplot_file, "\n")
  cat("🎉 ", timepoint, " KEGG 富集分析完成！\n", sep = "")
}


# ==========================================================
# 第 10 步：全部完成
# ==========================================================

cat("\n==================================================\n")
cat("🎉 所有时间点 KEGG 富集分析完成！\n")
cat("默认使用 organism = '", kegg_organism, "'，keyType = '", kegg_keytype, "'\n", sep = "")
cat("结果目录：", result_dir, "\n")
cat("==================================================\n")
```