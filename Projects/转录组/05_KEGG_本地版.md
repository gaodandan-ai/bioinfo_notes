适用于本地 KEGG 注释表质量更加可靠的情况
需要文件
	Cgl_kegg.csv

结构
```
geneid,ko_number,level2_pathway_name
Cgl2807,ko00010,Glycolysis / Gluconeogenesis
Cgl0222,ko00010,Glycolysis / Gluconeogenesis
Cgl0331,ko00010,Glycolysis / Gluconeogenesis
Cgl2911,ko00010,Glycolysis / Gluconeogenesis
```

```R
# ==========================================================
# 基于本地 KEGG 注释文件的富集分析（C. glutamicum）
#
# 使用文件：
# F:/thermal_transcriptome/library/Cgl_kegg.csv
#
# 输入：
#   output_xx/02_xx_significant_DEGs.csv
#
# 输出：
#   09_xx_KEGG_enrichment_results.csv
#   09_xx_KEGG_dotplot.png
#   10_xx_KEGG_barplot.png
# ==========================================================


# ==========================================================
# 第 1 步：加载包
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
# 第 2 步：路径设置
# ==========================================================

project_dir <- "F:/thermal_transcriptome"
result_dir  <- file.path(project_dir, "results", "multi_timepoint_analysis")

kegg_file <- "F:/thermal_transcriptome/library/Cgl_kegg.csv"


# ==========================================================
# 第 3 步：读取 KEGG 注释文件
# ==========================================================

cat("📂 读取本地 KEGG 注释文件...\n")

kegg_raw <- read.csv(
  kegg_file,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)

# 统一列名（防止大小写问题）
colnames(kegg_raw) <- c("GeneID", "PathwayID", "PathwayName")

# 去除空值
kegg_raw <- kegg_raw %>%
  filter(!is.na(GeneID), !is.na(PathwayID), !is.na(PathwayName))

cat("✅ 注释条目数：", nrow(kegg_raw), "\n")


# ==========================================================
# 第 4 步：构建 KEGG 映射关系
#
# TERM2GENE：PathwayID → GeneID
# TERM2NAME：PathwayID → PathwayName
# ==========================================================

gene2path <- kegg_raw %>%
  dplyr::select(PathwayID, GeneID) %>%
  distinct()

kegg_name <- kegg_raw %>%
  dplyr::select(PathwayID, PathwayName) %>%
  distinct()

# 背景基因（非常关键）
background_genes <- unique(gene2path$GeneID)

cat("✅ 背景基因数：", length(background_genes), "\n")


# ==========================================================
# 第 5 步：时间点设置
# ==========================================================

timepoints <- c("1h", "4h", "24h")


# ==========================================================
# 第 6 步：循环富集分析
# ==========================================================

for (timepoint in timepoints) {
  
  cat("\n==================================================\n")
  cat("开始进行 ", timepoint, " KEGG 富集分析\n", sep = "")
  cat("==================================================\n")
  
  # 输入 DEG 文件
  deg_file <- file.path(
    result_dir,
    paste0("output_", timepoint),
    paste0("02_", timepoint, "_significant_DEGs.csv")
  )
  
  if (!file.exists(deg_file)) {
    warning(paste0("⚠️ 未找到文件：", deg_file))
    next
  }
  
  deg_df <- read.csv(deg_file, stringsAsFactors = FALSE)
  
  # 提取 DEG 基因
  deg_genes <- unique(trimws(deg_df$Geneid))
  cat("✅ DEG 总数：", length(deg_genes), "\n")
  
  # 与背景匹配
  deg_genes_matched <- intersect(deg_genes, background_genes)
  cat("✅ 匹配 KEGG 注释的基因数：", length(deg_genes_matched), "\n")
  
  if (length(deg_genes_matched) < 5) {
    warning("⚠️ 可用基因太少，跳过")
    next
  }
  
  
  # ========================================================
  # KEGG 富集分析（核心）
  # ========================================================
  
  kegg_enrich <- enricher(
    gene = deg_genes_matched,
    universe = background_genes,
    TERM2GENE = gene2path,
    TERM2NAME = kegg_name,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  
  if (is.null(kegg_enrich) || nrow(kegg_enrich@result) == 0) {
    cat("⚠️ 无显著通路\n")
    next
  }
  
  kegg_df <- as.data.frame(kegg_enrich)
  
  cat("✅ 富集通路数：", nrow(kegg_df), "\n")
  
  
  # ========================================================
  # 保存结果
  # ========================================================
  
  output_dir <- file.path(result_dir, paste0("output_", timepoint))
  
  result_file <- file.path(
    output_dir,
    paste0("09_", timepoint, "_KEGG_enrichment_results.csv")
  )
  
  write.csv(kegg_df, result_file, row.names = FALSE)
  
  
  # ========================================================
  # dotplot
  # ========================================================
  
  p_dot <- dotplot(
    kegg_enrich,
    showCategory = min(15, nrow(kegg_df)),
    title = paste0(timepoint, " KEGG Enrichment")
  )
  
  ggsave(
    file.path(output_dir, paste0("09_", timepoint, "_KEGG_dotplot.png")),
    p_dot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  # ========================================================
  # barplot
  # ========================================================
  
  p_bar <- barplot(
    kegg_enrich,
    showCategory = min(15, nrow(kegg_df)),
    title = paste0(timepoint, " KEGG Enrichment")
  )
  
  ggsave(
    file.path(output_dir, paste0("10_", timepoint, "_KEGG_barplot.png")),
    p_bar,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  cat("🎉 ", timepoint, " 完成！\n", sep = "")
}


cat("\n🎉 所有时间点 KEGG 富集完成（本地注释版）\n")
```