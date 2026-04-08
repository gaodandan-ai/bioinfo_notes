```R
- `F:/thermal_transcriptome/library/cgl_kegg_anno_alt.txt`
- `F:/thermal_transcriptome/library/cgl_pathway_name.txt`
```

```R
# ==========================================================
# 基于本地 KEGG 注释的多时间点 GSEA 分析
#
# 项目目录：
#   F:/thermal_transcriptome
#
# 使用本地 KEGG 注释：
#   F:/thermal_transcriptome/library/cgl_kegg_anno_alt.txt
#   F:/thermal_transcriptome/library/cgl_pathway_name.txt
#
# 输入文件：
#   output_1h/01_1h_differential_result.csv
#   output_4h/01_4h_differential_result.csv
#   output_24h/01_24h_differential_result.csv
#
# 输出文件（以 1h 为例）：
#   16_1h_KEGG_GSEA_results.csv
#   16_1h_KEGG_GSEA_dotplot.png
#   17_1h_KEGG_GSEA_enrichment_top1.png
# ==========================================================


# ==========================================================
# 第 1 步：加载所需包
# ==========================================================

if (!require("BiocManager")) install.packages("BiocManager")

required_pkgs <- c("clusterProfiler", "dplyr", "ggplot2", "enrichplot", "BiocParallel")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}

BiocParallel::register(BiocParallel::SerialParam())


# ==========================================================
# 第 2 步：设置路径
# ==========================================================

project_dir <- "F:/thermal_transcriptome"
result_dir  <- file.path(project_dir, "results", "multi_timepoint_analysis")

anno_file    <- "F:/thermal_transcriptome/library/cgl_kegg_anno_alt.txt"
pathway_file <- "F:/thermal_transcriptome/library/cgl_pathway_name.txt"


# ==========================================================
# 第 3 步：读取本地 KEGG 注释
# 说明：
# - gene2path：TERM2GENE，第一列 pathway，第二列 gene
# - kegg_name：TERM2NAME，第一列 pathway，第二列 pathway name
# ==========================================================

gene2path_raw <- read.delim(
  anno_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

colnames(gene2path_raw) <- c("GeneID", "PathwayID")

gene2path <- gene2path_raw %>%
  dplyr::select(PathwayID, GeneID) %>%
  dplyr::distinct()

kegg_name_raw <- read.delim(
  pathway_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(kegg_name_raw) <- c("PathwayID", "PathwayName")

# 去掉通路名后面的物种后缀
kegg_name <- kegg_name_raw %>%
  dplyr::mutate(
    PathwayName = sub(
      " - Corynebacterium glutamicum ATCC 13032 \\(Kyowa Hakko\\)$",
      "",
      PathwayName
    )
  ) %>%
  dplyr::distinct()

cat("✅ KEGG 注释条目数：", nrow(gene2path), "\n")
cat("✅ KEGG 通路数：", nrow(kegg_name), "\n")


# ==========================================================
# 第 4 步：定义时间点
# ==========================================================

timepoints <- c("1h", "4h", "24h")


# ==========================================================
# 第 5 步：循环进行 GSEA
# ==========================================================

for (timepoint in timepoints) {
  
  cat("\n==================================================\n")
  cat("开始进行 ", timepoint, " KEGG GSEA 分析\n", sep = "")
  cat("==================================================\n")
  
  output_dir <- file.path(result_dir, paste0("output_", timepoint))
  
  diff_file <- file.path(
    output_dir,
    paste0("01_", timepoint, "_differential_result.csv")
  )
  
  if (!file.exists(diff_file)) {
    warning(paste0("⚠️ 未找到文件：", diff_file, "，跳过"))
    next
  }
  
  # --------------------------------------------------------
  # 5.1 读取全基因差异结果
  # --------------------------------------------------------
  diff_df <- read.csv(
    diff_file,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  
  # 检查关键列
  required_cols <- c("Geneid", "log2FC")
  missing_cols <- setdiff(required_cols, colnames(diff_df))
  
  if (length(missing_cols) > 0) {
    warning(
      paste0(
        "⚠️ ", timepoint, " 缺少关键列：",
        paste(missing_cols, collapse = ", "),
        "，跳过"
      )
    )
    next
  }
  
  # --------------------------------------------------------
  # 5.2 构建 GSEA 排序基因列表
  #
  # GSEA 要求：
  # - 一个命名数值向量
  # - names = GeneID
  # - values = ranking score
  # - 必须按降序排序
  #
  # 这里先用 log2FC 作为排序分数
  # --------------------------------------------------------
  gsea_input <- diff_df %>%
    dplyr::select(Geneid, log2FC, P.Value) %>%
    dplyr::filter(
      !is.na(Geneid), Geneid != "",
      !is.na(log2FC),
      !is.na(P.Value),
      P.Value > 0
    ) %>%
    dplyr::mutate(
      rank_score = sign(log2FC) * (-log10(P.Value))
    ) %>%
    dplyr::group_by(Geneid) %>%
    dplyr::summarise(rank_score = mean(rank_score), .groups = "drop")
  
  
  
  # 只保留在 KEGG 注释中出现的基因
  annotated_genes <- unique(gene2path$GeneID)
  gsea_input <- gsea_input %>%
    dplyr::filter(Geneid %in% annotated_genes)
  
  cat("✅ 可用于 GSEA 的已注释基因数：", nrow(gsea_input), "\n")
  
  if (nrow(gsea_input) < 10) {
    warning(paste0("⚠️ ", timepoint, " 可用于 GSEA 的基因太少，跳过"))
    next
  }
  
  geneList <- gsea_input$rank_score
  names(geneList) <- gsea_input$Geneid
  geneList <- sort(geneList, decreasing = TRUE)
  # --------------------------------------------------------
  # 5.3 运行 GSEA
  # 说明：
  # - GSEA 用的是全基因排序，不是显著基因子集
  # - minGSSize / maxGSSize 控制通路大小
  # - pvalueCutoff 这里先放 0.05
  # --------------------------------------------------------
  cat("🔧 开始运行 GSEA...\n")
  
  gsea_res <- tryCatch(
    {
      GSEA(
        geneList = geneList,
        TERM2GENE = gene2path,
        TERM2NAME = kegg_name,
        pvalueCutoff = 0.6,
        pAdjustMethod = "BH",
        minGSSize = 3,
        maxGSSize = 500,
        verbose = FALSE,
        seed = TRUE
      )
    },
    error = function(e) {
      warning(paste0("⚠️ ", timepoint, " GSEA 运行失败：", e$message))
      return(NULL)
    }
  )
  
  if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
    cat("⚠️ 无显著 GSEA 通路\n")
    next
  }
  
  gsea_df <- as.data.frame(gsea_res)
  cat("✅ 显著 GSEA 通路数：", nrow(gsea_df), "\n")
  
  # --------------------------------------------------------
  # 5.4 保存 GSEA 结果表
  # --------------------------------------------------------
  result_file <- file.path(
    output_dir,
    paste0("16_", timepoint, "_KEGG_GSEA_results.csv")
  )
  
  write.csv(gsea_df, result_file, row.names = FALSE)
  cat("✅ GSEA 结果已保存：", result_file, "\n")
  
  # --------------------------------------------------------
  # 5.5 绘制 GSEA dotplot
  # --------------------------------------------------------
  p_dot <- dotplot(
    gsea_res,
    showCategory = min(15, nrow(gsea_df)),
    title = paste0(timepoint, " KEGG GSEA")
  )
  
  dotplot_file <- file.path(
    output_dir,
    paste0("16_", timepoint, "_KEGG_GSEA_dotplot.png")
  )
  
  ggsave(
    filename = dotplot_file,
    plot = p_dot,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  cat("✅ GSEA dotplot 已保存：", dotplot_file, "\n")
  
  # --------------------------------------------------------
  # 5.6 绘制 top1 enrichment plot
  # 说明：
  # - gseaplot2 用于展示某条通路在排序基因列表中的富集轨迹
  # - 这里默认画排名第一的通路
  # --------------------------------------------------------
  top_pathway_id <- gsea_df$ID[1]
  
  p_gsea <- gseaplot2(
    gsea_res,
    geneSetID = top_pathway_id,
    title = gsea_df$Description[1]
  )
  
  top1_file <- file.path(
    output_dir,
    paste0("17_", timepoint, "_KEGG_GSEA_enrichment_top1.png")
  )
  
  ggsave(
    filename = top1_file,
    plot = p_gsea,
    width = 9,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  
  cat("✅ Top1 enrichment plot 已保存：", top1_file, "\n")
  cat("🎉 ", timepoint, " KEGG GSEA 完成！\n", sep = "")
}


cat("\n==================================================\n")
cat("🎉 所有时间点 KEGG GSEA 分析完成！\n")
cat("==================================================\n")
```