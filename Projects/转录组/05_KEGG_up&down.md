```R
# ==========================================================
# 批量进行 1h / 4h / 24h 的 KEGG 富集分析
# 分别对 上调基因 和 下调基因 做 KEGG
# ==========================================================

if (!require("BiocManager")) install.packages("BiocManager")

required_pkgs <- c("clusterProfiler", "dplyr", "ggplot2")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}

project_dir <- "F:/thermal_transcriptome"
result_dir  <- file.path(project_dir, "results", "multi_timepoint_analysis")

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

kegg_organism <- "cgl"
kegg_keytype  <- "kegg"

for (timepoint in names(timepoint_params)) {
  
  cat("\n==================================================\n")
  cat("开始处理 ", timepoint, " 的 KEGG 富集分析\n", sep = "")
  cat("==================================================\n")
  
  params <- timepoint_params[[timepoint]]
  current_output_dir <- file.path(result_dir, params$output_subdir)
  sig_deg_path <- file.path(current_output_dir, params$sig_deg_file)
  
  if (!file.exists(sig_deg_path)) {
    warning(paste0("⚠️ 未找到文件：", sig_deg_path, "，跳过 ", timepoint))
    next
  }
  
  cat("✅ 读取显著差异基因文件：", sig_deg_path, "\n")
  
  deg_df <- read.csv(
    sig_deg_path,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM",
    na.strings = ""
  )
  
  if (!"Geneid" %in% colnames(deg_df)) {
    warning(paste0("⚠️ ", timepoint, " 文件缺少 Geneid 列，跳过"))
    next
  }
  
  if (!"log2FC" %in% colnames(deg_df)) {
    warning(paste0("⚠️ ", timepoint, " 文件缺少 log2FC 列，跳过"))
    next
  }
  
  # 上调基因
  up_gene_list <- deg_df %>%
    filter(!is.na(Geneid), Geneid != "", !is.na(log2FC), log2FC > 0) %>%
    pull(Geneid) %>%
    unique()
  
  # 下调基因
  down_gene_list <- deg_df %>%
    filter(!is.na(Geneid), Geneid != "", !is.na(log2FC), log2FC < 0) %>%
    pull(Geneid) %>%
    unique()
  
  cat("✅ ", timepoint, " 上调基因数：", length(up_gene_list), "\n", sep = "")
  cat("✅ ", timepoint, " 下调基因数：", length(down_gene_list), "\n", sep = "")
  
  
  # ========================================================
  # 上调基因 KEGG
  # ========================================================
  if (length(up_gene_list) >= 5) {
    
    cat("🔧 开始进行 ", timepoint, " 上调基因 KEGG 富集分析...\n", sep = "")
    
    up_kegg_res <- tryCatch(
      {
        enrichKEGG(
          gene          = up_gene_list,
          organism      = kegg_organism,
          keyType       = kegg_keytype,
          pvalueCutoff  = 0.5,
          qvalueCutoff  = 0.2,
          pAdjustMethod = "BH",
          minGSSize     = 5,
          maxGSSize     = 500
        )
      },
      error = function(e) {
        warning(paste0("⚠️ ", timepoint, " 上调基因 KEGG 富集失败：", e$message))
        return(NULL)
      }
    )
    
    if (!is.null(up_kegg_res) && nrow(as.data.frame(up_kegg_res)) > 0) {
      
      up_kegg_df <- as.data.frame(up_kegg_res)
      
      write.csv(
        up_kegg_df,
        file.path(current_output_dir, paste0("09_", timepoint, "_up_KEGG_enrichment_results.csv")),
        row.names = FALSE,
        na = ""
      )
      
      p_up_dot <- dotplot(
        up_kegg_res,
        showCategory = min(15, nrow(up_kegg_df)),
        title = paste0(timepoint, " Up-regulated KEGG Enrichment")
      )
      
      ggsave(
        filename = file.path(current_output_dir, paste0("09_", timepoint, "_up_KEGG_dotplot.png")),
        plot = p_up_dot,
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white"
      )
      
      p_up_bar <- barplot(
        up_kegg_res,
        showCategory = min(15, nrow(up_kegg_df)),
        title = paste0(timepoint, " Up-regulated KEGG Enrichment")
      )
      
      ggsave(
        filename = file.path(current_output_dir, paste0("10_", timepoint, "_up_KEGG_barplot.png")),
        plot = p_up_bar,
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white"
      )
      
      cat("✅ ", timepoint, " 上调基因 KEGG 结果已保存\n", sep = "")
      
    } else {
      warning(paste0("⚠️ ", timepoint, " 上调基因未获得显著 KEGG 富集结果"))
    }
    
  } else {
    warning(paste0("⚠️ ", timepoint, " 上调基因数少于 5，跳过"))
  }
  
  
  # ========================================================
  # 下调基因 KEGG
  # ========================================================
  if (length(down_gene_list) >= 5) {
    
    cat("🔧 开始进行 ", timepoint, " 下调基因 KEGG 富集分析...\n", sep = "")
    
    down_kegg_res <- tryCatch(
      {
        enrichKEGG(
          gene          = down_gene_list,
          organism      = kegg_organism,
          keyType       = kegg_keytype,
          pvalueCutoff  = 0.5,
          qvalueCutoff  = 0.2,
          pAdjustMethod = "BH",
          minGSSize     = 5,
          maxGSSize     = 500
        )
      },
      error = function(e) {
        warning(paste0("⚠️ ", timepoint, " 下调基因 KEGG 富集失败：", e$message))
        return(NULL)
      }
    )
    
    if (!is.null(down_kegg_res) && nrow(as.data.frame(down_kegg_res)) > 0) {
      
      down_kegg_df <- as.data.frame(down_kegg_res)
      
      write.csv(
        down_kegg_df,
        file.path(current_output_dir, paste0("09_", timepoint, "_down_KEGG_enrichment_results.csv")),
        row.names = FALSE,
        na = ""
      )
      
      p_down_dot <- dotplot(
        down_kegg_res,
        showCategory = min(15, nrow(down_kegg_df)),
        title = paste0(timepoint, " Down-regulated KEGG Enrichment")
      )
      
      ggsave(
        filename = file.path(current_output_dir, paste0("09_", timepoint, "_down_KEGG_dotplot.png")),
        plot = p_down_dot,
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white"
      )
      
      p_down_bar <- barplot(
        down_kegg_res,
        showCategory = min(15, nrow(down_kegg_df)),
        title = paste0(timepoint, " Down-regulated KEGG Enrichment")
      )
      
      ggsave(
        filename = file.path(current_output_dir, paste0("10_", timepoint, "_down_KEGG_barplot.png")),
        plot = p_down_bar,
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white"
      )
      
      cat("✅ ", timepoint, " 下调基因 KEGG 结果已保存\n", sep = "")
      
    } else {
      warning(paste0("⚠️ ", timepoint, " 下调基因未获得显著 KEGG 富集结果"))
    }
    
  } else {
    warning(paste0("⚠️ ", timepoint, " 下调基因数少于 5，跳过"))
  }
}

cat("\n==================================================\n")
cat("🎉 所有时间点上调/下调基因 KEGG 富集分析完成！\n")
cat("结果目录：", result_dir, "\n")
cat("==================================================\n")
```