```R
# ==========================================================
# 基于本地 GO 注释文件的多时间点 GO 富集分析（最终版）
# 修改内容：
# 1. 画图风格参考第一段代码的学术风格
# 2. 上下调基因分开绘图
# 3. 每个时间点输出：
#    - 1个综合GO结果表
#    - 1张上调GO图
#    - 1张下调GO图
# ==========================================================


# ==========================================================
# 第 1 步：加载所需 R 包
# ==========================================================

if (!require("BiocManager")) install.packages("BiocManager")

required_pkgs <- c("clusterProfiler", "dplyr", "ggplot2", "stringr")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}


# ==========================================================
# 第 2 步：设置路径
# ==========================================================

project_dir <- "F:/thermal_transcriptome"
result_dir  <- file.path(project_dir, "results", "multi_timepoint_analysis")
go_file     <- "F:/thermal_transcriptome/library/Cgl_go.csv"

timepoints <- c("1h", "4h", "24h")


# ==========================================================
# 第 3 步：读取 GO 注释文件
# ==========================================================

cat("📂 正在读取本地 GO 注释文件...\n")

go_raw <- read.csv(
  go_file,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)

colnames(go_raw) <- c("GeneID", "GO_ID", "GO_Name", "GO_Class")

go_raw <- go_raw %>%
  dplyr::filter(
    !is.na(GeneID), GeneID != "",
    !is.na(GO_ID), GO_ID != "",
    !is.na(GO_Name), GO_Name != "",
    !is.na(GO_Class), GO_Class != ""
  ) %>%
  dplyr::distinct()

cat("✅ GO 注释总条目数：", nrow(go_raw), "\n")


# ==========================================================
# 第 4 步：标准化 GO 分类名称
# ==========================================================

go_raw <- go_raw %>%
  dplyr::mutate(
    GO_Class = dplyr::case_when(
      GO_Class == "biological_process" ~ "BP",
      GO_Class == "molecular_function" ~ "MF",
      GO_Class == "cellular_component" ~ "CC",
      TRUE ~ GO_Class
    )
  )

cat("✅ GO 分类统计：\n")
print(table(go_raw$GO_Class))


# ==========================================================
# 第 5 步：构建 GO 富集所需对象
# ==========================================================

background_genes <- unique(go_raw$GeneID)

go_bp <- go_raw %>% dplyr::filter(GO_Class == "BP")
go_mf <- go_raw %>% dplyr::filter(GO_Class == "MF")
go_cc <- go_raw %>% dplyr::filter(GO_Class == "CC")

bp_term2gene <- go_bp %>% dplyr::select(GO_ID, GeneID) %>% dplyr::distinct()
bp_term2name <- go_bp %>% dplyr::select(GO_ID, GO_Name) %>% dplyr::distinct()

mf_term2gene <- go_mf %>% dplyr::select(GO_ID, GeneID) %>% dplyr::distinct()
mf_term2name <- go_mf %>% dplyr::select(GO_ID, GO_Name) %>% dplyr::distinct()

cc_term2gene <- go_cc %>% dplyr::select(GO_ID, GeneID) %>% dplyr::distinct()
cc_term2name <- go_cc %>% dplyr::select(GO_ID, GO_Name) %>% dplyr::distinct()

cat("✅ 背景基因数：", length(background_genes), "\n")


# ==========================================================
# 第 6 步：定义 GO 富集函数
# ==========================================================

run_go_enrichment <- function(gene_vec, bg_genes, term2gene, term2name) {
  clusterProfiler::enricher(
    gene = gene_vec,
    universe = bg_genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.3,
    qvalueCutoff = 0.5,
    minGSSize = 3,
    maxGSSize = 5000
  )
}


# ==========================================================
# 第 7 步：提取富集结果并加标签
# direction: up / down
# go_class: BP / MF / CC
# top_n: 每类保留前几个 term
# ==========================================================

extract_go_df <- function(go_res, direction, go_class, top_n = 10) {
  
  if (is.null(go_res)) return(NULL)
  
  df <- as.data.frame(go_res)
  if (nrow(df) == 0) return(NULL)
  
  df <- df %>%
    dplyr::arrange(p.adjust, pvalue) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      Direction = direction,
      GO_Class = go_class
    )
  
  return(df)
}


# ==========================================================
# 第 8 步：把 GeneRatio 从 "a/b" 转成数值
# ==========================================================

parse_gene_ratio <- function(x) {
  sapply(x, function(v) {
    parts <- strsplit(v, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
}


# ==========================================================
# 第 9 步：学术风格主题（参考第一段代码）
# ==========================================================

theme_academic <- theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 10),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 15, unit = "pt")
    ),
    axis.title.x = element_text(
      size = 12,
      face = "bold",
      margin = margin(t = 10, unit = "pt")
    ),
    axis.title.y = element_text(
      size = 12,
      face = "bold",
      margin = margin(r = 10, unit = "pt")
    ),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 11, face = "bold", family = "Arial"),
    strip.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.key = element_rect(fill = "transparent"),
    panel.grid.major = element_line(size = 0.5, color = "gray80"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 0.8, color = "black"),
    plot.margin = margin(10, 10, 10, 10, unit = "pt")
  )


# ==========================================================
# 第 10 步：单独绘制 up 或 down 的 GO dotplot
# 风格参考第一段代码
# ==========================================================

plot_go_separate <- function(go_df, direction_filter, out_file, plot_title) {
  
  if (is.null(go_df) || nrow(go_df) == 0) {
    cat("⚠️ 无可绘制的 GO 结果\n")
    return(NULL)
  }
  
  plot_df <- go_df %>%
    dplyr::filter(Direction == direction_filter)
  
  if (nrow(plot_df) == 0) {
    cat("⚠️ ", direction_filter, " 无可绘制结果\n", sep = "")
    return(NULL)
  }
  
  plot_df$GeneRatio_num <- parse_gene_ratio(plot_df$GeneRatio)
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      GO_Class = factor(GO_Class, levels = c("BP", "MF", "CC")),
      GO_Term_Short = stringr::str_trunc(Description, width = 50, side = "right", ellipsis = "...")
    ) %>%
    dplyr::group_by(GO_Class) %>%
    dplyr::arrange(p.adjust, GeneRatio_num, .by_group = TRUE) %>%
    dplyr::mutate(
      GO_Term_Short = factor(GO_Term_Short, levels = rev(unique(GO_Term_Short)))
    ) %>%
    dplyr::ungroup()
  
  p <- ggplot(plot_df, aes(
    x = Count,
    y = GO_Term_Short,
    size = Count,
    color = -log10(p.adjust)
  )) +
    geom_point(alpha = 0.8, stroke = 0.3) +
    facet_wrap(~GO_Class, scales = "free_y", ncol = 1) +
    scale_size_continuous(
      range = c(3, 10),
      name = "Number of DEGs"
    ) +
    scale_color_gradient(
      low = "#2c7fb8",
      high = "#d7301f",
      name = "-log10(adj.P)"
    ) +
    labs(
      title = plot_title,
      x = "Number of DEGs in GO Term",
      y = "GO Term"
    ) +
    theme_academic
  
  ggsave(
    filename = out_file,
    plot = p,
    width = 12,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  
  cat("✅ ", direction_filter, " GO 图已保存：", out_file, "\n", sep = "")
}


# ==========================================================
# 第 11 步：对单个时间点进行分析
# 输出：
#   11_1h_GO_enrichment_combined.csv
#   11_1h_GO_dotplot_up.png
#   11_1h_GO_dotplot_down.png
# ==========================================================

run_one_timepoint_go <- function(timepoint, output_dir, background_genes) {
  
  deg_file <- file.path(output_dir, paste0("02_", timepoint, "_significant_DEGs.csv"))
  
  if (!file.exists(deg_file)) {
    warning(paste0("⚠️ 未找到文件：", deg_file, "，跳过"))
    return(NULL)
  }
  
  deg_df <- read.csv(
    deg_file,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  
  if (!"Geneid" %in% colnames(deg_df)) {
    warning(paste0("⚠️ ", timepoint, " DEG 文件缺少 Geneid 列，跳过"))
    return(NULL)
  }
  
  if (!"log2FC" %in% colnames(deg_df)) {
    warning(paste0("⚠️ ", timepoint, " DEG 文件缺少 log2FC 列，跳过"))
    return(NULL)
  }
  
  # --------------------------------------------------------
  # 分离上下调基因
  # --------------------------------------------------------
  up_genes <- deg_df %>%
    dplyr::filter(!is.na(Geneid), Geneid != "", !is.na(log2FC), log2FC > 0) %>%
    dplyr::pull(Geneid) %>%
    trimws() %>%
    unique()
  
  down_genes <- deg_df %>%
    dplyr::filter(!is.na(Geneid), Geneid != "", !is.na(log2FC), log2FC < 0) %>%
    dplyr::pull(Geneid) %>%
    trimws() %>%
    unique()
  
  up_genes_matched <- intersect(up_genes, background_genes)
  down_genes_matched <- intersect(down_genes, background_genes)
  
  cat("✅ ", timepoint, " 上调 DEG 数：", length(up_genes), "\n", sep = "")
  cat("✅ ", timepoint, " 下调 DEG 数：", length(down_genes), "\n", sep = "")
  cat("✅ ", timepoint, " 上调匹配 GO 注释数：", length(up_genes_matched), "\n", sep = "")
  cat("✅ ", timepoint, " 下调匹配 GO 注释数：", length(down_genes_matched), "\n", sep = "")
  
  all_res_list <- list()
  
  # --------------------------------------------------------
  # 上调基因：BP / MF / CC
  # --------------------------------------------------------
  if (length(up_genes_matched) >= 5) {
    
    cat("🔧 开始 ", timepoint, " 上调基因 GO 富集分析...\n", sep = "")
    
    up_bp_res <- run_go_enrichment(up_genes_matched, background_genes, bp_term2gene, bp_term2name)
    up_mf_res <- run_go_enrichment(up_genes_matched, background_genes, mf_term2gene, mf_term2name)
    up_cc_res <- run_go_enrichment(up_genes_matched, background_genes, cc_term2gene, cc_term2name)
    
    all_res_list[[length(all_res_list) + 1]] <- extract_go_df(up_bp_res, "up", "BP", top_n = 10)
    all_res_list[[length(all_res_list) + 1]] <- extract_go_df(up_mf_res, "up", "MF", top_n = 10)
    all_res_list[[length(all_res_list) + 1]] <- extract_go_df(up_cc_res, "up", "CC", top_n = 10)
    
  } else {
    cat("⚠️ ", timepoint, " 上调基因太少，跳过上调 GO 富集\n", sep = "")
  }
  
  # --------------------------------------------------------
  # 下调基因：BP / MF / CC
  # --------------------------------------------------------
  if (length(down_genes_matched) >= 5) {
    
    cat("🔧 开始 ", timepoint, " 下调基因 GO 富集分析...\n", sep = "")
    
    down_bp_res <- run_go_enrichment(down_genes_matched, background_genes, bp_term2gene, bp_term2name)
    down_mf_res <- run_go_enrichment(down_genes_matched, background_genes, mf_term2gene, mf_term2name)
    down_cc_res <- run_go_enrichment(down_genes_matched, background_genes, cc_term2gene, cc_term2name)
    
    all_res_list[[length(all_res_list) + 1]] <- extract_go_df(down_bp_res, "down", "BP", top_n = 10)
    all_res_list[[length(all_res_list) + 1]] <- extract_go_df(down_mf_res, "down", "MF", top_n = 10)
    all_res_list[[length(all_res_list) + 1]] <- extract_go_df(down_cc_res, "down", "CC", top_n = 10)
    
  } else {
    cat("⚠️ ", timepoint, " 下调基因太少，跳过下调 GO 富集\n", sep = "")
  }
  
  combined_go_df <- dplyr::bind_rows(all_res_list)
  
  if (is.null(combined_go_df) || nrow(combined_go_df) == 0) {
    cat("⚠️ ", timepoint, " 无可输出的 GO 富集结果\n", sep = "")
    return(NULL)
  }
  
  # 排序
  combined_go_df <- combined_go_df %>%
    dplyr::arrange(Direction, GO_Class, p.adjust, pvalue)
  
  # 保存综合结果表
  csv_file <- file.path(output_dir, paste0("11_", timepoint, "_GO_enrichment_combined.csv"))
  write.csv(combined_go_df, csv_file, row.names = FALSE)
  cat("✅ 综合 GO 结果表已保存：", csv_file, "\n")
  
  # 保存上调图
  plot_file_up <- file.path(output_dir, paste0("11_", timepoint, "_GO_dotplot_up.png"))
  plot_go_separate(
    combined_go_df,
    direction_filter = "up",
    out_file = plot_file_up,
    plot_title = paste0(timepoint, " Upregulated DEGs: GO Enrichment")
  )
  
  # 保存下调图
  plot_file_down <- file.path(output_dir, paste0("11_", timepoint, "_GO_dotplot_down.png"))
  plot_go_separate(
    combined_go_df,
    direction_filter = "down",
    out_file = plot_file_down,
    plot_title = paste0(timepoint, " Downregulated DEGs: GO Enrichment")
  )
  
  cat("🎉 ", timepoint, " GO 分析完成！\n\n", sep = "")
}


# ==========================================================
# 第 12 步：批量运行所有时间点
# ==========================================================

for (tp in timepoints) {
  cat("==================================================\n")
  cat("🚀 开始处理时间点：", tp, "\n")
  cat("==================================================\n")
  
  output_dir <- file.path(result_dir, paste0("output_", tp))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  run_one_timepoint_go(
    timepoint = tp,
    output_dir = output_dir,
    background_genes = background_genes
  )
}

cat("🎉 所有时间点 GO 富集分析完成！\n")
```