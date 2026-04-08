```R
# ==========================================================
# 批量绘制 1h / 4h / 24h 的火山图和 MA 图（增强版）
#
# 功能：
# 1) 每个时间点分别绘制火山图
# 2) 每个时间点分别绘制 MA 图
# 3) 火山图自动标注 top 10 上调 / top 10 下调基因
# 4) 导出火山图和 MA 图对应的作图数据 CSV
#
# 输入文件：
#   output_1h/01_1h_differential_result.csv
#   output_4h/01_4h_differential_result.csv
#   output_24h/01_24h_differential_result.csv
#
# 输出文件（以 1h 为例）：
#   output_1h/07_1h_volcano_plot.png
#   output_1h/07_1h_volcano_data.csv
#   output_1h/08_1h_MA_plot.png
#   output_1h/08_1h_MA_data.csv
# ==========================================================


# ==========================================================
# 第 1 步：加载需要的 R 包
# 说明：
# - dplyr：数据整理
# - ggplot2：画图
# - ggrepel：给火山图加不重叠标签
# ==========================================================

required_pkgs <- c("dplyr", "ggplot2", "ggrepel")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, update = FALSE)
    library(pkg, character.only = TRUE)
  }
}


# ==========================================================
# 第 2 步：设置项目路径
# ==========================================================

project_dir <- "F:/thermal_transcriptome"
result_dir  <- file.path(project_dir, "results", "multi_timepoint_analysis")


# ==========================================================
# 第 3 步：定义三个时间点参数
# ==========================================================

timepoint_params <- list(
  "1h" = list(
    output_subdir = "output_1h",
    contrast_name = "Treat1h_vs_Control1h",
    diff_file = "01_1h_differential_result.csv"
  ),
  "4h" = list(
    output_subdir = "output_4h",
    contrast_name = "Treat4h_vs_Control4h",
    diff_file = "01_4h_differential_result.csv"
  ),
  "24h" = list(
    output_subdir = "output_24h",
    contrast_name = "Treat24h_vs_Control24h",
    diff_file = "01_24h_differential_result.csv"
  )
)


# ==========================================================
# 第 4 步：循环处理每个时间点
# ==========================================================

for (timepoint in names(timepoint_params)) {
  
  cat("\n==================================================\n")
  cat("开始处理 ", timepoint, " 的火山图和 MA 图\n", sep = "")
  cat("==================================================\n")
  
  params <- timepoint_params[[timepoint]]
  current_output_dir <- file.path(result_dir, params$output_subdir)
  diff_file <- file.path(current_output_dir, params$diff_file)
  
  # --------------------------------------------------------
  # 4.1 检查输入文件是否存在
  # --------------------------------------------------------
  if (!file.exists(diff_file)) {
    warning(paste0("⚠️ 未找到文件：", diff_file, "，跳过 ", timepoint))
    next
  }
  
  cat("✅ 读取差异分析结果：", diff_file, "\n")
  
  # --------------------------------------------------------
  # 4.2 读取差异分析结果
  # --------------------------------------------------------
  out2 <- read.csv(
    diff_file,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM",
    na.strings = ""
  )
  
  # --------------------------------------------------------
  # 4.3 检查关键列
  # 说明：
  # - Geneid：标注基因名需要
  # - log2FC：横轴 / MA图纵轴
  # - adj.P.Val：显著性
  # - AveExpr：MA图横轴
  # --------------------------------------------------------
  required_cols <- c("Geneid", "log2FC", "adj.P.Val", "AveExpr")
  missing_cols <- setdiff(required_cols, colnames(out2))
  
  if (length(missing_cols) > 0) {
    warning(
      paste0(
        "⚠️ ", timepoint, " 缺少关键列：",
        paste(missing_cols, collapse = ", "),
        "，跳过该时间点"
      )
    )
    next
  }
  
  # --------------------------------------------------------
  # 4.4 过滤缺失值
  # --------------------------------------------------------
  out2 <- out2 %>%
    filter(
      !is.na(Geneid),
      !is.na(log2FC),
      !is.na(adj.P.Val),
      !is.na(AveExpr)
    )
  
  if (nrow(out2) == 0) {
    warning(paste0("⚠️ ", timepoint, " 没有可用于作图的数据，跳过"))
    next
  }
  
  
  # ========================================================
  # 第 5 步：构建火山图数据
  # 说明：
  # - group：上调 / 下调 / 不显著
  # - neg_log10_fdr：火山图纵轴
  # ========================================================
  
  volcano_data <- out2 %>%
    mutate(
      group = case_when(
        adj.P.Val < 0.05 & log2FC > 1  ~ "Upregulated",
        adj.P.Val < 0.05 & log2FC < -1 ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      neg_log10_fdr = -log10(adj.P.Val)
    )
  
  # 处理 adj.P.Val = 0 导致 Inf 的情况
  if (any(is.infinite(volcano_data$neg_log10_fdr))) {
    max_finite <- max(
      volcano_data$neg_log10_fdr[is.finite(volcano_data$neg_log10_fdr)],
      na.rm = TRUE
    )
    volcano_data$neg_log10_fdr[is.infinite(volcano_data$neg_log10_fdr)] <- max_finite + 1
  }
  
  
  # ========================================================
  # 第 6 步：选出要标注的 top 基因
  # 说明：
  # - 上调基因：log2FC > 1 且 adj.P.Val < 0.05
  # - 下调基因：log2FC < -1 且 adj.P.Val < 0.05
  #
  # 这里按显著性优先排序：
  # - 先取 adj.P.Val 最小的基因
  # - 如果你以后更想强调倍数变化，也可以改成按 abs(log2FC) 排序
  # ========================================================
  
  top_up_genes <- volcano_data %>%
    filter(group == "Upregulated") %>%
    arrange(adj.P.Val, desc(log2FC)) %>%
    slice_head(n = 10)
  
  top_down_genes <- volcano_data %>%
    filter(group == "Downregulated") %>%
    arrange(adj.P.Val, log2FC) %>%
    slice_head(n = 10)
  
  label_data <- bind_rows(top_up_genes, top_down_genes)
  
  cat("✅ 上调标注基因数：", nrow(top_up_genes), "\n")
  cat("✅ 下调标注基因数：", nrow(top_down_genes), "\n")
  
  
  # ========================================================
  # 第 7 步：保存火山图数据 CSV
  # 说明：
  # - 这个 CSV 就是火山图实际用到的数据表
  # - 你后面如果要复查、筛基因、重画图，都很方便
  # ========================================================
  
  volcano_csv_file <- file.path(
    current_output_dir,
    paste0("07_", timepoint, "_volcano_data.csv")
  )
  
  write.csv(volcano_data, volcano_csv_file, row.names = FALSE, na = "")
  cat("✅ 火山图数据已保存：", volcano_csv_file, "\n")
  
  
  # ========================================================
  # 第 8 步：绘制火山图
  # ========================================================
  
  cat("🔧 生成火山图...\n")
  
  volcano_plot <- ggplot(
    volcano_data,
    aes(x = log2FC, y = neg_log10_fdr, color = group)
  ) +
    geom_point(alpha = 0.6, size = 1) +
    
    # 阈值参考线
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.5
    ) +
    
    # 给 top 基因添加文字标签
    geom_text_repel(
      data = label_data,
      aes(label = Geneid),
      size = 3,
      box.padding = 0.4,
      point.padding = 0.2,
      segment.color = "gray50",
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    
    # 配色
    scale_color_manual(
      values = c(
        "Upregulated" = "#E64B35",
        "Downregulated" = "#3C5488",
        "Not significant" = "gray60"
      )
    ) +
    
    labs(
      x = paste0("log2(Fold Change) (", params$contrast_name, ")"),
      y = "-log10(Adjusted P-value)",
      title = paste0(timepoint, ": Volcano Plot"),
      color = "Group"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid = element_blank()
    )
  
  volcano_file <- file.path(
    current_output_dir,
    paste0("07_", timepoint, "_volcano_plot.png")
  )
  
  ggsave(
    filename = volcano_file,
    plot = volcano_plot,
    width = 7,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  cat("✅ 火山图已保存：", volcano_file, "\n")
  
  
  # ========================================================
  # 第 9 步：构建 MA 图数据
  # 说明：
  # - MA 图和火山图分组逻辑一致
  # - 只是横纵轴不同
  # ========================================================
  
  ma_data <- out2 %>%
    mutate(
      group = case_when(
        adj.P.Val < 0.05 & log2FC > 1  ~ "Upregulated",
        adj.P.Val < 0.05 & log2FC < -1 ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  
  # ========================================================
  # 第 10 步：保存 MA 图数据 CSV
  # ========================================================
  
  ma_csv_file <- file.path(
    current_output_dir,
    paste0("08_", timepoint, "_MA_data.csv")
  )
  
  write.csv(ma_data, ma_csv_file, row.names = FALSE, na = "")
  cat("✅ MA图数据已保存：", ma_csv_file, "\n")
  
  
  # ========================================================
  # 第 11 步：绘制 MA 图
  # ========================================================
  
  cat("🔧 生成 MA 图...\n")
  
  ma_plot <- ggplot(
    ma_data,
    aes(x = AveExpr, y = log2FC, color = group)
  ) +
    geom_point(alpha = 0.6, size = 1) +
    
    # 参考线
    geom_hline(
      yintercept = 0,
      linetype = "solid",
      color = "black",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = c(-1, 1),
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.5
    ) +
    
    # 配色
    scale_color_manual(
      values = c(
        "Upregulated" = "#E64B35",
        "Downregulated" = "#3C5488",
        "Not significant" = "gray60"
      )
    ) +
    
    labs(
      x = "Average Expression (log2)",
      y = paste0("log2(Fold Change) (", params$contrast_name, ")"),
      title = paste0(timepoint, ": MA Plot"),
      color = "Group"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid = element_blank()
    )
  
  ma_file <- file.path(
    current_output_dir,
    paste0("08_", timepoint, "_MA_plot.png")
  )
  
  ggsave(
    filename = ma_file,
    plot = ma_plot,
    width = 6,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
  cat("✅ MA 图已保存：", ma_file, "\n")
  cat("🎉 ", timepoint, " 处理完成！\n", sep = "")
}


# ==========================================================
# 第 12 步：全部完成
# ==========================================================

cat("\n==================================================\n")
cat("🎉 所有时间点的火山图、MA图和对应CSV都已生成完成！\n")
cat("结果目录：", result_dir, "\n")
cat("==================================================\n")
```