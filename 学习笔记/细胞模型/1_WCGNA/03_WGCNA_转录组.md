## 选 soft-threshold（软阈值，power）

WGCNA 不直接拿相关系数建网，而是把相关性做一个幂次变换：

```
a_ij = |cor(i,j)|^β    #加权   得到0~1的连接强度
```

这里的 `β` 就是 power。

作用是：

- 强相关会被强化
- 弱相关会被压低

所以第一步要选一个合适的 `β`
## 构建共表达网络

## 划分模块

## 模块和表型做相关

## 找hub gene


```R
rm(list = ls())

# ================================
# 0. 设置路径
# ================================
outdir <- "G:/thermal_wqz/WGCNA"
setwd(outdir)

# ================================
# 1. 读取输入文件
# ================================
datExpr <- read.csv("step8_datExpr.csv", row.names = 1, check.names = FALSE)
trait_df <- read.csv("step9_trait_aligned.csv", row.names = 1, check.names = FALSE)

dim(datExpr)
dim(trait_df)
head(datExpr[,1:5])
trait_df

# datExpr 必须是数值矩阵
datExpr <- as.data.frame(datExpr)
trait_df <- as.data.frame(trait_df)

# 只保留数值型性状用于相关分析
trait_num <- trait_df[, c("Treatment", "Time")]
trait_num <- as.data.frame(trait_num)

# 检查样本顺序
all(rownames(datExpr) == rownames(trait_num))

# ================================
# 2. 加载 WGCNA
# ================================

library(WGCNA)

# 允许多线程（Windows有时会忽略，不影响）
allowWGCNAThreads()

# 防止字符串自动变因子
options(stringsAsFactors = FALSE)

# ================================
# 3. 检查数据质量
# ================================
gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

if (!gsg$allOK) {
  # 如果有坏样本/坏基因，去掉
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  trait_num <- trait_num[rownames(datExpr), ]
}

# 保存检查后的矩阵
write.csv(datExpr, file.path(outdir, "step12_datExpr_checked.csv"))
write.csv(trait_num, file.path(outdir, "step12_trait_checked.csv"))

# ================================
# 4. 再做一次样本聚类（WGCNA标准步骤）
# ================================
sampleTree2 <- hclust(dist(datExpr), method = "average")

pdf(file.path(outdir, "step13_sampleTree_WGCNA.pdf"), width = 10, height = 6)
plot(sampleTree2, main = "Sample clustering before WGCNA",
     sub = "", xlab = "", cex.lab = 1.5)
dev.off()

# ================================
# 5. 选择 soft-threshold power
# ================================
powers <- c(1:20)

sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         verbose = 5,
                         networkType = "signed")

# 保存结果表
write.csv(sft$fitIndices, file.path(outdir, "step14_soft_threshold_table.csv"),
          row.names = FALSE)

# 画图1：Scale-free topology fit index
pdf(file.path(outdir, "step14_soft_threshold_plots.pdf"), width = 12, height = 5)

par(mfrow = c(1,2))

cex1 <- 0.9

plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers,
     cex = cex1)
abline(h = 0.80, col = "red")

# 画图2：Mean connectivity
plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels = powers,
     cex = cex1)

dev.off()

# ================================
# 6. 手动设置 power
# ================================
# 先默认设一个常用值；你看图后可以改
softPower <- 8

# 保存
writeLines(paste("Chosen softPower =", softPower),
           con = file.path(outdir, "step15_chosen_softPower.txt"))

# ================================
# 7. 构建网络并识别模块
# ================================
net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

# 提取模块颜色
moduleColors <- net$colors
table(moduleColors)

# 保存每个基因对应的模块
gene_module <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors
)
write.csv(gene_module, file.path(outdir, "step16_gene_module_assignment.csv"),
          row.names = FALSE)

# 画模块树
pdf(file.path(outdir, "step16_gene_dendrogram_modules.pdf"), width = 12, height = 8)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# ================================
# 8. 计算模块特征基因（Module Eigengenes）
# ================================
MEs0 <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# 保存模块特征基因矩阵
write.csv(MEs, file.path(outdir, "step17_module_eigengenes.csv"))

# ================================
# 9. 模块-性状相关分析
# ================================
moduleTraitCor <- cor(MEs, trait_num, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 保存结果
write.csv(moduleTraitCor, file.path(outdir, "step18_module_trait_cor.csv"))
write.csv(moduleTraitPvalue, file.path(outdir, "step18_module_trait_pvalue.csv"))

# 构造热图文字
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")",
                    sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

pdf(file.path(outdir, "step18_module_trait_heatmap.pdf"), width = 8, height = 10)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(trait_num),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = "Module-trait relationships")
dev.off()


# =========================================================
# 自动选前2个最重要模块
# 生成：
# 1) step18_Module-trait-relationship_boxplot.pdf
# 2) step18_gene-Module-trait-significance.pdf
# =========================================================

# -----------------------------
# 先选 Treatment 相关性最强的前2个模块
# -----------------------------
top_n <- 2
top_modules <- rownames(moduleTraitCor)[order(-abs(moduleTraitCor[, "Treatment"]))][1:top_n]

# 构建更好看的分组标签
trait_plot <- trait_num
trait_plot$Treatment_group <- factor(trait_plot$Treatment,
                                     levels = c(0, 1),
                                     labels = c("Control", "Treatment"))
trait_plot$Time_group <- factor(trait_plot$Time,
                                levels = c(1, 4, 24),
                                labels = c("1h", "4h", "24h"))

# -----------------------------
# 小函数：把 p 值转成显著性星号
# -----------------------------
p_to_star <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

# =========================================================
# 1. 模块-性状关系箱线图（升级版）
# =========================================================
pdf(file.path(outdir, "step18_Module-trait-relationship_boxplot.pdf"),
    width = 12, height = 5 * top_n)

par(mfrow = c(top_n, 2), mar = c(5, 5, 4, 1))

for (mod in top_modules) {
  me_values <- MEs[, mod]
  module_name <- sub("^ME", "", mod)
  
  # 提取热图中的相关性和p值
  cor_treat <- moduleTraitCor[mod, "Treatment"]
  p_treat   <- moduleTraitPvalue[mod, "Treatment"]
  cor_time  <- moduleTraitCor[mod, "Time"]
  p_time    <- moduleTraitPvalue[mod, "Time"]
  
  # -----------------------------
  # 左图：Treatment箱线图 + 显著性
  # -----------------------------
  boxplot(me_values ~ trait_plot$Treatment_group,
          main = paste0(mod,
                        "\nTreatment cor=", round(cor_treat, 3),
                        ", p=", signif(p_treat, 2)),
          xlab = "Treatment",
          ylab = "Module eigengene",
          outline = FALSE)
  
  stripchart(me_values ~ trait_plot$Treatment_group,
             vertical = TRUE,
             method = "jitter",
             pch = 19,
             add = TRUE)
  
  # 计算Treatment组间t检验
  p_box_treat <- tryCatch(
    t.test(me_values ~ trait_plot$Treatment_group)$p.value,
    error = function(e) NA
  )
  
  y_max <- max(me_values, na.rm = TRUE)
  y_min <- min(me_values, na.rm = TRUE)
  y_pos <- y_max + 0.12 * (y_max - y_min + 1e-6)
  
  segments(1, y_pos, 2, y_pos)
  segments(1, y_pos, 1, y_pos - 0.03 * (y_max - y_min + 1e-6))
  segments(2, y_pos, 2, y_pos - 0.03 * (y_max - y_min + 1e-6))
  text(1.5, y_pos + 0.04 * (y_max - y_min + 1e-6),
       labels = p_to_star(p_box_treat), cex = 1.2)
  
  # -----------------------------
  # 右图：Time箱线图 + ANOVA显著性
  # -----------------------------
  boxplot(me_values ~ trait_plot$Time_group,
          main = paste0(mod,
                        "\nTime cor=", round(cor_time, 3),
                        ", p=", signif(p_time, 2)),
          xlab = "Time",
          ylab = "Module eigengene",
          outline = FALSE)
  
  stripchart(me_values ~ trait_plot$Time_group,
             vertical = TRUE,
             method = "jitter",
             pch = 19,
             add = TRUE)
  
  # 计算Time的一元方差分析
  p_box_time <- tryCatch(
    summary(aov(me_values ~ trait_plot$Time_group))[[1]][["Pr(>F)"]][1],
    error = function(e) NA
  )
  
  y_max2 <- max(me_values, na.rm = TRUE)
  y_min2 <- min(me_values, na.rm = TRUE)
  y_pos2 <- y_max2 + 0.12 * (y_max2 - y_min2 + 1e-6)
  
  segments(1, y_pos2, 3, y_pos2)
  segments(1, y_pos2, 1, y_pos2 - 0.03 * (y_max2 - y_min2 + 1e-6))
  segments(3, y_pos2, 3, y_pos2 - 0.03 * (y_max2 - y_min2 + 1e-6))
  text(2, y_pos2 + 0.04 * (y_max2 - y_min2 + 1e-6),
       labels = p_to_star(p_box_time), cex = 1.2)
}

dev.off()


# =========================================================
# 2. Gene significance vs Module membership
# =========================================================
pdf(file.path(outdir, "step18_gene-Module-trait-significance.pdf"),
    width = 7, height = 6 * top_n)

par(mfrow = c(top_n, 1), mar = c(5, 5, 4, 2))

for (mod in top_modules) {
  module_name <- sub("^ME", "", mod)   # MEmagenta -> magenta
  mm_col <- mod                        # geneInfo里列名，如 MEmagenta
  
  module_genes <- geneInfo$Module == module_name
  
  if (sum(module_genes) > 0 && mm_col %in% colnames(geneInfo)) {
    
    cor_gs_mm <- cor(abs(geneInfo[module_genes, mm_col]),
                     abs(geneInfo[module_genes, "GS_Treatment"]),
                     use = "p")
    
    verboseScatterplot(
      abs(geneInfo[module_genes, mm_col]),
      abs(geneInfo[module_genes, "GS_Treatment"]),
      xlab = paste("Module Membership in", module_name, "module"),
      ylab = "Gene significance for Treatment",
      main = paste0(module_name,
                    " module: GS vs MM (r=",
                    round(cor_gs_mm, 3), ")"),
      cex.main = 1.2,
      cex.lab = 1.1,
      cex.axis = 1.0,
      col = module_name
    )
  }
}

dev.off()


# =========================================================
# 3. 保存前2个关键模块名单
# =========================================================
top_module_table <- data.frame(
  Module = top_modules,
  Treatment_cor = moduleTraitCor[top_modules, "Treatment"],
  Treatment_p = moduleTraitPvalue[top_modules, "Treatment"],
  Time_cor = moduleTraitCor[top_modules, "Time"],
  Time_p = moduleTraitPvalue[top_modules, "Time"]
)

write.csv(top_module_table,
          file.path(outdir, "step18_top2_modules_summary.csv"),
          row.names = FALSE)


# ================================
# 14. 保存工作空间
# ================================
save(datExpr, trait_num, net, moduleColors, MEs, moduleTraitCor,
     moduleTraitPvalue, geneInfo,
     file = file.path(outdir, "step22_WGCNA_main_results.RData"))
```



```R
# ================================
# 10. 找与性状最相关的模块
# ================================
# 以 Treatment 为例
treatment_cor <- moduleTraitCor[, "Treatment"]
treatment_p <- moduleTraitPvalue[, "Treatment"]

key_modules_treatment <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = treatment_cor,
  Pvalue = treatment_p
)
key_modules_treatment <- key_modules_treatment[order(-abs(key_modules_treatment$Correlation)), ]

write.csv(key_modules_treatment,
          file.path(outdir, "step19_key_modules_for_Treatment.csv"),
          row.names = FALSE)

# 以 Time 为例
time_cor <- moduleTraitCor[, "Time"]
time_p <- moduleTraitPvalue[, "Time"]

key_modules_time <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = time_cor,
  Pvalue = time_p
)
key_modules_time <- key_modules_time[order(-abs(key_modules_time$Correlation)), ]

write.csv(key_modules_time,
          file.path(outdir, "step19_key_modules_for_Time.csv"),
          row.names = FALSE)

# ================================
# 11. 计算 GS 和 MM
# ================================
# GS: gene significance
GS.Treatment <- as.numeric(cor(datExpr, trait_num$Treatment, use = "p"))
GS.Time <- as.numeric(cor(datExpr, trait_num$Time, use = "p"))

GS.p.Treatment <- corPvalueStudent(GS.Treatment, nrow(datExpr))
GS.p.Time <- corPvalueStudent(GS.Time, nrow(datExpr))

# MM: module membership
MM <- as.data.frame(cor(datExpr, MEs, use = "p"))
MM.pvalue <- as.data.frame(corPvalueStudent(as.matrix(MM), nrow(datExpr)))

# 整合成一个表
geneInfo <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors,
  GS_Treatment = GS.Treatment,
  GS_p_Treatment = GS.p.Treatment,
  GS_Time = GS.Time,
  GS_p_Time = GS.p.Time
)

geneInfo <- cbind(geneInfo, MM, MM.pvalue)

write.csv(geneInfo, file.path(outdir, "step20_gene_GS_MM_table.csv"),
          row.names = FALSE)


}
```
