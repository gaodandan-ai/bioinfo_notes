原始文件：PCA-思琪发-all-wqz-ori.csv  原始基因count
结构：
```
Geneid,C-1-1,C-1-2,C-1-3,P-1-1,P-1-2,P-1-3,C-4-1,C-4-2,C-4-3,P-4-1,P-4-2,P-4-3,C-24-1,C-24-2,C-24-3,P-24-1,P-24-2,P-24-3,C-0-1,C-0-2,C-0-3
Cgl0001,1024,1237,834,353,801,731,283,222,310,215,305,232,77,32,53,212,237,235,164,145,199
Cgl0002,18,8,9,11,10,1,5,19,11,13,11,7,4,0,6,28,16,10,14,26,19
Cgl0003,1120,1468,1274,688,1146,1052,1171,887,1359,1142,1079,1112,86,35,59,802,925,696,608,545,640
Cgl0004,563,540,519,255,541,448,284,212,377,327,315,296,13,8,17,316,380,314,157,135,177
```

### 1. 去除C-0-1,C-0-2,C-0-3（种子液样本，后续不分析）

```bash
awk -F',' '{
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19
}' OFS=',' PCA-思琪发-all-wqz-ori.csv > raw_count.csv

```

or
```bash
awk -F',' '{
    for(i=1;i<=19;i++){
        printf $i
        if(i<19) printf ","
    }
    printf "\n"
}' PCA-思琪发-all-wqz-ori.csv > raw_count.csv
```

or
```bash
cut -d',' -f1-19 PCA-思琪发-all-wqz-ori.csv > raw_count.csv
#`-d','`：按逗号分隔
#`-f1-19`：取前19列
```

```bash
head -1 raw_count.csv

Geneid,C-1-1,C-1-2,C-1-3,P-1-1,P-1-2,P-1-3,C-4-1,C-4-2,C-4-3,P-4-1,P-4-2,P-4-3,C-24-1,C-24-2,C-24-3,P-24-1,P-24-2,P-24-3
```


## 前处理 处理成WGCNA需要的格式

```R
rm(list = ls())

# ================================
# 0. 设置路径 & 创建输出文件夹
# ================================
setwd("G:/thermal_wqz")

outdir <- "G:/thermal_wqz/WGCNA"
dir.create(outdir, showWarnings = FALSE)


# ================================
# 1. 读取 raw count
# ================================
raw_counts <- read.csv("raw_count.csv", header = TRUE, check.names = FALSE)

# 保存原始表（备份）
write.csv(raw_counts, file.path(outdir, "step1_raw_counts_original.csv"), row.names = FALSE)


# ================================
# 2. 设置 Geneid 为行名
# ================================
rownames(raw_counts) <- raw_counts$Geneid
raw_counts <- raw_counts[, -1]

# 保存
write.csv(raw_counts, file.path(outdir, "step2_counts_with_rownames.csv"))


# ================================
# 3. 样本名 - → _
# ================================
colnames(raw_counts) <- gsub("-", "_", colnames(raw_counts))

# 保存
write.csv(raw_counts, file.path(outdir, "step3_counts_rename_samples.csv"))


# ================================
# 4. 构建 trait 表
# ================================
sample_names <- colnames(raw_counts)

trait_df <- data.frame(
  Sample = sample_names,
  Group = ifelse(grepl("^C_", sample_names), "Control", "Treatment"),
  Treatment = ifelse(grepl("^C_", sample_names), 0, 1),
  Time = ifelse(grepl("_1_", sample_names), 1,
                ifelse(grepl("_4_", sample_names), 4,
                       ifelse(grepl("_24_", sample_names), 24, NA)))
)

rownames(trait_df) <- trait_df$Sample

# 保存
write.csv(trait_df, file.path(outdir, "step4_trait_table.csv"))

# 查看 trait 表
trait_df

# ================================
# 5. 过滤低表达基因
# 这里先做一个比较常见的过滤：  
# 至少在一半样本中 count >= 10
# ================================
keep_genes <- rowSums(raw_counts >= 10) >= (ncol(raw_counts) / 2)
filtered_counts <- raw_counts[keep_genes, ]

# 看过滤后剩多少基因
dim(filtered_counts)

# 保存
write.csv(filtered_counts, file.path(outdir, "step5_filtered_counts.csv"))


# ================================
# 6. vst 标准化
# ================================
library(DESeq2)

coldata <- data.frame(row.names = colnames(filtered_counts))

dds <- DESeqDataSetFromMatrix(
  countData = round(filtered_counts),
  colData = coldata,
  design = ~ 1
)
# 由于这里只是为了做表达量变换，不做差异分析
# design 可以先设成 ~ 1

vsd <- vst(dds, blind = TRUE)
expr_matrix <- assay(vsd)

# 保存
write.csv(expr_matrix, file.path(outdir, "step6_vst_matrix.csv"))


# ================================
# 7. 过滤低变异基因
# WGCNA一般不建议拿所有基因直接跑
# 这里保留变异最大的前50%基因
# ================================
gene_var <- apply(expr_matrix, 1, var)
var_cutoff <- quantile(gene_var, 0.5)

expr_matrix_filt <- expr_matrix[gene_var > var_cutoff, ]

dim(expr_matrix_filt)
# 保存
write.csv(expr_matrix_filt, file.path(outdir, "step7_high_var_genes.csv"))


# ================================
# 8. 转置为WGCNA格式
# WGCNA要求：行=样本，列=基因
# ================================
datExpr <- t(expr_matrix_filt)

dim(datExpr)
# 保存
write.csv(datExpr, file.path(outdir, "step8_datExpr.csv"))


# ================================
# 9. 对齐 trait 表
# ================================
trait_df <- trait_df[rownames(datExpr), ]

# 保存
write.csv(trait_df, file.path(outdir, "step9_trait_aligned.csv"))


# ================================
# 10. 样本聚类
# ================================
sampleTree <- hclust(dist(datExpr), method = "average")

pdf(file.path(outdir, "step10_sample_clustering.pdf"), width = 10, height = 6)
plot(sampleTree, main = "Sample clustering",
     xlab = "", sub = "")
dev.off()


# ================================
# 11. PCA
# ================================
pca <- prcomp(datExpr, scale. = TRUE)

pdf(file.path(outdir, "step11_PCA.pdf"), width = 8, height = 6)
plot(pca$x[,1], pca$x[,2],
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 2), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 2), "%)"),
     main = "PCA",
     pch = 19)
text(pca$x[,1], pca$x[,2], labels = rownames(datExpr), pos = 3, cex = 0.7)
dev.off()


save(datExpr, trait_df, file = file.path(outdir, "WGCNA_input.RData"))
```