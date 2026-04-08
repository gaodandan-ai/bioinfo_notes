
## 样本聚类

直观展示样本相似性
![[Pasted image 20260324180302.png]]

## 因子聚类与模块划分

将因子按照表达模式聚类到不同功能模块

![[Pasted image 20260324180419.png]]

## 模块-表型关系

关键因子筛选基于模块成员度(MM)和因子显著性(GS)精准定位核心靶点
![[Pasted image 20260324180543.png]]

## 关键因子筛选

结合**模块成员度（MM）**（因子与模块关联强度）与**因子显著性（GS）**（因子与表型相关性），筛选高MM且高GS的核心靶点。

![[Pasted image 20260324195848.png]]


## 可视化网络

将模块里的基因关系变成一张网络图
	基因A 和 基因B 相关
	基因A 和 基因C 相关

```
A —— B —— C
 \         /
    D —— E
```

已经计算过
	相关性（correlation）
	加权网络（adjacency）
	TOM 拓扑重叠

最后有一个矩阵
![[Pasted image 20260324200246.png]]

随后设置一个阈值，只保留> n的连接

	Cytoscape
	VisANT

将这些连接关系矩阵画成图

#### hub gene

看哪些节点连接基因最多
与GS+MM找到的 hub gene 是两种独立的证据

### 生成Cytoscape的网络文件

```R
# =========================================================
# Cytoscape network export for magenta and pink modules
# 输出到: G:/thermal_wqz/WGCNA/Cytoscape
# =========================================================

library(WGCNA)

cyto_dir <- file.path(outdir, "Cytoscape")
dir.create(cyto_dir, showWarnings = FALSE)

# 你要导出的模块
target_modules <- c("magenta", "pink")

# 参数：建议先不要太低，不然边太多
edge_threshold <- 0.15

# 如果你想只保留最强的边，也可以限制每个模块最多保留多少条
max_edges <- 300

for (target_module in target_modules) {
  
  cat("Processing module:", target_module, "\n")
  
  # -----------------------------
  # 1. 提取模块基因
  # -----------------------------
  in_module <- moduleColors == target_module
  mod_genes <- colnames(datExpr)[in_module]
  
  if (length(mod_genes) < 2) {
    cat("Skip", target_module, ": too few genes\n")
    next
  }
  
  datExpr_mod <- datExpr[, mod_genes, drop = FALSE]
  
  # -----------------------------
  # 2. 计算 adjacency
  # Cytoscape通常用 adjacency 导出边
  # -----------------------------
  adj_mod <- adjacency(
    datExpr_mod,
    power = softPower,
    type = "signed"
  )
  
  # 行列名
  dimnames(adj_mod) <- list(mod_genes, mod_genes)
  
  # -----------------------------
  # 3. 导出 edge 表
  # 只保留上三角，避免重复
  # -----------------------------
  edge_idx <- which(adj_mod > edge_threshold, arr.ind = TRUE)
  edge_idx <- edge_idx[edge_idx[,1] < edge_idx[,2], , drop = FALSE]
  
  if (nrow(edge_idx) == 0) {
    cat("No edges above threshold for", target_module, "\n")
    next
  }
  
  edge_df <- data.frame(
    source = rownames(adj_mod)[edge_idx[,1]],
    target = colnames(adj_mod)[edge_idx[,2]],
    weight = adj_mod[edge_idx],
    module = target_module,
    stringsAsFactors = FALSE
  )
  
  # 按权重从大到小排序
  edge_df <- edge_df[order(-edge_df$weight), ]
  
  # 如果边太多，只取前 max_edges 条
  if (nrow(edge_df) > max_edges) {
    edge_df <- edge_df[1:max_edges, ]
  }
  
  write.csv(
    edge_df,
    file.path(cyto_dir, paste0("cytoscape_edges_", target_module, ".csv")),
    row.names = FALSE
  )
  
  # -----------------------------
  # 4. 导出 node 表
  # 合并 geneInfo 里的 GS / MM 信息
  # -----------------------------
  ME_col <- paste0("ME", target_module)
  
  node_df <- data.frame(
    id = mod_genes,
    module = target_module,
    stringsAsFactors = FALSE
  )
  
  # 加入 geneInfo 信息
  geneInfo_sub <- geneInfo[match(mod_genes, geneInfo$Gene), , drop = FALSE]
  
  node_df$GS_Treatment <- geneInfo_sub$GS_Treatment
  node_df$GS_Time <- geneInfo_sub$GS_Time
  
  if (ME_col %in% colnames(geneInfo_sub)) {
    node_df$MM <- geneInfo_sub[[ME_col]]
  } else {
    node_df$MM <- NA
  }
  
  # 方向标签
  node_df$Direction_Treatment <- ifelse(node_df$GS_Treatment > 0,
                                        "Up_in_Treatment",
                                        "Down_in_Treatment")
  
  # hub 分数：可用于 Cytoscape 排序/调大小
  node_df$HubScore <- abs(node_df$GS_Treatment) * abs(node_df$MM)
  
  # 节点连接度（基于导出的 edge）
  degree_table <- table(c(edge_df$source, edge_df$target))
  node_df$Degree <- as.numeric(degree_table[node_df$id])
  node_df$Degree[is.na(node_df$Degree)] <- 0
  
  # -----------------------------
  # 5. 如果你有注释表，就合并进去
  # -----------------------------
  anno_file <- file.path(outdir, "cgl_anno.csv")
  
  if (file.exists(anno_file)) {
    anno <- read.csv(
      anno_file,
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    node_df <- merge(
      node_df,
      anno[, c("cgl-id", "gene-name", "product")],
      by.x = "id",
      by.y = "cgl-id",
      all.x = TRUE
    )
  } else {
    node_df$`gene-name` <- NA
    node_df$product <- NA
  }
  
  # 排序
  node_df <- node_df[order(-node_df$HubScore, -node_df$Degree), ]
  
  write.csv(
    node_df,
    file.path(cyto_dir, paste0("cytoscape_nodes_", target_module, ".csv")),
    row.names = FALSE
  )
}

cat("Cytoscape files exported to:", cyto_dir, "\n")
```