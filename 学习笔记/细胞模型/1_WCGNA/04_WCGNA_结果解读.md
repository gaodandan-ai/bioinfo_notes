# 结果解读

## WGCNA回答的问题

	哪些基因协同变化
	哪些模块和热胁迫相关
	哪些基因是关键调控者

### step18_module_trait_heatmap.pdf
![[Pasted image 20260324143044.png]]


#### Treatment（处理 vs 对照）

- **正相关（>0）**
    - 处理组高表达
- **负相关（<0）**
    - 对照组高表达

---

####  Time（时间）

- 正相关：
    - 随时间上升（1h → 24h）
- 负相关：
    - 随时间下降


#### ⭐ 强相关模块

一般看：

- |cor| > 0.6（可以放宽到0.5）  # treatment  热胁迫响应模块， # Time 时间响应模块
- p value < 0.05

### step19_key_modules_for_Treatment.csv

排序好的结果，重点看排序好的前2~3个结果

![[Pasted image 20260324143701.png]]

### step20_gene_GS_MM_table.csv

找hub gene的核心表

两个关键指标

	GS (gene significance)
		基因 vs 表型的相关性，越高表示越和热胁迫先关
	MM(module membership)
		基因 vs 模块的相关性，越高表示越像模块的核心成员

hub gene通常定义

	|GS| > 0.6
	|MM| > 0.8

## 具体操作

	选择step19_key_modules_for_Treatment.csv前两个模块
	提取关键基因


```R
gene_module <- read.csv("step16_gene_module_assignment.csv")

magenta_genes <- gene_module$Gene[
  gene_module$Module == "magenta"
]

length(magenta_genes)

write.csv(magenta_genes, file.path(outdir, "step23_megenta_genes.csv"), row.names = FALSE)
```

	找hub gene

```R
geneInfo <- read.csv("step20_gene_GS_MM_table.csv")

hub_magenta <- geneInfo[
  geneInfo$Module == "magenta" &
  abs(geneInfo$GS_Treatment) > 0.6 &
  abs(geneInfo$MEmagenta) > 0.8
]

head(hub_magenta)
```

	与注释文件的基因信息合并

```R
#================读取注释文件=========
anno <- read.csv("G:/thermal_wqz/WGCNA/cgl_anno.csv",
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE)

head(anno)
colnames(anno)

#===========和hub gene合并=================

hub_annot <- merge(hub_magenta,
                   anno,
                   by.x = "Gene",
                   by.y = "cgl-id",
                   all.x = TRUE)

head(hub_annot)

write.csv(hub_annot,
          "G:/thermal_wqz/WGCNA/step25_hub_magenta_annotated.csv",
          row.names = FALSE)

#==========保留更重要的列（干净版）===========

hub_annot_small <- hub_annot[, c("Gene", "gene-name", "product", "Module", "GS_Treatment", "MEmagenta")]

write.csv(hub_annot_small,
          "G:/thermal_wqz/WGCNA/step25_hub_magenta_annotated_simple.csv",
          row.names = FALSE)

head(hub_annot_small)

```

step25_hub_magenta_annotated.csv

![[Pasted image 20260324153721.png]]

step25_hub_magenta_annotated_simple.csv
![[Pasted image 20260324153816.png]]