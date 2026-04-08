将两个文件的cgl-id 和Geneid 比对到一起，最后保留基因id，log2FC，P.Value.gene-name,product,生成新表

```bash
cd f:/thermal_wqz/multi_timepoint_analysis/output_24h/24h_DEG_product_match
```

第一种方法：bash
```bash
awk '
BEGIN{FS=OFS=","}
FNR==NR{
    if (NR>1) {
        gsub(/"/, "", $1)
        gsub(/"/, "", $3)
        gsub(/"/, "", $8)
        gene[$1]=$3
        product[$1]=$8
    }
    next
}
FNR==1{
    print "gene_id","log2FC","P.Value","gene_name","product"
    next
}
{
    gsub(/"/, "", $1)
    gsub(/"/, "", $2)
    gsub(/"/, "", $5)

    if ($5 in gene) {
        print $5,$1,$2,gene[$5],product[$5]
    }
}
' /f/thermal_wqz/multi_timepoint_analysis/output_24h/24h_DEG_product_match/cgl_anno.csv \
  /f/thermal_wqz/multi_timepoint_analysis/output_24h/24h_DEG_product_match/24h_significant_DEGs.csv \
  > /f/thermal_wqz/multi_timepoint_analysis/output_24h/24h_DEG_product_match/24h_DEG_annotated.csv

```

第二种方法：R语言
```R
library(readr)
library(dplyr)

# 读取文件
anno <- read_csv("F:/thermal_wqz/multi_timepoint_analysis/output_1h/1h_DEG_product_match/cgl_anno.csv")
deg  <- read_csv("F:/thermal_wqz/multi_timepoint_analysis/output_1h/1h_DEG_product_match/1h_significant_DEGs.csv")

# 查看列名（可选）
# colnames(anno)
# colnames(deg)

# 合并并筛选列
result <- deg %>%
  left_join(anno, by = c("Geneid" = "cgl-id")) %>%
  select(
    gene_id = Geneid,
    log2FC = log2FC_Treat1h_vs_Control1h,
    P.Value,
    gene_name = `gene-name`,
    product
  )

# 写出文件
write_csv(result, "F:/thermal_wqz/multi_timepoint_analysis/output_1h/1h_DEG_product_match/DEG_annotated.csv")
```

24h差异基因过多
再筛选log2FC>2的基因生成新表

```bash
awk -F',' '
NR==1 {print; next}
{
    gsub(/"/, "", $2)
    if ($2+0 > 2 || $2+0 < -2)
        print
}
' /f/thermal_wqz/multi_timepoint_analysis/output_24h/24h_DEG_product_match/24h_DEG_annotated.csv \
> /f/thermal_wqz/multi_timepoint_analysis/output_24h/24h_DEG_product_match/24h_DEG_annotated_log2FC_gt2.csv
```