需要两个文件
	G:/thermal_protemic/mapping_cgl_simple.csv
	G:/thermal_protemic/cgl_anno.csv

```R
# ====================== 1. 读取文件 ======================
map <- read.csv("G:/thermal_protemic/mapping_cgl_simple.csv", stringsAsFactors = FALSE)
deg1h <- read.csv("G:/thermal_protemic/1h/1h_DEG.csv", stringsAsFactors = FALSE, row.names = 1)
deg4h <- read.csv("G:/thermal_protemic/4h/4h_DEG.csv", stringsAsFactors = FALSE, row.names = 1)
deg24h <- read.csv("G:/thermal_protemic/24h/24h_DEG.csv", stringsAsFactors = FALSE, row.names = 1)

# ==================
==== 2. 给mapping文件改列名（适配你的格式） ======================
colnames(map) <- c("Protein", "Gene")

# ====================== 3. 批量匹配函数（蛋白 → 基因） ======================
map_deg <- function(deg_df, map_df){
  # 把蛋白ID变成一列
  deg_df$Protein <- rownames(deg_df)
  # 匹配基因
  deg_mapped <- merge(
    x = map_df,
    y = deg_df,
    by = "Protein",
    all.y = TRUE  # 保留所有DEG，只匹配有对应基因的
  )
  # 去掉没匹配到的
  deg_mapped <- deg_mapped[!is.na(deg_mapped$Gene), ]
  return(deg_mapped)
}

# ====================== 4. 批量处理三个时间点 ======================
deg1h_mapped <- map_deg(deg1h, map)
deg4h_mapped <- map_deg(deg4h, map)
deg24h_mapped <- map_deg(deg24h, map)

# ====================== 5. 输出结果CSV ======================
write.csv(deg1h_mapped, "G:/thermal_protemic/1h/1h_DEG_mapped.csv", row.names = FALSE, quote = FALSE)
write.csv(deg4h_mapped, "G:/thermal_protemic/4h/4h_DEG_mapped.csv", row.names = FALSE, quote = FALSE)
write.csv(deg24h_mapped, "G:/thermal_protemic/24h/24h_DEG_mapped.csv", row.names = FALSE, quote = FALSE)

# ====================== 6. 输出统计信息（方便你看） ======================
cat("===== 匹配完成 =====\n")
cat("1h 原始DEG数量：", nrow(deg1h), "  |  匹配到基因：", nrow(deg1h_mapped), "\n")
cat("4h 原始DEG数量：", nrow(deg4h), "  |  匹配到基因：", nrow(deg4h_mapped), "\n")
cat("24h 原始DEG数量：", nrow(deg24h), "  |  匹配到基因：", nrow(deg24h_mapped), "\n")
```