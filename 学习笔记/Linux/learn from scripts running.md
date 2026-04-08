为脚本设置可执行权限，并执行脚本
```
```bash
chmod +x test.sh
./text.sh
```

## 常见的技巧

写脚本时，开头最好加上

```bash
set -euo pipefail
```
	-e 有命令报错就退出
	-u 用了未定义的变量就报错
	pipefail:管道中出错能检测到

打印日志

```bash
echo "[INFO]开始处理 $sample"
```

跳过已经完成的样本,这样脚本中断后可以续跑

```bash
if [ -f "$bam" ]; then
	echo "$sample 已完成，跳过"
	continue
fi
```

nohup: no hang up(不挂断)
```bash
nohup
```
查看是否在跑
```
ps -ef | grep crispri_dual
```

结束进程

```bash
kill 进程号
```

查看日志
```
tail -f run.log
```

实例：
```
nohup python scripts/dual_gene/crispri_dual_gene_screening_pipeline.py \
  --rawdata-dir /data/gaodd/CRISPRi_thermal/raw_data/dual_gene \
  --reference-fasta /data/gaodd/CRISPRi_thermal/crispri_screening_pipeline/data/library3_slt.fasta \
  --anno-file /data/gaodd/CRISPRi_thermal/crispri_screening_pipeline/data/cgl-anno.csv \
  --output-dir /data/gaodd/CRISPRi_thermal/dual_gene_results \
  --threads 8 \
  --repeat-max-errors 3 \
> run.log 2>&1 &
```
nohup：忽略断线信号
> run.log  把输出写到日志
2>&1  错误信息也写到日志
&  后台运行


## 上调或者下调的基因写成新表

```bash
awk -F ',' 'NR==1 || $2>0' G:/thermal_protemic/24h/24h_DEG.csv > G:/thermal_protemic/24h/24h_up.csv
```
- `-F ','`：按逗号分割 CSV
- `NR==1`：保留表头
- `$2>0`：第 2 列是 logFC，大于 0 → 上调
- `$2<0`：第 2 列小于 0 → 下调
- `> 24h_up.csv`：输出新文件

```bash
awk -F ',' 'BEGIN{OFS=","} NR==1 || $2<0 {print $1,$2,$3,$10,$15}' G:/thermal_protemic/24h/24h_DEG.csv > G:/thermal_protemic/24h/24h_down_selected.csv
```
- `-F ','` 按逗号分隔
- `OFS=","` 输出也是逗号分隔
- `print $1,$2,$3,$10,$15` **只输出这 5 列**
- `NR==1` 保留表头
- `$2>0` 只留上调
- `$2<0` 只留下调


## `>>`：追加，`>`：覆盖