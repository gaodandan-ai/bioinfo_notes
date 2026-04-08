**加权基因共表达网络分析（Weighted Gene Co-expression Network Analysis, WGCNA）是一种强大的无监督系统生物学方法
## 本质/核心逻辑

将表达模式相似的基因聚在一块，找关键基因
样本数>=15

基因表达矩阵
```
gene  1h  4h  24h
geneA  10  20  30
geneB  12  22  44
```
表达模式相似的基因可能属于
	同一通路
	同一调控
	同一调控

## step1 用相关性得到相关性矩阵

```
cor(dnaK, groEL) ≈ 0.98（很像）
cor(dnaK, clpB) ≈ 0.2（不像）
```

## step2 变成网络

做一个强化连接的操作
aij​=∣cor(i,j)∣β

- cor 越大 → 关系越强
- β（soft threshold）：
    - 放大强关系
    - 压低弱关系

结果变成：基因网络

## step3 找基因群

每个模块=一组共表达基因
```
模块1（蓝色）：dnaK, groEL, hspR
模块2（黄色）：metA, metB, metC
```

## step4 找核心基因

在每个模块里找连接最多的基因