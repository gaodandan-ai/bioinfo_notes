## 什么是COBRApy

**COBRApy** 是一个用 Python 做 **代谢网络建模（metabolic modeling）** 的工具

	WGCNA：找“哪些基因一起变化”  
	COBRApy：算“这些基因/代谢通路到底怎么影响生长和代谢”

### 核心思想

	Flux Balance Analysis（FBA）：在所有约束条件下，计算“最优代谢流分配”

## 能做什么

#### 预测生长

```python
model.optimize()
```

#### 基因敲除分析

```python
model.genes.get_by_id("geneX").knock_out()
model.optimize()
```
