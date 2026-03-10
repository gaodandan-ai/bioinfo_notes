# Python 文件读取

## 学习目标

学会用 Python 读取文本文件。

---

## 代码

```python
with open("genes.txt") as f:
    for line in f:
        print(line.strip())
```

---

## 代码解释

```
open() 打开文件
for line in f 逐行读取
strip() 去掉换行符
```

---

## 练习

```python
with open("test.txt") as f:
    for line in f:
        print(line)
```

---

## 易错点

```
Python 缩进必须一致
文件路径要正确
```