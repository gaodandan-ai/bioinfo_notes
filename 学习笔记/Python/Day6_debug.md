## 核心原则

**报错不是问题，它是在“告诉你哪里错了”
	环境问题
	路径问题
	参数问题

### 第一步：只看最后1~3行

#### 1.环境问题

```bash
ModuleNotFoundError: No module named 'pandas'            
```

解决方案

```python
pip install pandas
```

#### 2.路径问题

```python
- `No such file or directory`
- `can't open file`
```

解决方案

```bash
pwd        # 看当前在哪
ls         # 看文件在不在
tree -L 2  # 看结构
```


#### 3.参数问题

```bash
the following arguments are required
--rawdata_dir ...
```

解决方案

```bash
python xxx.py --help
```