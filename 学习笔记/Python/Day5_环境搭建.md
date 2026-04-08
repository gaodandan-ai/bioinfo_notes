# 学习内容

	1.Anaconda/Miniconda
	2.基于SSH的远程开发

## Anaconda/Miniconda

自带丰富库，适合生物信息学分析
	Anaconda含有大量预装库，体量较大
	miniconda 按需安装

### 定期更新
```python
conda update -n base -c defaults conda
```

### 创建独立环境
```python
conda create -n biofino_py3 python=3.9
```

### 激活环境
```python
conda activate biofino_py3
```

# 思维模型

**数据 → 脚本 → pipeline → 项目**

### 一个完整的项目包含
	rawdata  原始数据
	scripts     逻辑
	results     输出
	run.sh      流程
	README  说明书
# README=说明书

### 介绍用途
```
This pipeline performs RNA-seq analysis...
```

### 文件结构
```
my_project/
├── rawdata/
├── scripts/
├── results/
├── README.md
└── run.sh
```

### 怎么运行
```python
bash run.sh
```

### 每一步在做什么
```
step1:
step2:
...
```

### 依赖环境
```
python 3.8
R 4.2
```

## 基于SSH的远程开发

	远程登录服务器
	远程执行命令
	安全传输数据


工作模式
```
本地电脑（Trae）
        ↓ SSH
远程服务器（Linux + Python + 数据）
        ↓
运行代码 / 处理数据 / 训练模型
```

服务器信息
```
- IP：`172.16.2.114`
- 端口：`22`
- 用户名：`gaodd`
```

### 第一步：配置SSH

打开文件夹
```
C:\Users\你的用户名\.ssh\config
```

加一段
```
Host lab-server
    HostName 172.16.2.114
    User gaodd
    Port 22
```

### 第二步：在trae里边连接

	1.安装插件remote-SSH
	2.ctrl+shift+P
	3.输入Remote-SSH: Connect to Host
	4.选择lab-server