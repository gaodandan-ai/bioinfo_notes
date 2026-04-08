# 主题

	1.终端提示符
	2.基础命令
---

# 一、学习目标

- 学会使用 Linux 基础命令
- 理解 pwd、ls、cd、mkdir 的作用
- 能在终端创建文件夹和进入目录

---

# 二、知识点总结

## 终端提示符

```
PC            → 用户名
@xxxx         → 电脑名
MINGW64       → Git Bash 环境
/g/bioinfo_notes → 当前目录
$             → 等待输入命令
```

## 1. 命令 / 函数

### pwd  查看当前目录

```bash
pwd
```

示例输出：

```bash
/data/gaodd
```

---

### ls 查看当前目录文件


```bash
ls
```

查看详细信息：

```bash
ls -l
```

---

### cd 进入某个目录


```bash
cd bioinfo_learning
```

返回上级：

```bash
cd ..
```

返回home:
```bash
cd ~
```

### mkdir  创建文件夹

创建文件夹
```bash
mkdir filename
```

### rm  删除

删除文件
```bash
rm filename.txt
```

删除文件夹
```bash
rm -r foldername
```

### cp 复制

```
cp file.txt
```

复制到文件夹

```
cp file.txt /Data
```

### mv 移动或者重命名

移动文件
```
rm file.txt /Data
```

重命名
```
rm oldname.txt newname.txt
```

### clear  清屏

清空终端

```
clear
```

---

# 三、代码示例

```bash
mkdir bioinfo_learning
cd bioinfo_learning
mkdir python
mkdir r
mkdir data
ls
```

---

# 四、我的理解

用自己的话写：

```
pwd 就是查看我现在在哪个目录
ls 就是查看当前目录有什么文件
cd 是进入目录
mkdir 是创建文件夹
```

---

# 五、实践练习

今天实际做过的操作：

```bash
mkdir test_linux
cd test_linux
touch test.txt
ls
```

---

# 六、易错点

```
cd 后面必须是目录
rm 删除文件后不能恢复
Linux 区分大小写
```

---

# 七、扩展知识

例如：

```
ls -lh 可以更容易看文件大小
```

---

# 八、以后可能用到的场景

例如：

```
生信分析时经常需要在服务器创建目录
存放 fastq / bam / vcf 数据
```