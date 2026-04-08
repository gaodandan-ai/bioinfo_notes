

# 一、学习目标

- less
- wc
- find
- 快速查看大文件
---

### less

基础语法
```bash
less file.txt
```
![[Pasted image 20260316211801.png]]

在文件中搜索
```
/关键词
n   #跳到下一个匹配

q返回
```

### wc  统计文件内容（行数  单词数   字符数）

Word count

```bash
wc file.txt
```

显示
```
行数  单词数 字符数  文件名
```

比如
```
PC@Komorebi MINGW64 /g/cg_sequencing1
$ wc beizhu.txt
  321   863 15234 beizhu.txt
```

最常用参数
```bash
wc -l file.txt   #统计行数
```

```
$ wc -l beizhu.txt
321 beizhu.txt
```


### find   查找文件

基础语法

```bash
find 目录 -name 文件名
```

示例
```bash
find . -name beizhu.txt   #. 表示在当前目录  -name 按照名字搜索
./beizhu.txt
./script/beizhu.txt
./02_scripts/beizhu.txt


find . -name "*.txt"    #查找所有的txt文件
```

# 总结

Linux处理文件

	cat 直接查看文件
	head 查看前10行
	tail 查看后10行
	less 分页查看
	wc 统计行数

```bash
less varient.vcf
grep -v "^#" varient.vcf | wc -l
#查看vcf文件，统计突变数量
```