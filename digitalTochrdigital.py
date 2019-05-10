#!/usr/bin/python
#-*-coding : utf-8-*-
#Copyright(c) 2019 - leimingJiang <leiming8886@163.com>
import os
import optparse
#FileName = raw_input('Please enter your a file name: ')
parse=optparse.OptionParser()
#directory
parse.add_option('-i','--inputdir',dest='filedir',action='store',metavar='directory of input gtf',help='enter your gtfs')
parse.add_option('-o','--outdir',dest='outdir',action='store',metavar='',help='please enter directory')
(options,args) = parse.parse_args()
inPutFileNames = options.filedir
newdir = options.outdir
#定义集合表示染色体1-23，和X，Y
chr=set('chr'+str(x) for x in list(range(1,24)))
chr.add('chrX')
chr.add('chrY')
digchr=set(str(x) for x in list(range(1,24)))
digchr.add('X')
digchr.add('Y')
if not os.path.exists(newdir):
    os.makedirs(newdir)
'''def digi_To_chrdigi(*path):
	for subpath in path:
		fr = open(subpath,'r')
		f = open(newdir+'/'+subpath, 'w')
		for line in fr.readlines():
			line=line.strip()
			array=line.split('\t')
			firs=array.pop(0)
			if firs not in chr:
				continue
			firs='chr'+firs
			array.insert(0,firs)
			f.write('\t'.join(array)+"\n")
	f.close()
	fr.close()
'''
#定义函数设置过滤后的文件类型 当然可以设置多个类型，后缀.gtf，返回文件名
#
def all_path_suffix(dirname,filter):
	result = []#所有的文件
	for maindir, subdir, file_name_list in os.walk(dirname):
		# print("1:",maindir) #当前主目录
        # print("2:",subdir) #当前主目录下的所有目录
        # print("3:",file_name_list)  #当前主目录下的所有文件
		for filename in file_name_list:
			#apath = os.path.join(maindir, filename)#合并成一个完整路径
			ext = os.path.splitext(filename)[1]# 获取文件后缀 [0]获取的是除了文件名以外的内容
			if ext == filter:
				result.append(filename)
	return result
filenames=[]
filenames=all_path_suffix(inPutFileNames,'.gtf')
print(filenames)
path=filenames
#digi_To_chrdigi(filenames)
for subpath in path:
	#subpath='K140_3rd.gtf' subpath
	fr=open(subpath,'r')
	f = open(newdir+'/'+subpath, 'w')
	for line in fr.readlines():
		line=line.strip()
		array=line.split('\t')
		firs=array.pop(0)
		if firs in chr:
			array.insert(0,firs)
			f.write('\t'.join(array)+"\n")
		if firs in digchr:
			firs='chr'+firs
			array.insert(0,firs)
			f.write('\t'.join(array)+"\n")
f.close()
fr.close()
