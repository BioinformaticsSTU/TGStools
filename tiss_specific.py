#!/usr/bin/python
#-*-coding : utf-8-*-
#Copyright(c) 2019 - leimingJiang <leiming8886@163.com>
import optparse
import pandas as pd
import re
import os
import datetime

'''
start_time = datetime.datetime.now()
#FileName = raw_input('Please enter your a file name: ')
parse=optparse.OptionParser()
parse.add_option('-i','--input',dest='files',action='store',metavar='input files',help='enter your transcript (contain the transcript ID)')
#,nargs="*"
parse.add_option('-o','--out',dest='outfile',action='store',metavar='output files',help='assign your output file')
parse.add_option('-t','--tissue',dest='tss',action='store',metavar='tissue name',help='please enter tissue name')
parse.add_option('-r','--reference',dest='hg',action='store',metavar='ref name',help='please enter hg38 or hg19')

(options,args) = parse.parse_args()
outPutFileName = options.outfile+'.gtf'
tss_name = options.tss
refgene = options.hg
'''

def tiss_specific(outfile, tss, hg, files):
	start_time = datetime.datetime.now()
	outPutFileName = outfile+'.gtf';
	tss_name = tss
	refgene = hg
	PATH=os.path.split(os.path.realpath(__file__))[0]
	inPutFileNames = list(files.split(","))
	df = pd.read_csv(PATH+"/"+"expression.txt",sep='\t')
	col_name_del_probes_De=list(df.columns)
	col_name_del_probes_De.remove("Probes")
	col_name_del_probes_De.remove("Description")

	#set dic key = tiss_name. value = transcript IDs
	dict_tiss={}
	# print(list(df.loc[df['adipose'] != 0,'Probes']))
	#dict_tiss dict key tiss,value probes
	for col_line in col_name_del_probes_De:
		col_line_Dedu=col_line.split('_')[0]
		dict_tiss[col_line_Dedu] = list(df.loc[df[col_line] != 0,'Probes'])
	all_pro=set()
	for x in dict_tiss:
		all_pro.update(dict_tiss[x])
	if refgene == 'hg38':
		ref_gtf=open(PATH+"/"+'lncRNA_hg38.gtf','r')
		print("hg38")
	if refgene == 'hg19':
		ref_gtf=open(PATH+"/"+'lncRNA_hg19.gtf','r')
		print("hg19")
	fr=ref_gtf.readlines()
	#set_lncRNA key probe,value temp_list[[chr,start,end,strand],[chr,start,end,strand]]
	set_lncRNA={}
	for line in fr:
		line=line.strip()
		arr_line=line.split('\t')
		chr=arr_line[0]
		start=int(arr_line[3])
		end=int(arr_line[4])
		strand=arr_line[6]
		temp_list=[chr,start,end,strand]
		last_arr_tem=arr_line[-1].split(';')[0]
		last_arr_tem=last_arr_tem.replace('\"','')
		last_arr_id=last_arr_tem.split(' ')[1]
		if last_arr_id not in set_lncRNA:
			set_lncRNA[last_arr_id]=[]
			set_lncRNA[last_arr_id].append(temp_list)
		else:
			if temp_list in set_lncRNA[last_arr_id]:
				continue
			else:
				set_lncRNA[last_arr_id].append(temp_list)
	ref_gtf.close()


	#open gtf to align
	#the first five line of gtf_union.gtf is test,3 overlop,2 return
	if len(inPutFileNames)==1:
		sam_gtf = open(inPutFileNames[0], 'r')
		sam_gtf_out = open(outPutFileName, 'w')
		if tss_name in dict_tiss:
			pick = tss_name
			Pick_pro = set(dict_tiss[pick])
			lncRNA_left_pro = set()
			for col_line in dict_tiss:
				if col_line == pick:
					continue
				for probe in dict_tiss[col_line]:
					lncRNA_left_pro.add(probe)
			for line in sam_gtf.readlines():
				line=line.strip()
				arr_line=line.split('\t')
				chr= arr_line[0]
				start=int(arr_line[3])
				end=int(arr_line[4])
				strand=arr_line[6]
				num_temp_list=0
				for key_id in all_pro:
					if key_id not in set_lncRNA:
						continue
					for list1 in set_lncRNA[key_id]:
						if list1[-1] == strand and list1[0] == chr:
							if (start - list1[1]) * (list1[2] - start) >= 0 or (end - list1[1]) * (list1[2] - end) >= 0:
								if key_id in Pick_pro:
									num_temp_list += 1
				if num_temp_list>0:
					sam_gtf_out.write(line + '\n')
		else:
			for line in sam_gtf.readlines():
				line=line.strip()
				arr_line=line.split('\t')
				chr= arr_line[0]
				start=int(arr_line[3])
				end=int(arr_line[4])
				strand=arr_line[6]
				num_temp_list=0
				for key_id in set_lncRNA:
					for list1 in set_lncRNA[key_id]:
						if list1[-1]==strand and list1[0] == chr:
							if (start- list1[1])*(list1[2]-start)>=0 or (end- list1[1])*(list1[2]-end)>=0:
								num_temp_list+=1
				if num_temp_list>0:
					continue
				else:
					sam_gtf_out.write(line+'\n')

		sam_gtf.close()
		sam_gtf_out.close()
	if len(inPutFileNames)==2:
		sam_gtf1 = open(inPutFileNames[0], 'r')
		sam_gtf2 = open(inPutFileNames[1], 'r')
		set_gtf2 = set()
		set_gtf1 = set()
		for line in sam_gtf1:
			line = line.strip()
			arr_line = line.split('\t')
			chr = arr_line[0]
			start = int(arr_line[3])
			end = int(arr_line[4])
			strand = arr_line[6]
			temp_list = [chr, start, end, strand]
			last_arr_tem = arr_line[-1].split(';')[0]
			last_arr_tem = last_arr_tem.replace('\"', '')
			last_arr_id = last_arr_tem.split(' ')[1]
			set_gtf1.add(last_arr_id)
		sam_gtf1.close()
		for line in sam_gtf2 :
			line = line.strip()
			arr_line = line.split('\t')
			chr = arr_line[0]
			start = int(arr_line[3])
			end = int(arr_line[4])
			strand = arr_line[6]
			temp_list = [chr, start, end, strand]
			last_arr_tem = arr_line[-1].split(';')[0]
			last_arr_tem = last_arr_tem.replace('\"', '')
			last_arr_id = last_arr_tem.split(' ')[1]
			set_gtf2.add(last_arr_id)
		sam_gtf2.close()
		sam = open(inPutFileNames[1], 'r')
		sam_gtf_out = open(outPutFileName, 'w')
		set_left=set_gtf2-set_gtf1
		for line in sam :
			line = line.strip()
			arr_line = line.split('\t')
			chr = arr_line[0]
			start = int(arr_line[3])
			end = int(arr_line[4])
			strand = arr_line[6]
			temp_list = [chr, start, end, strand]
			last_arr_tem = arr_line[-1].split(';')[0]
			last_arr_tem = last_arr_tem.replace('\"', '')
			last_arr_id = last_arr_tem.split(' ')[1]
			if last_arr_id in set_left:
				sam_gtf_out.write(line + '\n')
		sam.close()
		sam_gtf_out.close()
	end_time = datetime.datetime.now()
	print(end_time - start_time," second")


