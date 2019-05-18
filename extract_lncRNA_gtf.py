
#!/usr/bin/python
#-*-coding : utf-8-*-
#Copyright(c) 2019 - leimingJiang <leiming8886@163.com>
import optparse
import re
import os
######################### define input and output######################################
"""
parse=optparse.OptionParser()
parse.add_option('-i','--input',dest='input',action='store',metavar='input files',help='enter your transcript (contain the transcript ID)')
parse.add_option('-o','--out',dest='outfile',action='store',metavar='output files',help='assign your output file')
parse.add_option('-g','--gtf',dest='gtf',action='store',metavar='gtf file name',help='please enter your gtf files')

(options,args) = parse.parse_args()
inPutFileName = options.input
outPutFileName = options.outfile
FileType = options.gtf
"""

def extract_lncRNA_gtf(input, outfile, gtf):
	inPutFileName = input
	outPutFileName = outfile
	FileType = gtf
	try:
		if not os.path.exists(inPutFileName):
			print("Error: "+inPutFileName+" doesn't exist!")
	except:
		print("Error: "+inPutFileName+" doesn't exist!")
	try:
		if not outPutFileName:
			print("outfile doesn't exist!")
	except:
		print("outfile doesn't exist!")
	#set lncRNA transcript ID;
	fr_lnc= open(inPutFileName)
	lncRNA_set={}
	for line in fr_lnc.readlines():
		line.strip()
		arr=line.split('\t')
		if arr[0] == "transcript ID":
			continue
		lncRNA_set[arr[0]]=arr[1]
	fr_lnc.close()
	#open GTF
	fr_gtf = open(FileType)
	fr_out = open(outPutFileName+'.gtf','w')
	for line in fr_gtf.readlines():
		line.strip()
		arr=line.split('\t')
		array=arr[-1].split(';')
		del array[-1]
		for line1 in array:
			line1=line1.strip()
			line1=line1.replace('\"','')
			line1 = line1.strip()
			arr1=line1.split(' ')
			if arr1[0]=="transcript_id" and lncRNA_set.get(arr1[1]):
				fr_out.write(line)
	fr_gtf.close()
	fr_out.close()




