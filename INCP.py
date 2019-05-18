#!/usr/bin/python
#-*-coding : utf-8-*-
#Copyright(c) 2019 - leimingJiang <leiming8886@163.com>
#pip install matplotlib_venn
#pip install matplotlib
#pip install numpy
import optparse
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles,venn2,venn2_circles
######### define input and output ########

"""
#FileName = raw_input('Please enter your a file name: ')
parse=optparse.OptionParser()
parse.add_option('-i','--input',dest='file',action='store',metavar='input files',help='enter your transcript (sequence or gtf)')
parse.add_option('-p','--parallel',dest='parallel',action='store',metavar='prallel numbers',help='please enter your specified speed ratio')
parse.add_option('-g','--gtf',dest='gtf',action='store_true',metavar='gtf file name',help='please enter your gtf files')
parse.add_option('-r','--reference',dest='reference',action='store',metavar='',help='if your input file is gtf type please enter RefGenome directory')

(options,args) = parse.parse_args()
inPutFileName = options.file
Parallel = options.parallel
FileType = options.gtf
Directory = options.reference
"""

def sub_array(A,B):
    x=set(A)
    y=set(B)
    return list(x - y)

def intersect_array(A,B):
    x=set(A)
    y=set(B)
    return list(x & y)

def union_array(A,B):
    x=set(A)
    y=set(B)
    return list(x | y)

def TwoLineFasta (Seq_Array):
    Tmp_sequence_Arr = []
    Tmp_trans_str = ''
    for i in range(len(Seq_Array)):
        Seq_Array[i]=Seq_Array[i].strip()
        if '>' in Seq_Array[i]:
            if i == 0:
                Tmp_sequence_Arr.append(Seq_Array[i])
            else:
                Tmp_sequence_Arr.append(Tmp_trans_str)
                Tmp_sequence_Arr.append(Seq_Array[i])
                Tmp_trans_str = ''
        else:
            if i == len(Seq_Array) - 1:
                Tmp_trans_str = Tmp_trans_str + str(Seq_Array[i])
                Tmp_sequence_Arr.append(Tmp_trans_str)
            else:
                Tmp_trans_str = Tmp_trans_str + str(Seq_Array[i])
    return Tmp_sequence_Arr

def noTwolineFasta (Seq_Array,wid=50):
    Tmp_sequence_Arr = []
    width=int(wid)
    for i in range(len(Seq_Array)):
        Seq_Array[i]=Seq_Array[i].strip()
        len_fa=len(Seq_Array[i])
        start = 0
        end = start + width
        if '>' in Seq_Array[i]:
            Tmp_sequence_Arr.append(Seq_Array[i])
        else:
            while end < len_fa:
            #output=substr(seq_record.seq, i, 50)
                Tmp_sequence_Arr.append(str(Seq_Array[i][start:end]))
                start = start + width
                end += width
            if end >= len_fa:
                Tmp_sequence_Arr.append(str(Seq_Array[i][start:]))
    return Tmp_sequence_Arr

############
def INCP(file, parallel, gtf, reference):
	inPutFileName = file
	Parallel = parallel
	FileType = gtf
	Directory = reference

	PATH=os.path.split(os.path.realpath(__file__))[0]
	CNCIPATH=PATH+'/CNCI-master'
	PLEK=PATH+'/PLEK.1.2'
	outPutFileName=os.path.splitext(inPutFileName)[0]
	try:
		if not os.path.exists(inPutFileName):
			print("Error: "+inPutFileName+" doesn't exist!")
	except:
		print("Error: "+inPutFileName+" doesn't exist!")

	#CNCI and PLEK code
	if FileType:
		if not os.path.exists(Directory):
			print("please enter RefGenome directory of 2bit")
		os.system('python ' + CNCIPATH + '/CNCI.py -f '+inPutFileName+' -g -o '+outPutFileName+' -m ve -p '+Parallel+' -d ' +Directory)
	#fasta is not TwoLineFasta, fastaFiles = inPutFileName + '.fa', so need to convert format TwoLineFasta
		fastaFiles = outPutFileName + '.gtf.fa'
		fastaFiles_twoline=outPutFileName +'_plek'+'.fa'
		GtfInFiles = open(fastaFiles)
		inFilesArr = GtfInFiles.read()
		sequence_Arr = inFilesArr.split('\n')
		sLen = len(sequence_Arr) - 1#the last row is the null due to split (\n)
		del sequence_Arr[sLen]
		ARRAY =  TwoLineFasta(sequence_Arr)
		fr = open(fastaFiles_twoline,'w')
		for line in ARRAY:
			fr.write(line+"\n")
		fr.close()
		os.system('python ' + PLEK + '/PLEK.py '+' -fasta '+fastaFiles_twoline+' -out '+outPutFileName+'_PLEK'+' -thread 10 ')
	else:
		os.system('python ' + PLEK + '/PLEK.py '+' -fasta '+inPutFileName+' -out '+outPutFileName+'_PLEK'+' -thread 10 ')
		fastaFiles_notwoline=outPutFileName+"_notwoline"+'.fa'
		GtfInFiles = open(inPutFileName)
		inFilesArr = GtfInFiles.read()
		sequence_Arr = inFilesArr.split('\n')
		sLen = len(sequence_Arr) - 1#the last row is the null due to split (\n)
		del sequence_Arr[sLen]
		ARRAY =  noTwolineFasta(sequence_Arr,50)
		fr = open(fastaFiles_notwoline,'w')
		for line in ARRAY:
			fr.write(line+"\n")
		fr.close()
		os.system('python ' + CNCIPATH + '/CNCI.py -f '+fastaFiles_notwoline+' -o '+outPutFileName+' -m ve -p '+Parallel)

	#set hash of two output of the out_plek and out_cnci
	#[key for key in d]
	out_plek=outPutFileName+'_PLEK'
	out_cnci=outPutFileName+'/CNCI.index'
	#set hash of two output of the out_plek and out_cnci
	out_set_plek={}
	out_set_cnci={}
	pl_fr = open(out_plek)
	cn_fr = open(out_cnci)
	for line in  pl_fr.readlines():
		line=line.replace('>','')
		line1=line.strip()
		line_c=line1.split('\t')
		if line_c[0] == 'Coding':
			continue
		out_set_plek[line_c[2]]=line_c[0]+"\t"+line_c[1]
	gi_p=[key for key in out_set_plek]

	for line in  cn_fr.readlines():
		line1=line.strip()
		line_c=line1.split('\t')
		if line_c[1] ==  'coding' or line_c[1] ==  'index':
			continue
		out_set_cnci[line_c[0]]=line_c[1]+"\t"+line_c[2]
	gi_c=[key for key in out_set_cnci]

	#venny
	figure = plt.figure()
	venn2([set(gi_p), set(gi_c)], set_labels = ('PLEK', 'CNCI'),)
	plt.title("lncRNA predicted by PLEK and CNCI")
	figure.savefig('Venn_diagram_'+outPutFileName+'.pdf', bbox_inches='tight')
	plt.close()
	pl_fr.close()
	cn_fr.close()
	union=open(outPutFileName+"union_plek_cnci.txt",'w')
	inter=open(outPutFileName+"intersect_plek_cnci.txt",'w')
	union.write('transcript ID\tPLEK_index\tPLEK_score\tCNCI_index\tCNCI_score\n')
	inter.write('transcript ID\tPLEK_index\tPLEK_score\tCNCI_index\tCNCI_score\n')
	for key in union_array(gi_p,gi_c):
		if not out_set_cnci.get(key):
			out_set_cnci[key]="null\tnull"
		if not out_set_plek.get(key):
			out_set_plek[key]="null\tnull"
		union.write(key+'\t'+out_set_plek[key]+'\t'+out_set_cnci[key]+'\n')
	for key in intersect_array(gi_p,gi_c):
		inter.write(key+'\t'+out_set_plek[key]+'\t'+out_set_cnci[key]+'\n')
	union.close()
	inter.close()


