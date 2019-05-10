
import os
import re
import matplotlib.pyplot as plt

"""
path="histone";
gtf="../TEST/K140.gtf";
"""
# flag is useless? fantom5 and histone file are all bed format
def staDist(gtf, path, flag, prefix):
	def paintTSS(pdf_path,file_path):                                
		name_list = ['<50', '<100', '<500', '<1000', '>1000']
		with open(file_path, 'r') as f1:
			next(f1)
			list1 = f1.readlines()
			line = "".join(list1).split('\t')
			less50 = int(line[0])
			less100 = int(line[1])
			less500 = int(line[2])
			less1000 = int(line[3])
			more1000 = int(line[4].replace('\n', ''))
			num = less50 + less100 + less500 + less1000 + more1000
		num_list = [less50 / num, less100 / num, less500 / num, less1000 / num, more1000 / num]
		fig1 = plt.figure(1)
		plt.bar(range(len(num_list)), num_list, color=["#D62728", "#ECE812", "#9467BD", "#E377C2", "#2CA02C"], alpha=0.5, tick_label=name_list)
		fig1.savefig(pdf_path, dpi=150)
		plt.close()
	#######
	# get histone data
	chromosome1={
	"chr1":[],"chr2":[],"chr3":[],"chr4":[],"chr5":[],"chr6":[],"chr7":[],"chr8":[],"chr9":[],"chr10":[],"chr11":[],"chr12":[],"chr13":[],"chr14":[],"chr15":[],"chr16":[],"chr17":[],"chr18":[],"chr19":[],"chr20":[],"chr21":[],"chr22":[],"chrX":[],"chrY":[],"chrM":[]
	};
	chromosome2={
	"chr1":[],"chr2":[],"chr3":[],"chr4":[],"chr5":[],"chr6":[],"chr7":[],"chr8":[],"chr9":[],"chr10":[],"chr11":[],"chr12":[],"chr13":[],"chr14":[],"chr15":[],"chr16":[],"chr17":[],"chr18":[],"chr19":[],"chr20":[],"chr21":[],"chr22":[],"chrX":[],"chrY":[],"chrM":[]
	};
	chromosome=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"];
	INDEX={"-":chromosome1, "+":chromosome2};
	FILES=os.listdir(path);
	###########
	if(flag=='histone'):
		for FILE in FILES:
			fp=open(path+"/"+FILE, "r");
			line=fp.readline();
			while(line!=""):
				words=line.strip().split("\t");
				chr=words[0];
				strand=words[5];
				if(strand=="+"):
					point=int(words[1]);
				else:
					point=int(words[2]);
				if(chr not in chromosome):
					line=fp.readline();
					continue;
				INDEX[strand][chr].append([chr, strand, point]);
				line=fp.readline();
			fp.close();
	else:
		for FILE in FILES:
			fp=open(path+"/"+FILE, "r");
			line=fp.readline();
			while(line!=""):
				words=line.strip().split("\t");
				chr=words[0];
				strand=words[5];
				if(strand=="+"):
					point=int(words[1]);
				else:
					point=int(words[2]);
				if(chr not in chromosome):
					line=fp.readline();
					continue;
				INDEX[strand][chr].append([chr, strand, point]);
				line=fp.readline();
			fp.close();
	###########
	# get gtf data
	TRANS=[];
	file=gtf;
	fp=open(file, "r");
	line=fp.readline();
	while(line!=""):
		words=line.strip().split("\t");
		chr='chr'+words[0];
		start=words[3];
		end=words[4];
		strand=words[6];
		transcript_id=words[8].split(";")[1].strip().split(" ")[1].strip('"');
		if(strand=="+"):
			point=int(start);
		else:
			point=int(end);
		TRANS.append([chr, point, strand, transcript_id]);
		line=fp.readline();
	fp.close();
	trans_id=[];
	TSS={};
	for trans in TRANS:
		if(trans[3] not in trans_id):
			TSS.update({trans[3]:trans[0:3]});
			trans_id.append(trans[3]);
		else:
			if(trans[2]=="+"):
				if(trans[1]<TSS[trans[3]][1]):
					TSS[trans[3]][1]=trans[1];
			else:
				if(trans[1]>TSS[trans[3]][1]):
					TSS[trans[3]][1]=trans[1];
	#################
	# sta
	RESULT=[0, 0, 0, 0, 0];
	fp=open(prefix+"_DIST.txt", "w");
	fp.write('<=50'+'\t<=100'+'\t<=500'+'\t<=1000'+'\t>1000\n');
	for id in trans_id:
		trans=TSS[id];
		chr=trans[0];
		strand=trans[2];
		point=trans[1];
		name=id;
		score=[0, 0, 0, 0, 0];
		for index in INDEX[strand][chr]:
			if(abs(index[2]-point) < 1000):
				score[3]+=1;
				if(abs(index[2]-point) < 500):
					score[2]+=1;
					if(abs(index[2]-point) < 100):
						score[1]+=1;
						if(abs(index[2]-point) < 50):
							score[0]+=1;
							break;
		for i in range(4):
			if(score[i]>0):
				RESULT[i]+=1;
				break;
		else:
			RESULT[4]+=1;
	STRING="";
	mark=0;
	for i in RESULT:
		if(mark==0):
			mark=1;
			STRING+=str(i);
			continue;
		STRING+="\t"+str(i);
	fp.write(STRING+"\n");
	fp.close()
	paintTSS(prefix+"_DIST.pdf", prefix+"_DIST.txt")
