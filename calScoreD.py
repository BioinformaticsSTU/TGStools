

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
plt.switch_backend('agg');



########## calScoreD
def calScoreD(control, treated, prefix):
	SAMPLES=[];
	SAMPLES.append(control);
	treated=treated.replace(' ', '').split(",");
	for sample in treated:
		SAMPLES.append(sample);
	##########
	# get gene and transcripts from ioi file, and produce GT file
	for sample in SAMPLES:
		directory=sample+"_Events/";
		file=directory+sample+".ioi";
		if(not os.path.exists(file)):
			print("Error: "+file+" doesn't exist!");
			return();
		fp=open(file, "r");
		line=fp.readline();
		GENES_SET=[];
		file_w=directory+sample+"_GT.txt";
		fp_w=open(file_w, "w");
		while(line!=""):
			words=line.split("\t");
			gene_id=words[1];
			if(gene_id in GENES_SET):
				line=fp.readline();
				continue;
			total_transcripts=words[4];
			fp_w.write(gene_id+"\t"+total_transcripts);
			GENES_SET.append(gene_id);
			line=fp.readline();
		fp.close();
		fp_w.close();

	#################
	# get gene from GT file
	def get_gene(file):
		fp=open(file);
		line=fp.readline();
		line=fp.readline();
		GENE_SET=[];
		while(line!=""):
			gene=line.split("\t")[0];
			GENE_SET.append(gene);
			line=fp.readline();
		fp.close();
		return(GENE_SET);


	#################
	# get gene and transcripts
	def get_GT(GENE_COMM, file):
		GENE_tmp=[]; # all gene id from file
		TF_tmp=[]; # all transcripts from file
		TF=[]; # transcripts of GENE_COMM
		fp=open(file, "r");
		line=fp.readline();
		line=fp.readline();
		while(line!=""):
			words=line.split();
			gene_id=words[0];
			GENE_tmp.append(words[0]);
			TF_tmp.append(words[1]);
			line=fp.readline();
		fp.close();
		for gene in GENE_COMM:
			num=GENE_tmp.index(gene);
			TF.append(TF_tmp[num]);
		return(GENE_COMM, TF);

	# give two transcripts list, calculate score_D
	def get_score(TF_N, TF_C):
		tfs_n=TF_N.split(",");
		tfs_c=TF_C.split(",");
		comm=0;
		for tf_n in tfs_n:
			if(tf_n in tfs_c):
				comm=comm+1;
		all=len(tfs_n)+len(tfs_c)-comm;
		return(float(1-comm/all));

	# find common gene and use get_score to get score_D
	def compare(GENE_N, GENE_C, TF_N, TF_C):
		gene_D=[];
		score_D=[];
		TF_COMM_N=[];
		TF_COMM_C=[];
		i=0;
		for gene_n in GENE_N:
			j=0;
			for gene_c in GENE_C:
				if(gene_n==gene_c):
					score_D.append(get_score(TF_N[i], TF_C[j]));
					gene_D.append(gene_c);
					TF_COMM_N.append(TF_N[i]);
					TF_COMM_C.append(TF_C[j]);
					break;
				j=j+1;
			i=i+1;
		return(gene_D, TF_COMM_N, TF_COMM_C, score_D);

	###############
	# get common gene
	GENE_NOR=get_gene(control+"_Events/"+control+"_GT.txt");
	for sample in treated:
		file=sample+"_Events/"+sample+"_GT.txt";
		GENE_CAN=get_gene(file);
		GENE_NOR=[val for val in GENE_NOR if val in GENE_CAN];

	GENE_COMM=GENE_NOR;

	##########
	# get score_D and save data into file_compare
	directory=control+"_Events/";
	file_N=directory+control+"_GT.txt";
	[GENE_N, TF_N]=get_GT(GENE_COMM, file_N);

	# create a 0 score_D list
	final_score_D=[];
	for gene in GENE_COMM:
		final_score_D.append(0);

	# create a matrix
	score_D_mat=[];
	score_D_mat.append(TF_N);
	header_score_D=[];
	header_score_D.append(control+"_transcripts");

	for sample in treated:
		header_score_D.append(sample+"_transcripts");
		header_score_D.append(sample+"_score");
		file_compare=control+"_"+sample+".txt";
		directory=sample+"_Events/";
		file_C=directory+sample+"_GT.txt";
		[GENE_C, TF_C]=get_GT(GENE_COMM, file_C);
		[gene_D, TF_COMM_N, TF_COMM_C, score_D]=compare(GENE_N, GENE_C, TF_N, TF_C);
		score_D_mat.append(TF_COMM_C);
		score_D_mat.append(score_D);
		i=0;
		for score in final_score_D:
			final_score_D[i]=final_score_D[i]+score_D[i];
			i=i+1;

	header_score_D.append("score_D");
	score_D_mat.append(final_score_D);

	# save data into one file
	score_D_DF= pd.DataFrame(score_D_mat, index=header_score_D, columns=GENE_COMM).T;
	score_D_DF.to_csv(prefix+"_score_D.txt", sep="\t");

#############################



