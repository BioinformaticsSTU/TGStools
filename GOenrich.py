

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import copy
import gseapy
from gseapy.parser import Biomart
import pandas as pd
import math
from argparse import ArgumentParser, RawTextHelpFormatter
plt.switch_backend('agg');

"""
description = "";
parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter, add_help=True);
parser.add_argument("-i", "--input", required=True, help="score D result");
parser.add_argument("-t", "--threshold", required=True, help="threshold for adjusted p-value");
parser.add_argument("-n", "--number", required=True, help="number of top genes for analysis.");
parser.add_argument("-f", "--type", required=True, help="type of image, 'bar', 'scatter' and 'all' can be chosen.");
"""

###############################

def GOenrich(input, threshold, number, type, prefix):
	GO_name=prefix;
	p_threshold=float(threshold); # 0.05
	num=int(number); #500
	if(not os.path.exists(input)):
		print(input+" doesn't exist!");
		return();
	if(not type in ["bar", "scatter", "all"]):
		print("Error: arguement type wrong!");
		return();
	if(not isinstance(num, int)):
		print("Error: arguement number wrong!")
		return();
	if(not isinstance(p_threshold, float)):
		print("Error: arguement threshold wrong!")
		return();
	############
	# get GENE_COMM and final_score_D
	GENE_COMM=[];
	final_score_D=[];
	file_D=input;
	fp=open(file_D, "r");
	line=fp.readline();
	line=fp.readline();
	while(line!=""):
		words=line.split("\t");
		GENE_COMM.append(words[0]);
		final_score_D.append(words[-1].strip());
		line=fp.readline();

	fp.close();

	###########
	# sort data and get top 500 genes
	data = np.array([GENE_COMM, final_score_D]);
	idex=np.lexsort([data[1,:]]);
	data_sorted=data[:, idex[::-1]];
	GENES_top=data_sorted[0, :num];
	SCORES_top=data_sorted[1, :num];

	################
	## use BioMart and convert ensembl_id to gene_symbol
	bm = Biomart(verbose=False, host="asia.ensembl.org")
	marts = bm.get_marts()
	datasets = bm.get_datasets(mart='ENSEMBL_MART_ENSEMBL')
	attrs = bm.get_attributes(dataset='hsapiens_gene_ensembl')
	filters = bm.get_filters(dataset='hsapiens_gene_ensembl')
	results = bm.query(dataset='hsapiens_gene_ensembl', attributes= ["ensembl_gene_id", "hgnc_symbol"], filters={'ensembl_gene_id': GENES_top.tolist()}, filename="query.results.txt")
	gene_List=[];
	hgnc_symbol_list=results.hgnc_symbol.tolist();
	for gene in hgnc_symbol_list:
		if(isinstance(gene, str)):
			gene_List.append(gene);

	libs=["GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018"];
	gseapy.enrichr(gene_list=gene_List, description=GO_name, gene_sets=libs, outdir=GO_name);

	##########
	# fig for GO enrichment analysis. the threshold of adjusted p-value is less than 0.05, calculate -lg(adjusted p-value) and use gene num 
	TERMS=[];
	P_VALUE=[];
	GO_terms=[];
	GENES_NUM=[];
	file_GO=prefix+"_GO_reports.txt";
	fp_w=open(file_GO, "w");
	mark=0;

	for lib in libs:
		GO_file=GO_name+"/"+lib+"."+GO_name+".enrichr.reports.txt";
		fp=open(GO_file, "r");
		line=fp.readline();
		if(mark==0):
			fp_w.write(line);
			mark=1;
		line=fp.readline();
		while(line!=""):
			fp_w.write(line);
			words=line.split("\t");
			term=words[1];
			p_value=float(words[4]);
			genes_num=len(words[4].split(";"))
			if(p_value<p_threshold):
				tmp=term.split(" ");
				TERMS.append(tmp[-1].strip('(').strip(')'));
				del tmp[-1];
				GO_terms.append(" ".join(tmp));
				P_VALUE.append(p_value);
				GENES_NUM.append(genes_num);
			line=fp.readline();
		fp.close();
	fp_w.close();
	########
	GO_value=[];
	for value in P_VALUE:
		GO_value.append(round(-math.log(value, 10), 3));

	# MARK, sort GO_value
	MAT=frame = pd.DataFrame({"GO_terms":GO_terms,"TERMS":TERMS, "GO_value":GO_value});
	MAT_sorted=MAT.sort_values(by = 'GO_value',axis = 0,ascending = True);
	GO_terms=MAT_sorted['GO_terms'].tolist();
	TERMS=MAT_sorted['TERMS'].tolist();
	GO_value=MAT_sorted['GO_value'].tolist();
	##############
	if(type=="bar" or type=="all"):
		fig = plt.figure(figsize=(9, 6));
		matplotlib.rcParams['font.sans-serif'] = ['SimHei'];
		matplotlib.rcParams['axes.unicode_minus'] = False;
		# color='steelblue'
		plt.barh(range(len(GO_value)), GO_value, height=0.7, color='lightskyblue', alpha=0.8);
		plt.yticks(range(len(GO_value)), TERMS);
		plt.xlim(0, 10);
		plt.xlabel("-lg(P)");
		plt.title("GO enrichment");
		"""
		for x, y in enumerate(GO_value):
			plt.text(y + 0.2, x - 0.1, '%s' % y);
		"""
		i=0;
		for x, y in enumerate(GO_value):
			if(len(GO_terms[i])>70):
				plt.text(0.2, x-0.2, '%s' % GO_terms[i][0:70]+"...", fontsize=12);
			else:
				plt.text(0.2, x-0.2, '%s' % GO_terms[i], fontsize=12);
			#plt.text(0.2, x - 0.3, '%s' % y, fontsize=12);
			i=i+1;
		plt.show();
		fig.savefig(prefix+"_GO_barh.png");
	if(type=="scatter" or type=="all"):
		fig = plt.figure(figsize=(9, 6));
		cm = plt.cm.get_cmap('RdYlGn');
		NUM=[];
		for num in GENES_NUM:
			NUM.append(100*num);
		sc = plt.scatter(GO_value, TERMS, c=GO_value, vmin=0, s=NUM, cmap=cm);
		plt.colorbar(sc);
		plt.xlabel("-lg(P)");
		plt.title("GO enrichment");
		plt.show()
		fig.savefig(prefix+"_GO_scatter.png");

##########


