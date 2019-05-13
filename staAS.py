import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import copy
import pandas as pd
from  scipy.stats import chi2_contingency
from argparse import ArgumentParser, RawTextHelpFormatter
plt.switch_backend('agg');


class get_AS(object):
	def __init__(self, SAMPLES, prefix):
		SAMPLES_LIST=SAMPLES.replace(" ", "").split(",");
		LIST=[];
		for TMP in SAMPLES_LIST:
			tmp=TMP.split(".");
			del tmp[-1];
			TMP=".".join(tmp);
			LIST.append(TMP);
		self.SAMPLES = LIST;
		self.prefix= prefix;
		self.mat=[];

	def run(self):
		for sample in self.SAMPLES:
			gtf=sample+".gtf";
			if(not os.path.exists(gtf)):
				print(gtf+" doesn't exist!");
				return();
			directory=sample+"_Events/";
			if(not os.path.exists(directory)):
				os.mkdir(directory);
			header=directory+sample;
			os.system("python "+os.path.split(os.path.realpath(__file__))[0]+"/SUPPA-master/suppa.py generateEvents -i "+gtf+" -o "+header+" -e SE SS MX RI FL -f ioe");
			os.system("python "+os.path.split(os.path.realpath(__file__))[0]+"/SUPPA-master/suppa.py generateEvents -i "+gtf+" -o "+header+" -f ioi");

	def plot_bar(self, mat, prefix):
		df=self.mat;
		prefix=self.prefix;
		n_row=df.shape[0];
		n_col=df.shape[1];
		fig=plt.figure('barh_align'); 
		ax1=fig.add_subplot(1, 2, 1); 
		df.plot(kind='bar', figsize=(2*n_col, 2*n_row+5), color=["#E377C2", "#ECE812", "#9467BD", "#D62728", "#2CA02C", "#FF7F0E", "#1F77B4"], alpha=0.5);
		#df.plot(kind='bar', figsize=(2*n_col, 2*n_row+5), color=["#A6A8AA", "#EA8C8C", "#6EBDE4", "#8FD6A5", "#C6A7CF", "#FBC17A", "#FFF994"]);# 横柱形图 kind=barh
		plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., markerscale=10, fontsize=20);
		plt.xticks(fontsize=25, rotation=30);
		plt.yticks(fontsize=25);
		plt.ylabel("count", fontsize=30);
		# add a subplot with no frame
		ax2=fig.add_subplot(122, frameon=False);
		#delete ax2 from the figure
		fig.delaxes(ax2);
		plt.subplots_adjust(right=0.8);
		plt.savefig(prefix+'_bar.png', dpi=600, format='png');
		plt.show();

	def plot_barh_align(self, mat, prefix):
		df=self.mat;
		prefix=self.prefix;
		mat_proportion=copy.copy(df);
		sum_mat=df.sum(1).tolist();
		n_row=df.shape[0];
		n_col=df.shape[1];
		for i in range(0, n_row):
			for j in range(0, n_col):
				mat_proportion.iloc[i, j]=mat_proportion.iloc[i, j]/sum_mat[i];
		fig=plt.figure('barh_align'); 
		ax1=fig.add_subplot(1, 2, 1); 
		mat_proportion.plot(kind='bar', stacked=True, figsize=(1.2*n_col, 1*n_row+5), color=["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#ECE812", "#E377C2"], alpha=0.5);
		#mat_proportion.plot(kind='bar', stacked=True, figsize=(1.2*n_col, 1*n_row+5), color=["#A6A8AA", "#EA8C8C", "#6EBDE4", "#8FD6A5", "#C6A7CF", "#FBC17A", "#FFF994"]);
		plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., markerscale=10, fontsize=20);
		plt.xticks(fontsize=25, rotation=30);
		plt.yticks(fontsize=25);
		plt.ylabel("percentage", fontsize=30);
		# add a subplot with no frame
		ax2=fig.add_subplot(122, frameon=False);
		#delete ax2 from the figure
		fig.delaxes(ax2);
		plt.subplots_adjust(right=0.8);
		plt.show();
		#####
		plt.savefig(prefix+'_barh_align.png', dpi=600, format='png');
		plt.show();

	def plot_AS(self):
		data_DF=[];
		ASs=["A3", "A5", "AF", "AL", "MX", "RI", "SE"];
		for sample in self.SAMPLES:
			labels = [];
			fracs = [];
			for AS in ASs:
				directory=sample+"_Events/";
				IOE=directory+sample+"_"+AS+"_strict.ioe";
				count = len(open(IOE, 'r').readlines())-1;
				labels.append(AS);
				fracs.append(count);
			print(sample+":");
			print(labels);
			print(fracs);
			data_DF.append(fracs);
		self.mat=pd.DataFrame(data=data_DF, index=self.SAMPLES, columns=ASs);
		kf = chi2_contingency(self.mat)
		print('chisq-statistic=%.4f, p-value=%.4f, df=%i expected_frep=%s'%kf);
		if(kf[1]<0.05):
			print('Proportion of alternative splicing in each sample has not difference.');
		else:
			print('Proportion of alternative splicing in each sample has difference.');
		self.mat.to_csv(self.prefix+"_sta.txt", sep="\t");
		self.plot_bar(self.mat, self.prefix);
		self.plot_barh_align(self.mat, self.prefix);


def staAS(SAMPLES, prefix):
	obj=get_AS(SAMPLES, prefix);
	obj.run();
	obj.plot_AS();
	return;








