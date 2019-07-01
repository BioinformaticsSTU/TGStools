

def TransFilt(gtf, transDist, threshold, prefix):
	threshold=int(threshold);
	Trans=[];
	file_dist=transDist;
	fp=open(file_dist, "r");
	for line in fp.readlines():
		[trans, num]=line.strip().split("\t");
		if(int(num)<threshold):
			Trans.append(trans);
	fp.close();
	fp_w=open(prefix+"_trans.gtf", "w");
	file_gtf=gtf;
	fp=open(file_gtf, "r");
	for line in fp.readlines():
		#words=line.strip().split("\t");
		trans=line.strip().split("\t")[8].split(";")[1].strip().split(" ")[1].strip('"');
		if(trans in Trans):
			fp_w.write(line);
	fp.close();
	fp_w.close();

#####

#TransFilt("TEST/Sample1.gtf", "trans_dist.txt", 1000, "transFilt");












