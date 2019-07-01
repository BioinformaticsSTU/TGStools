

from argparse import ArgumentParser, RawTextHelpFormatter


def create_parser(subparsers, fun):
	if(fun=="TransDisp"):
		parser_geneDisplay = subparsers.add_parser('TransDisp', help='create a macroscopic image showing transcripts of queried gene');
		group_input = parser_geneDisplay.add_argument_group("Input files arguments");
		group_input.add_argument('-g','--gtf', required=True, help="gtf file");
		group_input.add_argument('-q','--quant', required=False, help="quantity of transcript");
		group_input.add_argument('-i','--id', required=True, help="gene id that you want to query");
		group_input.add_argument('-p','--path', required=True, help="directory which contain histone files or fatom5 files");
		return(parser_geneDisplay);
	if(fun=="StaDist"):
		parser_staDist = subparsers.add_parser('StaDist', help='statistics the distance of transcript-start-site to the closest histone or fantom5 data(bed format)');
		group_input = parser_staDist.add_argument_group("Input files arguments");
		group_input.add_argument('-g','--gtf', required=True, help="gtf file");
		group_input.add_argument('-f','--flag', required=True, choices=['histone', 'fantom5'], help="flag that tells programme the type of files in path");
		group_input.add_argument('-p','--path', required=True, help="directory which contain histone files or fatom5 files");
		group_input.add_argument('-r','--prefix', required=True, help="prefix of output files");
		return(parser_staDist);
	if(fun=="TransFilt"):
		parser_TransFilt = subparsers.add_parser('TransFilt', help='statistics the distance of transcript-start-site to the closest histone or fantom5 data(bed format)');
		group_input = parser_TransFilt.add_argument_group("Input files arguments");
		group_input.add_argument('-g','--gtf', required=True, help="gtf file");
		group_input.add_argument('-d','--distance', required=True, help="file of distance between transcript-start-site to the closest histone site");
		group_input.add_argument('-t','--threshold', required=True, help="threshold for distance between transcript-start-site of transcript to the closest histone site");
		group_input.add_argument('-p','--prefix', required=True, help="prefix of output files");
		return(parser_TransFilt);
	if(fun=="StaAS"):
		parser_staAS = subparsers.add_parser('StaAS', help='calculate the proportion of each alternative splicing event in different samples and create graphs');
		group_input = parser_staAS.add_argument_group("Input files arguments");
		group_input.add_argument("-n", "--names", required=True, help="names of samples");
		group_input.add_argument("-p", "--prefix", required=True, help="prefix of output file");
		return(parser_staAS);
	if(fun=="CalScoreD"):
		parser_calScoreD = subparsers.add_parser('CalScoreD', help='calculate score_D(its formula can be seen below) of each gene');
		group_input = parser_calScoreD.add_argument_group("Input files arguments");
		group_input.add_argument("-c", "--control", required=True, help="control sample");
		group_input.add_argument("-t", "--treated", required=True, help="treated samples");
		group_input.add_argument("-p", "--prefix", required=True, help="prefix of output file");
		return(parser_calScoreD);
	if(fun=="GOEnrich"):
		parser_GOenrich = subparsers.add_parser('GOEnrich', help='select top genes and do GO enrichment analysis');
		group_input = parser_GOenrich.add_argument_group("Input files arguments");
		group_input.add_argument("-i", "--input", required=True, help="score D result");
		group_input.add_argument("-t", "--threshold", type=float, required=True, help="threshold for adjusted p-value");
		group_input.add_argument("-n", "--number", type=int, required=True, help="number of top genes for analysis.");
		group_input.add_argument("-f", "--type", choices=['bar', 'scatter', 'all'],required=True, help="type of image, 'bar', 'scatter' and 'all' can be chosen.");
		group_input.add_argument("-p", "--prefix", required=True, help="prefix of output file");
		return(parser_GOenrich);
	if(fun=="LncPred"):
		parser_INCP = subparsers.add_parser('LncPred', help='an integration classification tool of CNCI and PLEK for identify coding or non-coding transcripts (fasta file and gtf file)');
		group_input = parser_INCP.add_argument_group("Input files arguments");
		group_input.add_argument('-i','--input', required=True, help="transcript file(sequence or gtf format)");
		group_input.add_argument('-p','--parallel', required=True, help="specified speed ratio");
		group_input.add_argument('-g','--gtf', required=False, dest='gtf',action='store_true', help="please enter your gtf files");
		group_input.add_argument('-r','--reference', required=False, help="if your input file is gtf format please enter RefGenome directory");
		return(parser_INCP);
	if(fun=="LncExt"):
		parser_extract_lncRNA_gtf = subparsers.add_parser('LncExt', help='extract lncRNA information of GTF format based on the tanscript ID of the candidate lncRNA');
		group_input = parser_extract_lncRNA_gtf.add_argument_group("input files arguments");
		group_input.add_argument('-i','--input', required=True, help="enter your transcript (contain the transcript ID)");
		group_input.add_argument('-o','--out', required=True, help="assign your output file");
		group_input.add_argument('-g','--gtf', required=True, help="please enter your gtf files");
		return(parser_extract_lncRNA_gtf);
	if(fun=="LncExtTiss"):
		parser_tiss_specific = subparsers.add_parser('LncExtTiss', help='extract cancer-specific lncRNA information of GTF format');
		group_input = parser_tiss_specific.add_argument_group("Input files arguments");
		group_input.add_argument('-i','--input', required=True, help="enter your transcript (contain the transcript ID)");
		group_input.add_argument('-o','--out', required=True, help="assign your output file");
		group_input.add_argument('-t','--tissue', required=False, help="please enter tissue name");
		group_input.add_argument('-r','--reference', choices=["hg38", "hg19"], required=True, help="please enter hg38 or hg19");
		return(parser_tiss_specific);

#######
def main():
	description = "";
	parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter, add_help=True);
	#subparsers = parser.add_subparsers(help='sub-command help');
	subparsers = parser.add_subparsers(dest='subcommand_name', help='sub-command help');

	# create the parser for the "TransDisp" command
	parser_geneDisplay=create_parser(subparsers, "TransDisp");
	# create the parser for the "StaDist" command
	parser_staDist=create_parser(subparsers, "StaDist");
	# create the parser for the "TransFilt" command
	parser_TransFilt=create_parser(subparsers, "TransFilt");
	# create the parser for the "StaAS" command
	parser_staAS=create_parser(subparsers, "StaAS");
	# create the parser for the "CalScoreD" command
	parser_calScoreD=create_parser(subparsers, "CalScoreD");
	# create the parser for the "GOEnrich" command
	parser_GOenrich=create_parser(subparsers, "GOEnrich");
	# create the parser for the "Incp" command
	parser_INCP=create_parser(subparsers, "Incp");
	# create the parser for the "LncExt" command
	parser_extract_lncRNA_gtf=create_parser(subparsers, "LncExt");
	# create the parser for the "LncExtTiss" command
	parser_tiss_specific=create_parser(subparsers, "LncExtTiss");

	#############
	args = parser.parse_args();
	subcommand = args.subcommand_name;
	############
	if subcommand == "TransDisp":
		from geneDisplay import geneDisplay
		if(args.quant):
			geneDisplay(args.gtf, args.id, args.quant, args.path);
		else:
			geneDisplay(args.gtf, args.id, "", args.path);

	if subcommand == "StaDist":
		from staDist import staDist
		staDist(args.gtf, args.path, args.flag, args.prefix);

	if subcommand == "TransFilt":
		from TransFilt import TransFilt
		TransFilt(args.gtf, args.distance, args.threshold, args.prefix);

	if subcommand == "StaAS":
		from staAS import staAS
		staAS(args.names, args.prefix);

	if subcommand == "CalScoreD":
		from calScoreD import calScoreD
		calScoreD(args.control, args.treated, args.prefix);

	if subcommand == "GOenrich":
		from GOenrich import GOenrich
		GOenrich(args.input, args.threshold, args.number, args.type, args.prefix);

	if subcommand == "LncPred":
		from INCP import INCP
		INCP(args.input, args.parallel, args.gtf, args.reference);

	if subcommand == "LncExt":
		from extract_lncRNA_gtf import extract_lncRNA_gtf
		extract_lncRNA_gtf(args.input, args.out, args.gtf);

	if subcommand == "LncExtTiss":
		from tiss_specific import tiss_specific
		tiss_specific(args.out, args.tissue, args.reference, args.input);

if __name__ == '__main__':
    main();





