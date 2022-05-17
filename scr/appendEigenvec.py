import pandas as pd
import glob
import re, sys, argparse

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", help="path to eigenvec file ",type=str,required=True)
	parser.add_argument("-met", help="method", type=str,required=True)
	parser.add_argument("-ass", help="assembly", type=str,required=True)
	parser.add_argument("-len", help="length", type=str,required=True)
	parser.add_argument("-o", help="path to output file ",type=str,required=True)
	args = parser.parse_args()

	df = pd.DataFrame()
	fileList=glob.glob(args.i)
	for f in fileList:
		tmp = pd.read_csv(f,skiprows=1, sep = "\t", header = None)
		tmp.columns = ['chrom', 'id1', 'id2', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
		tmp = tmp[['chrom', 'id1', 'PC1', 'PC2']]
		df = df.append(tmp)
	df.loc[:, 'method'] = args.met
	df.loc[:, 'assembly'] = args.ass
	df.loc[:, 'length'] = args.len

	pop=pd.read_csv("data/hprc_samples.tsv", sep = "\t")
	df2 = df.merge(pop, left_on = "id1", right_on = "Sample")
	df2 = df2[['chrom', 'Sample', 'PC1', 'PC2', 'Superpopulation', 'method', 'assembly', 'length']]
	df2.to_csv(args.o, sep = "\t", index=False)

if __name__ == "__main__":
	main()