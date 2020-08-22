import os,sys
from Bio import SeqIO
from Bio import Seq
import argparse
parser = argparse.ArgumentParser(description="Filter varscan output in vcf format by only keeping variants with the highest number of reads")
parser.add_argument("-o","--output",required=True,help="The varscn output file")
parser.add_argument("-i","--input",required=True,help="The varscan input file")
parser.add_argument("-g","--validateIndel",required=False,help="Validate indels by searching the flanking motif in the reads")
parser.add_argument("-1","--read1",required=False,help="First fastq file in reads, needed for indel validation")
parser.add_argument("-r","--reference",required=False,help="Reference fasta file, needed for indel validation")
parser.add_argument("-p","--installationDirectory",required=False,help="Full path to the GRACy directory")
parser.add_argument("-t","--threads",required=False,help="Number of threads")
parser.add_argument("-of","--outputFolder",required=False,help="Number of threads")


args = vars(parser.parse_args())
inputFile = args['input']
outputFile = args['output']
reference = args['reference']
numThreads = args['threads']
outputFolder = args['outputFolder']
if not args['validateIndel']:
	validate = 0
else:
	validate = int(args['validateIndel'])

read1 = args['read1']

installationDirectory = args['installationDirectory']
kmersInReads = []
if validate == 1:
	print("Counting kmers")
	#os.system(installationDirectory+"src/conda/bin/jellyfish count -m 50 -s 1G -t "+numThreads+" -C  -o kmerCount.jf "+read1)
	os.system(installationDirectory+"/src/conda/bin/kmc -cs20000 -fa  -k50 "+read1+" "+outputFolder+"/kmcDB " + outputFolder+"/ > "+outputFolder+"/null 2>&1")
	os.system(installationDirectory+"/src/conda/bin/kmc_dump  "+outputFolder+"/kmcDB "+ outputFolder+"/kmcDB_output > "+outputFolder+"/null 2>&1")

	kmerCount = {}
	infile = open(outputFolder+"/kmcDB_output")
	while True:
		line = infile.readline().rstrip()
		if not line:
			break
		fields = line.split("\t")
		if not fields[0] in kmerCount:
			kmerCount[fields[0]]= int(fields[1])

	refSeq = {}
	for seq_record in SeqIO.parse(reference,"fasta"):
		if not str(seq_record.id) in refSeq:
			refSeq[str(seq_record.id)] = str(seq_record.seq)
	




infile = open(inputFile)
outfile = open(outputFile,"w")

line = infile.readline().rstrip()
outfile.write(line+"\n")

while not "#CHROM" in line:
	line = infile.readline().rstrip()
	if not line:
		break
	outfile.write(line+"\n")

while True:
	line = infile.readline().rstrip()
	if not line:
		break
	fields = line.split("\t")
	info = fields[9].split(":")
	if validate == 1:
		numValidated = 0

		if len(fields[4])>1 and (int(info[5]) > 0.7*int(info[4])): #a deletion
			numValidated+=1
			if numValidated%100==0:
				print("Validated %d deletion" %numValidated)
			
			#print(refSeq[fields[0]][int(fields[1])-25:int(fields[1])-1]+"-"+fields[3]+"-"+refSeq[fields[0]][int(fields[1]):int(fields[1])+25])
			#print(refSeq[fields[0]][int(fields[1])-25:int(fields[1])-1]+"-"+fields[4]+"-"+refSeq[fields[0]][int(fields[1]):int(fields[1])+25-len(fields[4])+1])
			refAllele = refSeq[fields[0]][int(fields[1])-25:int(fields[1])-1]+fields[3]+refSeq[fields[0]][int(fields[1]):int(fields[1])+25]
			altAllele = refSeq[fields[0]][int(fields[1])-25:int(fields[1])-1]+fields[4]+refSeq[fields[0]][int(fields[1]):int(fields[1])+25-len(fields[4])+1]
			"""refReads = kmersInReads.count(refAllele) + kmersInReads.count(Seq.reverse_complement(refAllele))
			altReads = kmersInReads.count(altAllele) + kmersInReads.count(Seq.reverse_complement(altAllele))

			#print(refAllele)
			os.system(installationDirectory+"src/conda/bin/jellyfish query kmerCount.jf "+refAllele+" >count.txt")
			countFile = open("count.txt")
			countLine = countFile.readline().rstrip()
			if not countLine:
				refReads = 0
			else:
				countList = countLine.split(" ")
				refReads = int(countList[1])
			countFile.close()
			#print(refReads)

			#print(altAllele)
			os.system(installationDirectory+"src/conda/bin/jellyfish query kmerCount.jf "+altAllele+" >count.txt")
			countFile = open("count.txt")
			countLine = countFile.readline().rstrip()
			if not countLine:
				altReads=0
			else:
				countList = countLine.split(" ")
				altReads = int(countList[1])
			countFile.close()
			#print(altReads)"""
		

			if not altAllele in kmerCount:
				altAllele = Seq.reverse_complement(altAllele)
			
			if not refAllele in kmerCount:
				refAllele = Seq.reverse_complement(refAllele)

			if refAllele in kmerCount and altAllele in kmerCount: #a validation can be computed
				if kmerCount[altAllele]>kmerCount[refAllele]:
					outfile.write(line+"\n")
			else:
				if int(info[5]) > int(info[4]): #the alternate allele has a highe number of reads than the reference allele
					outfile.write(line+"\n")

		if len(fields[3])>1 and (int(info[5]) > 0.7*int(info[4])):
			#print("Validating deletion at position %s" %fields[1])
			#print(refSeq[int(fields[1])-25:int(fields[1])-1]+"-"+fields[3]+"-"+refSeq[int(fields[1])+len(fields[3])-1:int(fields[1])+25])
			#print(refSeq[int(fields[1])-25:int(fields[1])-1]+"-"+fields[4]+"-"+refSeq[int(fields[1])+len(fields[3])-1:int(fields[1])+25+len(fields[3])-1])
			refAllele = refSeq[fields[0]][int(fields[1])-25:int(fields[1])-1]+fields[3]+refSeq[fields[0]][int(fields[1])+len(fields[3])-1:int(fields[1])+25]
			altAllele = refSeq[fields[0]][int(fields[1])-25:int(fields[1])-1]+fields[4]+refSeq[fields[0]][int(fields[1])+len(fields[3])-1:int(fields[1])+25+len(fields[3])-1]
			"""refReads = kmersInReads.count(refAllele) + kmersInReads.count(Seq.reverse_complement(refAllele))
			altReads = kmersInReads.count(altAllele) + kmersInReads.count(Seq.reverse_complement(altAllele))

			#print(refAllele)
			os.system(installationDirectory+"src/conda/bin/jellyfish query kmerCount.jf "+refAllele+" >count.txt")
			countFile = open("count.txt")
			countLine = countFile.readline().rstrip()
			if not countLine:
				refReads = 0
			else:
				countList = countLine.split(" ")
				refReads = int(countList[1])
			countFile.close()
			#print(refReads)

			#print(altAllele)
			os.system(installationDirectory+"src/conda/bin/jellyfish query kmerCount.jf "+altAllele+" >count.txt")
			countFile = open("count.txt")
			countLine = countFile.readline().rstrip()
			if not countLine:
				altReads=0
			else:
				countList = countLine.split(" ")
				altReads = int(countList[1])
			countFile.close()
			#print(altReads)"""
		
			if not altAllele in kmerCount:
				altAllele = Seq.reverse_complement(altAllele)
			
			if not refAllele in kmerCount:
				refAllele = Seq.reverse_complement(refAllele)

			if refAllele in kmerCount and altAllele in kmerCount: #a validation can be computed
				if kmerCount[altAllele]>kmerCount[refAllele]:
					outfile.write(line+"\n")
			else:
				if int(info[5]) > int(info[4]): #the alternate allele has a highe number of reads than the reference allele
					outfile.write(line+"\n")


		if len(fields[3]) == 1 and len(fields[4]) == 1: #just a SNP
			if int(info[5]) > int(info[4]): #the alternate allele has a highe number of reads than the reference allele
				outfile.write(line+"\n")



		
	else:
		if int(info[5]) > int(info[4]): #the alternate allele has a highe number of reads than the reference allele
			outfile.write(line+"\n")


infile.close()
outfile.close()


	



