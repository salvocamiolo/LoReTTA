


import os,sys
from Bio import SeqIO
import numpy as np
import time
from datetime import datetime

import argparse

parser = argparse.ArgumentParser(description="Tool to perform reference guided de novo of PacBio reads")
parser.add_argument("-i","--inputReads",required=True, help="The Pacbio reads file in fastq format")
parser.add_argument('-ref',"--reference",required=True,help="Reference file in fasta format")
parser.add_argument("-o","--outputFolder",required=True,help="The output folder name")
parser.add_argument("-q","--quality",required=True, help="Phred quality threshold")
parser.add_argument("-wsize","--windowSize",required=True, help="Window size")
parser.add_argument("-wstep","--windowStep",required=True, help="Window step")
parser.add_argument("-t","--threads",required=True,help="Number of threads")
parser.add_argument("-p","--prefix",required=True,help="Prefix of output files")
parser.add_argument("-cr","--closingReads",required=True,help="Fasta file with reads to be used to close the gaps")
parser.add_argument("-ho","--homology",required=True,help="The homology the target and the reference genome share")

args = vars(parser.parse_args())
outputFolder = args['outputFolder']
inputReadsFile = args['inputReads']
refFile = args['reference']
threshold = int(args['quality'])
windowSize = int(args['windowSize'])
windowStep = int(args['windowStep'])
numThreads = args['threads']
prefix = args['prefix']
closingReads = args['closingReads']
homology = args['homology']
#Convert the input fastq file in fast format


def chopReads(inputFile):
	numSeq = 0
	outfile = open(inputFile+"_chopped.fasta","w")
	for seq_record in SeqIO.parse(inputFile,"fasta"):
		
		sequence = str(seq_record.seq)
		for a in range(0,len(sequence)-150,+150):
			numSeq+=1
			outfile.write(">ChoppedSeq_"+str(numSeq)+"\n"+sequence[a:a+150]+"\n")
	outfile.close()

logFile = open(outputFolder+"/"+prefix+"_assembly.log","w")
startTime = datetime.now()
current_time = startTime.strftime("%H:%M:%S")
logFile.write("De novo assembly started at "+str(current_time)+"\n\n")

#Load reference genome in memory

for seq_record in SeqIO.parse(refFile,"fasta"):
	refSeq = str(seq_record.seq)


reads = inputReadsFile
if ".fastq" in reads or ".fq" in reads:
	print("* * Converting input file from fastq to fasta....\n")
	
	
	os.system(installationDirectory+"/src/conda/bin/fq2fa "+inputReadsFile+" "+outputFolder+'/originalReads.fasta')

	reads = outputFolder+'/originalReads.fasta'


#Perform reference guided de novo assembly
if reads == "":

	exit()

else:
	
	readsSeq = {}

	print("* Reference guided de novo assembly")
	print("* * Loading reads in memory")
	

	#Loading high quality reads in memory
	for seq_record in SeqIO.parse(reads,"fasta"):
		if not str(seq_record.id) in readsSeq:
			readsSeq[str(seq_record.id)] = str(seq_record.seq)

	

	stage_a = open(outputFolder+"/local_assemblies.fasta","w")




	print("* * Assembly on sliding windows started")
	
	now = datetime.now()
	current_time = now.strftime("%H:%M:%S")
	logFile.write("Sliding window assembly started at "+str(current_time)+"\n\n")
	logFile.write("Range\tContig_size\n")
	lastWindow = False
	for a in range(0,len(refSeq),+windowStep):
		if lastWindow == True:
			break
	#for a in range(1):
		endPos = a+windowSize
		if endPos>len(refSeq):
			endPos=len(refSeq)
			windowSize = len(refSeq) - a
			lastWindow = True
		windowSuffix = str(a)+"-"+str(endPos) 
		print("* * * Assembling region "+str(a)+"-"+str(endPos))
		
		partSeq = refSeq[a:endPos]
		tempFasta = open(outputFolder+"/partReference.fasta","w")
		tempFasta.write(">partReference\n"+partSeq+"\n")
		tempFasta.close()

		os.system(installationDirectory+"/src/conda/bin/minimap2 -x map-pb -t "+numThreads+" "+outputFolder+"/partReference.fasta "+reads+" > "+outputFolder+"/outputMinimap_"+windowSuffix)

		os.system("awk '($11/$2)>"+homology+"' "+outputFolder+"/outputMinimap_"+windowSuffix+"  | sort -k2rn,2rn >  "+outputFolder+"/outputMinimap_filtered_"+windowSuffix)
		# awk '($10/$2)>0.5' |
		totalCollectedBases = 0
		infile = open(outputFolder+"/outputMinimap_filtered_"+windowSuffix)
		outfile = open(outputFolder+"/mapped.fasta_"+windowSuffix,"w")
		while True:
			line = infile.readline().rstrip()
			if not line:
				break
			fields = line.split("\t")
			outfile.write(">"+fields[0]+"\n"+readsSeq[fields[0]]+"\n")
			totalCollectedBases+=len(str(readsSeq[fields[0]]))
			if totalCollectedBases > windowSize*1000: #do not collect more than a number of reads leading to a 1000x coverage of the window
				break

		outfile.close()
		infile.close()
		os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/hqKmerAssembly.py -p "+installationDirectory+" -r "+outputFolder+"/mapped.fasta_"+windowSuffix+" -ref "+outputFolder+"/partReference.fasta -t "+numThreads+" -of "+outputFolder)


		if os.path.isfile(outputFolder+"/mapped.fasta_"+windowSuffix+"_localAssembly.fasta") == True:
			for seq_record in SeqIO.parse(outputFolder+"/mapped.fasta_"+windowSuffix+"_localAssembly.fasta","fasta"):

				maxScaffoldLength = len(str(seq_record.seq))
				longestContig = str(seq_record.seq)

			print("* * * Contig size: "+str(maxScaffoldLength))
			
			logFile.write("Range "+str(a)+"-"+str(a+windowSize)+"\t"+str(maxScaffoldLength)+"\n")

			stage_a.write(">Range_"+str(a)+"_"+str(endPos)+"\n"+longestContig+"\n")
			#print("Finito")
			#sys.stdin.read(1)
		else:
			print("No assembly")
			


	stage_a.close()


	#Joining contigs
	#Count contigs in local_assemblies
	numContigs = 0
	for seq_record in SeqIO.parse(outputFolder+"/local_assemblies.fasta","fasta"):
		numContigs+=1
	
	if numContigs>1:
		now = datetime.now()
		current_time = now.strftime("%H:%M:%S")
		logFile.write("Conting joining started at "+str(current_time)+"\n\n")

		print("* * Joining contigs.... ")
		
		os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/contigsJoiner.py -c "+outputFolder+"/local_assemblies.fasta -r "+refFile+" -p "+installationDirectory+" -o "+outputFolder)

		#Check final number of scaffold and attempt gap closure if > 1
		print("* * Attempting gap closure.... ")
		
		finalScaffols = SeqIO.to_dict(SeqIO.parse(outputFolder+"/scaffolds.fasta","fasta"))

		if len(finalScaffols)>1:
			print("\nAttempting gap closure.... ")
			
			now = datetime.now()
			current_time = now.strftime("%H:%M:%S")
			logFile.write("Gap closure started at "+str(current_time)+"\n\n")
			os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/oneReadContigsJoiner.py \
				-p "+installationDirectory+ " -c "+outputFolder+"/scaffolds.fasta -r "+refFile+" -x " + \
					outputFolder+" -s "+ closingReads +" -o scaffolds_gapClosed.fasta")
		else:
			os.system("cp "+outputFolder+"/scaffolds.fasta "+outputFolder+"/scaffolds_gapClosed.fasta")
	
	else:
		os.system("cp "+outputFolder+"/local_assemblies.fasta "+outputFolder+"/scaffolds_gapClosed.fasta")
	
	
	#Final alignment and consensus calling
	print("* * Calling consensus.... ")
	#print("* * * Subsampling.... ")
	
	now = datetime.now()
	current_time = now.strftime("%H:%M:%S")
	logFile.write("Consensus calling started at "+str(current_time)+"\n\n")
	"""outfile = open(outputFolder+"/subSample.fasta","w")
	totCoverage = 0
	for seq_record in SeqIO.parse(reads,"fasta"):
		totCoverage+=len(str(seq_record.seq))
		SeqIO.write(seq_record,outfile,"fasta")
		if totCoverage > 100*len(refSeq):
			break
	outfile.close()
"""
	print("* * * Assembly correction ")
	print("* * * Chopping reads.... ")
	
	chopReads(reads)
	print("* * * First consensus calling.... ")
	print("* * * Mapping original reads to the assembled sequence.... ")
	
	
	
	os.system(installationDirectory+"/src/conda/bin/minimap2 -a -x map-pb -t "+numThreads+" "+outputFolder+"/scaffolds_gapClosed.fasta "+reads+"_chopped.fasta"+" > "+outputFolder+"/alignment1.sam")
	print("* * * Converting sam to bam.... ")
	
	os.system(installationDirectory+"/src/conda/bin/samtools view -F 4 -bS -h "+outputFolder+"/alignment1.sam > "+outputFolder+"/alignment1.bam")
	print("* * * Sorting.... ")
	
	os.system(installationDirectory+"/src/conda/bin/samtools sort -o "+outputFolder+"/alignment_sorted1.bam "+outputFolder+"/alignment1.bam")
	print("* * * Indexing.... ")
	
	os.system(installationDirectory+"/src/conda/bin/samtools index "+outputFolder+"/alignment_sorted1.bam")

	print("* * * Creating pilleup.... ")
	
	os.system(installationDirectory+"/src/conda/bin/samtools mpileup -f "+outputFolder+ \
		"/scaffolds_gapClosed.fasta "+outputFolder +"/alignment_sorted1.bam > "+outputFolder+ \
			"/pileup1.txt")
	os.system("awk '$3!=\"N\"' "+outputFolder+"/pileup1.txt >"+outputFolder+"/pileup2.txt")
	os.system("mv "+outputFolder+"/pileup2.txt "+outputFolder+"/pileup1.txt")

	print("* * * Calling variants.... ")
	
	os.system(installationDirectory+"/src/conda/bin/varscan mpileup2cns "+outputFolder+"/pileup1.txt --variants --output-vcf --min-var-freq 0.5 --min-avg-qual 0 --strand-filter 0 --min-coverage 5   > "+outputFolder+"/output1.vcf")
	
	os.system(installationDirectory+"/src/conda/bin/bgzip -f -c "+outputFolder+"/output1.vcf > "+outputFolder+"/output1.vcf.gz")
	os.system(installationDirectory+"/src/conda/bin/tabix -f "+outputFolder+"/output1.vcf.gz")
	os.system("cat "+outputFolder+"/scaffolds_gapClosed.fasta | "+installationDirectory+"/src/conda/bin/bcftools consensus "+outputFolder+"/output1.vcf.gz > "+outputFolder+"/finalAssembly1.fasta")



	print("* * * second consensus calling.... ")
	print("* * * Mapping original reads to the assembled sequence.... ")
	
	
	
	os.system(installationDirectory+"/src/conda/bin/minimap2 -a -x map-pb -t "+numThreads+" "+outputFolder+"/finalAssembly1.fasta "+reads+"_chopped.fasta"+" > "+outputFolder+"/alignment2.sam")
	print("* * * Converting sam to bam.... ")
	
	os.system(installationDirectory+"/src/conda/bin/samtools view -F 4 -bS -h "+outputFolder+"/alignment2.sam > "+outputFolder+"/alignment2.bam")
	print("* * * Sorting.... ")
	
	os.system(installationDirectory+"/src/conda/bin/samtools sort -o "+outputFolder+"/alignment_sorted2.bam "+outputFolder+"/alignment2.bam")
	print("* * * Indexing.... ")
	
	os.system(installationDirectory+"/src/conda/bin/samtools index "+outputFolder+"/alignment_sorted2.bam")

	print("* * * Creating pilleup.... ")
	
	os.system(installationDirectory+"/src/conda/bin/samtools mpileup -f "+outputFolder+ \
		"/finalAssembly1.fasta "+outputFolder +"/alignment_sorted2.bam > "+outputFolder+ \
			"/pileup2.txt")

	os.system("awk '$3!=\"N\"' "+outputFolder+"/pileup2.txt >"+outputFolder+"/pileup3.txt")
	os.system("mv "+outputFolder+"/pileup3.txt "+outputFolder+"/pileup2.txt")

	print("* * * Calling variants.... ")
	
	os.system(installationDirectory+"/src/conda/bin/varscan mpileup2cns "+outputFolder+"/pileup2.txt --variants --output-vcf --min-avg-qual 0 --strand-filter 0 --min-coverage 5   > "+outputFolder+"/output2.vcf")
	os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/varscanFilter.py -i "+outputFolder+"/output2.vcf -o "+outputFolder+"/output_filtered2.vcf -1 "+reads+" -g 1  -r "+outputFolder+"/finalAssembly1.fasta -p "+installationDirectory +" -t "+numThreads+" -of "+outputFolder) 
	os.system(installationDirectory+"/src/conda/bin/bgzip -f -c "+outputFolder+"/output_filtered2.vcf > "+outputFolder+"/output_filtered2.vcf.gz")
	os.system(installationDirectory+"/src/conda/bin/tabix -f "+outputFolder+"/output_filtered2.vcf.gz")
	os.system("cat "+outputFolder+"/finalAssembly1.fasta | "+installationDirectory+"/src/conda/bin/bcftools consensus "+outputFolder+"/output_filtered2.vcf.gz > "+outputFolder+"/finalAssembly.fasta")





	print("\n\n De novo assembly finished!")
	
	endTime = datetime.now()
	current_time = endTime.strftime("%H:%M:%S")
	logFile.write("De novo assembly finished at "+str(current_time)+"\n\n")
	totalTime = (endTime -startTime).total_seconds()
	logFile.write("Total processing time: "+str(totalTime)+" seconds")
	logFile.close()
