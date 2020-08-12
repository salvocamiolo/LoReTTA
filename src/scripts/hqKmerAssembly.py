import os,sys
from Bio import SeqIO
from Bio import Seq
import argparse


parser = argparse.ArgumentParser(description="Attempt an assembly with highly frequent kmers")
parser.add_argument("-p","--installationDirectory",required=True,help="The installation directory")
parser.add_argument("-r","--reads",required=True,help="fasta file with reads mapped to a genome or a portion")
parser.add_argument("-ref","--reference",required=True,help="fasta file with a reference genome")
parser.add_argument("-t","--threads",required=True,help="Number of threads")
parser.add_argument("-of","--outputFolder",required=True,help="output folder")
args = vars(parser.parse_args())

installationDirectory = args['installationDirectory']
numThreads = args['threads']
reference  = args['reference']
reads = args['reads']
outputFolder = args['outputFolder']

for seq_record in SeqIO.parse(reference,"fasta"):
    refSequence = str(seq_record.seq)
    refLength = len(refSequence)

assembledSequence = ""
kmerSize = 61
kmerCoverage = 50
noAssembly = False
os.system("cp  "+reads+" "+reads+".fasta")
os.system("rm -rf "+outputFolder+"/kmerDB*")
while float(len(assembledSequence))/float(refLength) < 0.9:
    print("Trying kmer size %d / kmer coverage %f" %(kmerSize,kmerCoverage))
    #prepare kmer database 

    os.system(installationDirectory+"/src/conda/bin/kmc -fa  -k"+str(kmerSize)+" "+reads+ " " +outputFolder+"/"+(reads.split("."))[-1]+" "+ outputFolder+"/ > "+outputFolder+"/null 2>&1")

    os.system(installationDirectory+"/src/conda/bin/kmc_dump -ci"+str(int(kmerCoverage))+" "+outputFolder+"/"+(reads.split("."))[-1]+" "+ outputFolder+"/"+(reads.split("."))[-1]+"_output > "+outputFolder+"/null 2>&1")

    
    infile = open(outputFolder+"/"+(reads.split("."))[-1]+"_output")

    outfile = open(outputFolder+"/"+(reads.split("."))[-1]+"_output.fastq","w")
    numSeq = 0
    os.system("rm -rf "+outputFolder+"/outputSpades")
    while True:
        line = infile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        numSeq+=1
        outfile.write("@Sequence_"+str(numSeq)+"\n"+fields[0]+"\n+\n")
        for a in range(len(fields[0])):
            outfile.write("G")
        outfile.write("\n")
    outfile.close()
    if float(numSeq)>refLength*0.9:
        os.system(installationDirectory+"/src/conda/bin/spades.py -t "+numThreads+" -s "+outputFolder+"/"+(reads.split("."))[-1]+"_output.fastq --pacbio "+reads+".fasta --phred-offset 33 --careful -o "+outputFolder+"/"+(reads.split("."))[-1]+"_outputSpades > "+outputFolder+"/null 2>&1")
        print("Spades completed")


        maxScaffoldLength = 0

        if os.path.isfile(outputFolder+"/"+(reads.split("."))[-1]+"_outputSpades/scaffolds.fasta") == True:
            for seq_record in SeqIO.parse(outputFolder+"/"+(reads.split("."))[-1]+"_outputSpades/scaffolds.fasta","fasta"):
                if len(str(seq_record.seq)) > maxScaffoldLength:
                    maxScaffoldLength = len(str(seq_record.seq))
                    assembledSequence = str(seq_record.seq)

    if kmerCoverage >11:
        kmerCoverage = kmerCoverage - 10
    else:
        kmerCoverage = kmerCoverage - 2
    if kmerCoverage < 2:
        kmerCoverage = 50
        kmerSize = kmerSize - 10
        if kmerSize <31:
            noAssembly = True
            break

if noAssembly==False:
    outfile = open(reads+"_localAssembly.fasta","w")
    outfile.write(">local\n"+assembledSequence+"\n")
    outfile.close()
    

