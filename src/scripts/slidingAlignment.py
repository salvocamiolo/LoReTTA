import sys
from Bio import SeqIO
import os
import argparse


parser = argparse.ArgumentParser(description="Perform sliding window alignment")
parser.add_argument("-r","--reads",required=True, help="Reads to use in the alignment")
parser.add_argument("-f","--reference",required=True,help="Reference file to which mapping the reads")
parser.add_argument("-p","--installationDirectory",required=True,help="Conda environment installation directory")
parser.add_argument("-t","--threads",required=True,help="Number of threads")
parser.add_argument("-o","--outputFolder",required=True, help="Output file name")

args = vars(parser.parse_args())
refFile = args['reference']
reads = args['reads']
installationDirectory = args['installationDirectory']
numThreads = args['threads']
outputFolder = args['outputFolder']



for seq_record in SeqIO.parse(refFile,"fasta"):
    referenceSeq = str(seq_record.seq)
    seqID = str(seq_record.id)

outfile = open(outputFolder+"/sliding_alignment.sam","w")
outfile.write("@SQ\tSN:"+seqID+"\tLN:"+str(len(referenceSeq))+"\n")
outfile.close()
lastRound = False
for a in range(0,len(referenceSeq),+10000):
    if lastRound==True:
        break
    rangeToFill = 12000
    if (a+12000)>len(referenceSeq):
        rangeToFill = len(referenceSeq)-a
        lastRound = True

    print("Analyzing",a,a+rangeToFill)

    newSeq = ""
    for b in range(0,a):
        newSeq += "N"
    newSeq+=referenceSeq[a:a+rangeToFill]
    for b in range(a+rangeToFill,len(referenceSeq)):
        newSeq += "N"
    
    outfile = open("newSeq.fasta","w")
    outfile.write(">"+seqID+"\n"+newSeq+"\n")
    outfile.close()

    os.system(installationDirectory+"/src/conda/bin/minimap2 -t "+numThreads+" -a -x map-pb newSeq.fasta "+reads+" > temp.sam")
    os.system(installationDirectory+"/src/conda/bin/samtools view -F 4 temp.sam >>"+outputFolder+"/sliding_alignment.sam")
    os.system(installationDirectory+"/src/conda/bin/samtools view -bS -h "+outputFolder+"/sliding_alignment.sam >"+outputFolder+"/sliding_alignment.bam")
    os.system(installationDirectory+"/src/conda/bin/samtools sort -o "+outputFolder+"/sliding_alignment_sorted.bam "+outputFolder+"/sliding_alignment.bam")
    os.system(installationDirectory+"/src/conda/bin/samtools index "+outputFolder+"/sliding_alignment_sorted.bam ")
    os.system(installationDirectory+"/src/conda/bin/samtools mpileup -f "+refFile+" "+outputFolder+"/sliding_alignment_sorted.bam > "+outputFolder+"/pileup.txt")



