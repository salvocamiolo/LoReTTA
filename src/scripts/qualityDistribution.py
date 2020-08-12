import matplotlib.pyplot as plt 
import sys,os
from Bio import SeqIO

readFile = sys.argv[1]
outputFolder = sys.argv[2]

totQualityValues = []
for seq_record in SeqIO.parse(readFile,"fastq"):
    quality = seq_record.letter_annotations["phred_quality"]
    totQualityValues+=quality
    if len(totQualityValues)>1000000:
        break

fig = plt.figure(figsize=(10,10))
plt.hist(totQualityValues)
plt.xlabel("quality score")
fig.savefig("qualityScoreDistribution.png",dpi=300)



