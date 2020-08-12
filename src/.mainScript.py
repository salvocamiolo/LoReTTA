# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'GUI.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog 
from PyQt5.QtWidgets import QMessageBox

import os,sys
from Bio import SeqIO
import numpy as np
import time
import matplotlib.pyplot as plt
from datetime import datetime



class Ui_Form(object):
	def setupUi(self, Form,installationDirectory):
		Form.setObjectName("Form")
		Form.resize(846, 721)
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(10, 10, 101, 16))
		self.label.setObjectName("label")
		self.projectFolderLineEdit = QtWidgets.QLineEdit(Form)
		self.projectFolderLineEdit.setGeometry(QtCore.QRect(10, 30, 361, 21))
		self.projectFolderLineEdit.setObjectName("projectFolderLineEdit")
		self.projectFolderButton = QtWidgets.QPushButton(Form)
		self.projectFolderButton.setGeometry(QtCore.QRect(390, 30, 113, 21))
		self.projectFolderButton.setObjectName("projectFolderButton")
		self.label_2 = QtWidgets.QLabel(Form)
		self.label_2.setGeometry(QtCore.QRect(10, 90, 101, 16))
		self.label_2.setObjectName("label_2")
		self.referenceButton = QtWidgets.QPushButton(Form)
		self.referenceButton.setGeometry(QtCore.QRect(390, 370, 113, 21))
		self.referenceButton.setObjectName("referenceButton")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(10, 350, 101, 16))
		self.label_3.setObjectName("label_3")
		self.referenceLineEdit = QtWidgets.QLineEdit(Form)
		self.referenceLineEdit.setGeometry(QtCore.QRect(10, 370, 361, 21))
		self.referenceLineEdit.setObjectName("referenceLineEdit")
		self.label_13 = QtWidgets.QLabel(Form)
		self.label_13.setGeometry(QtCore.QRect(10, 440, 101, 16))
		self.label_13.setObjectName("label_13")
		self.windowSizeLineEdit = QtWidgets.QLineEdit(Form)
		self.windowSizeLineEdit.setGeometry(QtCore.QRect(10, 460, 101, 21))
		self.windowSizeLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.windowSizeLineEdit.setObjectName("windowSizeLineEdit")
		self.windowStepLineEdit = QtWidgets.QLineEdit(Form)
		self.windowStepLineEdit.setGeometry(QtCore.QRect(140, 460, 101, 21))
		self.windowStepLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.windowStepLineEdit.setObjectName("windowStepLineEdit")
		self.label_14 = QtWidgets.QLabel(Form)
		self.label_14.setGeometry(QtCore.QRect(140, 440, 101, 16))
		self.label_14.setObjectName("label_14")
		self.runButton = QtWidgets.QPushButton(Form)
		self.runButton.setGeometry(QtCore.QRect(600, 690, 113, 21))
		self.runButton.setObjectName("runButton")
		self.exitButton = QtWidgets.QPushButton(Form)
		self.exitButton.setGeometry(QtCore.QRect(720, 690, 113, 21))
		self.exitButton.setObjectName("exitButton")
		self.label_16 = QtWidgets.QLabel(Form)
		self.label_16.setGeometry(QtCore.QRect(20, 490, 111, 16))
		self.label_16.setObjectName("label_16")
		self.logTextEdit = QtWidgets.QTextEdit(Form)
		self.logTextEdit.setGeometry(QtCore.QRect(10, 510, 821, 161))
		self.logTextEdit.setObjectName("logTextEdit")
		self.numThreadsLineEdit = QtWidgets.QLineEdit(Form)
		self.numThreadsLineEdit.setGeometry(QtCore.QRect(380, 460, 101, 21))
		self.numThreadsLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.numThreadsLineEdit.setObjectName("numThreadsLineEdit")
		self.numThreadsLabel = QtWidgets.QLabel(Form)
		self.numThreadsLabel.setGeometry(QtCore.QRect(380, 440, 101, 16))
		self.numThreadsLabel.setObjectName("numThreadsLabel")
		self.label_4 = QtWidgets.QLabel(Form)
		self.label_4.setGeometry(QtCore.QRect(520, 10, 101, 16))
		self.label_4.setObjectName("label_4")
		self.frame = QtWidgets.QFrame(Form)
		self.frame.setGeometry(QtCore.QRect(520, 30, 311, 391))
		self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
		self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
		self.frame.setObjectName("frame")
		self.groupBox = QtWidgets.QGroupBox(self.frame)
		self.groupBox.setGeometry(QtCore.QRect(10, 10, 291, 171))
		self.groupBox.setObjectName("groupBox")
		self.o_numReadsLineEdit = QtWidgets.QLineEdit(self.groupBox)
		self.o_numReadsLineEdit.setEnabled(True)
		self.o_numReadsLineEdit.setGeometry(QtCore.QRect(10, 50, 113, 23))
		self.o_numReadsLineEdit.setReadOnly(True)
		self.o_numReadsLineEdit.setObjectName("o_numReadsLineEdit")
		self.label_5 = QtWidgets.QLabel(self.groupBox)
		self.label_5.setGeometry(QtCore.QRect(10, 30, 101, 16))
		self.label_5.setObjectName("label_5")
		self.o_avQualLineEdit = QtWidgets.QLineEdit(self.groupBox)
		self.o_avQualLineEdit.setEnabled(True)
		self.o_avQualLineEdit.setGeometry(QtCore.QRect(160, 50, 113, 23))
		self.o_avQualLineEdit.setReadOnly(True)
		self.o_avQualLineEdit.setObjectName("o_avQualLineEdit")
		self.label_6 = QtWidgets.QLabel(self.groupBox)
		self.label_6.setGeometry(QtCore.QRect(160, 30, 101, 16))
		self.label_6.setObjectName("label_6")
		self.o_estCoverageLineEdit = QtWidgets.QLineEdit(self.groupBox)
		self.o_estCoverageLineEdit.setEnabled(True)
		self.o_estCoverageLineEdit.setGeometry(QtCore.QRect(10, 120, 113, 23))
		self.o_estCoverageLineEdit.setReadOnly(True)
		self.o_estCoverageLineEdit.setObjectName("o_estCoverageLineEdit")
		self.label_7 = QtWidgets.QLabel(self.groupBox)
		self.label_7.setGeometry(QtCore.QRect(10, 100, 101, 16))
		self.label_7.setObjectName("label_7")
		self.groupBox_2 = QtWidgets.QGroupBox(self.frame)
		self.groupBox_2.setGeometry(QtCore.QRect(10, 210, 291, 171))
		self.groupBox_2.setObjectName("groupBox_2")
		self.r_numReadsLineEdit_3 = QtWidgets.QLineEdit(self.groupBox_2)
		self.r_numReadsLineEdit_3.setEnabled(True)
		self.r_numReadsLineEdit_3.setGeometry(QtCore.QRect(10, 50, 113, 23))
		self.r_numReadsLineEdit_3.setReadOnly(True)
		self.r_numReadsLineEdit_3.setObjectName("r_numReadsLineEdit_3")
		self.label_17 = QtWidgets.QLabel(self.groupBox_2)
		self.label_17.setGeometry(QtCore.QRect(10, 30, 101, 16))
		self.label_17.setObjectName("label_17")
		self.r_avQualLineEdit = QtWidgets.QLineEdit(self.groupBox_2)
		self.r_avQualLineEdit.setEnabled(True)
		self.r_avQualLineEdit.setGeometry(QtCore.QRect(160, 50, 113, 23))
		self.r_avQualLineEdit.setReadOnly(True)
		self.r_avQualLineEdit.setObjectName("r_avQualLineEdit")
		self.label_18 = QtWidgets.QLabel(self.groupBox_2)
		self.label_18.setGeometry(QtCore.QRect(160, 30, 101, 16))
		self.label_18.setObjectName("label_18")
		self.r_estCoverageLineEdit = QtWidgets.QLineEdit(self.groupBox_2)
		self.r_estCoverageLineEdit.setEnabled(True)
		self.r_estCoverageLineEdit.setGeometry(QtCore.QRect(10, 120, 113, 23))
		self.r_estCoverageLineEdit.setReadOnly(True)
		self.r_estCoverageLineEdit.setObjectName("r_estCoverageLineEdit")
		self.label_19 = QtWidgets.QLabel(self.groupBox_2)
		self.label_19.setGeometry(QtCore.QRect(10, 100, 101, 16))
		self.label_19.setObjectName("label_19")
		self.frame_2 = QtWidgets.QFrame(Form)
		self.frame_2.setGeometry(QtCore.QRect(10, 110, 501, 211))
		self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
		self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
		self.frame_2.setObjectName("frame_2")
		self.readsFileLineEdit = QtWidgets.QLineEdit(self.frame_2)
		self.readsFileLineEdit.setGeometry(QtCore.QRect(10, 30, 351, 21))
		self.readsFileLineEdit.setObjectName("readsFileLineEdit")
		self.readsFileButton = QtWidgets.QPushButton(self.frame_2)
		self.readsFileButton.setGeometry(QtCore.QRect(380, 30, 113, 21))
		self.readsFileButton.setObjectName("readsFileButton")
		self.qualityStatsButton = QtWidgets.QPushButton(self.frame_2)
		self.qualityStatsButton.setGeometry(QtCore.QRect(10, 70, 241, 21))
		self.qualityStatsButton.setObjectName("qualityStatsButton")
		self.filterButton = QtWidgets.QPushButton(self.frame_2)
		self.filterButton.setGeometry(QtCore.QRect(140, 120, 113, 21))
		self.filterButton.setObjectName("filterButton")
		self.qualityLineEdit = QtWidgets.QLineEdit(self.frame_2)
		self.qualityLineEdit.setGeometry(QtCore.QRect(10, 120, 101, 21))
		self.qualityLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.qualityLineEdit.setObjectName("qualityLineEdit")
		self.label_9 = QtWidgets.QLabel(self.frame_2)
		self.label_9.setGeometry(QtCore.QRect(10, 100, 101, 16))
		self.label_9.setObjectName("label_9")
		self.label_8 = QtWidgets.QLabel(self.frame_2)
		self.label_8.setGeometry(QtCore.QRect(10, 10, 161, 16))
		self.label_8.setObjectName("label_8")
		self.label_10 = QtWidgets.QLabel(self.frame_2)
		self.label_10.setGeometry(QtCore.QRect(10, 150, 161, 16))
		self.label_10.setObjectName("label_10")
		self.gapClosingReadsButton = QtWidgets.QPushButton(self.frame_2)
		self.gapClosingReadsButton.setGeometry(QtCore.QRect(380, 170, 113, 21))
		self.gapClosingReadsButton.setObjectName("gapClosingReadsButton")
		self.gapClosingReadsLineEdit = QtWidgets.QLineEdit(self.frame_2)
		self.gapClosingReadsLineEdit.setGeometry(QtCore.QRect(10, 170, 351, 21))
		self.gapClosingReadsLineEdit.setObjectName("gapClosingReadsLineEdit")
		self.numThreadsLabel_2 = QtWidgets.QLabel(Form)
		self.numThreadsLabel_2.setGeometry(QtCore.QRect(260, 440, 101, 16))
		self.numThreadsLabel_2.setObjectName("numThreadsLabel_2")
		self.homologyLineEdit = QtWidgets.QLineEdit(Form)
		self.homologyLineEdit.setGeometry(QtCore.QRect(260, 460, 101, 21))
		self.homologyLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.homologyLineEdit.setObjectName("homologyLineEdit")
		self.prefixLineEdit = QtWidgets.QLineEdit(Form)
		self.prefixLineEdit.setGeometry(QtCore.QRect(510, 460, 321, 21))
		self.prefixLineEdit.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.prefixLineEdit.setObjectName("prefixLineEdit")
		self.numThreadsLabel_3 = QtWidgets.QLabel(Form)
		self.numThreadsLabel_3.setGeometry(QtCore.QRect(510, 440, 171, 16))
		self.numThreadsLabel_3.setObjectName("numThreadsLabel_3")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)


		self.projectFolderButton.clicked.connect(self.selectProjectFolder)
		self.readsFileButton.clicked.connect(self.selectReadFile)
		self.referenceButton.clicked.connect(self.selectReferenceFile)
		self.runButton.clicked.connect(self.runAssembly)
		self.exitButton.clicked.connect(self.exitProgram)
		self.qualityStatsButton.clicked.connect(self.calculateStatistics)
		self.filterButton.clicked.connect(self.filterReads)
		self.gapClosingReadsButton.clicked.connect(self.selectGapClosingReadsFile)
		

	def filterReads(self):
		self.logTextEdit.append("* * Extracting high quality read portions....\n")
		self.logTextEdit.repaint()

		
		inputFile = self.readsFileLineEdit.text()
		threshold = self.qualityLineEdit.text()
		minLen = 150

		numSeq = 0
		totSequences = 0
		outputFolder = self.projectFolderLineEdit.text()
		logFile = open(outputFolder+"/"+self.prefixLineEdit.text()+"_readsFiltering.log","w")
		startTime = datetime.now()
		current_time = startTime.strftime("%H:%M:%S")
		logFile.write("Reads filtering started at "+str(current_time)+"\n\n")
		totNumBases = 0
		outfile = open(outputFolder+"/hq_reads.fasta","w")
		for seq_record in SeqIO.parse(inputFile,"fastq"):
			numSeq+=1
			totSequences+=1
			if totSequences%3000 == 0:
				self.logTextEdit.append("* * * "+str(totSequences)+" reads analyzed....")
				self.logTextEdit.repaint()
			sequence = str(seq_record.seq)
			totNumBases+=len(sequence)

			quality = seq_record.letter_annotations["phred_quality"]
			maskedSeq = ""
			qualityValues = []
			for a in range(len(quality)):
				qualityValues.append(float(quality[a]))
				if quality[a]>int(threshold):
					maskedSeq+=sequence[a]
				else:
					if len(maskedSeq)>int(minLen):
						outfile.write(">HQSequence_"+str(numSeq)+"\n"+maskedSeq+"\n")
					numSeq+=1
					maskedSeq=""
					qualityValues = []
			if len(qualityValues)==len(quality):
				outfile.write(">HQSequence_"+str(numSeq)+"\n"+maskedSeq+"\n")
				numSeq+=1
				maskedSeq=""

		outfile.close()
		self.logTextEdit.append("\n\nJob finished!")
		self.logTextEdit.repaint()
		endTime = datetime.now()
		current_time = endTime.strftime("%H:%M:%S")
		logFile.write("Reads filtering finished at "+str(current_time)+"\n\n")
		totalTime = (endTime -startTime).total_seconds()
		logFile.write("Total processing time: "+str(totalTime)+" seconds")
		logFile.close()


		

	
	def calculateStatistics(self):
		inputFile = self.readsFileLineEdit.text()
		threshold = self.qualityLineEdit.text()
		minLen = 150

		self.logTextEdit.append("Read statistics calculation started")
		self.logTextEdit.repaint()
		outputFolder = self.projectFolderLineEdit.text()

		logFile = open(outputFolder+"/"+self.prefixLineEdit.text()+"_reads_Quality_Statistics.log","w")
		startTime = datetime.now()
		current_time = startTime.strftime("%H:%M:%S")
		logFile.write("Reads quality stats elaboration started at "+str(current_time)+"\n\n")



		#Load reference genome in memory
		self.logTextEdit.append("Loading reference....")
		self.logTextEdit.repaint()
		refFile = self.referenceLineEdit.text()
		for seq_record in SeqIO.parse(refFile,"fasta"):
			refSeq = str(seq_record.seq)
		
		#Load the reads in memory
		self.logTextEdit.append("Calculating original read stats....")
		self.logTextEdit.repaint()
		
		numSeq = 0
		totSequences = 0
		
		totNumBases = 0
		inputSequences = {}
		
		qualityValues = []
		for seq_record in SeqIO.parse(inputFile,"fastq"):
			if not str(seq_record.id) in inputSequences:
				inputSequences[str(seq_record.id)] = seq_record
			numSeq+=1
			totSequences+=1
			
			if totSequences%5000 == 0:
				self.logTextEdit.append("* "+str(totSequences)+" reads analyzed....")
				self.logTextEdit.repaint()
			sequence = str(seq_record.seq)
			totNumBases+=len(sequence)

			quality = seq_record.letter_annotations["phred_quality"]
			if len(qualityValues)<1000000:
				qualityValues+=quality
			
				
		

		self.logTextEdit.append("Calculating stats....")
		self.logTextEdit.repaint()
		averageQuality = str(int(np.mean(qualityValues)))
		averageCoverage = str( int(float(totNumBases) / float(len(refSeq))  )  )+" X"
		self.o_numReadsLineEdit.setText(str(totSequences))
		self.o_avQualLineEdit.setText(averageQuality)
		self.o_estCoverageLineEdit.setText(averageCoverage)
		logFile.write("Original reads\n")
		logFile.write("Total number of reads: "+str(totSequences)+"\n")
		logFile.write("Average quality: "+averageQuality+"\n")
		logFile.write("Average coverage: "+averageCoverage+"\n\n")
		fig = plt.figure(figsize=(10,10))
		plt.hist(qualityValues,bins=200,label="All reads",color="red",alpha=0.5)
		
	


		self.logTextEdit.append("Calculating reference homologous read stats....")
		self.logTextEdit.append("* Aligning reads....")
		self.logTextEdit.repaint()

		numSeq = 0
		totSequences = 0
		qualityValues = []
		totNumBases = 0
		os.system(installationDirectory+"/src/conda/bin/minimap2 -x map-pb -t "+self.numThreadsLineEdit.text()+" "+refFile+" "+inputFile+" > "+outputFolder+"/outputMinimap")


		outfile = open(outputFolder+"/refSpecificReads.fastq","w")
		infile = open(outputFolder+"/outputMinimap")
		while True:
			line = infile.readline().rstrip()
			if not line:
				break
			fields = line.split("\t")
			if (float(fields[3])-float(fields[2]))/float(fields[1]) >float(self.homologyLineEdit.text()): 
				SeqIO.write(inputSequences[fields[0]],outfile,"fastq")
		outfile.close()

		inputFile = outputFolder+"/refSpecificReads.fastq"
		for seq_record in SeqIO.parse(inputFile,"fastq"):

			numSeq+=1
			totSequences+=1
			
			if totSequences%5000 == 0:
				self.logTextEdit.append("* "+str(totSequences)+" reads analyzed....")
				self.logTextEdit.repaint()
			sequence = str(seq_record.seq)
			totNumBases+=len(sequence)

			quality = seq_record.letter_annotations["phred_quality"]
			if len(qualityValues)<1000000:
				qualityValues+=quality
			
				
		outfile.close()
		averageQuality = str(int(np.mean(qualityValues)))
		averageCoverage = str( int(float(totNumBases) / float(len(refSeq))  )  )+" X"
		self.logTextEdit.append("Calculating stats....")
		self.logTextEdit.repaint()
		self.r_numReadsLineEdit_3.setText(str(totSequences))
		self.r_avQualLineEdit.setText(averageQuality)
		self.r_estCoverageLineEdit.setText(averageCoverage)
		plt.hist(qualityValues,bins=200,label="Reference mapping reads",color="blue",alpha=0.5)
		plt.xlabel("Quality phred score")
		plt.legend()
		fileName = (self.readsFileLineEdit.text().split("/"))[-1]
		fig.savefig(outputFolder+"/"+fileName+"_qualityDist.png",dpi=300)
		logFile.write("Reference homologues reads\n")
		logFile.write("Total number of reads: "+str(totSequences)+"\n")
		logFile.write("Average quality: "+averageQuality+"\n")
		logFile.write("Average coverage: "+averageCoverage+"\n")
		self.logTextEdit.append("\n\nJob finished!")
		self.logTextEdit.repaint()

		endTime = datetime.now()
		current_time = endTime.strftime("%H:%M:%S")
		logFile.write("Reads quality stats elaboration finished at "+str(current_time)+"\n\n")
		totalTime = (endTime -startTime).total_seconds()
		logFile.write("Total processing time: "+str(totalTime)+" seconds")
		logFile.close()




	def selectProjectFolder(self):
		folderName = QFileDialog.getExistingDirectory(None, "Select project folder","./")
		self.projectFolderLineEdit.setText(folderName)

		if os.path.isdir(folderName) == False:
			os.system("mkdir "+folderName)

	def selectGapClosingReadsFile(self):
		filename, __ = QFileDialog.getOpenFileName(None,"Select read file","./")
		self.gapClosingReadsLineEdit.setText(filename)


	def selectReadFile(self):
		filename, __ = QFileDialog.getOpenFileName(None,"Select read file","./")
		self.readsFileLineEdit.setText(filename)


	def selectReferenceFile(self):
		filename, __ = QFileDialog.getOpenFileName(None,"Select reference file","./")
		self.referenceLineEdit.setText(filename)

	def exitProgram(self):
		exit(1)

	def chopReads(self,inputFile):
		numSeq = 0
		outfile = open(inputFile+"_chopped.fasta","w")
		for seq_record in SeqIO.parse(inputFile,"fasta"):
			
			sequence = str(seq_record.seq)
			for a in range(0,len(sequence)-150,+150):
				numSeq+=1
				outfile.write(">ChoppedSeq_"+str(numSeq)+"\n"+sequence[a:a+150]+"\n")
		outfile.close()



	# ****************************************************
	# ***************** Main algorithm *******************
	# ****************************************************
	
	#Check required fields
	def runAssembly(self):
		
		
		self.logTextEdit.append("Reference guided de novo assembly started\n")
		self.logTextEdit.repaint()

		if os.path.isdir(self.projectFolderLineEdit.text()) == "":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid project folder was not selected")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should select a folder where all the intermediary file will be stored.\n ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		if os.path.isfile(self.readsFileLineEdit.text()) == False:
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid read file was not selected")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should select a fastq file with the PacBio reads.\n ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		if os.path.isfile(self.referenceLineEdit.text()) == False:
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid reference file was not selected")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should select a fasta formatted file reference genome.\n ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		outputFolder = self.projectFolderLineEdit.text()
		logFile = open(outputFolder+"/"+self.prefixLineEdit.text()+"_assembly.log","w")
		startTime = datetime.now()
		current_time = startTime.strftime("%H:%M:%S")
		logFile.write("De novo assembly started at "+str(current_time)+"\n\n")

		#Load reference genome in memory
		refFile = self.referenceLineEdit.text()
		for seq_record in SeqIO.parse(refFile,"fasta"):
			refSeq = str(seq_record.seq)

		
		reads = self.readsFileLineEdit.text()
		if ".fastq" in reads or ".fq" in reads:
			self.logTextEdit.append("* * Converting input file from fastq to fasta....\n")
			self.logTextEdit.repaint()
			
			os.system(installationDirectory+"/src/conda/bin/fq2fa "+self.readsFileLineEdit.text()+" "+outputFolder+'/originalReads.fasta')

			reads = outputFolder+'/originalReads.fasta'
		
		
		#Perform reference guided de novo assembly
		if reads == "":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("Something did not work with the quality filtering step!")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("The quality filtering step did not produce the expected hq_reads.fasta file in the output folder\n ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		else:
			
			readsSeq = {}
			windowSize = int(self.windowSizeLineEdit.text())
			windowStep = int(self.windowStepLineEdit.text())

			self.logTextEdit.append("* Reference guided de novo assembly")
			self.logTextEdit.append("* * Loading reads in memory")
			self.logTextEdit.repaint()

			#Loading high quality reads in memory
			for seq_record in SeqIO.parse(reads,"fasta"):
				if not str(seq_record.id) in readsSeq:
					readsSeq[str(seq_record.id)] = str(seq_record.seq)

			

			stage_a = open(outputFolder+"/local_assemblies.fasta","w")




			self.logTextEdit.append("* * Assembly on sliding windows started")
			self.logTextEdit.repaint()
			now = datetime.now()
			current_time = now.strftime("%H:%M:%S")
			logFile.write("Sliding window assembly started at "+str(current_time)+"\n\n")
			logFile.write("Range\tContig_size\n")
			lastWindow = False
			for a in range(0,len(refSeq),+windowStep):
				if lastWindow==True:
					break
			#for a in range(1):
				endPos = a+windowSize
				if endPos>len(refSeq):
					endPos=len(refSeq)
					windowSize = len(refSeq) - a
					lastWindow = True
				windowSuffix = str(a)+"-"+str(endPos) 
				self.logTextEdit.append("* * * Assembling region "+str(a)+"-"+str(endPos))
				self.logTextEdit.repaint()
				partSeq = refSeq[a:endPos]
				tempFasta = open(outputFolder+"/partReference.fasta","w")
				tempFasta.write(">partReference\n"+partSeq+"\n")
				tempFasta.close()

				os.system(installationDirectory+"/src/conda/bin/minimap2 -x map-pb -t "+self.numThreadsLineEdit.text()+" "+outputFolder+"/partReference.fasta "+reads+" > "+outputFolder+"/outputMinimap_"+windowSuffix)

				homology = self.homologyLineEdit.text()
				os.system("awk '($11/$2)>"+homology+"' "+outputFolder+"/outputMinimap_"+windowSuffix+" | sort -k2rn,2rn >  "+outputFolder+"/outputMinimap_filtered_"+windowSuffix)

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
				os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/hqKmerAssembly.py -p "+installationDirectory+" -r "+outputFolder+"/mapped.fasta_"+windowSuffix+" -ref "+outputFolder+"/partReference.fasta -t "+self.numThreadsLineEdit.text()+" -of "+outputFolder)

				
					
				

				if os.path.isfile(outputFolder+"/mapped.fasta_"+windowSuffix+"_localAssembly.fasta") == True:
					for seq_record in SeqIO.parse(outputFolder+"/mapped.fasta_"+windowSuffix+"_localAssembly.fasta","fasta"):

						maxScaffoldLength = len(str(seq_record.seq))
						longestContig = str(seq_record.seq)

					self.logTextEdit.append("* * * Contig size: "+str(maxScaffoldLength))
					self.logTextEdit.repaint()
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

				self.logTextEdit.append("* * Joining contigs.... ")
				self.logTextEdit.repaint()
				os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/contigsJoiner.py -c "+outputFolder+"/local_assemblies.fasta -r "+refFile+" -p "+installationDirectory+" -o "+outputFolder)

				#Check final number of scaffold and attempt gap closure if > 1
				self.logTextEdit.append("* * Attempting gap closure.... ")
				self.logTextEdit.repaint()
				finalScaffols = SeqIO.to_dict(SeqIO.parse(outputFolder+"/scaffolds.fasta","fasta"))

				if len(finalScaffols)>1:
					self.logTextEdit.append("\nAttempting gap closure.... ")
					self.logTextEdit.repaint()
					now = datetime.now()
					current_time = now.strftime("%H:%M:%S")
					logFile.write("Gap closure started at "+str(current_time)+"\n\n")
					os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/oneReadContigsJoiner.py \
						-p "+installationDirectory+ " -c "+outputFolder+"/scaffolds.fasta -r "+refFile+" -x " + \
							outputFolder+" -s "+ self.gapClosingReadsLineEdit.text() +" -o scaffolds_gapClosed.fasta")
				else:
					os.system("cp "+outputFolder+"/scaffolds.fasta "+outputFolder+"/scaffolds_gapClosed.fasta")
			
			else:
				os.system("cp "+outputFolder+"/local_assemblies.fasta "+outputFolder+"/scaffolds_gapClosed.fasta")
			
			
			#Final alignment and consensus calling
			self.logTextEdit.append("* * Calling consensus.... ")
			#self.logTextEdit.append("* * * Subsampling.... ")
			self.logTextEdit.repaint()
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
			outfile.close()"""

			self.logTextEdit.append("* * * Assembly correction ")
			self.logTextEdit.append("* * * Chopping reads.... ")
			self.logTextEdit.repaint()
			self.chopReads(reads)
			self.logTextEdit.append("* * * First assembly correction.... ")
			self.logTextEdit.append("* * * Mapping original reads to the assembled sequence.... ")
			self.logTextEdit.repaint()
			
			
			os.system(installationDirectory+"/src/conda/bin/minimap2 -a -x map-pb -t "+self.numThreadsLineEdit.text()+" "+outputFolder+"/scaffolds_gapClosed.fasta "+reads+"_chopped.fasta"+" > "+outputFolder+"/alignment1.sam")
			self.logTextEdit.append("* * * Converting sam to bam.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/samtools view -F 4 -bS -h "+outputFolder+"/alignment1.sam > "+outputFolder+"/alignment1.bam")
			self.logTextEdit.append("* * * Sorting.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/samtools sort -o "+outputFolder+"/alignment_sorted1.bam "+outputFolder+"/alignment1.bam")
			self.logTextEdit.append("* * * Indexing.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/samtools index "+outputFolder+"/alignment_sorted1.bam")

			self.logTextEdit.append("* * * Creating pilleup.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/samtools mpileup -f "+outputFolder+ \
				"/scaffolds_gapClosed.fasta "+outputFolder +"/alignment_sorted1.bam > "+outputFolder+ \
					"/pileup1.txt")
			os.system("awk '$3!=\"N\"' "+outputFolder+"/pileup1.txt >"+outputFolder+"/pileup2.txt")
			os.system("mv "+outputFolder+"/pileup2.txt "+outputFolder+"/pileup1.txt")

			self.logTextEdit.append("* * * Calling variants.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/varscan mpileup2cns "+outputFolder+"/pileup1.txt --variants --output-vcf --min-var-freq 0.5 --min-avg-qual 0 --strand-filter 0 --min-coverage 5   > "+outputFolder+"/output1.vcf")
			
			os.system(installationDirectory+"/src/conda/bin/bgzip -f -c "+outputFolder+"/output1.vcf > "+outputFolder+"/output1.vcf.gz")
			os.system(installationDirectory+"/src/conda/bin/tabix -f "+outputFolder+"/output1.vcf.gz")
			os.system("cat "+outputFolder+"/scaffolds_gapClosed.fasta | "+installationDirectory+"/src/conda/bin/bcftools consensus "+outputFolder+"/output1.vcf.gz > "+outputFolder+"/finalAssembly1.fasta")



			self.logTextEdit.append("* * * Second assembly correction.... ")
			self.logTextEdit.append("* * * Mapping original reads to the assembled sequence.... ")
			self.logTextEdit.repaint()
			
			
			os.system(installationDirectory+"/src/conda/bin/minimap2 -a -x map-pb -t "+self.numThreadsLineEdit.text()+" "+outputFolder+"/finalAssembly1.fasta "+reads+"_chopped.fasta"+" > "+outputFolder+"/alignment2.sam")
			self.logTextEdit.append("* * * Converting sam to bam.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/samtools view -F 4 -bS -h "+outputFolder+"/alignment2.sam > "+outputFolder+"/alignment2.bam")
			self.logTextEdit.append("* * * Sorting.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/samtools sort -o "+outputFolder+"/alignment_sorted2.bam "+outputFolder+"/alignment2.bam")
			self.logTextEdit.append("* * * Indexing.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/samtools index "+outputFolder+"/alignment_sorted2.bam")

			self.logTextEdit.append("* * * Creating pilleup.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/samtools mpileup -f "+outputFolder+ \
				"/finalAssembly1.fasta "+outputFolder +"/alignment_sorted2.bam > "+outputFolder+ \
					"/pileup2.txt")

			os.system("awk '$3!=\"N\"' "+outputFolder+"/pileup2.txt >"+outputFolder+"/pileup3.txt")
			os.system("mv "+outputFolder+"/pileup3.txt "+outputFolder+"/pileup2.txt")

			self.logTextEdit.append("* * * Calling variants.... ")
			self.logTextEdit.repaint()
			os.system(installationDirectory+"/src/conda/bin/varscan mpileup2cns "+outputFolder+"/pileup2.txt --variants --output-vcf --min-avg-qual 0 --strand-filter 0 --min-coverage 5   > "+outputFolder+"/output2.vcf")
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/varscanFilter.py -i "+outputFolder+"/output2.vcf -o "+outputFolder+"/output_filtered2.vcf -1 "+reads+" -g 1  -r "+outputFolder+"/finalAssembly1.fasta -p "+installationDirectory +" -of "+outputFolder+" -t "+self.numThreadsLineEdit.text()) 
			os.system(installationDirectory+"/src/conda/bin/bgzip -f -c "+outputFolder+"/output_filtered2.vcf > "+outputFolder+"/output_filtered2.vcf.gz")
			os.system(installationDirectory+"/src/conda/bin/tabix -f "+outputFolder+"/output_filtered2.vcf.gz")
			os.system("cat "+outputFolder+"/finalAssembly1.fasta | "+installationDirectory+"/src/conda/bin/bcftools consensus "+outputFolder+"/output_filtered2.vcf.gz > "+outputFolder+"/finalAssembly.fasta")



			

			
			




			self.logTextEdit.append("\n\n De novo assembly finished!")
			self.logTextEdit.repaint()
			endTime = datetime.now()
			current_time = endTime.strftime("%H:%M:%S")
			logFile.write("De novo assembly finished at "+str(current_time)+"\n\n")
			totalTime = (endTime -startTime).total_seconds()
			logFile.write("Total processing time: "+str(totalTime)+" seconds")
			logFile.close()




	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Form"))
		self.label.setText(_translate("Form", "Project folder"))
		self.projectFolderButton.setText(_translate("Form", "Select"))
		self.label_2.setText(_translate("Form", "Read files"))
		self.referenceButton.setText(_translate("Form", "Select"))
		self.label_3.setText(_translate("Form", "Reference file"))
		self.label_13.setText(_translate("Form", "Window size"))
		self.windowSizeLineEdit.setText(_translate("Form", "20000"))
		self.windowStepLineEdit.setText(_translate("Form", "10000"))
		self.label_14.setText(_translate("Form", "Window step"))
		self.runButton.setText(_translate("Form", "Run"))
		self.exitButton.setText(_translate("Form", "Exit"))
		self.label_16.setText(_translate("Form", "Log "))
		self.numThreadsLineEdit.setText(_translate("Form", "8"))
		self.numThreadsLabel.setText(_translate("Form", "Num. threads"))
		self.label_4.setText(_translate("Form", "Read statistics"))
		self.groupBox.setTitle(_translate("Form", "Original"))
		self.o_numReadsLineEdit.setText(_translate("Form", "--"))
		self.label_5.setText(_translate("Form", "Read number"))
		self.o_avQualLineEdit.setText(_translate("Form", "--"))
		self.label_6.setText(_translate("Form", "Average  quality"))
		self.o_estCoverageLineEdit.setText(_translate("Form", "--"))
		self.label_7.setText(_translate("Form", "Coverage"))
		self.groupBox_2.setTitle(_translate("Form", "Reference specific"))
		self.r_numReadsLineEdit_3.setText(_translate("Form", "--"))
		self.label_17.setText(_translate("Form", "Read number"))
		self.r_avQualLineEdit.setText(_translate("Form", "--"))
		self.label_18.setText(_translate("Form", "Average  quality"))
		self.r_estCoverageLineEdit.setText(_translate("Form", "--"))
		self.label_19.setText(_translate("Form", "Coverage"))
		self.readsFileButton.setText(_translate("Form", "Select"))
		self.qualityStatsButton.setText(_translate("Form", "Plot quality score distribution"))
		self.filterButton.setText(_translate("Form", "Filter"))
		self.qualityLineEdit.setText(_translate("Form", "30"))
		self.label_9.setText(_translate("Form", "Min. quality"))
		self.label_8.setText(_translate("Form", "Assembly"))
		self.label_10.setText(_translate("Form", "Gap closing"))
		self.gapClosingReadsButton.setText(_translate("Form", "Select"))
		self.numThreadsLabel_2.setText(_translate("Form", "Homology"))
		self.homologyLineEdit.setText(_translate("Form", "0.7"))
		self.prefixLineEdit.setText(_translate("Form", "output"))
		self.numThreadsLabel_3.setText(_translate("Form", "Output files prefix"))
	

	



if __name__ == "__main__":
	import sys
	app = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	ui = Ui_Form()
	ui.setupUi(Form,installationDirectory)
	Form.show()
	sys.exit(app.exec_())