installationPath=$(pwd)
echo "#!"$installationPath"/src/conda/bin/python" >LoReTTA
echo "installationDirectory = \""$installationPath"/\"" >> LoReTTA
echo " " >> LoReTTA
cat LoReTTA ./src/.mainScript.py >temp ; mv temp LoReTTA
chmod +x LoReTTA

installationPath=$(pwd)
echo "#!"$installationPath"/src/conda/bin/python" >LoReTTA_cl
echo "installationDirectory = \""$installationPath"/\"" >> LoReTTA_cl
echo " " >> LoReTTA_cl
cat LoReTTA_cl ./src/.clScript.py >temp ; mv temp LoReTTA_cl
chmod +x LoReTTA_cl



if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 already installed"
else
	echo "Installing miniconda 3"
	cd src
	bash Miniconda3-latest-Linux-x86_64.sh -b -p ./conda >> installation.log 2>&1
	cd ../
	if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 successfully installed"
	fi
fi



./src/conda/bin/conda config --set notify_outdated_conda false

echo ""
echo "Checking lastz...."
./src/conda/bin/conda list > condaList
if grep -Fq lastz condaList; then
		echo "lastz is already installed"
else
	echo "Installing lastz. Please wait...."
    ./src/conda/bin/conda install -c bioconda -y  lastz=1.0.4 >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq lastz condaList; then
		echo "lastz was successfully installed"
	else
		echo "lastz was not installed. Please check file installation.log for details"
	fi
fi

echo ""
echo "Checking blast...."
./src/conda/bin/conda list > condaList
if grep -Fq blast condaList; then
		echo "blast is already installed"
else
	echo "Installing blast. Please wait...."
    ./src/conda/bin/conda install -c bioconda -y  blast=2.9.0 >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq blast condaList; then
		echo "blast was successfully installed"
	else
		echo "blast was not installed. Please check file installation.log for details"
	fi
fi

echo ""
echo "Checking samtools...."
./src/conda/bin/conda list > condaList
if grep -Fq samtools condaList; then
		echo "samtools is already installed"
else
	echo "Installing samtools. Please wait...."
    ./src/conda/bin/conda install -c bioconda -y  samtools=1.3.1 >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq samtools condaList; then
		echo "samtools was successfully installed"
	else
		echo "samtools was not installed. Please check file installation.log for details"
	fi
fi

echo ""
echo "Checking pyqt...."
./src/conda/bin/conda list > condaList
if grep -Fq pyqt condaList; then
		echo "pyqt is already installed"
else
	echo "Installing pyqt. Please wait...."
    ./src/conda/bin/conda install -c anaconda -y pyqt=5.9.2 >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq pyqt condaList; then
		echo "pyqt was successfully installed"
	else
		echo "pyqt was not installed. Please check file installation.log for details"
	fi
fi

echo ""
echo "Checking biopython...."
./src/conda/bin/conda list > condaList
if grep -Fq biopython condaList; then
		echo "biopython is already installed"
else
	echo "Installing biopython. Please wait...."
    ./src/conda/bin/conda install -c bioconda -y  biopython=1.76 >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq biopython condaList; then
		echo "biopython was successfully installed"
	else
		echo "biopython was not installed. Please check file installation.log for details"
	fi
fi

echo ""
echo "Checking bowtie2...."
./src/conda/bin/conda list > condaList
if grep -Fq bowtie2 condaList; then
		echo "bowtie2 is already installed"
else
	echo "Installing bowtie2. Please wait...."
    ./src/conda/bin/conda install -c bioconda -y bowtie2 >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq bowtie2 condaList; then
		echo "bowtie2 was successfully installed"
	else
		echo "bowtie2 was not installed. Please check file installation.log for details"
	fi
fi	
	
echo ""
echo "Checking minimap2...."
./src/conda/bin/conda list > condaList
if grep -Fq minimap2 condaList; then
		echo "minimap2 is already installed"
else
	echo "Installing minimap2. Please wait...."
    ./src/conda/bin/conda install -c bioconda -y minimap2 >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq minimap2 condaList; then
		echo "minimap2 was successfully installed"
	else
		echo "minimap2 was not installed. Please check file installation.log for details"
	fi
fi	

echo ""
echo "Checking varscan...."
./src/conda/bin/conda list > condaList
if grep -Fq varscan condaList; then
		echo "varscan is already installed"
else
	echo "Installing varscan. Please wait...."
    ./src/conda/bin/conda install -c bioconda -y varscan >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq varscan condaList; then
		echo "varscan was successfully installed"
	else
		echo "varscan was not installed. Please check file installation.log for details"
	fi
fi	

echo ""
echo "Checking tabix...."
./src/conda/bin/conda list > condaList
if grep -Fq tabix condaList; then
		echo "tabix is already installed"
else
	echo "Installing tabix. Please wait...."
    ./src/conda/bin/conda install -c bioconda -y tabix >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq tabix condaList; then
		echo "tabix was successfully installed"
	else
		echo "tabix was not installed. Please check file installation.log for details"
	fi
fi	
	
echo ""
echo "Checking bcftools...."
./src/conda/bin/conda list > condaList
if grep -Fq bcftools condaList; then
		echo "bcftools is already installed"
else
	echo "Installing bcftools. Please wait...."
    ./src/conda/bin/conda install -c bioconda -y bcftools >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq bcftools condaList; then
		echo "bcftools was successfully installed"
	else
		echo "bcftools was not installed. Please check file installation.log for details"
	fi
fi	
	

echo ""
echo "Checking matplotlib...."
./src/conda/bin/conda list > condaList
if grep -Fq matplotlib condaList; then
		echo "matplotlib is already installed"
else
	echo "Installing matplotlib. Please wait...."
    ./src/conda/bin/conda install -y matplotlib >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq matplotlib condaList; then
		echo "matplotlib was successfully installed"
	else
		echo "matplotlib was not installed. Please check file installation.log for details"
	fi
fi	

echo ""
echo "Checking bam2fastq...."
./src/conda/bin/conda list > condaList
if grep -Fq bam2fastq condaList; then
		echo "bam2fastq is already installed"
else
	echo "Installing bam2fastq. Please wait...."
    ./src/conda/bin/conda install -c yuxiang -y bam2fastq >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq bam2fastq condaList; then
		echo "bam2fastq was successfully installed"
	else
		echo "bam2fastq was not installed. Please check file installation.log for details"
	fi
fi	
	
echo ""
echo "Checking spades...."
./src/conda/bin/conda list > condaList
if grep -Fq spades condaList; then
		echo "spades is already installed"
else
	echo "Installing spades. Please wait...."
    ./src/conda/bin/conda install -c bioconda -y spades=3.11 >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq spades condaList; then
		echo "spades was successfully installed"
	else
		echo "spades was not installed. Please check file installation.log for details"
	fi
fi	

echo ""
echo "Checking kmc...."
./src/conda/bin/conda list > condaList
if grep -Fq kmc condaList; then
		echo "kmc is already installed"
else
	echo "Installing kmc. Please wait...."
    ./src/conda/bin/conda install -c bioconda -y kmc >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq kmc condaList; then
		echo "kmc was successfully installed"
	else
		echo "kmc was not installed. Please check file installation.log for details"
	fi
fi	

echo ""
echo "Checking idba...."
./src/conda/bin/conda list > condaList
if grep -Fq idba condaList; then
		echo "idba is already installed"
else
	echo "Installing idba. Please wait...."
    ./src/conda/bin/conda install -c bioconda -y idba >> installation.log
	./src/conda/bin/conda list > condaList
	if grep -Fq idba condaList; then
		echo "idba was successfully installed"
	else
		echo "idba was not installed. Please check file installation.log for details"
	fi
fi	
	


echo ""
if  test -f "./src/conda2/bin/conda"; then
	echo "Miniconda2 already installed"
else
	cd src
	echo "Installing miniconda 2"
	bash Miniconda2-latest-Linux-x86_64.sh -b -p ./conda2 >> installation.log 2>&1
	cd ../
	if  test -f "./src/conda2/bin/conda"; then
	echo "Miniconda2 successfully installed"
	fi
fi

echo ""
echo "Checking spades_3.10...."
./src/conda2/bin/conda list > condaList
if grep -Fq spades condaList; then
		echo "spades_3.10 is already installed"
else
	echo "Installing spades_3.10. Please wait...."
    ./src/conda2/bin/conda install -c bioconda -y  spades=3.10.1 >> installation.log
	./src/conda2/bin/conda list > condaList
	if grep -Fq spades condaList; then
		echo "spades_3.10 was successfully installed"
	else
		echo "spades_3.10 was not installed. Please check file installation.log for details"
	fi
fi



	
	




