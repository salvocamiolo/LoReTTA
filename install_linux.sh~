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


	cd src
	bash Miniconda3-latest-Linux-x86_64.sh -b -p ./conda
	cd ../
	if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 successfully installed"
	./src/conda/bin/conda install -c bioconda -y  lastz=1.0.4
	./src/conda/bin/conda install -c bioconda -y  blast=2.9.0
	./src/conda/bin/conda install -c bioconda -y  samtools=1.3.1
	./src/conda/bin/conda install -c anaconda -y pyqt=5.9.2
	./src/conda/bin/conda install -c bioconda -y  biopython=1.76
	./src/conda/bin/conda install -c bioconda -y bowtie2
	./src/conda/bin/conda install -c bioconda -y minimap2
	./src/conda/bin/conda install -c bioconda -y varscan
	./src/conda/bin/conda install -c bioconda -y tabix
	./src/conda/bin/conda install -c bioconda -y bcftools
	./src/conda/bin/conda install -y matplotlib
	./src/conda/bin/conda install -c yuxiang -y bam2fastq
	./src/conda/bin/conda install -c bioconda -y spades=3.11
	./src/conda/bin/conda install -c bioconda -y kmc
	./src/conda/bin/conda install -c bioconda -y idba
	
	fi
fi


if  test -f "./src/conda/bin/conda2"; then
	echo "Miniconda2 already installed"
else


	cd src
	bash Miniconda2-latest-Linux-x86_64.sh -b -p ./conda2
	cd ../
	if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 successfully installed"
	./src/conda2/bin/conda install -c bioconda -y  spades=3.10.1
	
	
	fi
fi



