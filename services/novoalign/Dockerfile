FROM ubuntu:16.04

ENV NOVO_FILE novocraftV3.09.00.Linux3.10.0.tar.gz
ENV NOVO_LICENSE novoalign.lic

ENV NOVO_PATH /opt/novocraft
RUN mkdir $NOVO_PATH
WORKDIR $NOVO_PATH

# Copy Novocraft suite and NovoAlign license if present
COPY ./assets/$NOVO_FILE ./assets/$NOVO_LICENSE* $NOVO_PATH/

# Extract Novocraft suite to $NOVO_PATH
RUN tar -xzf $NOVO_FILE -C /opt
WORKDIR /

# Add to PATH
ENV PATH $NOVO_PATH:$PATH

# Install USeq in order to build the reference genome from annotations, see
# http://www.novocraft.com/documentation/novoalign-2/novoalign-user-guide/rnaseq-analysis-mrna-and-the-spliceosome/

# Install dependencies
RUN apt-get update
RUN apt-get install -y wget unzip default-jre

# Install R 3.5
RUN apt-get install -y apt-transport-https
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN apt-get update
RUN apt-get install -y r-base

# Install required R packages
RUN R -e "install.packages(c('DESeq', 'BiocManager'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('qvalue')"

# Download USeq
ENV DST_PATH /opt/useq
ENV USEQ_VESION 9.1.9
ENV USEQ_ZIP_NAME USeq_$USEQ_VESION.zip
RUN wget https://github.com/HuntsmanCancerInstitute/USeq/releases/download/USeq_$USEQ_VESION/$USEQ_ZIP_NAME -P $DST_PATH
RUN unzip $DST_PATH/$USEQ_ZIP_NAME -d $DST_PATH
RUN rm $DST_PATH/$USEQ_ZIP_NAME
RUN mv $DST_PATH/USeq_$USEQ_VESION/* $DST_PATH
RUN rm -r $DST_PATH/USeq_$USEQ_VESION

# Download FASTA references per chromosome
RUN wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/* -P references/hg19
RUN wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/* -P references/hg38

# Copy annotations
COPY ./assets/refFlat.hg19.txt.gz ./annotations/hg19.txt.gz
RUN gunzip ./annotations/hg19.txt.gz
COPY ./assets/refFlat.hg38.txt.gz ./annotations/hg38.txt.gz
RUN gunzip ./annotations/hg38.txt.gz

COPY ./assets/build_genome.sh ./build_genome.sh
RUN chmod +x ./build_genome.sh
