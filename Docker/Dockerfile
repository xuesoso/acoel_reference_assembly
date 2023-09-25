FROM ubuntu:20.04

#### Install dependency packages ####
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y build-essential git zlib1g-dev libssl-dev libbz2-dev liblzma-dev cmake wget


#### Set tool directory ####
WORKDIR /tools

#### Download RATTLE ####
#### We are checking out the specific version of "3c1d21bb46" as it was the version that generated our transcriptome
RUN git clone --recurse-submodules https://github.com/comprna/RATTLE.git && cd /tools/RATTLE/ && git checkout 3c1d21bb46
RUN cd /tools/RATTLE/ && bash ./build.sh && mv rattle /usr/bin/

#### Download and install seqtk ####
RUN git clone https://github.com/lh3/seqtk.git && cd seqtk && git checkout c91a6afd5d9fa51c5069315a3a89c3d4b246fc22 && make && mv seqtk /usr/bin/

#### Install minimap version 2.17-r941 ####
# RUN cd /tools/ && git clone https://github.com/lh3/seqtk.git && cd seqtk && git checkout c91a6afd5d9fa51c5069315a3a89c3d4b246fc22 && make && mv seqtk /usr/bin/
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 && tar -xvf minimap2-2.17_x64-linux.tar.bz2 && \
rm minimap2-2.17_x64-linux.tar.bz2
RUN cd /tools/minimap2-2.17_x64-linux && mv minimap2 /usr/bin/

#### Install samtools version 1.10 ####
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && tar -xvf samtools-1.10.tar.bz2 && rm samtools-1.10.tar.bz2
RUN cd /tools/samtools-1.10 && ./configure --without-curses && make && make install && mv samtools /usr/bin/

#### Install pilon version 1.23 ####
RUN apt-get install -y default-jre
RUN wget https://github.com/broadinstitute/pilon/releases/download/v1.23/pilon-1.23.jar && mv pilon-1.23.jar /usr/bin/

#### Install TransDeocoder version 5.5.0 ####
RUN apt-get install -y libdb-dev r-base
RUN wget -q -O - https://cpanmin.us | perl - App::cpanminus && cpanm install DB_File && cpanm install URI::Escape
RUN wget https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.1.tar.gz && tar -xvf TransDecoder-v5.7.1.tar.gz && rm TransDecoder-v5.7.1.tar.gz
RUN cd /tools/TransDecoder-TransDecoder-v5.7.1 && mv * /usr/bin/
#### We need to install a different version of `get_longest_ORF_per_transcript.pl` otherwise it throws a parse header error
# RUN wget https://github.com/TransDecoder/TransDecoder/blob/devel/util/get_longest_ORF_per_transcript.pl -O /usr/bin/util/get_longest_ORF_per_transcript.pl

#### Install Salmon version 1.3.0 ####
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.3.0/salmon-1.3.0_linux_x86_64.tar.gz && tar -xvf salmon-1.3.0_linux_x86_64.tar.gz && rm salmon-1.3.0_linux_x86_64.tar.gz
RUN cd /tools/salmon-latest_linux_x86_64 && mv * /usr/bin/ && ln -s /usr/bin/bin/salmon /usr/bin/salmon

#### Install miniconda with python 3.9 ####
ENV CONDA_DIR /tools/miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.5.2-0-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh -b -p /tools/miniconda && rm miniconda.sh
#### Put conda in path so we can use conda activate ####
ENV PATH=$CONDA_DIR/bin:$PATH

#### Install pysam
RUN pip install pysam

#### Install NanoFilt ####
RUN wget https://github.com/wdecoster/nanofilt/archive/refs/tags/v2.6.0.tar.gz -O nanofilt.tar.gz && tar -xvf nanofilt.tar.gz && rm nanofilt.tar.gz
RUN cd nanofilt-2.6.0 && python setup.py install

#### Install Trimmomatic-0.39 ####
RUN apt-get install -y unzip
RUN wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip && unzip Trimmomatic-0.39.zip && rm Trimmomatic-0.39.zip

#### Set working directory
WORKDIR /tmp

#### Configure permission settings ####
RUN chmod ugo+rwx /tmp
# RUN mkdir /work && chmod ugo+rwx /work

#### Create a user with the same UID as the host user ####
ARG USER_ID
ARG GROUP_ID
RUN addgroup --gid $GROUP_ID user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user
USER user

#### Run Jupyter Notebook backend ####
# CMD jupyter notebook --port=8887 --no-browser --ip=0.0.0.0 --allow-root --NotebookApp.password="" --NotebookApp.token=""
