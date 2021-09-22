FROM ubuntu:18.04

LABEL maintainer="evobio721@gmail.com" \
      version="1.0" \
      github="https://github.com/btmartin721/ClineHelpR" \
      description="Docker image for ClineHelpR that includes all required dependencies" \
      code_source="https://towardsdatascience.com/conda-pip-and-docker-ftw-d64fe638dc45"

# Specify bash as the default shell
SHELL [ "/bin/bash", "--login", "-c" ]

ARG username=user
ARG uid=1000
ARG gid=100
ARG IMAGE_NAME=clinehelpr
ARG IMAGE_TAG=1.0
ENV USER $username
ENV UID $uid
ENV GID $gid
ENV HOME /home/$USER
ENV IMAGE_NAME $IMAGE_NAME
ENV IMAGE_TAG $IMAGE_TAG

# Add a user so that they aren't executing as root
RUN adduser --disabled-password \
    --gecos "Non-root user" \
    --uid $UID \
    --home $HOME \
    $USER

###############################################
### Install tools not available on conda cloud
###############################################

# Install dependencies for BGC and tools for downloading files
RUN apt-get -y update && \
    apt-get -y install build-essential && \
    apt-get -y install libgsl-dev && \
    apt-get -y install libhdf5-serial-dev && \
    apt-get -y install wget && \
    apt-get -y install git

##############################################
### Install miniconda
##############################################

# Move conda environment to /tmp and changer ownership to user
COPY environment.yml /tmp/
RUN chown $UID /tmp/environment.yml

# Copy entrypoint.sh to /usr/loca/bin and give user privelidges 
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chown $UID /usr/local/bin/entrypoint.sh && \
    chmod u+x /usr/local/bin/entrypoint.sh

# Change from root to user
USER $USER

# install miniconda
ENV CONDA_DIR $HOME/miniconda3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

# make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH

# make conda activate command available from /bin/bash --login shells
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile

# make conda activate command available from /bin/bash --interative shells
RUN conda init bash

# Create a project directory inside user home
ENV PROJECT_DIR $HOME/app
RUN mkdir $PROJECT_DIR
WORKDIR $PROJECT_DIR

# Build the conda environment
ENV ENV_PREFIX $PROJECT_DIR/env
RUN conda update --name base --channel defaults conda && \
    conda env create --prefix $ENV_PREFIX --file /tmp/environment.yml --force && \
    conda clean --all --yes

# Create entry point to make sure conda is activated at runtime
ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]

# Create bin and src directories
RUN mkdir src && mkdir bin
ENV SRC_DIR $PROJECT_DIR/src
ENV BIN_DIR $PROJECT_DIR/bin
WORKDIR $SRC_DIR

##########################################
### Install BGC
##########################################

# Copy BGC source code to image
COPY bgcdist1.03.tar.gz .
RUN tar -xzvf bgcdist1.03.tar.gz
WORKDIR bgcdist/

# Compile BGC
RUN h5c++ -Wall -O2 -o bgc bgc_main.C bgc_func_readdata.C bgc_func_initialize.C \
    bgc_func_mcmc.C bgc_func_write.C bgc_func_linkage.C bgc_func_ngs.C \
    bgc_func_hdf5.C mvrandist.c -lgsl -lgslcblas

# Compile estpost
RUN h5c++ -Wall -O3 -o estpost estpost_h5.c -lgsl -lgslcblas

# Make binaries executable
RUN chmod u+x bgc && \
    chmod u+x estpost

# Add everything to the path
ENV PATH="${PROJECT_DIR}/bin:${PROJECT_DIR}/env/bin:${PATH}"

#####################################
### Install R dependencies
#####################################

# Install INTROGRESS R package and setup R for use with Jupyter notebooks
RUN R -e "install.packages('introgress', dependencies=TRUE, repos='http://cran.us.r-project.org')" && \
    R -e "install.packages('IRkernel', repos='http://cran.us.r-project.org')" && \
    R -e "IRkernel::installspec(user = TRUE)"

# Install ClineHelpR API
RUN R -e "devtools::install_github('btmartin721/ClineHelpR')" && \
    R -e "download.file(url = 'https://github.com/btmartin721/ClineHelpR/archive/refs/heads/master.zip', destfile = 'ClineHelpR-master.zip')" && \
    R -e "unzip(zipfile = 'ClineHelpR-master.zip')"

# Copy ClineHelpR scripts to container path
COPY scripts $SRC_DIR/scripts

RUN mkdir $SRC_DIR/tkchafin

WORKDIR $SRC_DIR/tkchafin

# Get supporting scripts from tkchafin/scripts repo
RUN git clone --depth 1 https://github.com/tkchafin/scripts.git

# Copy supporting scripts to scripts directory
RUN cp scripts/phylip2bgc.pl $SRC_DIR/scripts/ && \
    cp scripts/phylip2introgress.pl $SRC_DIR/scripts/ && \
    cp scripts/nremover.pl $SRC_DIR/scripts/

# Put it in the path
ENV PATH=${SRC_DIR}/scripts:${PATH}

WORKDIR $PROJECT_DIR

# default command will launch JupyterLab server for development
CMD [ "jupyter", "lab", "--no-browser", "--ip", "0.0.0.0" ]

