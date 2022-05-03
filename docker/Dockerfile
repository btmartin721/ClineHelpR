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
    
# Set project directory
ENV PROJECT_DIR $HOME/app

# Set other frequently used environmental variables
ENV ENV_PREFIX $PROJECT_DIR/env
ENV SRC_DIR $PROJECT_DIR/src
ENV BIN_DIR $PROJECT_DIR/bin

# Add a user so that they aren't executing as root
# Then install dependencies for BGC and tools for downloadin files using apt-get
RUN adduser --disabled-password \
    --gecos "Non-root user" \
    --uid $UID \
    --home $HOME \
    $USER && \
    apt-get -y update && \
    apt-get -y install build-essential \
    libgsl-dev \
    libhdf5-serial-dev \
    wget \
    git

##############################################
### Install miniconda
##############################################

# Move conda environment to /tmp and changer ownership to user
COPY environment.yml /tmp/

# Copy entrypoint.sh to /usr/local/bin and give user privelidges
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh

RUN chown $UID /tmp/environment.yml && \
    chown $UID /usr/local/bin/entrypoint.sh && \
    chmod u+x /usr/local/bin/entrypoint.sh

# Change from root to user
USER $USER

# install miniconda
ENV CONDA_DIR $HOME/miniconda3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh && \
    mkdir $HOME/app

# make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH

# Create a project directory inside user home
WORKDIR $PROJECT_DIR

# make conda activate command available from /bin/bash --login shells
# make conda activate command available from /bin/bash --interative shells
# Then build the conda environment
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile && \
    conda init bash && \
    conda update --name base --channel defaults conda && \
    conda env create --prefix $ENV_PREFIX --file /tmp/environment.yml --force && \
    conda clean --all --yes && \
    mkdir src && \
    mkdir bin && \
    mkdir src/tkchafin

WORKDIR $SRC_DIR

# Create entry point to make sure conda is activated at runtime
ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]

##########################################
### Install BGC
##########################################

# Copy BGC source code to image
COPY bgcdist1.03.tar.gz .
RUN tar -xzvf bgcdist1.03.tar.gz
WORKDIR bgcdist/

# Compile BGC and estpost
RUN h5c++ -Wall -O2 -o bgc bgc_main.C bgc_func_readdata.C bgc_func_initialize.C \
    bgc_func_mcmc.C bgc_func_write.C bgc_func_linkage.C bgc_func_ngs.C \
    bgc_func_hdf5.C mvrandist.c -lgsl -lgslcblas && \
    h5c++ -Wall -O3 -o estpost estpost_h5.c -lgsl -lgslcblas && \
    chmod u+x bgc && \
    chmod u+x estpost

# Add everything to the path
ENV PATH="${PROJECT_DIR}/bin:${PROJECT_DIR}/env/bin:${PATH}"

#####################################
### Install R dependencies
#####################################

# Install INTROGRESS R package and ClineHelpR API, and setup R for use with Jupyter notebooks
RUN echo "Updated genind2introgress functio" && \
    R -e "install.packages('introgress', dependencies=TRUE, repos='http://cran.us.r-project.org')" && \
    R -e "install.packages('IRkernel', repos='http://cran.us.r-project.org')" && \
    R -e "IRkernel::installspec(user = TRUE)" && \
    R -e "devtools::install_github('btmartin721/ClineHelpR')" && \
    R -e "download.file(url = 'https://github.com/btmartin721/ClineHelpR/archive/refs/heads/master.zip', destfile = 'ClineHelpR-master.zip')" && \
    R -e "unzip(zipfile = 'ClineHelpR-master.zip')"

# Copy ClineHelpR scripts to container path
COPY scripts $SRC_DIR/scripts

WORKDIR $SRC_DIR/tkchafin

USER root

# Get supporting scripts from tkchafin/scripts repo
# Copy supporting scripts to scripts directory
RUN git clone --depth 1 https://github.com/tkchafin/scripts.git && \
    cp scripts/phylip2bgc.pl $SRC_DIR/scripts/ && \
    cp scripts/phylip2introgress.pl $SRC_DIR/scripts/ && \
    cp scripts/nremover.pl $SRC_DIR/scripts/

USER $USER

# Put it in the path
ENV PATH=${SRC_DIR}/scripts:${PATH}

WORKDIR $PROJECT_DIR

# default command will launch JupyterLab server for development
CMD [ "jupyter", "lab", "--no-browser", "--ip", "0.0.0.0" ]

