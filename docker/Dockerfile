FROM ubuntu:18.04

LABEL maintainer="evobio721@gmail.com" \
      version="1.0" \
      github="https://github.com/btmartin721/ClineHelpR" \
      description="Docker image for ClineHelpR that includes all required dependencies"

# Specify bash as the default shell
SHELL [ "/bin/bash", "--login", "-c" ]

ARG username=current_user
ARG uid=1000
ARG gid=100
ARG IMAGE_NAME=ClineHelpR
ARG IMAGE_TAG=1.0
ENV USER $username
ENV UID $uid
ENV GID $gid
ENV HOME /home/$USER
ENV IMAGE_NAME $IMAGE_NAME
ENV IMAGE_TAG $IMAGE_TAG

RUN apt-get update && \
    apt-get install -y wget && \
    rm -rf /var/lib/apt/lists/*

RUN adduser --disabled-password \
    --gecos "Non-root user" \
    --uid $UID \
    --gid $GID \
    --home $HOME \
    $USER

COPY environment.yml /tmp/
RUN chown $UID:$GID /tmp/environment.yml

COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chown $UID:$GID /usr/local/bin/entrypoint.sh && \
    chmod u+x /usr/local/bin/entrypoint.sh

USER $USER

# install miniconda
ENV MINICONDA_VERSION 4.10.2
ENV CONDA_DIR $HOME/miniconda3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-$MINICONDA_VERSION-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

# make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH

# make conda activate command available from /bin/bash --login shells
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile

# make conda activate command available from /bin/bash --interative shells
RUN conda init bash

# create a project directory inside user home
ENV PROJECT_DIR $HOME/app
RUN mkdir $PROJECT_DIR
WORKDIR $PROJECT_DIR

# build the conda environment
ENV ENV_PREFIX $PWD/env
RUN conda update --name base --channel defaults conda && \
    conda env create --prefix $ENV_PREFIX --file /tmp/environment.yml --force && \
    conda clean --all --yes

# Create entry point to make sure conda is activated at runtime
ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]

# install other tools not available on conda cloud
RUN apt-get -y update && \
    apt-get -y install build-essential && \
    apt-get -y install libgsl-dev && \
    apt-get -y install libhdf5-serial-dev

# Compile BGC
RUN h5c++ -Wall -O2 -o bgc bgc_main.C bgc_func_readdata.C bgc_func_initialize.C bgc_func_mcmc.C bgc_func_write.C bgc_fu$

# Install INTROGRESS R package
RUN R -e "install.packages('introgress', dependencies=TRUE, repos='http://cran.us.r-project.org')"

# default command will launch JupyterLab server for development
CMD [ "jupyter", "lab", "--no-browser", "--ip", "0.0.0.0" ]

