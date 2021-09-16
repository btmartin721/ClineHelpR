FROM continuumio/miniconda
LABEL maintainer="evobio721@gmail.com"

COPY clinehelpr_env.yml .
RUN \
    conda env update -n root -f clinehelpr_env.yml \
    && conda clean -a

# install other tools not available on conda cloud
RUN apt-get update && \
    apt-get -y install build-essential && \
    apt-get -y install libgsl-dev && \
    apt-get -y install libhdf5-serial-dev
