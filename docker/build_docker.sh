#!/bin/bash

sudo docker image build -t clinehelpr:1.0 \
	--build-arg uid=$UID   \
	--build-arg gid=$(id -g)   \
	--file Dockerfile \
	../

