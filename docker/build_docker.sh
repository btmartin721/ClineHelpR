#!/bin/bash

sudo docker image build --no-cache -t btmartin721/clinehelpr:latest \
	--build-arg uid=$UID   \
	--build-arg gid=$(id -g)   \
	--file Dockerfile \
	../

