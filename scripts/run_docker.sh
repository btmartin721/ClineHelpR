#!/bin/bash

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "Run ClineHelpR Docker image"
   echo
   echo "Usage: run_docker.sh <0 or 1> <PATH_TO_PROJECT_DIRECTORY>"
   echo "Options:"
   echo "-h     Print this Help."
   echo "-s     Shell to run docker in; 0 -> BASH shell, 1 -> Jupyter Notebook"
   echo "-p     Path to project directory (optional); If not specified, uses current working directory"
   echo
}

# Get the options
while getopts "hs:p:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      s) # Shell to use
	 METHOD=${OPTARG}
	 if [ ! "$METHOD" -ne 0 ] && [ ! "$METHOD" -ne 1 ]; then
		 echo "Command-line argument must equal 0 or 1!"
                 echo "Usage to run in BASH shell: run_docker.sh 0 <PROJECT_DIRECTORY_PATH>"
                 echo "Usage to run in Jupyter Notebook: run_docker.sh 1 <PROJECT_DIRECTORY_PATH>"
                 Help
                 exit 1;
         fi 
	 ;;
      p) # Path to project directory
	 PROJECT_DIR=${OPTARG}
	 ;;
      \?) # Incorrect arguments
	  echo "Error: Invalid option"
	  Help
	  exit;;
   esac
done

if [ -z "$METHOD" ]; then
        echo "Error: -s option is required!";
        Help
        exit 1;
fi

if [ -z "$PROJECT_DIR" ]; then
	echo "Project directory not specified; using current working directory"
	PROJECT_DIR=${PWD};
fi

if [ "$PROJECT_DIR" = '.' ]; then
	PROJECT_DIR=${PWD};
fi

if [ "$PROJECT_DIR" = './' ]; then
        PROJECT_DIR=${PWD};
fi

if [ ! -d "$PROJECT_DIR" ]; then
	echo "Specified project directory does not exist!"
	echo "Usage: run_docker.sh -s <0 or 1> -p <PROJECT_DIRECTORY_PATH>"
	Help
	exit 2;
fi

if [ ! -d "$PROJECT_DIR/data" ]; then
	mkdir ${PROJECT_DIR}/data;
fi

if [ ! -d "$PROJECT_DIR/notebooks" ]; then
        mkdir ${PROJECT_DIR}/notebooks;
fi

if [ ! -d "$PROJECT_DIR/results" ]; then
        mkdir ${PROJECT_DIR}/results;
fi

if [ "$METHOD" -eq 0 ]; then
	sudo docker container run --rm -it \
		--volume ${PROJECT_DIR}/data:/home/user/app/data \
		--volume ${PROJECT_DIR}/notebooks:/home/user/app/notebooks \
		--volume ${PROJECT_DIR}/results:/home/user/app/results \
		--publish 8888:8888 \
		btmartin721/clinehelpr:latest \
		/bin/bash;
else
	sudo docker container run --rm --tty \
                --volume ${PROJECT_DIR}/data:/home/user/app/data \
                --volume ${PROJECT_DIR}/notebooks:/home/user/app/notebooks \
                --volume ${PROJECT_DIR}/results:/home/user/app/results \
                --publish 8888:8888 \
                btmartin721/clinehelpr:latest;
fi

