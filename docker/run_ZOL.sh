#!/bin/bash

# function
get_abs_filename() {
  # $1 : relative filename
  filename=$1
  parentdir=$(dirname "${filename}")

  if [ -d "${filename}" ]; then
      echo "$(cd "${filename}" && pwd)"
  elif [ -d "${parentdir}" ]; then
    echo "$(cd "${parentdir}" && pwd)/$(basename "${filename}")"
  fi
}

## The following is also adapted from analogous file from BIG-SCAPE
if [[ $# -eq 0 || $1 == "-h" || $1 == "--help" || $1 == "-v" || $1 == "--version" ]]; then
	docker pull raufs/zol:latest
	docker run \
	--detach=false \
	--rm \
	--user=$(id -u):$(id -g) \
	raufs/zol:latest \
elif [[ $1 == 'prepTG' ]]; then
	set -o errexit
	set -o nounset

  ORIG_ARGS=$@
  INPUT_DIR="NA"
  INPUT_PARENT_DIR="NA"
  OUTPUT_DIR="NA"
  OUTPUT_PARENT_DIR="NA"
  REFERENCE_PROTEOME_DIR="NA"
  REFERENCE_PROTEOME_PARENT_DIR="NA"

  for var in "$@"
  do
    if [[ "$var" == '-i' || "$var" == '--input_dir' ]]; then
      shift
      ABS_VALUE=(get_abs_filename $1)
    elif [[ "$var" == '-o' || "$var" == '--output_dir' ]]; then
      shift
    elif [[ "$var" == '-rp' || "$var" == '--reference_proteome' ]]; then
      shift
    else
      shift
    fi
  done

	# get input directory and parent directory


  # get output directory and parent directory

	if [[ $# -eq 0 ]]; then
		echo You must specify an output dir
		exit 1;
	else
		readonly OUTPUT_FILE=$(basename $1)
		echo output ${OUTPUT_FILE}
		readonly OUTPUT_DIR=$(dirname $(realpath $1))
		shift
	fi

	if [ ! -d ${OUTPUT_DIR} ]; then
   		mkdir ${OUTPUT_DIR}
	fi

	# Links within the container
	readonly CONTAINER_INPUT_DIR=/home/prepTG_input
	readonly CONTAINER_REFPROTEOME_DIR=/home/prepTG_refproteome
	readonly CONTAINER_OUTPUT_DIR=/home/prepTG_output

  # run individual programs
	docker pull raufs/zol:latest
	docker run \
	    --volume ${INPUT_DIR}:${CONTAINER_SRC_DIR}:ro \
	    --volume ${OUTPUT_DIR}:${CONTAINER_DST_DIR}:rw \
	    --detach=false \
	    --rm \
	    --user=$(id -u):$(id -g) \
	     raufs/zol:latest \
	    -i ${CONTAINER_SRC_DIR}/${INPUT_FILE} \
	    -o ${CONTAINER_DST_DIR}/${OUTPUT_FILE}\
    	$@
elif [[ ]]; then

elif [[ ]]; then

fi
