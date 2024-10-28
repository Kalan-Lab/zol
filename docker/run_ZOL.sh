#!/bin/bash

# AUTHOR: Rauf Salamzade
# AFFILIATION: Kalan Lab, UW-Madison
# run_ZOL.sh - wrapper to run 3 major programs of zol using Docker image.

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



elif [[ $1 == 'cgc' ]]; then
        set -o errexit
        set -o nounset

        # Links within the container
        readonly CONTAINER_INPUT_DIR=/home/cgc_input
        readonly CONTAINER_OUTPUT_DIR=/home/cgc_output

  # variables for updating/input paths to account for Docker mounting

  DOCKER_VOLUME_ARGS=""
  CGC_ARGS=""
  OUTPUT_PARENT_DIR="NA"
    while [[ ! $# -eq 0 ]]; do
      if [[ "$1" == '-i' || "$1" == '--zol-results-dir' ]]; then
        shift
        ABS_VALUE=$(get_abs_filename $1)
        INPUT_DIR=$(basename $ABS_VALUE)
        INPUT_PARENT_DIR=$(dirname $ABS_VALUE)
        CGC_ARGS+="-i $CONTAINER_INPUT_DIR/$INPUT_DIR "
        DOCKER_VOLUME_ARGS+="--volume $INPUT_PARENT_DIR:$CONTAINER_INPUT_DIR:ro "
        shift
      elif [[ "$1" == '-o' || "$1" == '--output-dir' ]]; then
        shift
        ABS_VALUE=$(get_abs_filename $1)
        OUTPUT_DIR=$(basename $ABS_VALUE)
        OUTPUT_PARENT_DIR=$(dirname $ABS_VALUE)
        CGC_ARGS+="-o $CONTAINER_OUTPUT_DIR/$OUTPUT_DIR "
        DOCKER_VOLUME_ARGS+="--volume $OUTPUT_PARENT_DIR:$CONTAINER_OUTPUT_DIR:rw "
        shift
      else
        CGC_ARGS+="$1 "
        shift
      fi
    done

        if [[ ! -d ${OUTPUT_PARENT_DIR} && $OUTPUT_PARENT_DIR != "NA" ]]; then
                mkdir ${OUTPUT_PARENT_DIR}
        fi

  # run cgc
  docker pull raufs/zol:latest
  docker run ${DOCKER_VOLUME_ARGS} --detach=false --rm --user=$(id -u):$(id -g) raufs/zol:latest ${CGC_ARGS}

elif [[ $1 == 'cgcg' ]]; then
  set -o errexit
  set -o nounset

  # Links within the container
  readonly CONTAINER_INPUT_DIR=/home/cgcg_input
  readonly CONTAINER_OUTPUT_DIR=/home/cgcg_output

  # variables for updating/input paths to account for Docker mounting

  DOCKER_VOLUME_ARGS=""
  CGCG_ARGS=""
  OUTPUT_PARENT_DIR="NA"
    while [[ ! $# -eq 0 ]]; do
      if [[ "$1" == '-i' || "$1" == '--zol-results-dir' ]]; then
        shift
        ABS_VALUE=$(get_abs_filename $1)
        INPUT_DIR=$(basename $ABS_VALUE)
        INPUT_PARENT_DIR=$(dirname $ABS_VALUE)
        CGCG_ARGS+="-i $CONTAINER_INPUT_DIR/$INPUT_DIR "
        DOCKER_VOLUME_ARGS+="--volume $INPUT_PARENT_DIR:$CONTAINER_INPUT_DIR:ro "
        shift
      elif [[ "$1" == '-o' || "$1" == '--output-dir' ]]; then
        shift
        ABS_VALUE=$(get_abs_filename $1)
        OUTPUT_DIR=$(basename $ABS_VALUE)
        OUTPUT_PARENT_DIR=$(dirname $ABS_VALUE)
        CGCG_ARGS+="-o $CONTAINER_OUTPUT_DIR/$OUTPUT_DIR "
        DOCKER_VOLUME_ARGS+="--volume $OUTPUT_PARENT_DIR:$CONTAINER_OUTPUT_DIR:rw "
        shift
      else
        CGCG_ARGS+="$1 "
        shift
      fi
    done

        if [[ ! -d ${OUTPUT_PARENT_DIR} && $OUTPUT_PARENT_DIR != "NA" ]]; then
                mkdir ${OUTPUT_PARENT_DIR}
        fi

  # run cgcg
  docker pull raufs/zol:latest
  docker run ${DOCKER_VOLUME_ARGS} --detach=false --rm --user=$(id -u):$(id -g) raufs/zol:latest ${CGCG_ARGS}

elif [[ $1 == 'prepTG' ]]; then
	set -o errexit
	set -o nounset

	# Links within the container
	readonly CONTAINER_INPUT_DIR=/home/prepTG_input
	readonly CONTAINER_REFPROTEOME_DIR=/home/prepTG_refproteome
	readonly CONTAINER_OUTPUT_DIR=/home/prepTG_output

  # variables for updating/input paths to account for Docker mounting

  DOCKER_VOLUME_ARGS=""
  PREPTG_ARGS=""
  OUTPUT_PARENT_DIR="NA"
    while [[ ! $# -eq 0 ]]; do
      if [[ "$1" == '-i' || "$1" == '--input-dir' ]]; then
        shift
        ABS_VALUE=$(get_abs_filename $1)
        INPUT_DIR=$(basename $ABS_VALUE)
        INPUT_PARENT_DIR=$(dirname $ABS_VALUE)
        PREPTG_ARGS+="-i $CONTAINER_INPUT_DIR/$INPUT_DIR "
        DOCKER_VOLUME_ARGS+="--volume $INPUT_PARENT_DIR:$CONTAINER_INPUT_DIR:ro "
        shift
      elif [[ "$1" == '-o' || "$1" == '--output-dir' ]]; then
        shift
        ABS_VALUE=$(get_abs_filename $1)
        OUTPUT_DIR=$(basename $ABS_VALUE)
        OUTPUT_PARENT_DIR=$(dirname $ABS_VALUE)
        PREPTG_ARGS+="-o $CONTAINER_OUTPUT_DIR/$OUTPUT_DIR "
        DOCKER_VOLUME_ARGS+="--volume $OUTPUT_PARENT_DIR:$CONTAINER_OUTPUT_DIR:rw "
        shift
      elif [[ "$1" == '-rp' || "$1" == '--reference-proteome' ]]; then
        shift
        ABS_VALUE=$(get_abs_filename $1)
        REFERENCE_PROTEOME=$(basename $ABS_VALUE)
        REFERENCE_PROTEOME_PARENT_DIR=$(dirname $ABS_VALUE)
        PREPTG_ARGS+="-rp $CONTAINER_REFPROTEOME_DIR/$REFERENCE_PROTEOME "
        DOCKER_VOLUME_ARGS+="--volume $REFERENCE_PROTEOME_PARENT_DIR:$CONTAINER_REFPROTEOME_DIR:ro "
        shift
      else
        PREPTG_ARGS+="$1 "
        shift
      fi
    done

	if [[ ! -d ${OUTPUT_PARENT_DIR} && $OUTPUT_PARENT_DIR != "NA" ]]; then
   		mkdir ${OUTPUT_PARENT_DIR}
	fi

  # run prepTG
  docker pull raufs/zol:latest
  docker run ${DOCKER_VOLUME_ARGS} --detach=false --rm --user=$(id -u):$(id -g) raufs/zol:latest ${PREPTG_ARGS}

elif [[ $1 == 'fai' ]]; then
  set -o errexit
  set -o nounset
  
  # Links within the container
  readonly CONTAINER_INPUT_DIR=/home/fai_input
  readonly CONTAINER_PROTQUE_DIR=/home/fai_protque
  readonly CONTAINER_REFGENOME_DIR=/home/fai_refgenome
  readonly CONTAINER_KEYPROT_DIR=/home/fai_keyprot
  readonly CONTAINER_TARGETGEN_DIR=/home/fai_targetgen
  readonly CONTAINER_OUTPUT_DIR=/home/fai_output

  # variables for updating/input paths to account for Docker mounting
  DOCKER_VOLUME_ARGS=""
  FAI_ARGS=""
  OUTPUT_PARENT_DIR="NA"
  while [[ ! $# -eq 0 ]]; do
    if [[ "$1" == '-i' || "$1" == '--input-dir' ]]; then
      shift
      ABS_VALUE=$(get_abs_filename $1)
      INPUT_DIR=$(basename $ABS_VALUE)
      INPUT_PARENT_DIR=$(dirname $ABS_VALUE)
      FAI_ARGS+="-i $CONTAINER_INPUT_DIR/$INPUT_DIR "
      DOCKER_VOLUME_ARGS+="--volume $INPUT_PARENT_DIR:$CONTAINER_INPUT_DIR:ro "
      shift
    elif [[ "$1" == '-r' || "$1" == '--reference-genome' ]]; then
      shift
      ABS_VALUE=$(get_abs_filename $1)
      REFGENOME=$(basename $ABS_VALUE)
      REFGENOME_PARENT_DIR=$(dirname $ABS_VALUE)
      FAI_ARGS+="-r $CONTAINER_REFGENOME_DIR/$REFGENOME "
      DOCKER_VOLUME_ARGS+="--volume $REFGENOME_PARENT_DIR:$CONTAINER_REFGENOME_DIR:ro "
      shift
    elif [[ "$1" == '-pq' || "$1" == '--protein-queries' ]]; then
      shift
      ABS_VALUE=$(get_abs_filename $1)
      PROTEIN_QUERY=$(basename $ABS_VALUE)
      PROTEIN_QUERY_PARENT_DIR=$(dirname $ABS_VALUE)
      FAI_ARGS+="-pq $CONTAINER_PROTQUE_DIR/$PROTEIN_QUERY "
      DOCKER_VOLUME_ARGS+="--volume $PROTEIN_QUERY_PARENT_DIR:$CONTAINER_PROTQUE_DIR:ro "
      shift
    elif [[ "$1" == '-kp' || "$1" == '--key-proteins' ]]; then
      shift
      ABS_VALUE=$(get_abs_filename $1)
      KEY_QUERY_PROTEINS=$(basename $ABS_VALUE)
      KEY_QUERY_PROTEINS_PARENT_DIR=$(dirname $ABS_VALUE)
      FAI_ARGS+="-kp $CONTAINER_KEYPROT_DIR/$KEY_QUERY_PROTEINS "
      DOCKER_VOLUME_ARGS+="--volume $KEY_QUERY_PROTEINS_PARENT_DIR:$CONTAINER_KEYPROT_DIR:ro "
      shift
    elif [[ "$1" == '-tg' || "$1" == '--target-genomes' ]]; then
      shift
      ABS_VALUE=$(get_abs_filename $1)
      TARGET_GENOMES=$(basename $ABS_VALUE)
      TARGET_GENOMES_PARENT_DIR=$(dirname $ABS_VALUE)
      FAI_ARGS+="-tg $CONTAINER_TARGETGEN_DIR/$TARGET_GENOMES "
      DOCKER_VOLUME_ARGS+="--volume $TARGET_GENOMES_PARENT_DIR:$CONTAINER_TARGETGEN_DIR:ro "
      shift
    elif [[ "$1" == '-o' || "$1" == '--output-dir' ]]; then
      shift
      ABS_VALUE=$(get_abs_filename $1)
      OUTPUT_DIR=$(basename $ABS_VALUE)
      OUTPUT_PARENT_DIR=$(dirname $ABS_VALUE)
      FAI_ARGS+="-o $CONTAINER_OUTPUT_DIR/$OUTPUT_DIR "
      DOCKER_VOLUME_ARGS+="--volume $OUTPUT_PARENT_DIR:$CONTAINER_OUTPUT_DIR:rw "
      shift
    else
      FAI_ARGS+="$1 "
      shift
    fi
  done
	if [[ ! -d ${OUTPUT_PARENT_DIR} && $OUTPUT_PARENT_DIR != "NA" ]]; then
   		mkdir ${OUTPUT_PARENT_DIR}
	fi

  # run fai
  docker pull raufs/zol:latest
  docker run ${DOCKER_VOLUME_ARGS} --detach=false --rm --user=$(id -u):$(id -g) raufs/zol:latest ${FAI_ARGS}

elif [[ $1 == 'zol' ]]; then
  set -o errexit
  set -o nounset

  # Links within the container
  readonly CONTAINER_INPUT_DIR=/home/fai_input
  readonly CONTAINER_CUSTOMDB_DIR=/home/zol_customdb
  readonly CONTAINER_FOCLIST_DIR=/home/zol_foclisting
  readonly CONTAINER_COMLIST_DIR=/home/zol_complisting
  readonly CONTAINER_OUTPUT_DIR=/home/zol_output

  # variables for updating/input paths to account for Docker mounting
  DOCKER_VOLUME_ARGS=""
  ZOL_ARGS=""
  OUTPUT_PARENT_DIR="NA"
  while [[ ! $# -eq 0 ]]; do
    if [[ "$1" == '-i' || "$1" == '--input-dir' ]]; then
      shift
      ABS_VALUE=$(get_abs_filename $1)
      INPUT_DIR=$(basename $ABS_VALUE)
      INPUT_PARENT_DIR=$(dirname $ABS_VALUE)
      ZOL_ARGS+="-i $CONTAINER_INPUT_DIR/$INPUT_DIR "
      DOCKER_VOLUME_ARGS+="--volume $INPUT_PARENT_DIR:$CONTAINER_INPUT_DIR:ro "
      shift
    elif [[ "$1" == '-f' || "$1" == '--focal-genbanks' ]]; then
      shift
      ABS_VALUE=$(get_abs_filename $1)
      FOCLISTING=$(basename $ABS_VALUE)
      FOCLISTING_PARENT_DIR=$(dirname $ABS_VALUE)
      ZOL_ARGS+="-f $CONTAINER_FOCLIST_DIR/$FOCLISTING "
      DOCKER_VOLUME_ARGS+="--volume $FOCLISTING_PARENT_DIR:$CONTAINER_FOCLIST_DIR:ro "
      shift
    elif [[ "$1" == '-fc' || "$1" == '--comparator-genbanks' ]]; then
      shift
      ABS_VALUE=$(get_abs_filename $1)
      COMLISTING=$(basename $ABS_VALUE)
      COMLISTING_PARENT_DIR=$(dirname $ABS_VALUE)
      ZOL_ARGS+="-fc $CONTAINER_COMLIST_DIR/$COMLISTING "
      DOCKER_VOLUME_ARGS+="--volume $COMLISTING_PARENT_DIR:$CONTAINER_COMLIST_DIR:ro "
      shift
    elif [[ "$1" == '-cd' || "$1" == '--custom-database' ]]; then
      shift
      ABS_VALUE=$(get_abs_filename $1)
      CUSTOM_DB=$(basename $ABS_VALUE)
      CUSTOM_DB_PARENT_DIR=$(dirname $ABS_VALUE)
      ZOL_ARGS+="-cd $CONTAINER_CUSTOMDB_DIR/$CUSTOM_DB "
      DOCKER_VOLUME_ARGS+="--volume $CUSTOM_DB_PARENT_DIR:$CONTAINER_CUSTOMDB_DIR:ro "
      shift
    elif [[ "$1" == '-o' || "$1" == '--output-dir' ]]; then
      shift
      ABS_VALUE=$(get_abs_filename $1)
      OUTPUT_DIR=$(basename $ABS_VALUE)
      OUTPUT_PARENT_DIR=$(dirname $ABS_VALUE)
      ZOL_ARGS+="-o $CONTAINER_OUTPUT_DIR/$OUTPUT_DIR "
      DOCKER_VOLUME_ARGS+="--volume $OUTPUT_PARENT_DIR:$CONTAINER_OUTPUT_DIR:rw "
      shift
    else
      ZOL_ARGS+="$1 "
      shift
    fi
  done
	if [[ ! -d ${OUTPUT_PARENT_DIR} && $OUTPUT_PARENT_DIR != "NA" ]]; then
   		mkdir ${OUTPUT_PARENT_DIR}
	fi

  # run zol
	docker pull raufs/zol:latest
	docker run ${DOCKER_VOLUME_ARGS} --detach=false --rm --user=$(id -u):$(id -g) raufs/zol:latest ${ZOL_ARGS}

else
  docker pull raufs/zol:latest
	docker run \
	--detach=false \
	--rm \
	--user=$(id -u):$(id -g) \
	raufs/zol:latest \

fi
