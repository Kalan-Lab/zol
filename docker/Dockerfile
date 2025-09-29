# This code is adapted from BiG-SCAPE's Dockerfile - if you are using zol - you should definitely check out
# BiG-SCAPE/CORASON incase it suites your needs better (https://www.nature.com/articles/s41589-019-0400-9) - eg.
# you are interested in clustering diverse gene-clusters based on protein domain simlarity or are interested
# in investigating the variability of contexts for a single reference gene!

FROM continuumio/miniconda3
LABEL maintainer="Rauf Salamzade - Kalan Lab, UW-Madison"

WORKDIR /usr/src
SHELL ["/bin/bash", "-c"]

# create and activate conda environment
RUN apt-get update && apt-get install -y git wget curl gcc
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba
RUN conda create -p /usr/src/zol_conda_env/ -c conda-forge -c bioconda zol

# Install zol
WORKDIR /usr/src/zol/ 
ENV PATH=/usr/src/zol_conda_env/bin:$PATH
ENV CONDA_PREFIX=/usr/src/zol_conda_env
ENV ZOL_EXEC_PATH=/usr/src/zol_conda_env/bin/
ENV ZOL_DATA_PATH=/usr/src/zol_conda_env/share/zol/db/
RUN mkdir -p "$ZOL_DATA_PATH" \
    && setup_annotation_dbs.py -m -y \
    && chmod -R 555 /usr/src/zol/ \
    && chmod 777 /home

USER 1000:1000
RUN mkdir /home/prepTG_input /home/prepTG_refproteome /home/prepTG_output && \
    mkdir /home/fai_input /home/fai_protque /home/fai_refgenome /home/fai_keyprot /home/fai_targetgen /home/fai_output && \
    mkdir /home/zol_input /home/zol_customdb /home/zol_foclisting /home/zol_complisting /home/zol_output /home/cgc_input && \
    mkdir /home/cgc_output /home/cgcg_input /home/cgcg_output

WORKDIR /home
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

ENTRYPOINT ["ZOL"]
CMD ["--help"]
