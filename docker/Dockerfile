# This code is adapted from BiG-SCAPE's Dockerfile - if you are using zol - you should definitely check out
# BiG-SCAPE/CORASON incase it suites your needs better (https://www.nature.com/articles/s41589-019-0400-9) - eg.
# you are interested in clustering diverse gene-clusters based on protein domain simlarity or are interested
# in investigating the variability of contexts for a single reference gene!

FROM continuumio/miniconda3
LABEL maintainer="Rauf Salamzade - Kalan Lab, UW-Madison"

WORKDIR /usr/src
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y git wget

# Clone zol github repo and create conda environment
RUN git clone https://github.com/Kalan-Lab/zol && rm -rf zol/test_case.tar.gz

# Create and activate conda environment
RUN conda install -n base conda-libmamba-solver && conda config --set solver libmamba && \
    conda env create -f /usr/src/zol/zol_env.yml -p /usr/src/zol_conda_env/ && conda clean --all -y && \
    echo "source activate /usr/src/zol_conda_env/" > ~/.bashrc && source ~/.bashrc && \
    apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

# Install zol
WORKDIR /usr/src/zol/ 
ENV PATH /usr/src/zol_conda_env/bin:$PATH
RUN python setup.py install && pip install -e . && setup_annotation_dbs.py -m && chmod -R 555 /usr/src/zol/ && chmod 777 /home

USER 1000:1000
RUN mkdir /home/prepTG_input /home/prepTG_refproteome /home/prepTG_output && \
    mkdir /home/fai_input /home/fai_protque /home/fai_refgenome /home/fai_keyprot /home/fai_targetgen /home/fai_output && \
    mkdir /home/zol_input /home/zol_customdb /home/zol_foclisting /home/zol_complisting /home/zol_output

WORKDIR /home
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

ENTRYPOINT ["ZOL"]
CMD ["--help"]
