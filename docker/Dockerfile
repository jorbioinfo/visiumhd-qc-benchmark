FROM bioconductor/bioconductor_docker:RELEASE_3_23

ENV DEBIAN_FRONTEND=noninteractive
ENV R_REMOTES_UPGRADE=never

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3-pip python3-venv libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

ENV VENV=/opt/pyenv
RUN python3 -m venv $VENV \
 && $VENV/bin/pip install --no-cache-dir --upgrade pip \
 && $VENV/bin/pip install --no-cache-dir \
        numpy==2.1.3 pandas==2.2.3 scipy==1.14.1 pyarrow==18.1.0 \
        matplotlib==3.9.3 \
 && $VENV/bin/python -c "import numpy,scipy,pyarrow,pandas,matplotlib; print('python deps OK')"
ENV PATH=$VENV/bin:$PATH

RUN R -q -e 'BiocManager::install(c( \
        "SpotSweeper", "SpatialExperiment", "SingleCellExperiment", \
        "scuttle", "scater", "BiocNeighbors", "DropletUtils", "escheR" \
    ), ask = FALSE, update = FALSE)' \
 && R -q -e 'stopifnot(all(c("SpotSweeper","SpatialExperiment","scuttle","DropletUtils") \
        %in% rownames(installed.packages())))'

RUN R -q -e 'install.packages(c("arrow","dplyr","tidyr","stringr","jsonlite", \
        "ggplot2","ggvenn","magrittr","Matrix","data.table"), Ncpus = 4)' \
 && R -q -e 'stopifnot("arrow" %in% rownames(installed.packages()))'

RUN R -q -e 'BiocManager::install("SpatialArtifacts", ask = FALSE, update = FALSE)' \
 && R -q -e 'stopifnot("SpatialArtifacts" %in% rownames(installed.packages())); \
        cat("SpatialArtifacts", as.character(packageVersion("SpatialArtifacts")), "OK\n")'

RUN R -q -e 'BiocManager::install("scrapper", ask = FALSE, update = FALSE)' \
 && R -q -e 'stopifnot("scrapper" %in% rownames(installed.packages())); \
        cat("scrapper", as.character(packageVersion("scrapper")), "OK\n")'

RUN R -q -e 'BiocManager::install(c("miQC", "flexmix"), ask = FALSE, update = FALSE)' \
 && R -q -e 'library(miQC); library(flexmix); \
        cat("miQC", as.character(packageVersion("miQC")), \
            "flexmix", as.character(packageVersion("flexmix")), "OK\n")'

RUN /opt/pyenv/bin/pip install --no-cache-dir \
        "scanpy==1.11.*" "anndata==0.11.*" "leidenalg==0.10.*" \
        "igraph==0.11.*" "scikit-learn==1.5.*" \
 && /opt/pyenv/bin/python -c "import scanpy, leidenalg, sklearn; \
        print('scanpy', scanpy.__version__, 'sklearn', sklearn.__version__, 'OK')"

WORKDIR /work
CMD ["bash"]
