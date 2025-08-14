FROM rocker/r-ver:4.3.3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    pandoc \
    pandoc-citeproc \
    libgit2-dev \
    libssh2-1-dev \
    libpq-dev \
    libbz2-dev \
    liblzma-dev \
    libicu-dev \
    make \
    gcc \
    g++ \
    gfortran \
    python3 \
    python3-pip \
    curl \
    wget \
    vim \
    nano \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# Install R packages needed for development
RUN install2.r --error --deps TRUE \
    remotes \
    devtools \
    testthat \
    knitr \
    rmarkdown \
    roxygen2 \
    usethis \
    pkgdown \
    covr \
    styler \
    lintr \
    && install2.r --error --deps TRUE \
    igraph \
    ggplot2 \
    dplyr \
    tidyr \
    tibble \
    readr \
    stringr \
    purrr \
    rlang \
    reshape2 \
    pheatmap \
    RColorBrewer \
    circlize \
    && install2.r --error --deps TRUE \
    BiocManager \
    && R -e "BiocManager::install(c('STRINGdb', 'clusterProfiler', 'GSEABase', 'GSVA', 'Biobase', 'org.Hs.eg.db', 'AnnotationDbi'), ask=FALSE, update=FALSE)"

# Install radian for better R terminal
RUN pip3 install radian

# Create non-root user for VS Code
RUN groupadd --gid 1000 vscode \
    && useradd --uid 1000 --gid vscode --shell /bin/bash --create-home vscode \
    && echo "vscode ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/vscode \
    && chmod 0440 /etc/sudoers.d/vscode

# Set working directory
WORKDIR /home/vscode/workspace

# Switch to non-root user
USER vscode

# Install additional R packages in user's library
RUN R -e "install.packages('languageserver', repos='https://cran.rstudio.com/')"

# Set default shell
ENV SHELL /bin/bash