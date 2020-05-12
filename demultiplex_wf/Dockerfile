FROM ubuntu:18.04

LABEL maintainer="John Sundh" email=john.sundh@nbis.se
LABEL description="Image to run demultiplexing + concatenation of ITS data"

# Use bash as shell
SHELL ["/bin/bash", "-c"]

WORKDIR /analysis

# Install necessary tools
RUN apt-get update && \
    apt-get install -y --no-install-recommends bzip2 ca-certificates \
        curl unzip wget language-pack-en && apt-get clean

# Install Miniconda and add to PATH
RUN curl https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh -O && \
    bash Miniconda3-py37_4.8.2-Linux-x86_64.sh -bf -p /usr/miniconda3/ && \
    rm Miniconda3-py37_4.8.2-Linux-x86_64.sh && \
    /usr/miniconda3/bin/conda clean --all && \
    ln -s /usr/miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /usr/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Add conda to PATH and set locale
ENV PATH="/usr/miniconda3/bin:${PATH}"
ENV LC_ALL en_US.UTF-8
ENV LC_LANG en_US.UTF-8

# Add conda environment file
COPY environment.yml .

# Install into base environment
RUN conda env update -f environment.yml -n base && conda clean --all

# Add workflow
COPY bin/deML /usr/miniconda3/bin/deML
COPY Snakefile .
COPY common.py .
COPY config.yml .
COPY scripts scripts
COPY info info

CMD ["/bin/bash"]