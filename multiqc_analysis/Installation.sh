# Install miniconda
# Create a new environment with Python version 3.7
conda create --name py3.7 python=3.7

# Activate the environment
conda activate py3.7

# Install multiqc
conda install -c bioconda -c conda-forge multiqc

# Install fastqc
conda install -c bioconda fastqc
