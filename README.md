Overview

This repo contains python code to generate a dot plot that identifies the top 5 marker genes for each cluster identified in our TBI and Sham data.

Data

The data was generated initially from a read counts matrix that was created from a single-cell RNA sequencing experiments using TBI and Sham mouse models. Using this matrix,
we created an AnnData object used in Scanpy to perform quality control and normalization on the single-cell data. Then, we used PCA to reduce the dimensionality and 
performed Leiden clustering to identify cell clusters in the data. Then, we ranked the genes in the each cluster to identify marker genes of each cluster. The dot plot shows
the marker genes for each cluster. We will then use this data to identify the cell type that each cluster represents.

Folder structure

The data is located in the tbi_matrix folder, with the file name DropSeqTBI.digital_expression.txt. All of the other files are just in the repository.

Installation

The TBI_proj_dendrogram.py file can be ran using Python 3.12. The required packages are located in the requirements.txt file, which is also pasted below:

anndata==0.10.6
numpy==1.26.4
pandas==2.2.1
pooch==1.8.1
scanpy==1.10.0
skimage==0.0
