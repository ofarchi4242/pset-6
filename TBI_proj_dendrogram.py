#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd


# In[2]:


df = pd.read_csv('tbi_matrix\DropSeqTBI.digital_expression.txt', sep='\t')


# In[3]:


df_T = df.T


# In[4]:


# Core scverse libraries
import scanpy as sc
import anndata as ad
import skimage

# Data retrieval
import pooch


# In[5]:


sc.settings.set_figure_params(dpi=50, facecolor="white")


# In[6]:


adata = sc.AnnData(df_T)


# In[7]:


# Loop through matrix rows (cells) and assign cells to sample IDs based on cell names
for cell_name in adata.obs.index:
    if 'Sham' in cell_name:
        adata.obs.at[cell_name, 'sample_id'] = 'Sham'
    elif 'TBI' in cell_name:
        adata.obs.at[cell_name, 'sample_id'] = 'TBI'


# In[8]:


# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("mt-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")


# In[9]:


sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)


# In[10]:


sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)


# In[11]:


# Step 2: Filter cells based on the threshold for mitochondrial gene percentage
threshold = 10  # Set your desired threshold (e.g., less than 10% mitochondrial genes)

# Create a mask for cells that meet the threshold
mito_filter = adata.obs['pct_counts_mt'] < threshold

# Filter cells in the Annotated Data object
adata_filtered = adata[mito_filter]


# In[12]:


# Saving count data
adata_filtered.layers["counts"] = adata_filtered.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata_filtered)
# Logarithmize the data
sc.pp.log1p(adata_filtered)


# In[13]:


sc.pp.highly_variable_genes(adata_filtered, n_top_genes=2000, batch_key="sample_id")


# In[14]:


sc.tl.pca(adata_filtered)


# In[15]:


sc.pp.neighbors(adata_filtered)


# In[16]:


sc.tl.umap(adata_filtered)


# In[17]:


# Perform Principal Component Analysis (PCA)
sc.pp.pca(adata_filtered, n_comps=50, svd_solver='arpack')



# In[18]:


# Perform Leiden clustering using the reduced-dimensional space from PCA
sc.tl.leiden(adata_filtered, resolution=0.5)  # Adjust resolution parameter as needed


# In[19]:


for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(
        adata_filtered, key_added=f"leiden_res_{res:4.2f}", resolution=res
    )


# In[20]:


sc.tl.rank_genes_groups(adata_filtered, groupby="leiden_res_0.50", method="wilcoxon")


# In[21]:


sc.pl.rank_genes_groups_dotplot(
    adata_filtered, groupby="leiden_res_0.50", standard_scale="var", n_genes=5
)


# In[ ]:




