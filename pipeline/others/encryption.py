#!/usr/bin/env python
# coding: utf-8

# In[30]:


import importlib.util
import sys
import os 
import pandas as pd
import numpy as np
from scipy.io import mmread, mmwrite

BUILDER_PATH = '/home/biodb/data/abio_database_pipeline/datasetBuilder.py'

builder_spec = importlib.util.spec_from_file_location('builder', BUILDER_PATH)
builder = importlib.util.module_from_spec(builder_spec)
builder_spec.loader.exec_module(builder)


# In[31]:


class Encryption(builder.DatasetBuilder):
    
    def add_matrix(self):
        
        tpm = self.read_template_tsv('expressionMatrix_TPM.tsv')
        tpm['Analytical_Biosciences'] = 0
        self.save_template_tsv(tpm, 'expressionMatrix_TPM.tsv')
        added_matrix = ['TPM']
        
        if os.path.getsize(self.norm_dir) < 100:
            pass
        else:
            norm = self.read_template_tsv('expressionMatrix_normalized.tsv')
            norm['Analytical_Biosciences'] = 0
            self.save_template_tsv(norm, 'expressionMatrix_normalized.tsv')
            added_matrix.append('normalized')
            
        if os.path.getsize(self.raw_dir) < 100:
            pass
        else:
            raw = self.read_template_tsv('expressionMatrix_rawCounts.tsv')
            raw['Analytical_Biosciences'] = 0
            self.save_template_tsv(raw, 'expressionMatrix_rawCounts.tsv')
            added_matrix.append('rawCounts')
        return added_matrix
    
    def add_gene_anno(self):
        gene = self.read_template_tsv('geneAnnotation.tsv')
        if 'Analytical_Biosciences' not in gene['geneSymbol'].tolist():
            gene = gene.append({'geneSymbol':'Analytical_Biosciences','ensemblID':'notAvailable','ensemblIDVersion':'','geneID':''}, ignore_index=True)
            self.save_template_tsv(gene, 'geneAnnotation.tsv')
        return 'geneAnno'
    
    def add_meta(self):
        uns = self.read_template_json('unstructuredData.json')
        meta = uns['metadata']
        meta['copyright'] = 'Analytical Biosciences'
        uns['metadata'] = meta
        self.save_template_json(uns, 'unstructuredData.json')
        return 'metadata'
    
    def encryption(self):
        if 'expressionMatrix_TPM.mtx' in os.listdir(self.process_dir):
            pass
        else:
            self.add_matrix()
            self.add_gene_anno()
        self.add_meta()


# In[32]:


path = '/home/biodb/data/dataset_collection/datasets/3_standard_dataset/'
for i in os.listdir(path):
    if i.startswith('No_'):
        if i != 'No_1':
            dataset_path = path + i
            for j in os.listdir(dataset_path):
                if j.startswith('part_'):
                    starting_dir = path + i + '/' + j
                    try:
                        x = Encryption(starting_dir)
                        x.add_gene_anno()
                    except Exception as e:
                        print([i,j,e])


# In[ ]:




