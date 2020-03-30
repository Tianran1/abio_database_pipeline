import os
import numpy as np
import pandas as pd
from typing import Optional, Union, Dict, Any, List,Tuple
from fbpca import pca
from geosketch import gs
#import algorithm

from pipeline.datasets_curation.datasetBuilder import IOMethod

class Downsample(IOMethod):

    def __init__(self, starting_dir = os.path.abspath(os.path.dirname(os.getcwd()))):
        super().__init__(starting_dir)
        self.starting_dir = starting_dir

    def tpm_downsample(self):

        if os.path.exists(self.process_dir+"/expressionMatrix_TPM.mtx"):
            TPM_npy = self.read_template_mtx("expressionMatrix_TPM.mtx")
            U,s,Vt = pca(TPM_npy,k=100)
            TPM_dimred = U[:,:100] * s[:100]
            N = 4000
            sketch_index = gs(TPM_dimred, N, replace = False)
            print(sketch_index)
            TPM_downsample =pd.DataFrame(TPM_npy.toarray()[sketch_index,:])
            
        else:
            TPM = self.read_template_tsv("expressionMatrix_TPM.tsv")
            # 考虑去除用于加密的基因
            # 
            TPM_data = TPM.iloc[:,2:]
            TPM_npy = TPM_data.values
            U,s,Vt = pca(TPM_npy,k=100)
            TPM_dimred = U[:,:100] * s[:100]
            N = 4000
            sketch_index = gs(TPM_dimred, N, replace = False)
            TPM_downsample = TPM_data.loc[sketch_index,:]

        if not os.path.exists(self.starting_dir + "/downsample_data"):
            os.makedirs(self.starting_dir + "/downsample_data")
        TPM_downsample.to_csv(self.starting_dir+"/downsample_data/cut_expressionMatrix_TPM.tsv",sep = "\t",index = False)

        return sketch_index
    
    # execute
    def downsample(self, tpm_downsampled = False):
        cellAnnotation = self.read_template_tsv("cellAnnotation.tsv")
    
    
        if cellAnnotation.shape[0] <= 4000:
            return_message = "no need for downsample"
            return return_message
        else:
            if not tpm_downsampled:
                sketch_index = self.tpm_downsample()
                cellAnnotation_downsample = cellAnnotation.loc[sketch_index,:]
            else:
                cellAnnotation = cellAnnotation.set_index('cellID',drop = False)
                cellID_downsample = pd.read_csv(self.starting_dir + "/downsample_data/cut_cellAnnotation.tsv",sep='\t')['cellID'].tolist()
                cellAnnotation_downsample = cellAnnotation.loc[cellID_downsample,:]
            
            cellAnnotation_downsample.to_csv(self.starting_dir+"/downsample_data/cut_cellAnnotation.tsv",sep = "\t",index = False)
                    