# abio_database_pipeline
## how to use
import importlib.util  
import sys  
import os   
import pandas as pd   
import numpy as np  
from scipy.io import mmread, mmwrite  

sys.path.append('/home/biodb/data/personal_folder/ztr/abio_database_pipeline/')  
from pipeline.datasets_curation import datasetBuilder  
from pipeline.datasets_curation import inspector  

my_builder = datasetBuilder.DatasetBuilder()  
my_inspector = inspector.Inspector()  
