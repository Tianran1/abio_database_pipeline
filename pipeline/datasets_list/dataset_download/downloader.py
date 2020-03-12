# Copyright 2019-present Analytical Biosciences.
#
# Author: Zedao Liu
# Date: 20190829
# 
#
# use: 
# my_GEO_downloader = GEODownloader(gse_name='GSE102580', working_dir=os.getcwd())
# this build the following directory hierarchy under working_dir
# GSE102580
# ├── GSE102580
# │   ├── downloaded_data
# │   │   ├── GPL_info.json
# │   │   ├── GSE_metadata.json
# │   │   ├── GSE_supplementary_data
# │   │   └── GSM_info.tsv
# │   └── processed_data
# │   │   ├── cellAnnotation.tsv
# │   │   ├── expressionMatrix_normalized.tsv
# │   │   ├── expressionMatrix_rawCounts.tsv
# │   │   ├── expressionMatrix_TPM.tsv
# │   │   ├── geneAnnotation.tsv
# │   │   ├── README.json
# │   │   └── unstructuredData.json
# │   └── code
# │   │   └── script.ipynb
# │   └── SOFT_file
# │       └── GSE102580_family.soft.gz
# └── SOFT_file
#     └── GSE102580_family.soft.gz
#
#
# Maintainence: 
# download target dir might need refinement.
#
#


import GEOparse

import os
import re
import json
from pprint import pprint

import wget

import pandas as pd
import numpy as np
from pathlib import Path

from shutil import copyfile
from os.path import dirname

script_template_dir = "/home/biodb/data/tutorial/template/code"


class GEODownloader(object):
    """top level wrapper"""

    def __init__(self,
                 gse_name: str = None,
                 working_dir: str = os.getcwd(),
                 ):
        self.gse_name = gse_name
        self.working_dir = working_dir

    def download(self):
        root_gse_dir = RootGSEDir(gse_name=self.gse_name, working_dir=self.working_dir)
        root_gse_dir.build_GSE()


class BaseGSEDir(object):
    """
    some GSE series in GEO has subseries,
    so we build an abstract BaseGSEDir class, to hold basic GSE operations,
    including Dir tree construction self._build_dir(),
    GEO downloading self._get_GSE(),
    wrapped into a single command self.build_GSE().
    """

    def __init__(self,
                 gse_name: str = None,
                 dataset_builder: object = None,
                 working_dir: str = os.getcwd(),
                 ):
        """
        The original design aims to include the DatasetBuilder to aid in
        building the Dir tree, but later adopted another strategy.
        For now, just leave the dataset_builder parameter as default None
        """

        # key parameters
        self.builder = dataset_builder

        # gse_name & GSEObject
        self.gse_name = gse_name
        self.GSEObject = None

        # dir
        self.working_dir = working_dir
        self.gse_dir = working_dir + '/' + self.gse_name
        self.SOFT_file_dir = self.gse_dir + '/SOFT_file'

        # info
        self.pubmed_id:int = 0

    def build_GSE(self):

        self._build_dir()
        self.GSEObject = self._get_GSE(gse=self.gse_name, destdir=self.SOFT_file_dir)
        self.pubmed_id = int(self.GSEObject.metadata['pubmed_id'][0])

    def _build_dir(self):
        """
        build a directory starting with GSE within working dir
        """

        if not os.path.exists(self.gse_dir):
            os.makedirs(self.gse_dir)

        if not os.path.exists(self.SOFT_file_dir):
            os.makedirs(self.SOFT_file_dir)

    # wrapper of GEOparse.get_GEO() to use unified parameters
    def _get_GSE(self, gse: str = None, destdir: str = None):
        """download SOFT file which contains GSE information"""

        GSEObject = GEOparse.get_GEO(geo=gse, destdir=destdir, annotate_gpl=True, silent=True)

        return GSEObject


class RootGSEDir(BaseGSEDir):
    """
    the root GSE dir class, representing a concrete GSE dir,
    contains SubGSEDir attributes to initiate tree-like traversal downlaod.
    """

    def __init__(self,
                 gse_name: str = None,
                 dataset_builder: object = None,
                 working_dir: str = os.getcwd(),
                 ):
        """
        extends the parent calss __init__()
        adding additional class attributes sub_gse_name_list & SubGSEDir_list to hold subGSE info
        """
        super().__init__(gse_name=gse_name, dataset_builder=dataset_builder, working_dir=working_dir)

        self.sub_gse_name_list = []
        self.SubGSEDir_list = []

    def build_GSE(self):
        """
        extends the parent build by retrieving sub gse series information
        """

        super().build_GSE()
        self.sub_gse_name_list = self._get_sub_gse_name_list()
        self.SubGSEDir_list = self._get_SubGSEDir_list()

        self._build_sub_gse()

    def _get_sub_gse_name_list(self):
        """
        sub GSE information is available within the GSE objects provided by the GEOParse package.
        GSEObject.metadata['relation']
        """

        # parse root-subseries relationships
        relationships = self.GSEObject.metadata['relation']
        sub_gse_name_list = [self.gse_name]

        print('Sub series info:')
        # subGSE match using regular expression, may need refinement.
        if len(relationships) > 1:
            print('relationship is %s' % str(relationships))
            relationships = sub_gse_name_list + [i for i in relationships if i.startswith('SuperSeries')]
            sub_gse_name_list = [re.search('GSE[0-9]+', i).group(0) for i in relationships]
            if len(sub_gse_name_list) > 1:
                try:
                    sub_gse_name_list.remove(self.gse_name)
                except ValueError:
                    pass

        if len(sub_gse_name_list) > 1:
            print('\t', 'this GSE series contains the following subseries %s' % str(sub_gse_name_list))
        else:
            print('\t', 'this GSE series contains no subseries, only %s' % str(sub_gse_name_list))

        return sub_gse_name_list

    def _get_SubGSEDir_list(self):
        """
        build a list of SubGSEDir objects to help traversal download.
        """

        SubGSEDir_list = []

        SubGSEDir_list = [SubGSEDir(gse_name=i, working_dir=self.gse_dir) for i in self.sub_gse_name_list]

        return SubGSEDir_list

    def _build_sub_gse(self):
        """
        for each SubGSEDir object within the SubGSEDir list, build the SubGSEDir.
        this is the core function link GSEDir & SubGSEDir, and core to the traversal.
        """

        for SubGSEDir in self.SubGSEDir_list:
            print('\n', 'building SubGSEDir: ', SubGSEDir)
            SubGSEDir.build_GSE()


class SubGSEDir(BaseGSEDir):
    """
    SubGSEDir concrete class, representing a SubGSEDir directory
    """

    def __init__(self,
                 gse_name: str = None,
                 dataset_builder: object = None,
                 working_dir: str = os.getcwd(),
                 ):
        """
        extends from the root init

        *** The initial design aims to integrate partioning when download.
        *** But now it seems better to first download the data, and do the partioning later.
        *** The key of this API interface is that:
        *** The GEODownloader API should extract ALL information from GEO,
        *** from which, partioning & curation should confidently start.

        """
        super().__init__(gse_name=gse_name, dataset_builder=dataset_builder, working_dir=working_dir)

        self.parts = []

        self.downloaded_data_dir = self.gse_dir + '/downloaded_data'
        self.GSE_metadata_file = self.downloaded_data_dir + '/GSE_metadata.json'
        self.GSE_supplementary_file_dir = self.downloaded_data_dir + '/GSE_supplementary_data'
        self.GSM_file = self.downloaded_data_dir + '/GSM_info.tsv'
        self.GPL_file = self.downloaded_data_dir + '/GPL_info.json'
        self.supplementary_file_dict = dict()

    def __str__(self):
        return self.gse_name

    def build_GSE(self):

        super().build_GSE()

        self._download_data()

        self.build_dataset()

    def build_dataset(self):

        """ the 2nd parameter: directory of template file 'script.ipynb'
            changes are needed if file path changes """
        ds_builder = Dataset_builder(self.gse_dir, script_template_dir)
        ds_builder.buildDataset()

    def _download_data(self):

        if not os.path.exists(self.downloaded_data_dir):
            os.makedirs(self.downloaded_data_dir)
        print('\n', 'downloading data:')

        # GSE metadata
        metadata = self.GSEObject.metadata
        if not os.path.exists(self.GSE_metadata_file):
            Path(self.GSE_metadata_file).touch()
        with open(self.GSE_metadata_file, 'w') as f:
            json.dump(self.GSEObject.metadata, f)
        print('\t' * 2, 'finished downloading GSE metadata, saved in %s' % self.GSE_metadata_file)

        # GPL
        GPLs = self.GSEObject.gpls
        GPL_dict = dict()
        for gpl in GPLs.keys():
            GPL_dict[gpl] = GPLs[gpl].metadata
        if not os.path.exists(self.GPL_file):
            Path(self.GPL_file).touch()
        with open(self.GPL_file, 'w') as f:
            json.dump(GPL_dict, f)
        print('\t' * 2, 'finished downloading GPL info, saved in %s' % self.GPL_file)

        # GSM sample information
        df_GSM_info = self.GSEObject.phenotype_data
        column_names = df_GSM_info.columns.tolist()
        column_names = [i for i in column_names if not i.startswith('contact')]
        df_GSM_info = df_GSM_info.loc[:, column_names]
        df_GSM_info.insert(loc=0, column='GSM', value=df_GSM_info.index.tolist())
        if not os.path.exists(self.GSM_file):
            Path(self.GSM_file).touch()
        df_GSM_info.to_csv(self.GSM_file, index=False, sep='\t')
        print('\t', 'finished downloading GSM metadata, saved in %s' % self.GSM_file)

        # supplementary data (matrices)
        supplementary_file_list = self.GSEObject.metadata['supplementary_file']
        self.supplementary_file_dict = {i.split('/')[-1]: i for i in supplementary_file_list}
        print('\nsupplementary_file_list:')
        pprint(self.supplementary_file_dict)
        if not os.path.exists(self.GSE_supplementary_file_dir):
            os.makedirs(self.GSE_supplementary_file_dir)
        for file in self.supplementary_file_dict.keys():
            for i in range(5):
                try:
                    print('\n downloading %s:' % file)
                    print(str(i+1))
                    wget.download(url=self.supplementary_file_dict[file], out=self.GSE_supplementary_file_dir)
                    break
                except Exception as e:
                    print('downloading error of %s,' % file, 'error: %s' % e)
                    continue

        ########!!!!!!
        # may need more code to ensure complete download
        ########
        os.system('rm -rf '+self.gse_dir+'/SOFT_file')
        father_gse_dir = os.path.dirname(self.gse_dir)
        os.system('rm -rf '+father_gse_dir+'/SOFT_file')


class Dataset_builder(object):
    """build dir tree"""

    def __init__(self,
                 starting_dir: str = os.getcwd(),
                 script_template_dir: str = None
                 ):
        """  script_template_dir: directory of template file 'script.ipynb' """

        self.starting_dir = starting_dir
        self.script_template_dir = script_template_dir
        self.process_dir = starting_dir + '/processed_data'
        self.code_dir = starting_dir + '/code'
        self.script_dir = self.code_dir + '/script.ipynb'
        self.cellanno_dir = self.process_dir + '/cellAnnotation.tsv'
        self.geneanno_dir = self.process_dir + '/geneAnnotation.tsv'
        self.tpm_dir = self.process_dir + '/expressionMatrix_TPM.tsv'
        self.norm_dir = self.process_dir + '/expressionMatrix_normalized.tsv'
        self.raw_dir = self.process_dir + '/expressionMatrix_rawCounts.tsv'
        self.unstruc_dir = self.process_dir + '/unstructuredData.json'
        self.read_dir = self.process_dir + '/README.json'

    def buildDataset(self):

        def _build_dataset_hierarchy():
            if not os.path.exists(self.process_dir):
                os.makedirs(self.process_dir)
            if not os.path.exists(self.code_dir):
                os.makedirs(self.code_dir)
            return_message = 'successfully generated dir: [1]processed_data [2]code'
            return return_message

        def _copy_script():
            try:
                os.system('cp /home/biodb/data/tutorial/template/code/script.ipynb '+self.script_dir)
                #copyfile(self.script_template_dir, self.script_dir)
            except IOError as e:
                print("'script.ipynb' replication failed")
            except:
                print("Unexpected error occured during 'script.ipynb' replication:", sys.exc_info())

        def _generate_cell_annotation():
            if not os.path.exists(self.cellanno_dir):
                colNames = np.array(
                    ['cellID', 'clusterID', 'clusterName', 'sampleID', 'cellOntologyName', 'cellOntologyID',
                     'meta_sourceName', 'meta_tissue', 'meta_description',
                     'FACSMarker',
                     'tSNE1', 'tSNE2', 'UMAP1', 'UMAP2'])
                dfCellAnnotation = pd.DataFrame(columns=colNames)
                dfCellAnnotation.to_csv(self.cellanno_dir, sep='\t', index=False)
            return_message = 'successfully generated cellAnnotation.tsv'
            return return_message

        def _generate_gene_annotation():
            if not os.path.exists(self.geneanno_dir):
                colNames = np.array(['geneSymbol', 'ensemblID', 'ensemblIDVersion', 'geneID'])
                dfGeneAnnotation = pd.DataFrame(columns=colNames)
                dfGeneAnnotation.to_csv(self.geneanno_dir, sep='\t', index=False)
            return_message = 'successfully generated geneAnnotation.tsv'
            return return_message

        def _generate_expression_matrix():
            if not os.path.exists(self.raw_dir):
                colNames = np.array(['cellID', 'geneSymbol1', 'geneSymbol2'])
                dfExpressionMatrix_rawCounts = pd.DataFrame(columns=colNames)
                dfExpressionMatrix_rawCounts.to_csv(self.raw_dir, sep='\t', index=False)

            if not os.path.exists(self.norm_dir):
                colNames = np.array(['normalizationMethod', 'cellID', 'geneSymbol1', 'geneSymbol2'])
                dfExpressionMatrix_normalized = pd.DataFrame(columns=colNames)
                dfExpressionMatrix_normalized.to_csv(self.norm_dir, sep='\t', index=False)

            if not os.path.exists(self.tpm_dir):
                colNames = np.array(['normalizationMethod', 'cellID', 'geneSymbol1', 'geneSymbol2'])
                dfExpressionMatrix_TPM = pd.DataFrame(columns=colNames)
                dfExpressionMatrix_TPM.to_csv(self.tpm_dir, sep='\t', index=False)

            return_message = 'successfully generated file: [1] expressionMatrix_rawCounts.tsv, [2] expressionMatrix_normalized.tsv'
            return return_message

        def _generate_unstructured_data():
            if not os.path.exists(self.unstruc_dir):
                unstructuredData = {
                    'metadata': {
						"datasetID": "", 
						"subDataset": "", 
						"description": "", 
						"title": "", 
						"authors": "", 
						"accessionNumber": "", 
						"pubmedID": 0, 
						"abstract": "", 
						"sourceID": "", 
						"numberOfCells": 0, 
						"libraryPreparationMethod": "", 
						"sequencingPlatform": "", 
						"clusteringMethod": "", 
						"biomarkerDerivationMethod": "", 
						"fastqURL": "", 
						"figureURL": "", 
						"isFigurePublic": "", 
						"taxonomyID": 0, 
						"genomeBuild": "", 
						"annotation": "", 
						"journal": "", 
						"publicationDate": "", 
						"citation": 0, 
						"tissue": "", 
						"clusterAvailability": "", 
						"disease": "", 
						"methodology": "", 
						"nonAdult": "", 
						"cancer": "", 
						"neuroscience": "", 
						"developmentalBiology": "", 
						"immunology": "", 
						"cellAtlas": "", 
						"isBadtSNE": ""
					},
                    'markerGenes': {
                        'pancreaticACell': {
                            'geneSymbol': ['a', 'b', 'c'],
                            'ensemblID': ['ENSG001', 'ENSG002', 'ENSG003'],
                            'pValue': [0.001, 0.0002, 0.0004],
                            'statistics': [0.99, 0.98, 0.78],
                            'statisticsType': ['FScore', 'FScore', 'FScore']
                        },
                        'pancreaticBCell': {
                            'geneSymbol': ['d', 'e', 'f'],
                            'ensemblID': ['ENSG001', 'ENSG002', 'ENSG003'],
                            'pValue': [0.001, 0.0002, 0.0004],
                            'statistics': [0.99, 0.98, 0.78],
                            'statisticsType': ['FScore', 'FScore', 'FScore']
                        }
                    }
                }
                with open(self.unstruc_dir, "w") as f:
                    json.dump(unstructuredData, f)

            return_message = 'successfully generated unstructuredData.json'
            return return_message

        def _generate_README():
            if not os.path.exists(self.read_dir):
                README = {
                    'qualityControl': 'notPassed',
                    'author': 'AuthorName',
                    'date': '20190724',
					'modificationDate': '20190724',
                    'unfinishedParts': ['markerGenes', 'tSNE'],
                    'authorComments': 'this dataset is part of a larger collection (Tabular Muris)',
                    'otherComments': {}
                }
                with open(self.read_dir, "w") as f:
                    json.dump(README, f)
            return_message = 'successfully generated README.json'
            return return_message

        _build_dataset_hierarchy()
        _copy_script()
        _generate_cell_annotation()
        _generate_gene_annotation()
        _generate_expression_matrix()
        _generate_unstructured_data()
        _generate_README()

        return_message = 'successfully build templates'
        return return_message