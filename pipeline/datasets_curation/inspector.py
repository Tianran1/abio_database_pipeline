import os
import json
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import sparse, io
import sys

from pipeline.datasets_curation.datasetBuilder import IOMethod

class VariableControl():
    
    standard_dataset_files = ['processed_data', 'downloaded_data', 'code']
    standard_processed_data_files = ["cellAnnotation.tsv", 
                                    "geneAnnotation.tsv",
                                    "expressionMatrix_normalized.tsv",
                                    "expressionMatrix_rawCounts.tsv",
                                    "expressionMatrix_TPM.tsv",
                                    "unstructuredData.json",
                                    "README.json"]
    optional_processed_data_files = ["scibet.tsv","Diff_genes.json","paga.json"]
    ignored_files = ["SOFT_file",".ipynb_checkpoints","downsample_data","stampTime.json"]
    standard_cellanno_variables = ['cellID', 'clusterID', 'clusterName', 'sampleID', 'cellOntologyName', 
                                   'cellOntologyID', 'FACSMarker', 'tSNE1', 'tSNE2', 'UMAP1', 'UMAP2']
    standard_geneanno_variables = ['geneSymbol','ensemblID','ensemblIDVersion','geneID']
    standard_metadata_variables = ['datasetID', 'subDataset','description','title', 'accessionNumber', 
    'abstract','authors', 'journal', 'sourceID', 'numberOfCells','libraryPreparationMethod', 
    'sequencingPlatform', 'pubmedID', 'clusteringMethod', 'biomarkerDerivationMethod', 'fastqURL', 
    'figureURL', 'taxonomyID','genomeBuild','annotation','publicationDate','citation','tissue',
    'clusterAvailability','disease','methodology','nonAdult','cancer','neuroscience','developmentalBiology',
    'immunology','cellAtlas','isFigurePublic','isBadtSNE']
    libmethod_keywords = ['10x chrominum','drop-seq','microwell-seq','C1 Fluidigm','inDrops',
                              'Smart-seq2','Smart-seq','CEL-seq']
    journal_keywords = ['Cancer cell','Cancer discovery','Cell','Cell stem cell',
                        'Cell system','Immmunity','Molecular cell','Nature',
                        'Nature biomedical engineering','Nature cell biology',
                        'Nature communications','Nature genetics',
                        'Nature medicine','Nature methods',
                        'Nature neuroscience','Neuron','Science','Science immunology',
                        'Science translational medicine','eLife']
    tissue_keywords = ['aorta','bladder','bone marrow','brain','diaphragm','fat','heart','kidney',
                       'large intestine','limb muscle','liver','lung,mammary gland',
                       'pancreas','skin','spleen','thymus','tongue','trachea']

class Inspector(IOMethod, VariableControl):
    
    def __init__(self, starting_dir = os.path.abspath(os.path.dirname(os.getcwd()))):
        super().__init__(starting_dir)
        self.starting_dir = starting_dir
        self.inspection_results = dict()
        self.print_line = dict()

    def inspect_expType(self):
        if 'expressionMatrix_TPM.mtx' in os.listdir(self.process_dir):
            self.exp_type = 'mtx'
        else:
            self.exp_type = 'tsv'
        return self.exp_type

    def inspect(self):
        
        self.inspect_directory_hierarchy()
        self.inspect_cell_annotation()
        self.inspect_gene_annotation()
        self.inspect_unstructured_data()
        self.inspect_expression_matrix()
        self.inspect_README()
        self._print_line()
        
        return self.inspection_results
    
    def inspect_directory_hierarchy(self):
        
        directory_hierarchy_results = dict() 
        
        #============================== 
        # inspect dataset files
        #============================== 
        directory_hierarchy_results['present directory'] = os.getcwd()
        current_dataset_files = os.listdir(self.starting_dir)
        mismatched_file = self._compare_name_sequences(current_dataset_files, self.standard_dataset_files)
        if mismatched_file == None:
            pass
        else:
            for i in self.ignored_files:
                if i in mismatched_file:
                     mismatched_file.remove(i)
            if mismatched_file == set():
                pass
            else:
                dataset_file_result = 'having mismatched file %s' %mismatched_file
                directory_hierarchy_results['dataset_file_result'] = dataset_file_result
        
        #============================== 
        # inspect processed_data names
        #==============================
        processed_data_files = os.listdir(self.process_dir)
        mismatched_file = self._compare_name_sequences(processed_data_files, self.standard_processed_data_files)
        if mismatched_file == None:
            pass
        else:
            for i in self.ignored_files + self.optional_processed_data_files:
                if i in mismatched_file:
                    mismatched_file.remove(i)
            if mismatched_file == set():
                pass
            else:
                processed_data_files_result = 'having mismatched file %s' %mismatched_file
                directory_hierarchy_results['processed_data_files_result'] = processed_data_files_result

        # inspect optional processed_data names
        df_cell = self.read_template_tsv('cellAnnotation.tsv')
        if set(df_cell['clusterName']) != {'notAvailable'}:
            if set(self.optional_processed_data_files) < set(processed_data_files):
                pass
            else:
                directory_hierarchy_results['optional_files_result'] = 'ERROR!should run auto_calculation!'

        self.inspection_results['directory_hierarchy_results'] = directory_hierarchy_results
        return self.inspection_results
   
    def inspect_cell_annotation(self):
        
        cell_annotation_result = dict()
        
        df_cell = self.read_template_tsv('cellAnnotation.tsv')
        
        # inspect mismatched keywords
        mismatched_keywords = self._compare_name_sequences(df_cell.columns.tolist(), self.standard_cellanno_variables)
        if mismatched_keywords != None:
            cell_meta_keywords = set()
            for i in mismatched_keywords:
                if i.startswith('meta'):
                    cell_meta_keywords.add(i)
            mismatched_keywords = mismatched_keywords - cell_meta_keywords
            if mismatched_keywords != set():
                cell_annotation_result['conserved_keyword_result'] = 'having mismatched keyword %s' %mismatched_keywords
            cell_annotation_result['cell_metadata_result'] = cell_meta_keywords
        
        # inpsect na values:
        if df_cell.isna().any().any():
            na_columns = df_gene.columns[df_gene.isna().any()]
            na_column_result = 'ERROR!NA value in the following columns %s' % na_columns
            cell_annotation_result['na_column_result'] = na_column_result

        # inspect cell ontology
        cellOntology_result = set()
        for i in df_cell['cellOntologyID']:
            if i.startswith('CL:') or i == 'notAvailable':
                pass
            else:
                cellOntology_result.add("cellOntologyID doesn't begin with 'CL:'")
        if len(set(df_cell['cellOntologyID'])) != len(set(df_cell['cellOntologyName'])):
            cellOntology_result.add("ERROR!cellOntologyID not matched with cellOntologyName")
        if cellOntology_result == set():
            pass
        else:
            cell_annotation_result['cellOntology'] = cellOntology_result

        # inspect cluster
        cluster_result = set()
        if len(df_cell['clusterName'].unique()) == 1:
            cluster_result.add('Only one cluster: %s!' %df_cell['clusterName'][0])
        if len(df_cell['clusterID'].unique()) != len(df_cell['clusterName'].unique()):
            cluster_result.add("ERROR!clusterID not matched with clusterName")
        for i in set(df_cell['clusterName']):
            if not isinstance(i, str):
                cluster_result.add("clusterName should be str!")
                break
        if cluster_result == set():
            pass
        else:
            cell_annotation_result['cluster'] = cluster_result
            
        # inspect sample meta
        sample_meta = set()
        for i in ['meta_sourceName','meta_tissue','meta_description','PCA1','PCA2']:
            if i in df_cell.columns.tolist():
                if set(df_cell[i]) == {'notAvailable'} or df_cell[i].isna().all():
                    sample_meta.add(i)
        if sample_meta == set():
            pass
        else:
            cell_annotation_result['meta_info'] = 'ERROR!please delete %s'%sample_meta
            
        # inspect missing dim_red
        dim_red_names = ['tSNE1', 'tSNE2', 'UMAP1', 'UMAP2']
        for i in dim_red_names:
            if df_cell[i].isna().any():
                dim_red_result = 'lack dim red'
                cell_annotation_result['dim_red_result'] = dim_red_result
                break
        
        
        self.inspection_results['cell_annotation_result'] = cell_annotation_result
        self.print_line['cell_annotation'] = df_cell[['cellID', 'clusterID', 'clusterName', 'sampleID', 'cellOntologyName', 'cellOntologyID']].iloc[0,:].tolist()
        
        return self.inspection_results
    
    def inspect_gene_annotation(self):
        
        gene_annotation_result = dict()

        # inspect mismatched keywords
        df_gene = self.read_template_tsv('geneAnnotation.tsv')
        mismatched_keywords = self._compare_name_sequences(df_gene.columns.tolist(), self.standard_geneanno_variables)
        if mismatched_keywords != None:
            conserved_keyword_result = 'having mismatched keyword %s' %mismatched_keywords
            gene_annotation_result['conserved_keyword_result'] = conserved_keyword_result
        
        # inspect geneSymbol & ensemblID
        if df_gene.iloc[:,:1].isna().any().any():
            na_columns = df_gene.columns[df_gene.isna().any()]
            na_column_result = 'ERROR!NA value in the following columns %s' % na_columns
            gene_annotation_result['na_column_result'] = na_column_result
        
        null_columns = []
        if set(df_gene['geneSymbol']) == {'notAvailable'}:
                null_columns.append('geneSymbol') 
        if set(df_gene['ensemblID']) == {'notAvailable'}:
                null_columns.append('ensemblID')
        if null_columns != []:
            gene_annotation_result['null_column_result'] = 'ERROR!Null value in %s' %null_columns
      
        self.inspection_results['gene_annotation_result'] = gene_annotation_result
        self.print_line['gene_annotation'] = df_gene.iloc[0,:].tolist()

        return self.inspection_results
    
    def inspect_expression_matrix(self):
        
        raw_expression_matrix_result = dict()
        normalized_expression_matrix_result = dict()
        TPM_expression_matrix_result = dict()
        
        #============================== 
        # inspect raw expression matrix
        #==============================
        with open(self.raw_dir, 'r') as file:
            column_names = file.readline()
            first_line = file.readline()
        try:
            assert column_names.split()[0] == 'cellID'
        except AssertionError:
            raw_expression_matrix_result['cellID_keyword_error'] = 'ERROR!rawCounts: wrong cellID keyword'
        try:
             for i in first_line.split('\t')[1:]:
                j = int(i)
        except ValueError:
            raw_expression_matrix_result['data_type_error'] = 'ERROR!rawCounts: not integer format'
                
        self.inspection_results['raw_expression_matrix_result'] = raw_expression_matrix_result

        #===================================== 
        # inspect normalized expression matrix
        #===================================== 
        with open(self.norm_dir, 'r') as file:
            column_names = file.readline()
            first_line = file.readline()
            genes_norm = column_names.split('\t')[2:]
            if first_line != '':
                normalized_cell_length = 1
                while True:
                    x = file.readline()
                    if x == '':
                        break
                    else: 
                        normalized_cell_length = normalized_cell_length + 1
                try:
                    assert column_names.split()[0] == 'normalizationMethod'
                except AssertionError:
                    normalized_expression_matrix_result['normalized_keyword_error'] = 'ERROR!normalized: wrong normalization keyword' 
                try:
                    assert column_names.split()[1] == 'cellID'
                except AssertionError:
                    normalized_expression_matrix_result['normalized_keyword_error'] = 'ERROR!normalized: wrong cellID keyword'
                try:
                    for i in first_line.split('\t')[2:]:
                        type(eval(i)) == float
                except:
                    normalized_expression_matrix_result['data_type_error'] = 'ERROR!normalized: not float format'     
                if first_line != '':
                    normalized_expression_matrix_result['normalizationMethod'] = first_line.split('\t')[0]
                        
        self.inspection_results['normalized_expression_matrix_result'] = normalized_expression_matrix_result
        
        #============================== 
        # inspect TPM expression matrix
        #============================== 
        with open(self.tpm_dir, 'r') as file:
            column_names = file.readline()
            first_line = file.readline()
            genes_TPM = column_names.split('\t')[2:]
            if first_line != '':
                TPM_cell_length = 1
                TPM_cellID = [first_line.split('\t')[1]]
                while True:
                    x = file.readline()
                    if x == '':
                        break
                    else: 
                        TPM_cell_length = TPM_cell_length + 1
                        TPM_cellID.append(x.split('\t')[1])
                try:
                    assert column_names.split()[0] == 'normalizationMethod'
                except AssertionError:
                    TPM_expression_matrix_result['tpm_keyword_error'] = 'ERROR!TPM: wrong normalization keyword' 
                try:
                    assert column_names.split()[1] == 'cellID'
                except AssertionError:
                    TPM_expression_matrix_result['tpm_keyword_error'] = 'ERROR!TPM: wrong cellID keyword'    
                try:
                    for i in first_line.split('\t')[2:]:
                        j = float(i)
                except ValueError:
                    TPM_expression_matrix_result['data_type_error'] = 'ERROR!TPM: not float format'  
                try:
                    if round(sum([float(x) for x in first_line.split('\t')[2:]])) == 1000000:
                        pass
                    else:
                        TPM_expression_matrix_result['number_error'] = 'ERROR!TPM: rowSums not equal to 1000000'
                except:
                    TPM_expression_matrix_result['number_error'] = 'ERROR!TPM: Wrong format!'
                TPM_expression_matrix_result['normalizationMethod'] = first_line.split('\t')[0]

        self.inspection_results['TPM_expression_matrix_result'] = TPM_expression_matrix_result
    
        #=================== 
        # inspect dimension
        #===================
        dimension_result = dict()

        df_cell = self.read_template_tsv('cellAnnotation.tsv')
        cell_length = df_cell.shape[0]
        df_gene = self.read_template_tsv('geneAnnotation.tsv')
        gene_length = df_gene.shape[0]
        
        # normalized
        if os.path.getsize(self.norm_dir) > 100:
            if normalized_cell_length != cell_length:
                dimension_result['normalized_cell_length_result'] = 'ERROR!not matched with cellAnnotation'
            if len(genes_norm) != gene_length:
                dimension_result['normalized_gene_length_result'] = 'ERROR!not matched with geneAnnotation'
        
        # TPM
        if os.path.getsize(self.tpm_dir) > 100:
            if TPM_cell_length != cell_length:
                dimension_result['TPM_cell_length_result'] = 'ERROR!not matched with cellAnnotation'
            if TPM_cellID != df_cell['cellID'].tolist():
                dimension_result['TPM_cellID_result'] = 'ERROR!not matched with cellAnnotation'
            if len(genes_TPM) != gene_length:
                dimension_result['TPM_gene_length_result'] = 'ERROR!not matched with geneAnnotation' 
        
        self.inspection_results['dimension_result'] = dimension_result
        return self.inspection_results        
    
    def inspect_unstructured_data(self):
        
        unstructured_data_result = dict()
        
        unstructured_data = self.read_template_json('unstructuredData.json')
            
        #============================== 
        # inspect metadata
        #==============================

        #========checking standard keys for mismatch ========
        metadata = unstructured_data['metadata']
        current_keys = [*metadata]
        mismatched_names = self._compare_name_sequences(current_keys, self.standard_metadata_variables)

        metadata_result = dict()
        if mismatched_names != None:
            metadata_result['mismatched_keyword'] = 'having mismatched keyword %s' %mismatched_names

        #========checking standard keys for datatype & content ========
        keyword_content = []
        integer_keywords = ['pubmedID', 'numberOfCells', 'taxonomyID','citation']
        string_keywords = list(set(self.standard_metadata_variables[1:]) - set(integer_keywords))
        url_keywords = ['sourceID', 'fastqURL', 'figureURL']

        try:
            for i in integer_keywords:
                if type(metadata[i]) != int:
                    keyword_content.append('ERROR!%s should be integer not string' %i)
                if metadata[i] == 0:
                    keyword_content.append('no content in %s' %i)
            for i in string_keywords:
                if metadata[i] == '':
                    keyword_content.append('no content in %s' %i)
            for i in url_keywords:
                if not metadata[i].startswith('http'):
                    keyword_content.append('%s keywords should startwith http' %i)
            if not metadata['libraryPreparationMethod'] in libmethod_keywords:
                keyword_content.append('libraryPreparationMethod:%s should be in the keywords list' %metadata['libraryPreparationMethod'])   
            if not metadata['journal'] in journal_keywords:
                keyword_content.append('journal:%s should be in the keywords list' %metadata['journal'])
            if not metadata['tissue'] in tissue_keywords:
                keyword_content.append('tissue:%s should be in the keywords list' %metadata['tissue'])
            if not isinstance(metadata['sequencingPlatform'],str):
                keyword_content.append('sequencingPlatform should be str')
        except:
            keyword_content.append('might have mismatched keywords')

        metadata_result['keyword_content'] = keyword_content

        unstructured_data_result['metadata_result'] = metadata_result

        #============================== 
        # inspect marker_genes
        #==============================

        marker_genes = unstructured_data['markerGenes']
        df_cell = self.read_template_tsv('cellAnnotation.tsv')
        if marker_genes == {}:
            unstructured_data_result['marker_gene_contents'] = 'no markers genes!'
        else:
            standard_keys = []
            cell_types = [*marker_genes]
            mismatched_Names = None

            for cell_type in cell_types:
                markers = marker_genes[cell_type]
                current_keys = [*markers]
                standard_keys = ['geneSymbol', 'ensemblID', 'pValue', 'statistics', 'statisticsType']
                mismatched_names = self._compare_name_sequences(current_keys, standard_keys)
                for key in current_keys:
                    try:
                        if type(markers[key]) != list or len(markers[key]) < 5:
                            empty_marker_genes = True
                            empty_marker_genes_error.append('problem in %s %s, this generally occurs as it should be a list(vector in R) here but not a character'
                                                %(cell_type, key))
                    except Exceptions as e:
                        empty_marker_genes.append(e)

            if mismatched_names != None:
                marker_gene_keywords = 'having mismatched keyword %s' %mismatched_names
                unstructured_data_result['marker_gene_keywords'] = marker_gene_keywords    
            eid = marker_genes[cell_types[0]]['ensemblID']
            try:
                if pd.unique(eid) == ['notAvailable']:
                    unstructured_data_result['ensemblID'] = 'ERROR!lack ensemblID'
            except:
                pass

        self.inspection_results['unstructured_data'] = unstructured_data_result

        return self.inspection_results
        
    def inspect_README(self):
        
        README_result = dict()
        
        with open(self.read_dir, 'r') as file:
            README = json.load(file)
            
            if README['qualityControl'] != 'notPassed':
                README_result['qualityControl'] = 'if you are not admin, do not rewrite the notPassed stamp.'
        
        README_result['unfinished'] = README['unfinishedParts']
        README_result['comments'] = README['authorComments']
        self.inspection_results['README_result'] = README_result
        
        return self.inspection_results

    def _compare_name_sequences(self, current_sequence, standard_sequence):
       
        """
        Parameters:
        ----------
            current_sequence, standard_sequence: list()
        Returns:
        ----------
            mismatched_names: set()in
        """
        
        if set(standard_sequence) == set(current_sequence):
            pass
        else: 
            mismatched_names = set(current_sequence).symmetric_difference(set(standard_sequence))
            return mismatched_names
        
    def _print_line(self):
        
        self.inspection_results['print_line'] = self.print_line 
        
        return self.inspection_results

    def _red_style(self, string: str = ""):
        
        return '\033[1;31m' + string + '\033[0m' 