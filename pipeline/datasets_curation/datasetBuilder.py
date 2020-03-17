import os
import json
import pandas as pd
import numpy as np
from pathlib import Path
import scanpy as sc
import anndata as an
import GEOparse
from typing import Optional, Union, Dict, Any, List,Tuple
from collections import Counter
from pymed import PubMed
from wos import WosClient
import wos.utils
import xml.dom.minidom
from scipy import sparse, io
import seaborn as sns
import sys
import louvain
import umap

ref_dir = os.path.realpath(sys.modules['pipeline'].__file__).split('__')[0] + 'resources/'
gene_ref_dir = ref_dir + 'gene_references/'

class IOMethod(object):

    def __init__(self, starting_dir: str = ""):
        self.down_dir = starting_dir + '/downloaded_data'
        self.process_dir = starting_dir + '/processed_data'
        self.code_dir = starting_dir + '/code'
        self.cellanno_dir = self.process_dir + '/cellAnnotation.tsv'
        self.geneanno_dir = self.process_dir + '/geneAnnotation.tsv'
        self.tpm_dir = self.process_dir + '/expressionMatrix_TPM.tsv'
        self.norm_dir = self.process_dir + '/expressionMatrix_normalized.tsv'
        self.raw_dir = self.process_dir + '/expressionMatrix_rawCounts.tsv'
        self.unstruc_dir = self.process_dir + '/unstructuredData.json'
        self.read_dir = self.process_dir + '/README.json'
        self.downloaded_metadata_dir = self.down_dir + "/GSE_metadata.json"

    def read_template_tsv(self, tsv_name=None):
        """
        Parameters:
        -----------
            tsv_name: str()
                template tsv file to be read.

        Returns:
        -----------
            df: pd.DataFrame()
        """

        df = pd.read_csv(self.process_dir + '/' + tsv_name, sep='\t')

        return df

    def save_template_tsv(self, dataframe=None, tsv_name=None):
        """
        Parameters:
        -----------
            dataframe: pandas.DataFrame()
            tsv_name: str()

        Returns:
        -----------
            return_message: str()
        """

        dataframe.to_csv(self.process_dir + '/' + tsv_name, sep='\t', index=False)

        return_message = 'overwritten %s' % tsv_name
        return return_message

    # use json.load()
    def read_template_json(self, json_name=None):
        """
        Parameters:
        -----------
            json_name: str()

        Returns:
        -----------
            my_json: dict()

        """
        with open(self.process_dir + '/' + json_name, 'r') as file:
            my_json = json.load(file)

        return my_json

    # use json.dump()
    def save_template_json(self, json_dict=None, json_name=None):
        """
        Parameters:
        -----------
            json_dict: dict()
            json_name: str()

        Returns:
        -----------
            return_message: str()

        """

        with open(self.process_dir + '/' + json_name, 'w') as file:
            json.dump(json_dict, file)

        return_message = 'overwritten json'
        return return_message
    
    # sparce mtx format
    def read_template_mtx(self, mtx_name=None):
        """
        Parameters:
        -----------
            mtx_name: str()

        Returns:
        -----------
            return_message: str()

        """
        mtx = io.mmread(self.process_dir + '/' + mtx_name)

        return mtx
    
    def save_template_mtx(self, mtx=None, mtx_name=None, genes=None, cellID=None):
        """
        Parameters:
        -----------
            mtx_name: str()
            genes: list()
            cellID: list()

        Returns:
        -----------
            return_message: str()

        """

        if not sparse.isspmatrix_coo(mtx):
            return 'not coo sparse matrix'
        try:
            if len(genes) != mtx.shape[1]: 
                return 'genes not matched with matrix'
            elif len(cellID) != mtx.shape[0]:
                return 'cellID not matched with matrix'
        except:
            pass
        if mtx_name == 'expressionMatrix_rawCounts.mtx':
            mtx_type = 'rawCounts'
            cell_raw = pd.DataFrame()
            cell_raw['cellID'] = cellID
            self.save_template_tsv(cell_raw,'cellID_' + mtx_type + '.tsv')
            gene_raw = pd.DataFrame()
            gene_raw['genes'] = genes
            self.save_template_tsv(gene_raw,'genes_' + mtx_type + '.tsv')

        if mtx_name == 'expressionMatrix_TPM.mtx':
            mtx_type = 'TPM'    

        io.mmwrite(self.process_dir + '/expressionMatrix_' + mtx_type + '.mtx', mtx)
            
        return_message = 'sparse matrix saved'
        return return_message  

class BaseDatasetBuilder(IOMethod):

    def buildDataset(self):

        def _build_dataset_hierarchy():
            if not os.path.exists(self.down_dir):
                os.makedirs(self.down_dir)
            if not os.path.exists(self.process_dir):
                os.makedirs(self.process_dir)
            if not os.path.exists(self.code_dir):
                os.makedirs(self.code_dir)
            return_message = 'successfully generated dir: [1] downloaded_data, [2] processed_data, [3] code'
            return return_message

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
                        'datasetID': '',  # No_1
                        'subDataset': '',
                        'description': '',
                        'correspondingFigure': '',
                        'title': '',
                        'authors':'',
                        'accessionNumber': '',  # GSE11111
                        'pubmedID': 0,  # 123456,
                        'keywords':[],
                        'abstract': '',  # 
                        'sourceID': '',  # 'https://doi.org/10.1016/j.celrep.2018.11.056',
                        'numberOfCells': 0,  # 300,
                        'libraryPreparationMethod': '',
                        # 10x chrominum, drop-seq, microwell-seq, C1 Fluidigm, inDrops, Smart-seq2, CEL-seq
                        'sequencingPlatform': '',  # 'Illumina HiSeq 2500'
                        'clusteringMethod': '',  # 'louvain',
                        'biomarkerDerivationMethod': '',  # 有cluster分类的填: wilcoxon，没有的填: notAvailable
                        'fastqURL': '',  # 'https://www.ebi.ac.uk/ega/datasets/EGAD00001003910'
                        'figureURL': '',  # ''
                        'isFigurePublic': '', # ''
                        'taxonomyID': 0,  #
                        'genomeBuild': '',  # 'hg38'
                        'annotation': '',  # 'ensembl_v93'
                        'journal': '',
                        'publicationDate':'',
                        'citation': 0,
                        'tissue':[''],
                        'tissueOntology':[''],
                        'clusterAvailability':'',
                        'disease':'',
                        'methodology':'',
                        'cancer':'',
                        'neuroscience':'',
                        'developmentalBiology':'',
                        'immunology':'',
                        'cellAtlas':'',
                        'tSNEAvailability':'',
                        'isBadtSNE':''
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

                """
                auto-generation of unstructuredData:
                      title
                      accessionNumber
                      abstract
                      pubmedID
                      taxonomyID
                """
                
                try: 
                    with open(self.downloaded_metadata_dir, 'r') as f:
                        downloaded_metadata = json.load(f)
                        unstructuredData["metadata"]["title"] = downloaded_metadata["title"][0]
                        unstructuredData["metadata"]["accessionNumber"] = downloaded_metadata["geo_accession"][0]
                        unstructuredData["metadata"]["abstract"] = downloaded_metadata["summary"][0]
                        unstructuredData["metadata"]["pubmedID"] = int(downloaded_metadata["pubmed_id"][0])
                        unstructuredData["metadata"]["taxonomyID"] = int(downloaded_metadata["sample_taxid"][0])
                except:
                    pass

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
                    'modificationDate':'20190724',
                    'unfinishedParts': ['markerGenes', 'tSNE'],
                    'authorComments': 'this dataset is part of a larger collection (Tabular Muris)',
                    'otherComments': {}
                }
                with open(self.read_dir, "w") as f:
                    json.dump(README, f)
            return_message = 'successfully generated README.json'
            return return_message

        _build_dataset_hierarchy()
        _generate_cell_annotation()
        _generate_gene_annotation()
        _generate_expression_matrix()
        _generate_unstructured_data()
        _generate_README()

        return_message = 'successfully build templates'
        return return_message

    # should consider, pd.read_csv(index_col=0)

class ScibetClassifier:
    
    def __init__(self):
        self.reference_core = None
        self.reference_genes = None
        self.reference_cell_types = None
        
    def calculate_core(self, TPM, genes, cell_types):
        """
        Parameters:
        ----------
            TPM: np.array()
            genes: np.array()
            cell_types: np.array()
            
        Returns:
        ----------
            reference_core: np.array()
                log1p transformed probability
            reference_genes: np.array()
                selected genes, now 2000
            reference_cell_types: np.array()
                labels
        """
        # use pd.Dataframe for easier calculations
        df_TPM = pd.DataFrame(TPM)
        df_TPM.index = cell_types
        df_TPM.columns = genes
        
        # averaging
        df_TPM_averaged = df_TPM.groupby(by=df_TPM.index).mean()
        self.df_TPM_averaged = df_TPM_averaged
        
        # gene filtering
        selected_genes = []
        for gene in df_TPM_averaged.columns:
            if gene.startswith('m') or gene.startswith('RP'):
                select_gene = False
            else:
                select_gene = True
            selected_genes.append(select_gene)
        
        df_TPM_filtered_mt_rp = df_TPM_averaged.loc[:, selected_genes]
        df_TPM_filtered = df_TPM_filtered_mt_rp.loc[:, df_TPM_filtered_mt_rp.mean(0) < 2000]
        informative_genes, informative_t_scores, informative_TPM, _ =             select_informative_genes(df_TPM_filtered.to_numpy(), df_TPM_filtered.columns.to_numpy(), number_of_genes=2000)
        df_TPM_informative_genes = pd.DataFrame(informative_TPM, columns=informative_genes, index=df_TPM_filtered.index)
        self.df_TPM_informative_genes = df_TPM_informative_genes
        df_TPM_informative_genes += 1

        # calculate probability and log transform
        df_prob = df_TPM_informative_genes.divide(df_TPM_informative_genes.sum(1), axis=0)
        df_prob_log = np.log1p(df_prob)
        
        self.reference_core = df_prob_log.to_numpy()
        self.reference_genes = df_prob_log.columns.to_numpy()
        self.reference_cell_types = df_prob_log.index.to_numpy()

        return self.reference_core, self.reference_genes, self.reference_cell_types
    
    def predict_cell_type(self, test_TPM, test_genes):
        """
        Parameters:
        -----------
            test_TPM: np.array()
            test_genes: np.array()
        Returns:
        -----------
            preddicted_cell_types: np.array()
        """
        df_test = pd.DataFrame(test_TPM, columns=test_genes)
        df_test_selected_genes = df_test.reindex(columns=self.reference_genes, fill_value=0)
        df_test_logged = np.log1p(df_test_selected_genes)
        df_pred = pd.DataFrame(np.dot(df_test_logged, self.reference_core.T), columns=self.reference_cell_types)
        predicted_cell_types = df_pred.T.idxmax(0).to_numpy()
        self.predicted_cell_types = predicted_cell_types
        
        return self.predicted_cell_types

    def generate_scibet_cell_types(self,
        TPM : np.ndarray = None,
        genes: list = None,
        taxonomy_id: int = None):
        df_reference = pd.read_csv(ref_dir + 'GSE11111_scibet_core.csv', index_col=0)
        reference_core = df_reference.to_numpy()
        mouse_reference_genes = df_reference.columns.tolist()
        reference_cell_types = df_reference.index.tolist()
        
        df_gene_conversion = pd.read_csv(ref_dir + 'homomuris.csv')
        mouse_to_human_dict = {df_gene_conversion['mgi'][i]:df_gene_conversion['symbol'][i] for i in range(df_gene_conversion.shape[0])}
        human_reference_genes = []
        for i in mouse_reference_genes:
            try: 
                human_reference_genes.append(mouse_to_human_dict[i])
            except:
                human_reference_genes.append('notAvailable')
        
        self.reference_core = reference_core
        self.reference_cell_types = reference_cell_types
        
        if taxonomy_id == 10090:
            self.reference_genes = mouse_reference_genes
            predicted_cell_types = self.predict_cell_type(TPM, genes)
            return predicted_cell_types

        elif taxonomy_id == 9606:
            self.reference_genes = human_reference_genes
            predicted_cell_types = self.predict_cell_type(TPM, genes)
            return predicted_cell_types

        else:
            pass

    def generate_scibet_HCL(self,
        TPM : np.ndarray = None,
        genes: list = None,
        taxonomy_id: int = None):

        if taxonomy_id == 9606:
            df_reference = pd.read_csv(ref_dir + 'HCL_model.csv', index_col=0)
            reference_core = df_reference.to_numpy()
            reference_genes = df_reference.columns.tolist()
            reference_cell_types = df_reference.index.tolist()
            self.reference_core = reference_core
            self.reference_cell_types = reference_cell_types
            self.reference_genes = reference_genes
            predicted_cell_types = self.predict_cell_type(TPM, genes)
            return predicted_cell_types
        else:
            pass

class Auto_calculation():
    def __init__(self):
        pass

    def calculate_diff_genes(self, TPM, genes, cell_types):
        """
        Parameters:
            TPM: np.array()
            genes: list()
            cell_types: list()
        """
        ann = an.AnnData(np.log1p(TPM))
        ann.var_names = genes
        ann.obs['cluster'] = cell_types
        diffGenes = dict()
        for each in ['t-test','wilcoxon']:
            markers = dict()
            sc.tl.rank_genes_groups(ann,'cluster',method = each,n_genes=1000,rankby_abs=True)
            for x in set(cell_types):
                markers[x] = {
                    'geneSymbol': ann.uns['rank_genes_groups']['names'][x].tolist(),
                    'logFC': ann.uns['rank_genes_groups']['logfoldchanges'][x].tolist(),
                    'statistics': ann.uns['rank_genes_groups']['scores'][x].tolist(),
                    'pValue': ann.uns['rank_genes_groups']['pvals'][x].tolist(),
                    'qValue': ann.uns['rank_genes_groups']['pvals_adj'][x].tolist()}
            diffGenes[each] = markers
        
        return diffGenes

    def select_informative_genes(self,
        TPM : Optional[Union[np.ndarray]] = None,
        genes : Optional[Union[np.ndarray]] = None,
        number_of_genes : Optional[int] = 0,
        labels : Optional[list] = []):

        """
        Returns:
        informative_genes: np.array()
        informative_t_scores: np.array()
        informative_TPM: np.array()
        p_values: np.array()
        """

        # unsupervised gene selection
        if len(labels) == 0:
            labeled_TPM = TPM

        # supervised gene selection
        elif len(labels) != 0:

            groups = set(labels)
            group_TPMs = []

            for group in groups:
                group_indexes = [index for index, x in enumerate(labels) if x == group]
                group_TPM = TPM[group_indexes, :]
                group_TPM = np.mean(group_TPM, axis=0)
                group_TPMs.append(group_TPM)

            labeled_TPM = np.vstack(group_TPMs)
            # print('\n\n\n\n\n\n', labeled_TPM.shape, '\n\n\n\n')

        # calculate t_scores    
        TPM_1 = labeled_TPM + 1    
        log_E = np.log2(np.mean(TPM_1, axis=0))
        E_log = np.mean(np.log2(TPM_1), axis=0)
        t_scores = log_E - E_log

        # return values
        descending_indexes = t_scores.argsort()[::-1][:number_of_genes]
        informative_genes = genes[descending_indexes]
        informative_t_scores = t_scores[descending_indexes]
        informative_t_scores = np.around(informative_t_scores, decimals=1)

        informative_genes_indexes = []
        for i in informative_genes.tolist():
            informative_genes_indexes.append(genes.tolist().index(i))
        informative_TPM = TPM[:, informative_genes_indexes]

        p_vals = []
        if len(labels) != 0:
            label_sets = set(labels)
            for i in range(number_of_genes):
                expression_vector = informative_TPM[:,i]
                grouped_expressions = []
                for label in label_sets:
                    grouped_expression = [expression_vector[i] for i in range(len(expression_vector)) if labels[i]==label]
                    grouped_expressions.append(grouped_expression)
                p_vals.append(scipy.stats.f_oneway(*grouped_expressions)[1])

        p_vals = np.array(p_vals)
    
        return informative_genes, informative_t_scores, informative_TPM, p_vals
    
    def calculate_scibet(self, TPM, genes, cell_types):
        """
        Parameters:
            TPM: np.array()
            genes: list()
            cell_types: list()
        """
        # use pd.Dataframe for easier calculations
        df_TPM = pd.DataFrame(TPM)
        df_TPM.index = cell_types
        df_TPM.columns = genes

        # averaging
        df_TPM_averaged = df_TPM.groupby(by=df_TPM.index).mean()

        # gene filtering
        selected_genes = []
        for gene in df_TPM_averaged.columns:
            if gene.startswith('m') or gene.startswith('RP'):
                select_gene = False
            else:
                select_gene = True
            selected_genes.append(select_gene)

        df_TPM_filtered_mt_rp = df_TPM_averaged.loc[:, selected_genes]
        df_TPM_filtered = df_TPM_filtered_mt_rp.loc[:, df_TPM_filtered_mt_rp.mean(0) < 2000]
        informative_genes, informative_t_scores, informative_TPM, _ = self.select_informative_genes(
            TPM = df_TPM_filtered.to_numpy(), genes = df_TPM_filtered.columns.to_numpy(), number_of_genes=2000)
        df_TPM_informative_genes = pd.DataFrame(informative_TPM, columns=informative_genes, index=df_TPM_filtered.index)

        df_TPM_informative_genes += 1

        # calculate probability and log transform
        df_prob = df_TPM_informative_genes.divide(df_TPM_informative_genes.sum(1), axis=0)       
        df_prob_log = np.log1p(df_prob)

        reference_core = df_prob_log.to_numpy()
        reference_genes = df_prob_log.columns.to_numpy()
        reference_cell_types = df_prob_log.index.to_numpy()
        
        scibet = pd.DataFrame(reference_core)
        scibet.columns = reference_genes
        scibet.index = reference_cell_types
        
        return scibet

    def calculate_paga(self, TPM, genes, cell_types):

        df_genes = genes
        informative_genes, informative_t_scores, informative_TPM, p_vals = self.select_informative_genes(TPM,genes,2000)
        df_paga = pd.DataFrame(informative_TPM,cell_types,informative_genes)
        adata = an.AnnData(df_paga.to_numpy(), obs={"cell_types": df_paga.index}, var={"genes": df_paga.columns})
        n = 50
        sc.tl.pca(adata, svd_solver='arpack', n_comps=n)
        sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
        sc.tl.diffmap(adata)
        sc.pp.neighbors(adata, n_neighbors=5, use_rep='X_diffmap')

        my_annotations = []
        cell_types_list = adata.obs['cell_types'].tolist()
        for i in range(0, len(cell_types_list)):
            my_annotations.append(cell_types_list[i])

        adata.obs["my_annotations"] = my_annotations
        adata.obs['my_annotations'] = adata.obs['my_annotations'].astype('category')

        sc.tl.paga(adata, groups='my_annotations')
        paga_info = adata.uns['paga']
        paga_nodes = adata.obs["my_annotations"].cat.categories.tolist()

        paga_connectivities = paga_info["connectivities"].toarray()
        my_connectivities = paga_connectivities
        unique = []
        for i in range(0, my_connectivities.shape[0]):
            for j in range(0, my_connectivities.shape[1]):
                if my_connectivities[i][j] != 0:
                    if my_connectivities[i][j] not in unique:
                        unique.append(my_connectivities[i][j])
                    else:
                        my_connectivities[i][j] = 0


        node_name = [node for node in paga_nodes]
        counts = Counter(my_annotations)
        node_size = [counts[node] for node in paga_nodes]

        #my_connectivities = np.around(my_connectivities, 3).tolist()

        connectivities = []
        for each in range(len(my_connectivities)):
            for num in range(len(my_connectivities[each])):
                if my_connectivities[each][num] != 0:
                    connectivities.append([each,num,my_connectivities[each][num]])

        paga_json = dict()
        paga_json['node_name'] = node_name
        paga_json['node_size'] = node_size
        paga_json['connectivities'] = connectivities
        
        return paga_json
    
    def calculate_scibetHCL(self, TPM, genes, taxonomy_id):
        """
        TPM: np.array()
        genes: np.array()
        """
        my_classifier = ScibetClassifier()
        return my_classifier.generate_scibet_HCL(TPM, genes, taxonomy_id)

class DerivationMethod(IOMethod):

    def get_metadata(self, pubmed_id: str = ""):
        """
        Parameters: 
            pubmed_id : str() "28152281"
        Returns: 
            journal: str() "Nature"
            publication_date: str() "2019-10-08"
            keywords:list
            authors: [
                {'lastname': 'Ryskov',
                'firstname': 'A P',
                'initials': 'AP'
                },
                ...
            ]
        """

        # get journal and publication_date via PubMed API
        pubmed = PubMed(tool="MyTool", email="my@email.address")
        articles = pubmed.query(str(pubmed_id), max_results=1)
        for i in articles:
            article = i
        keywords = article.keywords
        journal = article.journal
        if journal == "Science (New York, N.Y.)":
            journal = "Science"

        publication_date = str(article.publication_date)
        authors = article.authors
        for i in range(len(authors)):
            authors[i].pop('affiliation')

        return journal, publication_date, authors, keywords
    
    def gene_ref(self):
        
        species = self.read_template_json('unstructuredData.json')['metadata']['taxonomyID']
        if species == 9606:
            df = pd.read_csv(gene_ref_dir + 'human_unique_id_length.tsv', sep='\t')
        elif species == 10090:
            df = pd.read_csv(gene_ref_dir + 'mouse_unique_id_length.tsv', sep='\t')
        elif species == 7955:
            df = pd.read_csv(gene_ref_dir + 'zebra_fish_unique_id_length.tsv', sep='\t')
        elif species == 7227:
            df = pd.read_csv(gene_ref_dir + 'fruitfly_unique_id_length.tsv', sep='\t')
        elif species == 9541:
            df = pd.read_csv(gene_ref_dir + 'macaque_unique_id_length.tsv', sep='\t')
        elif species == 9598:
            df = pd.read_csv(gene_ref_dir + 'chimpanzee_unique_id_length.tsv', sep='\t')
        elif species == 9031:
            df = pd.read_csv(gene_ref_dir + 'chicken_unique_id_length.tsv', sep='\t')
        elif species == 6239:
            df = pd.read_csv(gene_ref_dir + 'Caenorhabditis_unique_id_length.tsv', sep='\t')   
        elif species == 9823:
            df = pd.read_csv(gene_ref_dir + 'pig_unique_id_length.tsv', sep='\t') 
        elif species == 10116:
            df = pd.read_csv(gene_ref_dir + 'rat_unique_id_length.tsv', sep='\t') 
        return df

    def toTPM(self, 
              gene_len=None, 
              exp_type='tsv', 
              ensemblID=False):
        
        raw = self.read_template_tsv('expressionMatrix_rawCounts.tsv')
        cellid = raw['cellID'].tolist()
        names = raw.columns.tolist()[1:]
        raw = raw.iloc[:, 1:].to_numpy()
        raw = sparse.csr_matrix(raw)
        
        if cellid == []:
            exp_type = 'mtx'
            raw = self.read_template_mtx('expressionMatrix_rawCounts.mtx')
            raw = raw.tocsr()
            cellid = self.read_template_tsv('cellID_rawCounts.tsv')['cellID'].tolist()
            names = self.read_template_tsv('genes_rawCounts.tsv')['genes'].tolist()
        
        names_index = np.where([not i.startswith('ERCC-') for i in names])[0]
        names = [i for i in names if not i.startswith('ERCC-')]
        raw = raw[:, names_index]
        
        libmethod = self.read_template_json('unstructuredData.json')['metadata']['libraryPreparationMethod']
        libmethod_keywords = ['10x chromium','drop-seq','microwell-seq','C1 Fluidigm','inDrops',
                          'Smart-seq2','Smart-seq','CEL-seq','CEL-seq2','MARS-seq','msSCRB-seq','SCRB-seq']
        full_len_keywords = ['C1 Fluidigm', 'Smart-seq2','Smart-seq']
        if libmethod in libmethod_keywords:
            if libmethod in full_len_keywords:
                gene_len = True
            else:
                gene_len = False
        
        ref_dict = {}
        if gene_len:
            df = self.gene_ref()
            if names[round(len(names)/2)].startswith('ENS', 0) or ensemblID:
                for i, x in enumerate(df['Gene_ID'].tolist()):
                    ref_dict[x] = df['Gene_length'].tolist()[i]
            else:
                for i, x in enumerate(df['Gene_name'].tolist()):
                    ref_dict[x] = df['Gene_length'].tolist()[i]

            leng = []
            for i, name in enumerate(names):
                try:
                    leng.append(int(ref_dict[name]))
                except:
                    leng.append(0)

            leng = np.array(leng)
            import statistics as st
            leng[np.where(leng == 0)] = st.median(np.array(leng[leng > 0]))
            length = [leng[i] for i in raw.indices]
            raw.data = raw.data/np.array(length)
        
        raw = raw.tocsc()
        rsum = raw.sum(axis=1).tolist()
        rsums = [rsum[i][0] for i in raw.indices]
        raw.data = raw.data/np.array(rsums)
        tpm = raw * 1e6
        tpm.data[np.where(tpm.data == np.inf)] = 0
        #tpm = tpm.fillna(0)
        
        if exp_type == 'tsv':
            tpm = tpm.todense()
            tpm = pd.DataFrame(tpm)
            tpm.columns = names
            tpm.insert(0, "normalizationMethod", "TPM from raw data")
            tpm.insert(1, "cellID", cellid)
            self.save_template_tsv(tpm, 'expressionMatrix_TPM.tsv')
        else:
            tpm = tpm.tocoo()
            self.save_template_mtx(mtx = tpm, mtx_name = 'expressionMatrix_TPM.mtx')
            self.generate_geneAnno(genes = names)
            
        return 'TPM generated'

    def calculate_ensemblID(self, gene=None):

        ref_dict = {}
        df = self.gene_ref()
        for i, symbol in enumerate(df['Gene_name'].tolist()):
            ref_dict[symbol] = df['Gene_ID'].tolist()[i]

        IDs = []
        for i, name in enumerate(gene):
            try:
                IDs.append(ref_dict[name])
            except:
                IDs.append('notAvailable')

        return IDs

    def calculate_geneSymbol(self, gene=None):

        ref_dict = {}
        df = self.gene_ref()
        for i, ID in enumerate(df['Gene_ID'].tolist()):
            ref_dict[ID] = df['Gene_name'].tolist()[i]

        names = []
        for i, ID in enumerate(gene):
            try:
                names.append(ref_dict[ID])
            except:
                names.append('notAvailable')

        return names

    def generate_geneAnno(self,ensemblID = None,genes = None):

        df_gene_anno = pd.DataFrame()
        if genes == None:
            with open(self.tpm_dir, 'r') as file:
                x = file.readline()[:-1]
            genes = x.split('\t')[2:]

        if genes[round(len(genes)/2)].startswith('ENS', 0):
            ensemblID = True
        if ensemblID: 
            names = self.calculate_geneSymbol(gene=genes)
            df_gene_anno['geneSymbol'] = names
            df_gene_anno['ensemblID'] = genes
            df_gene_anno['ensemblIDVersion'] = ''
            df_gene_anno['geneID'] = ''
        else:
            IDs = self.calculate_ensemblID(gene=genes)
            df_gene_anno['geneSymbol'] = genes
            df_gene_anno['ensemblID'] = IDs
            df_gene_anno['ensemblIDVersion'] = 'V90'
            df_gene_anno['geneID'] = ''

        df_gene_anno = df_gene_anno.fillna('notAvailable')
        self.save_template_tsv(df_gene_anno, 'geneAnnotation.tsv')
        
    def calculate_dim_red(self, tSNE=False, UMAP=False, fast = False):
        
        self.results = dict()
        df_TPM = self.read_template_tsv('expressionMatrix_TPM.tsv')
        if df_TPM['cellID'].tolist() == []:
            TPM = self.read_template_mtx('expressionMatrix_TPM.mtx')
            TPM = TPM.tocsr()
            TPM = TPM.log1p()
            fast = True
        else:
            TPM = df_TPM.iloc[:, 2:].to_numpy()
            TPM = np.log(TPM + 1) 
         
        #genes = df_TPM.columns.to_numpy()[2:]
        df_cell_anno = self.read_template_tsv('cellAnnotation.tsv')
        if TPM.shape[0] < 50:
            n = df_cell_anno.shape[0]
        else:
            n = 50
        print('TPM loaded')

        ann = an.AnnData(TPM)
        sc.pp.highly_variable_genes(ann, n_top_genes=500)
        sc.pp.pca(ann, n_comps=n, zero_center=True)
        X_pca = ann.obsm['X_pca']
        tSNE_init = X_pca[:, :2]
        print('feature selection and PCA compression finished ')

        if UMAP:
            reducer = umap.UMAP()
            X_embedded = reducer.fit_transform(X_pca)
            df_cell_anno['UMAP1'] = X_embedded[:, 0].tolist()
            df_cell_anno['UMAP2'] = X_embedded[:, 1].tolist()
            self.results['UMAP1'] = X_embedded[:, 0].tolist()
            self.results['UMAP2'] = X_embedded[:, 1].tolist()
            print('UMAP finished')

        if tSNE:
            if not fast:
                from sklearn.manifold import TSNE
                X_embedded = TSNE(n_components=2, metric='cosine', init=tSNE_init).fit_transform(X_pca)
                df_cell_anno['tSNE1'] = X_embedded[:, 0].tolist()
                df_cell_anno['tSNE2'] = X_embedded[:, 1].tolist()
                self.results['tSNE1'] = X_embedded[:, 0].tolist()
                self.results['tSNE2'] = X_embedded[:, 1].tolist()
                print('tSNE finished')
            else:
                sys.path.append(os.path.realpath(sys.modules['pipeline'].__file__).split('__')[0] + 'resources/software/FIt-SNE')
                from fast_tsne import fast_tsne
                X_embedded = fast_tsne(X_pca, perplexity=50, initialization=tSNE_init)
                df_cell_anno['tSNE1'] = X_embedded[:, 0].tolist()
                df_cell_anno['tSNE2'] = X_embedded[:, 1].tolist()
                print('tSNE finished')

        self.save_template_tsv(df_cell_anno, 'cellAnnotation.tsv')
    
    def tSNEplot(self):
        """
        judge good or bad tSNE
        """
        
        df_cell = self.read_template_tsv('cellAnnotation.tsv')
        if len(df_cell['clusterName'].unique().tolist()) == 1:
            print('cluster not available') 
        else:
            sns.scatterplot(x = df_cell['tSNE1'], y = df_cell['tSNE2'], hue = df_cell['clusterName'], legend = False)
    
    def auto_calculation(self,diff_genes=True,exp_type='tsv'):
        """
        some functions applied directly to the data automatically: markerGenes, scibet, paga
        using TPM,cellAnno and geneAnno data
        """ 
        
        my_calculation = Auto_calculation()
        taxonomy_id = self.read_template_json('unstructuredData.json')['metadata']['taxonomyID']
        df_cell = self.read_template_tsv('cellAnnotation.tsv')

        if len(df_cell['clusterName'].unique().tolist()) == 1:
            message = 'clusterName ERROR!'
            return message
        else:
            
            TPM = self.read_template_tsv('expressionMatrix_TPM.tsv')
            df_gene = self.read_template_tsv('geneAnnotation.tsv')
            y = df_gene.index[np.logical_or(df_gene['geneSymbol'] == 'notAvailable',df_gene['geneSymbol'] == 'Analytical_Biosciences')].tolist()
            df_gene = df_gene.drop(y,axis=0)
            genes = df_gene['geneSymbol'].to_list()
        
            if TPM['cellID'].tolist() == []:
                exp_type = 'mtx'
                TPM = self.read_template_mtx('expressionMatrix_TPM.mtx')
                TPM = TPM.toarray()
            else:
                TPM = TPM.iloc[:, 2:].to_numpy()
            TPM = np.delete(TPM, y, axis=1)

            # calculate scibet HCL
            if taxonomy_id == 9606:
                df_cell['meta_scibetHCL'] = my_calculation.calculate_scibetHCL(TPM, genes, taxonomy_id)
                self.save_template_tsv(df_cell,'cellAnnotation.tsv')

            # remove cluster with only one cell and notAvailable clusters or geneSymbols
            count = pd.value_counts(df_cell['clusterName'])
            try:
                t = count.index[count == 1].tolist()
                y1 = []
                for i in t:
                    y1.append(df_cell.index[df_cell['clusterName'] == i].tolist()[0])
                    
                TPM = np.delete(TPM, y1, axis=0)
                df_cell = df_cell.drop(y1, axis=0)
                df_cell = df_cell.reset_index(drop=True)
            except:
                pass
            y2 = df_cell.index[df_cell['clusterName'] == 'notAvailable'].tolist()
            TPM = np.delete(TPM, y2, axis=0)
            df_cell = df_cell.drop(y2, axis=0)
            
            df_gene = df_gene.fillna('notAvailable')
            genes = df_gene['geneSymbol'].to_numpy()
            cell_types = df_cell['clusterName'].astype(str).to_numpy()
            
            if diff_genes:
                diff = my_calculation.calculate_diff_genes(TPM, genes, cell_types)
                self.save_template_json(diff,'Diff_genes.json')
            
            scibet = my_calculation.calculate_scibet(TPM, genes, cell_types)
            scibet.to_csv(self.process_dir + '/scibet.tsv',sep = '\t')
            paga = my_calculation.calculate_paga(TPM, genes, cell_types)
            self.save_template_json(paga,'paga.json')
    
    def sample_info(self, GSE=None):
        gse = GEOparse.get_GEO(GSE)
        info = gse.phenotype_data
        return (info)

    def calculate_cluster(self, RUN=None):
        
        my_classifier = ScibetClassifier()
        df_cell = self.read_template_tsv('cellAnnotation.tsv')
        taxonomy_id = self.read_template_json('unstructuredData.json')['metadata']['taxonomyID']

        if set(df_cell['clusterName']) == {'notAvailable'} or RUN == True:
            df_gene = self.read_template_tsv('geneAnnotation.tsv')
            y = df_gene.index[df_gene['geneSymbol'] == 'notAvailable'].tolist()
            genes = df_gene.drop(y,axis=0)['geneSymbol'].to_list()
            if os.path.exists(self.process_dir+"/expressionMatrix_TPM.mtx"):
                if df_cell.shape[0] > 500000:
                    print('Too many cells !')
                else:
                    TPM = self.read_template_mtx('expressionMatrix_TPM.mtx')
                    TPM = TPM.todense()
                    TPM = np.delete(TPM, y, axis=1)
                    ann1 = an.AnnData(np.log1p(TPM))
                    sc.pp.highly_variable_genes(ann1, n_top_genes=500)
                    sc.pp.pca(ann1, n_comps=50, zero_center=True)
                    sc.pp.neighbors(ann1, n_neighbors=15, n_pcs=40)
                    sc.tl.louvain(ann1)
                    cluster = ann1.obs['louvain'].tolist()
                    cluster = [int(i)+1 for i in cluster]
                    df_cell['clusterID'] = cluster
                    df_cell['clusterName'] = ['cluster_'+str(i) for i in cluster]
                    df_cell['clusteringMethod'] = 'louvain'
                    if taxonomy_id == 10090 or taxonomy_id == 9606:
                        df_cell['clusterName_scibet'] = my_classifier.generate_scibet_cell_types(TPM, genes, taxonomy_id)
                    if taxonomy_id == 9606:
                        df_cell['meta_scibetHCL'] = my_classifier.generate_scibet_HCL(TPM, genes, taxonomy_id)
                    self.save_template_tsv(df_cell,'cellAnnotation.tsv')
            else:
                TPM = self.read_template_tsv('expressionMatrix_TPM.tsv')
                TPM = TPM.iloc[:, 2:].to_numpy()
                TPM = np.delete(TPM, y, axis=1)
                ann1 = an.AnnData(np.log1p(TPM))
                sc.pp.highly_variable_genes(ann1, n_top_genes=500)
                sc.pp.pca(ann1, n_comps=50, zero_center=True)
                sc.pp.neighbors(ann1, n_neighbors=15, n_pcs=40)
                sc.tl.louvain(ann1)
                cluster = ann1.obs['louvain'].tolist()
                cluster = [int(i)+1 for i in cluster]
                df_cell['clusterID'] = cluster
                df_cell['clusterName'] = ['cluster_'+str(i) for i in cluster]
                df_cell['clusteringMethod'] = 'louvain'
                if taxonomy_id == 10090 or taxonomy_id == 9606:
                        df_cell['clusterName_scibet'] = my_classifier.generate_scibet_cell_types(TPM, genes, taxonomy_id)
                if taxonomy_id == 9606:
                    df_cell['meta_scibetHCL'] = my_classifier.generate_scibet_HCL(TPM, genes, taxonomy_id)
                self.save_template_tsv(df_cell,'cellAnnotation.tsv')

class DatasetBuilder(BaseDatasetBuilder, DerivationMethod):
    def __init__(self, starting_dir = os.path.abspath(os.path.dirname(os.getcwd()))):
        super().__init__(starting_dir)
