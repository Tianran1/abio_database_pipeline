from pymed import PubMed
import os
import pandas as pd
from Bio import Entrez
import time

'''
将全局变量DATE_FROM改成比上次更新时间提前几天，以确保检索不会遗漏，格式"YYYY/MM/DD"。例如，上次10月1日更新，则一个月后11月1日运行该notebook，更改DATE_FROM="2019/09/28"。

更改旧的数据集list的路径OLD_LIST_DIR。
'''

OLD_LIST_DIR="pubmed_dataset_list_unique_new.csv"

NEW_LIST_DIR="pubmed_dataset_list_updated_"+time.strftime("%Y%m%d", time.localtime())+".csv"

DATE_FROM="2019/10/01"

ALLOWED_JOURNALS = [
    'Cancer Cell',
    'Cell',
    'Cell Stem Cell',
    'Cell Syst',
    'Elife',
    'Immunity',
    'Mol. Cell',
    'Nat Biomed Eng',
    'Nat Commun',
    'Nat. Cell Biol.',
    'Nat. Genet.',
    'Nat. Immunol.',
    'Nat. Med.',
    'Nat. Methods',
    'Nat. Neurosci.',
    'Nature',
    'Neuron',
    'Science',
    'Sci Immunol',
    'Sci Transl Med',
    'Cancer Discov'
]

class PubMedAPI():
    """
    wrapper of PubMed API,
    based on the pymed package
    """
    def __init__(self):
        self.pubmed = PubMed(tool="MyTool", email="my@email.address")

    def get_multiple_article_info(self,
                                  query_string: str = "",
                                  max_results: int = 20,
                                  ):
        """
        :param query_string: string, query sentence; example:"single-cell[Title]"
        :param max_results: int, maximum results
        """
        def _construct_query_string_by_journal_and_time():
            """
            pre-construct query_string
            filter: allowed journals and time(DATE_FROM-present)
            """
            query_string_by_journal = '"[Journal] OR "'.join(ALLOWED_JOURNALS)
            query_string_by_journal = '("' + query_string_by_journal + '"[Journal]) AND '
            query_string_by_journal = query_string_by_journal + '(("'+DATE_FROM+'"[Date - Publication] : "3000"[Date - Publication])) AND '
            return query_string_by_journal

        # construct query string (add journal filter)
        query_string = _construct_query_string_by_journal_and_time() + query_string
        results = self.pubmed.query(query_string, max_results)

        title = list()
        journal = list()
        pubmed_id = list()

        for article in results:
            try:
                journal.append(article.journal)
            except AttributeError:
                continue
            # title[:-1:] can delete "." in the end of every title
            title.append(article.title[:-1:])
            # "splitlines()[0]" can deal with occasional format anomaly
            pubmed_id.append(article.pubmed_id.splitlines()[0])

        df = pd.DataFrame(data = {
            "title": title,
            "journal": journal,
            "pubmed_id": pubmed_id
            })

        return df
		

def getGSEFromPubmedID(pubmed_id:str="") ->str:
    try:
        handle =Entrez.elink(dbfrom = 'pubmed', id = pubmed_id, db = 'gds')
        link_results = Entrez.read(handle)
        id_list=list()
        for item in link_results[0]["LinkSetDb"][0]["Link"]:
            id_list.append("GSE"+item['Id'][3::])
        id_list_total = "_".join(id_list)
        return id_list_total
    except:
        return "notAvailable"		


add_list=PubMedAPI().get_multiple_article_info("(single-cell[Title/Abstract] OR scRNA[Title/Abstract])",3000)

gseID=list()
for pubmed_id in add_list["pubmed_id"]:
    gseID.append(getGSEFromPubmedID(pubmed_id))
add_list["gseID"]=gseID

old_list = pd.read_csv(OLD_LIST_DIR,  index_col=0)

#过滤掉已使用的GSE号，如果一篇文章的所有GSE号均已使用，将整篇文章过滤掉。
used_gseID=set()
for each in old_list["gseID"]:
    used_gseID=used_gseID.union(set(each.split("_")))
used_gseID.discard("notAvailable")

add_gse=add_list["gseID"]
processed_add_gse=list()
for each in add_gse:
    if each=="notAvailable":
        processed_add_gse.append(each)
    else:
        processed_each=[i for i in list(each.split("_")) if not i in used_gseID]
        used_gseID=used_gseID.union(set(processed_each))  # update used_gseID
        processed_add_gse.append("_".join(processed_each))

add_list["gseID"]=processed_add_gse
add_list=add_list[add_list["gseID"]!=""]

new_list=old_list.append(add_list)
new_list.drop_duplicates(subset="pubmed_id", keep='first', inplace=True)
new_list.index=range(0,len(new_list))
new_list.to_csv(NEW_LIST_DIR, header='column_names', sep=',')

