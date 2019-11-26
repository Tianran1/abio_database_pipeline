import importlib.util
import sys
downloader_path = '/home/ztr/abio_database_pipeline/pipeline/datasets_list/dataset_download/downloader.py'
downloader_spec = importlib.util.spec_from_file_location('downloader', downloader_path)
downloader = importlib.util.module_from_spec(downloader_spec)
downloader_spec.loader.exec_module(downloader)

GSE_CODE = ['GSE112004','GSE95601','GSE103275','GSE98131','GSE123135','GSE112294']

downloader.batch_download(gse_names = GSE_CODE, working_dir = '/home/biodb/data/user_admin/unallocated_data')
