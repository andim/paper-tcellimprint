import os.path
import re
import urllib.request
import pandas as pd

from lib.config import *


path_urls = [# databases of tcrs with annotated specificities
             ('vdjdb.zip', r'https://github.com/antigenomics/vdjdb-db/releases/download/2019-08-08/vdjdb-2019-08-08.zip'),
             # metadata on Britanova dataset
             ('metadata-britanova.txt', 'https://zenodo.org/record/826447/files/metadata.txt'),
             # metadata on Emerson dataset
             ('metadata-emerson.xlsx', 'https://static-content.springer.com/esm/art%3A10.1038%2Fng.3822/MediaObjects/41588_2017_BFng3822_MOESM47_ESM.xlsx'),
             # naive frequencies from Britanova dataset
             ('naive-percentage-britanova.xlsx', 'https://doi.org/10.1371/journal.pcbi.1005572.s016'),
             # full dataset from emerson Nature genetics 2017
             ('emerson-2017-natgen.zip', 'https://s3-us-west-2.amazonaws.com/publishedproject-supplements/emerson-2017-natgen/emerson-2017-natgen.zip')
            ]
for path, url in path_urls:
    if not os.path.exists(data_directory+path):
        print(path)
        urllib.request.urlretrieve(url, data_directory+path)


# download data from Britanova study
zenodobase = 'https://zenodo.org/record/826447/files/'
df = pd.read_csv(data_directory+'metadata-britanova.txt', sep='\t')
for f in df['file_name']:
    path = britanova_rawdata_directory+f+'.gz'
    if not os.path.exists(path):
        print(path)
        url = zenodobase+f+'.gz'
        urllib.request.urlretrieve(url, path)
