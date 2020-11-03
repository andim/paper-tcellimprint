import os
import pandas as pd
import numpy as np
import scipy.stats
import seaborn as sns
import palettable

from lib import *

files = glob.glob(chu_rawdata_directory + '*.tsv.gz')
subjects, dates, phenotypes = [], [], []
for f in files:
    subject, date, phenotype  = re.split("[-_/.]", f)[-5:-2]
    subjects.append(subject)
    dates.append(date)
    phenotypes.append(phenotype)
meta = pd.DataFrame.from_dict(dict(subject=[s[-1] for s in subjects],
                                   date=pd.to_datetime(dates, yearfirst=True),
                                   phenotype=phenotypes,
                                   path=[os.path.split(f)[-1] for f in files]))
meta.to_csv(data_directory + 'metadata-chu.csv', index=False)
