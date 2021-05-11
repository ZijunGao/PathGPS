import sys
import pandas as pd
import numpy as np

fields = ["variant","pval"]
data_ref = pd.read_csv(sys.argv[1], sep = "\t", skipinitialspace=True, usecols=fields) # nrows = 10;
data_add = pd.read_csv(sys.argv[2], sep = "\t", skipinitialspace=True, usecols=fields) # nrows = 10;

data_merge = pd.merge(left=data_ref, right=data_add, how = "inner", left_on='variant', right_on='variant')

data_merge['pval'] = np.minimum(data_merge['pval_x'], data_merge['pval_y'])
del data_merge['pval_x']
del data_merge['pval_y']

data_merge.to_csv(sys.argv[3], header = True, sep = "\t")
