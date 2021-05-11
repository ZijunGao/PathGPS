import sys
import pandas as pd
import numpy as np

data_beta = pd.read_csv(sys.argv[1], sep = "\t")

fields = ["variant","beta","se","pval"]
data_add = pd.read_csv(sys.argv[2], sep = "\t", skipinitialspace=True, usecols=fields) # nrows = 10

data_merge = pd.merge(left=data_beta, right=data_add, how = "inner", left_on='variant', right_on='variant')

betaNew = "".join(('beta.',sys.argv[4]))
seNew = "".join(('se.',sys.argv[4]))
pvalNew = "".join(('pval.',sys.argv[4]))
data_merge2 = data_merge.rename(columns={'beta': betaNew, 'se': seNew, 'pval': pvalNew})
data_merge2.to_csv(sys.argv[3], header = True, index = False, sep = "\t")
