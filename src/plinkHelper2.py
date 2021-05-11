import sys
import pandas as pd
import numpy as np

data_beta = pd.read_csv(sys.argv[1], sep = "\t") # nrows = 100;  skipinitialspace=True,
print("finished reading data_beta")

fields = ["SNP","P", "TOTAL"]
data = pd.read_csv(sys.argv[2], sep = "\s+", usecols=fields) # nrows = 100;
data['SNP'] = data['SNP'].astype(str)
print("finished reading clumped data")

data_beta = pd.merge(left=data, right= data_beta, how = "inner", left_on='SNP', right_on='ID')
del data_beta['ID']
print("finished merging")

data_beta.to_csv(sys.argv[3], header = True, sep = "\t", index = False)
