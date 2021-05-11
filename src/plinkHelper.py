import sys
import pandas as pd
import numpy as np

fields = ["variant","pval"]
data = pd.read_csv(sys.argv[1], sep = "\t", skipinitialspace=True, usecols=fields) # nrows = 100/Users/zijungao/Desktop/biobankTemp
print("finished reading in minimal p-values")
data_variant = data["variant"]
data_variant = data_variant.str.split(pat = ":", expand = True)
print("finished str.split")
data_variant.columns = ["chromosome", "position", "ref", "alt"] # chr:pos:ref:alt
data["chromosome.position"] = data_variant["chromosome"] + ":" + data_variant[ "position"]
print("finished joining")
print(data[:10])

# reference ID
fields = ["chromosome","position", "ID"]
ref_ID = pd.read_csv(sys.argv[2], skipinitialspace=True, usecols=fields, sep = "\s+") # nrows = 100
print(ref_ID[:10])
ref_ID["chromosome.position"] = ref_ID["chromosome"].astype(str) + ":" + ref_ID[ "position"].astype(str)
print("finished processing ref_ID")

data = pd.merge(left=data, right= ref_ID, how = "inner", left_on='chromosome.position', right_on='chromosome.position')
data = data.dropna() # remove NA
print("finished merging")

data_plink = data[["ID", "pval"]]
data_plink.columns = ["SNP", "P"]
data_plink.to_csv(sys.argv[3], header = True, sep = " ", index = False)
print("finished writing plink file")

data = data[["ID", "chromosome.position", "variant"]]
data.to_csv(sys.argv[4], header = True, sep = "\t", index = False)
