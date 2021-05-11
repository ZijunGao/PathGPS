# preprocessing of the UK Biobank data

ref = read.table("./data/biobankSummary.tsv", sep = "\t", fill = TRUE, header = TRUE)
# 4359 in total
# minimal sample size: 819
# minimal case number: 100

# remove continuous_raw
indexContinuous_raw = which(ref$variable_type == "continuous_raw")
ref = ref[-indexContinuous_raw,] # 4084 left

# remove famale
ref = ref[-4084,] # 4087 left

# remove phenotypes with missing female or male data
femaleDictionary = maleDictionary = numeric(0)
fullAddress = read.table("./data/fullAddress.txt", sep = "\t", header = FALSE)
for(i in 1:dim(fullAddress)[1]){
  name = strsplit(gsub(".*-O ","",fullAddress[i,]),split = "[.]")[[1]]
  if(name[1] %in% ref$phenotype){
    if(name[4] == "female"){
      femaleDictionary = c(femaleDictionary, name[1])
    } else if(name[4] == "male")
      maleDictionary = c(maleDictionary, name[1])
  }
}
nameDictionary = intersect(femaleDictionary, maleDictionary)
indexFemaleMale = which(ref$phenotype %in% nameDictionary)
ref = ref[indexFemaleMale,] # 2719 left

# extract addresses for female and male
# remove redundant v2 phenotypes
femaleAddress = maleAddress = numeric(0)
for(i in 1:dim(fullAddress)[1]){
  name = strsplit(gsub(".*-O ","",fullAddress[i,]),split = "[.]")[[1]]
  if(name[1] %in% ref$phenotype){
    if(name[5] == "v2"){next} # exclude redundant v2
    if(name[4] == "female"){
      femaleAddress = c(femaleAddress, fullAddress[i,])
    } else if(name[4] == "male"){
      maleAddress = c(maleAddress, fullAddress[i,])
    }
  }
}


directory = "/Users/zijungao/Desktop/biobankTemp"
# compute minimal p-values
# can be run in parallel
step = 50
indexStart = seq(1,length(femaleAddress), by = step)
indexEnd = pmin(indexStart + step-1, length(femaleAddress))
terminalFile = paste(directory, "/SNPSelection.txt", sep = "")
pythonFile = "SNPSelectionHelper.py"
write("", file = terminalFile)
for(l in 1:length(indexStart)){ # length(indexStart)
  tempGzFile = paste("temp", l, ".tsv.gz", sep = "")
  tempFile = paste("temp", l, ".tsv", sep = "")
  refFileStart = paste(directory, "/data_ref", ".tsv", sep = "")
  refFile = paste(directory, "/data_ref", l, ".tsv", sep = "")
  addFile = paste(directory, "/", tempFile, sep = "")
  for (i in seq(indexStart[l], indexEnd[l])){ # seq(indexStart[l], indexEnd[l])
    address = substr(femaleAddress[i], start = 6, stop = 10000)
    # download
    write(paste("wget", address, sep = " "), file = terminalFile, append = TRUE)
    # unzip and remove.gz file
    name = gsub(".*-O ","",address)
    newName = paste(substr(name, start = 1, stop = nchar(name)-3),"gz", sep="")
    write(paste("mv", name, tempGzFile, sep = " "), file = terminalFile, append = TRUE)
    write(paste("gunzip", tempGzFile, sep = " "), file = terminalFile, append = TRUE)
    # run the python file on the temp data
    if(i == indexStart[l]){
      write(paste("python", pythonFile, refFileStart, addFile, refFile, sep = " "), file = terminalFile, append = TRUE)
    } else {
      write(paste("python", pythonFile, refFile, addFile, refFile, sep = " "), file = terminalFile, append = TRUE)}
    # remove the temp data
    write(paste("rm", tempFile, sep = " "), file = terminalFile, append = TRUE)
  }
}

# excute the terminalFile in terminal file
# chmod +x ./XXX.txt
#  ./XXX.txt

# merge minimal p-values
terminalFile = "/Users/zijungao/Desktop/biobankTemp/SNPSelectionMerge.txt"
write("", file = terminalFile)
directory = "/Users/zijungao/Desktop/biobankTemp"
refFileStart = paste(directory, "/data_ref", ".tsv", sep = "")
refFileAll = paste(directory, "/data_ref_all", ".tsv", sep = "")
for(i in 1:2){ # length(dataRefName)
  addFile = paste(directory, "/data_ref", i, ".tsv", sep = "")
  if(i == 1){
    write(paste("python", pythonFile, refFileStart, addFile, refFileAll, sep = " "), file = terminalFile, append = TRUE)
  } else {
    write(paste("python", pythonFile, refFileAll, addFile, refFileAll, sep = " "), file = terminalFile, append = TRUE)}
}

# plink
plinkDirectory = "/Users/zijungao/Desktop/Resources/plink_mac_20201019"
terminalFile = paste(directory, "/plink", ".txt", sep = "")
minimalPValueFile = paste(directory, "/data_ref", ".tsv", sep = "")
refIDFile = paste(directory, "/refID", ".txt", sep = "")
plinkFile = paste(plinkDirectory, "/clumpData", sep = "")
pythonFile = paste(directory, "/plinkHelper.py", sep = "")
betaFullFile = paste(directory, "/data_beta_full.tsv", sep = "")

write("", file = terminalFile)
write(paste("python", pythonFile, minimalPValueFile, refIDFile, plinkFile, betaFullFile, sep = " "), file = terminalFile, append = TRUE)
# ~10 min to generate a .clumped file in terminal
write(paste(plinkDirectory, "/plink --bfile refData --clump clumpData --clump-p1 1 --clump-p2 1 --clump-r2 0.001 --clump-kb 10000 --out biobankAll", sep = ""), append = TRUE)

# post processing
terminalFile = paste(directory, "/plink2", ".txt", sep = "")
clumpFile = paste(plinkDirectory, "/biobankAll.clumped" , sep = "")
betaFile = paste(directory, "/data_beta.tsv", sep = "")
write(paste("python", "plinkHelper2.py", betaFullFile, clumpFile, betaFile, sep = " "), file = terminalFile)


# index SNP beta, se, p-value extraction
step = 50
indexStart = seq(1,length(maleAddress), by = step)
indexEnd = pmin(indexStart + step-1, length(maleAddress))
terminalFile = "/Users/zijungao/Desktop/biobankTemp/indexSNPExtraction.txt"
pythonFile = "indexSNPExtractionHelper.py"
write("", file = terminalFile)
for(l in 1:1){ # length(indexStart)
  tempGzFile = paste("temp", l, ".tsv.gz", sep = "")
  tempFile = paste("temp", l, ".tsv", sep = "")
  refFileStart = paste(directory, "/data_beta", ".tsv", sep = "")
  refFile = paste(directory, "/data_beta", l, ".tsv", sep = "")
  addFile = paste(directory, "/", tempFile, sep = "")
  for (i in 1:2){ # seq(indexStart[l], indexEnd[l])
    address = substr(maleAddress[i], start = 6, stop = 10000)
    # download
    write(paste("wget", address, sep = " "), file = terminalFile, append = TRUE)
    # unzip and remove.gz file
    name = gsub(".*-O ","",address)
    newName = paste(substr(name, start = 1, stop = nchar(name)-3),"gz", sep="")
    write(paste("mv", name, tempGzFile, sep = " "), file = terminalFile, append = TRUE)
    write(paste("gunzip", tempGzFile, sep = " "), file = terminalFile, append = TRUE)
    # run the python file on the temp data
    if(i == indexStart[l]){
      write(paste("python", pythonFile, refFileStart, addFile, refFile, nameDictionary[i], sep = " "), file = terminalFile, append = TRUE)
    } else {
      write(paste("python", pythonFile, refFile, addFile, refFile,  nameDictionary[i], sep = " "), file = terminalFile, append = TRUE)}
    # remove the temp data
    write(paste("rm", tempFile, sep = " "), file = terminalFile, append = TRUE)
  }
}


