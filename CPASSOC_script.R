library(data.table)
library(dplyr)

source("C:\\MY_FOLDERS\\Genomics Project\\CPASSOC\\CPASSOC\\CPASSOC\\FunctionSet.R")

in_file_1 <- "path to first qc filtered SIAGE summary stats file "
in_file_2 <- "path to second qc filtered SIAGE summary stats file "

# this could work for BOLT too, but one needs to use the proper column headers
# in this file I use SAIGE headers, not BOLT
# if you are using a mix of SAIGE and BOLT files, then maybe make sure allele columns match!
# that would be the case if SAIGE and BOLT operate on UKBiobank files (the same dataset)
# If you have used our qc filter to filter your files, then p.value has been replaced by p_value
# there is a new qc filter with improved qc filtering script
# at all three servers, these are 
# on cedar and graham: 
# QC-Filter: /project/6048803/UKBiobank/QC_filters/EUR_QC_filters/EUR_QC_filter.tsv
# QC-Filtering script: /project/6048803/UKBiobank/QC_filters/EUR_QC_filters/EUR_QC_filtering_script.sh
# on beluga: 
# QC-Filter: lustre03/project/6048803/UKBiobank/QC_filters/EUR_QC_filters/EUR_QC_filter.tsv
# QC-Filtering script: lustre03/project/6048803/UKBiobank/QC_filters/EUR_QC_filters/EUR_QC_filtering_script.sh


f_1 <- fread(in_file_1)
N_1 <- f_1[1, N]
f_1 <- f_1[,.(rsid, CHR, POS, Allele1, Allele2, BETA, SE, p_value)]


gc()
f_2 <- fread(in_file_2)
N_2 <- f_2[1, N]
f_2 <- f_2[,.(rsid, Allele1, Allele2, BETA, SE, p_value)]


gc()

f_1[, Z:=BETA/SE]
f_2[, Z:=BETA/SE]

gc()

f_12 <- as.data.table(inner_join(f_1, f_2, by = c("rsid", "Allele1", "Allele2")))
# rm(f_1)
# rm(f_2)

gc()


### Preprocessing setup:

s <- c(N_1, N_2)

X <- f_12[,.(Z.x, Z.y)]
CorrMA <- cor(X)


### Shet calculation

start_time <- Sys.time()

### the main important (longest) calculation:

para <- EstimateGamma(N = 1E4, SampleSize = s, CorrMatrix = CorrMA)
j <- SHet(X = X, SampleSize = s, CorrMatrix = CorrMA)
p_Shet <- pgamma(q = j-para[3], shape = para[1], scale = para[2], lower.tail = F);

###

end_time <- Sys.time()
print(end_time - start_time) #53.13861 min


### Add the p_Shet column to the merged data table f_12

f_12 <- data.table(f_12 , p_Shet=p_Shet)

###

gc()

### Shom calculation:

start_time <- Sys.time()
a <- SHom(X = X, SampleSize = s, CorrMatrix = CorrMA)
p_Shom <- pchisq(a, df = 1, ncp = 0, lower.tail = F)

end_time <- Sys.time()
print(end_time - start_time)

### Add the p_Shet column to the merged data table f_12

f_12 <- data.table(f_12 , p_Shom=p_Shom)

###

### rename the columns 
colnames(f_12) <- c("rsid","CHR","POS","Allele1","Allele2","BETA_1","SE_1","p_value_1",
                    "Z_1","BETA_2","SE_2","p_value_2","Z_2","p_Shet","p_Shom")

### generate output:

write.table(f_12, 
            "path to the output", 
            append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)
