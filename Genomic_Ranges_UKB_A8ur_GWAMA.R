library(data.table)
library(dplyr)
library(GenomicRanges)

##### A function that an be used to calculate novel regions, not used in analysis, but possibly useful in the future

find_novel <- function(novel, background){
  bck <- makeGRangesFromDataFrame(background)
  bck <- reduce(bck)
  nov <- makeGRangesFromDataFrame(novel)
  nov <- reduce(nov)
  
  nov_bck <- setdiff(nov, bck)
  nov_UNI_bck <- union(nov, bck)
  nov_bck <- as.data.table(nov_bck)
  nov_UNI_bck <- as.data.table(nov_UNI_bck)
  
  return( dplyr::intersect(nov_bck, nov_UNI_bck) )
}

#######

# A8 is post FUMA asthma 8 european all
# A8_u is the new updated post FUMA ashtma 8 european urban
# A8_ur is the post GWAMA and post FUMA results

# import:

A8 <- fread("F:\\NIKO\\Genomics\\UK_Biobank\\FUMA_outputs\\Asthma_age_groups\\FUMA_job104779\\GenomicRiskLoci.txt")

A8u <- fread("F:\\NIKO\\Genomics\\UK_Biobank\\FUMA_outputs\\Asthma8_urban_rural\\FUMA_job98394\\GenomicRiskLoci.txt")

A8_ur <- fread('F:\\NIKO\\Genomics\\UK_Biobank\\FUMA_outputs\\Asthma8_urban_rural\\FUMA_job112561\\GenomicRiskLoci.txt')

# anything new in A8_ur that is not in A8
# anything new in A8u that is not in A8
# anything new in A8_ur GWAMA metanalysis that is not in A8 and not in A8u

A8 <- A8[,.(chr,start,end)]
A8u <- A8u[,.(chr,start,end)]
A8_ur <- A8_ur[,.(chr,start,end)]

##### A8u not in A8:

A81 <- makeGRangesFromDataFrame(A8)
A81 <- reduce(A81)

A8u1 <- makeGRangesFromDataFrame(A8u)
A8u1 <- reduce(A8u1)

A8u_A8 <- setdiff(A8u1, A81)
A8_UNI_A8u <- union(A8u1, A81)

A8u_A8 <- as.data.table(A8u_A8)
A8_UNI_A8u <- as.data.table(A8_UNI_A8u)

A8u_novel <- dplyr::intersect(A8u_A8, A8_UNI_A8u)

rm(A81)
rm(A8u1)

# nothing novel in Asthma 8 urban Eur, sinlge matrix, new correct filter
# compared to Asthma 8 all Eur, mulitmatrix, new correct filter

##### A8_ur not in A8: 
A82 <- makeGRangesFromDataFrame(A8)
A82 <- reduce(A82)

A8_ur2 <- makeGRangesFromDataFrame(A8_ur)
A8_ur2 <- reduce(A8_ur2)

A8ur_A8 <- setdiff(A8_ur2, A82)
A8_UNI_A8ur <- union(A8_ur2, A82)

A8ur_A8 <- as.data.table(A8ur_A8)
A8_UNI_A8ur <- as.data.table(A8_UNI_A8ur)

A8ur_novel <- dplyr::intersect(A8ur_A8, A8_UNI_A8ur)

# nothing novel in GWAMA Asthma 8 urban and rural Eur, sinlge matrix, new correct filter
# compared to Asthma 8 all Eur, mulitmatrix, new correct filter

rm(A82)
rm(A8_ur2)

##### A8_ur not in A8 union A8u (this turned out to be redundant as A8u is contained in A8 and is thus contained in the union of A8 with A8u):

A83 <- makeGRangesFromDataFrame(A8)
A83 <- reduce(A83)

A8u3 <- makeGRangesFromDataFrame(A8u)
A8u3 <- reduce(A8u3)

A8ur3 <- makeGRangesFromDataFrame(A8_ur)
A8ur3 <- reduce(A8ur3)

A8_UNI_A8u <- union(A83, A8u3)
A8ur_minus_A8_UNI_A8u <- setdiff(A8ur3, A8_UNI_A8u)
A8_UNI_A8u_UNI_A8ur <- union(A8ur3, A8_UNI_A8u)

A8ur_minus_A8_UNI_A8u <- as.data.table(A8ur_minus_A8_UNI_A8u) 
A8_UNI_A8u_UNI_A8ur <- as.data.table(A8_UNI_A8u_UNI_A8ur)

A8ur_novel3 <- dplyr::intersect(A8ur_minus_A8_UNI_A8u, A8_UNI_A8u_UNI_A8ur)

# nothing novel in GWAMA Asthma 8 urban and rural Eur, sinlge matrix, new correct filter 
# comapred to Asthma 8 urban, single matrix, and Asthma 8 all Eur, mulitmatrix, both new correct filter


##### A8_u not in A8_ur: 
A8_ur4 <- makeGRangesFromDataFrame(A8_ur)
A8_ur4 <- reduce(A8_ur4)

A8u4 <- makeGRangesFromDataFrame(A8u)
A8u4 <- reduce(A8u4)

A8u_A8ur <- setdiff(A8u4, A8_ur4)
A8u_UNI_A8ur <- union(A8u4, A8_ur4)

A8u_A8ur<- as.data.table(A8u_A8ur)
A8u_UNI_A8ur <- as.data.table(A8u_UNI_A8ur)

A8u_novel <- dplyr::intersect(A8u_A8ur, A8u_UNI_A8ur)

A8u_not_in_A8ur <- find_novel(A8u, A8_ur)

# nothing novel in GWAMA Asthma 8 urban and rural Eur, sinlge matrix, new correct filter
# compared to Asthma 8 all Eur, mulitmatrix, new correct filter

rm(A8u4)
rm(A8_ur4)
