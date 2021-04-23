library(data.table)
library(dplyr)

# chronic pain summary stats table, qc filtered, for EUR from bgen sample
mcp <- fread('F:\\NIKO\\Genomics\\UK_Biobank\\Sumstats\\qc_filtered\\Chronic_Pain\\chronic_pain_BOLT_qc_FUMA.txt')
# genomic rsik loci file from FUMA
grl <- fread('F:\\NIKO\\Genomics\\UK_Biobank\\FUMA_outputs\\chronic_pain_EUR\\FUMA_job98963\\GenomicRiskLoci.txt')

View(mcp[c(1:100),])

names(grl)[4] <- 'CHR'

# extract only the snp realted columns from mcp, ignoring stats:
smcp <- mcp[,.(CHR, BP, SNP, ALLELE0, ALLELE1)]
View(smcp[c(1:100),])

names(smcp)[c(2,3,4,5)] <- c('POS', 'rsID', 'Allele1', 'Allele2') 

# set up initial empty data table:
loci <- data.table(CHR=NULL,POS=NULL,rsID=NULL,Allele1=NULL,Allele2=NULL,GenomicLocus=NULL)

# loop over the genomic risk loci = number of rows of grl:
for(j in c(1:64)){
  chrom = grl[j,CHR]
  gl = grl[j,GenomicLocus]
  I <- grl[j, .(start,end)]
  s <- smcp[CHR==chrom,][ (I[,start] <= POS) & (POS <= I[,end]), ]
  s[,GenomicLocus:=gl]
  loci <- bind_rows(loci,s)
  rm(s)
}

# export the loci table:
write.table(loci, 'F:\\NIKO\\Genomics\\UK_Biobank\\FUMA_outputs\\chronic_pain_EUR\\snps_in_gen_risk_loci\\snps_in_risk_loci.tsv',append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)

# export chromosome indexed set of files with only the rsids of the snps from the loci table, without repetitions:
path <- 'F:\\NIKO\\Genomics\\UK_Biobank\\FUMA_outputs\\chronic_pain_EUR\\snps_in_gen_risk_loci\\risk_loci_snp_POS_chr'

loci <- loci[!duplicated(loci, by=c('CHR', 'POS')), ]

for(i in c(1:20)){
  path1 <- paste(path, i,'.txt',sep='')
  write.table(loci[CHR==i, .(POS)], 
              path1, 
              append = FALSE, sep = "\t", quote = FALSE, col.names=FALSE, row.names=FALSE)
}

