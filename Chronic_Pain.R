library(data.table)
library(dplyr)

############# YOUR INPUT GOES HERE: ############################################################

path_to_pain_duration <- "E:\\Genomics\\Resources\\Output_data_ph2\\pain_duration.txt" 
 
path_to_bgen_sample <- "E:\\Genomics\\Resources\\UKB_sample_file_newest\\ukb_bgen_sample_new.txt"

path_to_eur_sample <- "E:\\Genomics\\Resources\\FID_IID_EUR_bgen_sample\\IDs_UKB_EUR_3PCs_0920.txt"

path_to_covariates <- "E:\\Genomics\\Resources\\Output_data_ph2\\demogr_geno.txt"

path_to_output_table <- "E:\\Genomics\\Resources\\Chronic_Pain\\Pheno_bgen_sample_eur\\Pheno_chrn_pain_bgen_eur.txt" 

path_to_output_table1 <- "E:\\Genomics\\Resources\\Chronic_Pain\\Pheno_bgen_sample_eur_women\\Pheno_chrn_pain_bgen_eur_women.txt" 


############# NOW RUN THE WHOLE SCRIPT: ########################################################

pain <- fread(path_to_pain_duration) # full table with all genotyped and non-genotyped samples
sample_bgen <- fread(path_to_bgen_sample) # column of sample with imputed geontypes, in bgen files
eur <- fread(path_to_eur_sample) # column of sample with imputed genotypes of european origin
cov <- fread(path_to_covariates) # covariates
names(pain)[1] <- 'IID'
names(sample_bgen) <- 'IID'
names(eur) <- 'IID'
names(cov)[1] <- 'IID'

sample_bgen <- sample_bgen[order(sample_bgen[,IID]),]
sample_bgen <- sample_bgen[c(15:nrow(sample_bgen)),]


############# The function calc_chrn_pain calculates the number of chronic pain sites:

# This is an aolder version of the new version calc_chrn_pain_updt()
# same result different implementation: in uses data frames dplyr and the other more data tables
calc_chrn_pain <- function(pain_in){
  # select only columns that matter:
  pain_ <- select(pain_in, IID, Prf_no_Ansr_v0, Non_Abve_v0 , ends_with("m_v0"))
  # rename them, to avoid clashing symbols, i.e. "+" in the headers' names
  names(pain_) <- c("IID", "Prf_no_Ansr_v0", "Non_Abve_v0", "All_over_3m_v0", 
                    "Neck_Shldr_pn_3m_v0", "Hip_pn_3m_v0",        
                    "Back_pn_3m_v0", "Stom_Abdmn_pn_3m_v0",
                    "Knee_pn_3m_v0", "Headch_3m_v0", "Face_pn_3m_v0")
  # extract only the people who gave an answer to the pain questions: 
  pain_ <- filter(pain_, Prf_no_Ansr_v0 == 0 )
  # convert all NA into 0
  pain_[is.na(pain_)] <- 0
  # the function that calcualtes the number of chronic pain sites: 
  num_sites <- function(row){
    r <- row[indx]
    if(r[1] == 1){ return(8) }
    else if(max(r) == 1){ return(sum(r == 1)) }
    else if(max(r) == 0 & min(r) == 0){ return(0) }
    else { return(NA) }
  }
  # starting the construction of the chronic pain phenotype:
  # take only the IID from the pain_ table
  chrn_pain <- select(pain_, IID)
  # add the column of count of chronic pain siates:
  chrn_pain <- mutate(chrn_pain, count_pain_sites_3m_v0 = NA)
  # indx are the indexes of the headers that refer to 3 month pain:
  indx <- which(grepl("_3m_", names(pain_)))
  # calculate the number of chronic pain sites for each sample:
  chrn_pain[, 'count_pain_sites_3m_v0'] <- as.data.table(apply(pain_, 1, num_sites))
  # Add an FID column
  pain_ <- select(chrn_pain, IID)
  names(pain_) <- 'FID'
  chrn_pain <- bind_cols(pain_, chrn_pain)
  #names(chrn_pain)[c(1,2)] <- c('FID','IID')
  return(chrn_pain)
} 


calc_chrn_pain_updt <- function(pain_in){
  # select only columns that matter:
  # pain_ <- select(pain_in, IID, Prf_no_Ansr_v0, Non_Abve_v0 , ends_with("m_v0"))
  col <- c('IID', 'Prf_no_Ansr_v0', 'Non_Abve_v0', names(pain_in)[which(grepl("m_v0", names(pain_in)))])
  pain_in <- pain_in[,..col]
  # rename them, to avoid clashing symbols, i.e. "+" in the headers' names
  names(pain_in) <- c("IID", "Prf_no_Ansr_v0", "Non_Abve_v0", "All_over_3m_v0", 
                    "Neck_Shldr_pn_3m_v0", "Hip_pn_3m_v0",        
                    "Back_pn_3m_v0", "Stom_Abdmn_pn_3m_v0",
                    "Knee_pn_3m_v0", "Headch_3m_v0", "Face_pn_3m_v0")
  # extract only the people who gave an answer to the pain questions: 
  pain_in <- pain_in[Prf_no_Ansr_v0 == 0,]
  # convert all NA into 0
  pain_in[is.na(pain_in)] <- 0
  # the function that calcualtes the number of chronic pain sites: 
  num_sites <- function(row){
    r <- row[indx]
    if(r[1] == 1){ return(8) }
    else if(max(r) == 1){ return(sum(r == 1)) }
    else if(max(r) == 0 & min(r) == 0){ return(0) }
    else { return(NA) }
  }
  # starting the construction of the chronic pain phenotype:
  # take only the IID from the pain_ table
  chrn_pain <-pain_in[, .(IID)]
  # add the column of count of chronic pain siates:
  indx <- which(grepl("_3m_", names(pain_in)))
  chrn_pain <- data.table(chrn_pain, count_pain_sites_3m_v0=apply(pain_in, 1, num_sites))
  # indx are the indexes of the headers that refer to 3 month pain:
  
  # calculate the number of chronic pain sites for each sample:
  # chrn_pain[, 'count_pain_sites_3m_v0'] <- as.data.table(apply(pain_, 1, num_sites))
  # Add an FID column
  pain_in <- select(chrn_pain, IID)
  names(pain_in) <- 'FID'
  chrn_pain <- bind_cols(pain_in, chrn_pain)
  #names(chrn_pain)[c(1,2)] <- c('FID','IID')
  return(chrn_pain)
} 

#############################################################################################


########################## Master table, complete UKB sample chronic pain ###################

pheno_pain <- calc_chrn_pain_updt(pain)
pheno_pain_old <- calc_chrn_pain(pain)
pheno_pain <- pheno_pain[!is.na(count_pain_sites_3m_v0), ]
cov <- cov[,c(1,2,6,9,10,c(12:51))]
pheno_pain <- left_join(pheno_pain, cov, by=c('IID'))

#############################################################################################

########################## Table, full bgen sample chronic pain ######################

pain_sample_bgen <- inner_join(pheno_pain, sample_bgen, by='IID')

#############################################################################################

##########################  Table,  EUR bgen sample chronic pain ######################

pain_EUR_bgen <- inner_join(pheno_pain, eur, by='IID')

#############################################################################################

########################## Table, EUR bgen sample women chronic pain ######################
# pain_EUR_bgen_women <- pain_EUR_bgen[Genetic_sex == 0, ]
pain_EUR_bgen_women <- pain_EUR_bgen[Sex == 0, ]

# 0 women, 1 men
#############################################################################################

write.table(pain_EUR_bgen, 
            path_to_output_table,
            append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)

write.table(pain_EUR_bgen_women, 
            path_to_output_table1,
            append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)
