# R coloc_ GGE vs Short Sleeps
# R 4.4.0
# 20240618

# Library
library(coloc)
library(dplyr)
library(tidyverse)
library(ieugwasr) #To get LD matrix for Susie
library(vroom) # Faster than freader and read.table
library(genetics.binaRies)
# If genetics.binaRies is not available (not available in CRAN), download from:
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("MRCIEU/genetics.binaRies")



# Data Path
GGE.path <-  " " # Directory to trait 1 sumstat
short.sleep.path <-  " " # directory to trait 2 sumstat
export.path <-  " "
coloc.result.path <-  " "


# Part A: Cleaning GWAS sumstats file
# a/ Selecting only necessary columns for Coloc analysis (refer to tutorial Github).
gge.gwas.sumstats <-  vroom(file = file.path(GGE.path, "reformatted_ILAE3_Caucasian_GGE_final.txt"), 
                         delim = "\t", show_col_types = F)

gge.gwas.sumstats$Var <-  gge.gwas.sumstats$SE^2 # Calculate Var from Standard Errors

gge.trait <-  gge.gwas.sumstats %>% 
  select(CHR, BP, Allele1, Allele2, Effective_N, MarkerName, Beta, Var) 

# b/ Investigating and Cleaning missing data (NA) from the dataset
gge.trait %>% is.na() %>% table() # 5291 SNPs (only on Chr23) has NA beta value => cannot fix => Filter out
gge.trait <-  gge.trait %>% na.omit()
gge.trait %>% is.na() %>% table() # No NA left

# c/ Reformatting Short-sleep GWAS sumstats to the same column name as GGE sumstat (to ease the use of extraction function later)

shortsleep.gwas.sumstat  <-  vroom(file = file.path(short.sleep.path, "shortsumstats.txt"), 
                                delim = "\t", show_col_types = F)

shortsleep.gwas.sumstat$Var <-  (shortsleep.gwas.sumstat$SE_SHORTSLEEP)^2

sleep.trait <-  shortsleep.gwas.sumstat %>% select(CHR, BP, SNP, ALLELE1, ALLELE0, BETA_SHORTSLEEP, Var)

# Changing column name to match with GGE sumstat file
colnames(sleep.trait)[which(colnames(sleep.trait) %in% c("SNP","ALLELE1", "ALLELE0","BETA_SHORTSLEEP"))] = c("MarkerName","Allele1", "Allele2", "Beta")

sleep.trait %>% is.na() %>% table() #No missing data



# Part B: Extracting scanning regions from GGE traits.

### Step 1: Identify Lead SNPs for GGE, using leadSNPs file from FUMA 
gge.lead.snps= vroom(file = " ", # path to FUMA_leadsnps.txt file 
                     delim ="\t", show_col_types = F)

table(gge.lead.snps$p <5e-8) # All lead SNPs passed GWAS significant threshold as expected from FUMA setting


### Step 2: Define regions around lead SNPs of GGE, using 250,000 kb window
region_size <-  250000 # Can change the scanning window

# Function take lead SNPs and return the a table (leadSNPs; chr; start and stop of the region)  
get_region <-  function(snp_row) {
  start <-  snp_row$pos - region_size
  end <-  snp_row$pos + region_size
  return(data.frame(snp= snp_row$rsID, chr= snp_row$chr, start= start, end=end))
}

# Define the scanning regions around GGE lead SNPs
regions <-  do.call(rbind, lapply(1:nrow(gge.lead.snps), function(i) get_region(gge.lead.snps[i, ])))


### Step 3: Extract all SNPs within each defined region 

# Create a SNPs extraction function
extract_snps_in_region <-  function(region, gwas_data) {
  snps_in_region= gwas_data[gwas_data$CHR == region$chr & 
                              gwas_data$BP >= region$start &
                              gwas_data$BP <= region$end, ]
  snps_in_region$region <- paste(region$chr, region$start, region$end, sep = ":")
  return(snps_in_region)
}

# Run SNPs extraction function for GGE
gge.snps_in_region <-  do.call(rbind, lapply(1:nrow(regions), function(i) extract_snps_in_region(regions[i, ], gge.trait)))

# Run SNPs extraction function for Short-sleep
sleep.snps_in_region <-  do.call(rbind, lapply(1:nrow(regions), function(i) extract_snps_in_region(regions[i, ], sleep.trait)))


# Part C: Running Coloc Analysis

# Pre-analysis: Output the scanning regions before running coloc 
unique(sleep.snps_in_region$region)

# Coloc Analysis: Plug and Play
# Paste the region name (see above) to input_region, and run the script.
input_region=  " " # For example: 2:57958035:58458035

# Checking for any duplicated SNPs as it can cause issue with coloc analysis

# #For GGE
gge <-   gge.snps_in_region[which(gge.snps_in_region$region == input_region), ]
gge <-  gge[!duplicated(gge$MarkerName), ] #Remove any duplicate

# For Short Sleep
shortsleep <-  sleep.snps_in_region[which(sleep.snps_in_region$region == input_region), ] 
shortsleep <-  shortsleep[!duplicated(shortsleep$MarkerName), ] #Remove any duplicate

# Filter the same SNPs in the region (Between GGE and Short sleep) (for Coloc-Susie later)
share_snps <-  intersect(gge$MarkerName, shortsleep$MarkerName)
gge= gge %>% filter(MarkerName %in% share_snps)
shortsleep <-  shortsleep %>% filter(MarkerName %in% share_snps)


dataset1 <- list(
  beta = gge$Beta,
  varbeta= gge$Var,
  snp = gge$MarkerName,
  position= gge$BP,
  type = "cc"
) # For GGE (trait 1)


dataset2 <- list(
  beta = shortsleep$Beta,
  varbeta= shortsleep$Var,
  snp = shortsleep$MarkerName,
  position= shortsleep$BP,
  type = "cc"
) # For Short-sleep (trait 2)

# Coloc sensitivity test
check_dataset(dataset1)
check_dataset(dataset1) 

# Plotting the coloc sensitivity tests
par(mfrow= c(2,1))
plot_dataset(dataset1)
plot_dataset(dataset2)

# Coloc results
my.res= coloc.abf(dataset1, dataset2)
my.res$summary

sensitivity(my.res, "H4 > 0.80")

# Extract the ranking of SNPs based on the probability 
subset(my.res$results, SNP.PP.H4 > 0.8)

#Need to change the export file path 
vroom_write(x = my.res$summary, 
            file = file.path(coloc.result.path, paste0(input_region, "_250kb_result_GGE_to_ShortSleep.txt")), 
            delim = "\t")


# Part D: Coloc-Susie for overcome the assumption of normal Coloc

# a/ Creating LD matrix for Susie: 2:57738194:58238194

# Step 1: Extract SNPS in 2:57738194:58238194 for both GGE and ShortSleep

# Filter to include only these common snps
filtered_gge <-  gge[gge$MarkerName %in% share_snps, ]
filtered_shortslep <-  shortsleep[shortsleep$MarkerName %in% share_snps, ]

# Align the SNPs in order for both traits, else, it will create error for LD matrix
filtered_gge <-  filtered_gge[match(share_snps, filtered_gge$MarkerName), ]
filtered_shortslep <-  filtered_shortslep[match(share_snps, filtered_shortslep$MarkerName), ]

# Check if the rsID in the same order between 2 traits
identical(filtered_shortslep$MarkerName, filtered_shortslep$MarkerName)

# Use ld_matrix_local of ieugwasr to generate LD matrix, given the 1000G_phase3 reference bim file
plink_path <-  " " #PLINK path

# Creating matrix
ld.matrix <-  ieugwasr::ld_matrix(
  variants =  share_snps,
  with_alleles = FALSE, #Set this to FALSE to keep the original rsID
  pop = "EUR",
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = " " # EUR bim file
) 

length(colnames(ld.matrix)) # Not all share_snps in LD.matrix => Need filter again
length(share_snps)  


# 2nd round of filter
snp_in_LD.matrix <-  colnames(ld.matrix)
share_snps <- intersect(snp_in_LD.matrix, share_snps)
# share_snps %>% length() 

# Then redo the matrix
ld.matrix <-  ieugwasr::ld_matrix(
  variants =  share_snps, 
  with_alleles = FALSE, #Set this to FALSE to keep the original rsID
  pop = "EUR",
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = " " # EUR bim file
) 


# At this step, I cannot perform Susie analysis because of the issue with Allele flipping 
# => Need to keep the allele side according to the reference panel used to generated LD matrix 

# b/ Flipping Allele (Extremely crucial steps for Susie as need to pass the sensitivity test to run the analysis)

# Also need to re-filter the snps in both dataset
# Checking alignment between LD matrix and the two traits
filtered_gge <-  filtered_gge[match(share_snps, filtered_gge$MarkerName), ]
filtered_shortslep <-  shortsleep[shortsleep$MarkerName %in% share_snps, ]

identical(filtered_gge$MarkerName, filtered_shortslep$MarkerName) #Matching order


# Import EUR bim file (for alignment)
EUR.bim <-  vroom(" ", delim = "\t", show_col_types = F) # EUR bim file
colnames(EUR.bim) <-  c("CHR", "MarkerName", "Genetic Distance" ,"BP", "Allele2", "Allele1") #Allele_2 = Reference Allele of bim file
EUR.bim = EUR.bim %>% select(CHR, MarkerName, BP, Allele1, Allele2) # Select only important columns


# Flip_1: Flipping 1 trait accroding to EUR.bim (Using ShortSleep)
#Using ShortSleep to flip according to EUR.bim
filtered_EUR <-  EUR.bim %>% filter(MarkerName %in% share_snps)
merge_1st_df <-  merge(filtered_shortslep, filtered_EUR, by= "MarkerName", suffixes = c("_ShortSleep", "_EUR"))


merge_1st_df <-  merge_1st_df %>% 
  mutate(
    Allele1_ShortSleep = tolower(Allele1_ShortSleep),
    Allele2_ShortSleep = tolower(Allele2_ShortSleep),
    Allele1_EUR = tolower(Allele1_EUR),
    Allele2_EUR = tolower(Allele2_EUR)
  )

flip_function_1st <-  function(data_frame) {
  flip_require = FALSE
  if(data_frame$Allele1_ShortSleep == data_frame$Allele2_EUR &&
     data_frame$Allele2_ShortSleep == data_frame$Allele1_EUR) {
    
    #Swap alleles
    temp_allele1 = data_frame$Allele1_ShortSleep
    data_frame$Allele1_ShortSleep = data_frame$Allele2_ShortSleep
    data_frame$Allele2_ShortSleep = temp_allele1
    data_frame$BP_ShortSleep = -data_frame$BP_ShortSleep # Flip the beta sign
    flip_require = TRUE #If flipping is required --> TRUE
  }
  data_frame$flip_require = flip_require
  return(data_frame)  
}

aligned_data_1st = merge_1st_df %>%
  rowwise() %>%
  mutate(data = list(flip_function_1st(cur_data()))) %>%
  unnest(cols= c(data), names_repair = "unique")

aligned_data_1st$flip_require %>% table() # No need to flip Short Sleep


# Since Short Sleep is aligned perfectly witn the EUR bim file=> I used shortsleep as reference to flip GGE

# Flip_2: Flipping trait 2 according to trait 1
merge_2nd_df <-  merge(filtered_shortslep, filtered_gge, by= "MarkerName", suffixes = c("_ShortSleep", "_GGE"))

# Turn all the alleles to the same cases (using lower cases)
merge_2nd_df <-  merge_2nd_df %>% 
  mutate(
    Allele1_ShortSleep = tolower(Allele1_ShortSleep),
    Allele2_ShortSleep = tolower(Allele2_ShortSleep),
    Allele1_GGE = tolower(Allele1_GGE),
    Allele2_GGE = tolower(Allele2_GGE)
  )


flip_function_2nd <-  function(data_frame) {
  flip_require = FALSE
  if(data_frame$Allele1_ShortSleep == data_frame$Allele2_GGE &&
     data_frame$Allele2_ShortSleep == data_frame$Allele1_GGE) {
    
    #Swap alleles
    temp_allele1 = data_frame$Allele1_GGE
    data_frame$Allele1_GGE = data_frame$Allele2_GGE
    data_frame$Allele2_GGE = temp_allele1
    data_frame$Beta_GGE = -data_frame$Beta_GGE # Flip the beta sign
    flip_require = TRUE #If flipping is required --> TRUE
  }
  data_frame$flip_require = flip_require
  return(data_frame)  
}

aligned_data_2nd <-  merge_2nd_df %>%
  rowwise() %>%
  mutate(data = list(flip_function_2nd(cur_data()))) %>%
  unnest(cols= c(data), names_repair = "unique")

# aligned_data$flip_require %>% table() # 
aligned_data_2nd$flip_require %>% table() 


# Now merge the data back to the traits
# Arrange the dataset based on marker names
filtered_gge <-  filtered_gge %>% arrange(MarkerName)
filtered_shortslep  <-  filtered_shortslep %>% arrange(MarkerName)
aligned_data_2nd <-  aligned_data_2nd %>% arrange(MarkerName...1)

# Filter out the dataset for Susie (beta, var(beta), snp, position)
susie_GGE  <-  aligned_data_2nd %>% 
  select(MarkerName...1, Beta_GGE...30, Var_GGE...31, BP_GGE...26) # Careful with the chosen column

susie_ShortSleep <-  aligned_data_2nd %>% 
  select(MarkerName...1, Beta_ShortSleep...22, Var_ShortSleep...23, BP_ShortSleep...19)


# Check for ordering of SNPs between Sleep vs GEE vs the LD matrix
identical(susie_GGE$MarkerName...1, susie_ShortSleep$MarkerName...1 ) # Marker name of 2 traits in same order
identical(susie_GGE$MarkerName...1, ld.matrix %>% colnames()) # Not in same order with LD matrix => Arrange LD matrix

# Arrange the LD matrix
ordered_ld.matrix <-  ld.matrix[susie_GGE$MarkerName...1, susie_GGE$MarkerName...1]

# Check for order of SNPs 
identical(susie_GGE$MarkerName...1, susie_ShortSleep$MarkerName...1 ) # Marker name of 2 traits in same order
identical(susie_GGE$MarkerName...1, ordered_ld.matrix %>% colnames()) # Marker name in some order as the LD matrix


# c/ Prepare data for Susie
coloc_dataset1_GGE <- list(
  beta = susie_GGE$Beta_GGE...30,
  varbeta= susie_GGE$Var_GGE...31,
  snp = susie_GGE$MarkerName...1,
  position= susie_GGE$BP_GGE...26,
  type = "cc",
  LD= ordered_ld.matrix,
  N=as.integer("49388")
)

coloc_dataset2_ShortSleep <- list(
  beta = susie_ShortSleep$Beta_ShortSleep...22,
  varbeta= susie_ShortSleep$Var_ShortSleep...23,
  snp = susie_ShortSleep$MarkerName...1,
  position= susie_ShortSleep$BP_ShortSleep...19,
  type = "cc",
  LD= ordered_ld.matrix,
  N=as.integer("446118") 
)

# Checking the dataset
check_dataset(coloc_dataset1_GGE, req="LD")
check_dataset(coloc_dataset2_ShortSleep, req="LD") 

# par(mar= c(4,4,1,2))
check_alignment(coloc_dataset1_GGE) # 79.3% alignment for GGE
text(-100, 40000, "GGE", col="blue", cex= 2)

check_alignment(coloc_dataset2_ShortSleep) # 78.6% alignment for ShortSleep
text(-50, 40000, "ShortSleep", col="blue", cex= 2)



# Running SuSie
# Testing if Susie run
S1 <-  runsusie(coloc_dataset1_GGE) # converged: TRUE
summary(S1)

S2 <-  runsusie(coloc_dataset2_ShortSleep) # converged: TRUE
summary(S2)

susie.res <-  coloc.susie(S1,S2)
susie.res$summary

#Need to change the export file path 
vroom_write(x= susie.res$summary, 
            file= file.path(coloc.result.path, paste0(input_region, "_Susie_result_GGE_to_ShortSleep.txt")), 
            delim = "\t")


sensitivity(susie.res, "H4 > 0.7", row=1, dataset1 = coloc_dataset1_GGE, dataset2 = coloc_dataset2_ShortSleep)
sensitivity(susie.res, "H4 > 0.7", row=2, dataset1 = coloc_dataset1_GGE, dataset2 = coloc_dataset2_ShortSleep)

