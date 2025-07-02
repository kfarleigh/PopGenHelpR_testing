library(mmod)
library(adegenet)
library(hierfstat)
library(StAMPP)
library(stringr)
library(vcfR)
library(PopGenHelpR)
library(tidyr)


#######################################
##### Farleigh et al. (2021) data #####
#######################################

data("HornedLizard_VCF")
data("HornedLizard_Pop")

# Make horned lizard data into a genind
Genind <- vcfR2genind(HornedLizard_VCF)
Genind@pop <- as.factor(HornedLizard_Pop$Population)
ploidy(Genind) <- 2

# Make horned lizard data into a genlight
Glight <- vcfR2genlight(HornedLizard_VCF)
Glight@pop <- as.factor(HornedLizard_Pop$Population)
ploidy(Glight) <- 2

# Calculate differentiation measures
PGH_dif_HL <- Differentiation(HornedLizard_VCF, pops = HornedLizard_Pop, statistic = "all")

PGH_dif_HL_mat <- as.matrix(PGH_dif_HL$JostsD)


# Jost's D with mmod for the horned lizard data 
mmod <- pairwise_D(Genind)

mmod_mat <- as.matrix(mmod)


## Fst and Nei's D with StAMPP for the horned lizard data 
# Fst
Stmp_fst <- stamppFst(Glight, nboots = 0)

# Population Nei's D
Stmp_popND <- stamppNeisD(Glight)

# Individual Nei's D
Stmp_indND <- stamppNeisD(Glight, pop = FALSE)

#######################################
##### Weins & Collela (2025) data #####
#######################################

wiens_files <- list.files(pattern = "weins*", full.names = T)

wiens_pops <- read.csv(wiens_files[1])
wiens_vcf <- read.vcfR(wiens_files[2]) 

### Differentiation 
Genind2 <- vcfR2genind(wiens_vcf)
Genind2@pop <- as.factor(wiens_pops$Population)
ploidy(Genind2) <- 2


Glight2 <- vcfR2genlight(wiens_vcf)
Glight2@pop <- as.factor(wiens_pops$Population)
ploidy(Glight2) <- 2

# VCF
PGH_wiens_dif_vcf <- Differentiation(data = wiens_files[2], pops = wiens_pops, statistic = "all")

# vcfR object
PGH_wiens_dif_vcfR <- Differentiation(data = wiens_vcf, pops = wiens_pops, statistic = "all")

# geno file 
PGH_wiens_dif_geno <- Differentiation(data = wiens_files[3], pops = wiens_pops, statistic = "all", missing_value = '-1') # Stopped here


### Let's make sure that PGH gives the same estimates no matter what was used as input 
# Fst
table(PGH_wiens_dif_geno$Fst == PGH_wiens_dif_vcf$Fst)
table(PGH_wiens_dif_geno$Fst == PGH_wiens_dif_vcfR$Fst)

# Population Nei's D

table(PGH_wiens_dif_geno$NeisD_pop == PGH_wiens_dif_vcf$NeisD_pop)
table(PGH_wiens_dif_geno$NeisD_pop == PGH_wiens_dif_vcfR$NeisD_pop)

# Individual Nei's D

table(PGH_wiens_dif_geno$NeisD_ind == PGH_wiens_dif_vcf$NeisD_ind)
table(PGH_wiens_dif_geno$NeisD_ind == PGH_wiens_dif_vcfR$NeisD_ind)

# Jost's D

table(PGH_wiens_dif_geno$JostsD == PGH_wiens_dif_vcf$JostsD)
table(PGH_wiens_dif_geno$JostsD == PGH_wiens_dif_vcfR$JostsD)


# Jost's D with mmod 
mmod_wiens_JD <- pairwise_D(Genind2)


## Fst and Nei's D with StAMPP
# Fst
Stmp_fst2 <- stamppFst(Glight2, nboots = 0)

# Population Nei's D
Stmp_popND2 <- stamppNeisD(Glight2)


# Individual Nei's D
Stmp_indND2 <- stamppNeisD(Glight2, pop = FALSE)


#####################################
##### Murphy et al. (2025) data #####
#####################################

murphy_vcf2 <- read.vcfR('Murphy_CF_final.recode.vcf')
murphy_pops <- read.delim("Murphy_CF_popmap_final.txt", header = F)
murphy_pops$V2[murphy_pops$V2 == "5R"] <- "fiveR"

colnames(murphy_pops) <- c("Inds", "population")


Genind4 <- vcfR2genind(murphy_vcf2)
Genind4@pop <- as.factor(murphy_pops$population)
ploidy(Genind4) <- 2

Glight4 <- vcfR2genlight(murphy_vcf2)
Glight4@pop <- as.factor(murphy_pops$population)
ploidy(Glight4) <- 2


# Test differentiation with PGH
PGH_murphy_dif2_JDwNA <- Differentiation(data = murphy_vcf2, pops = murphy_pops, statistic = "all")

PGH_murphy_wNA_JD <- as.matrix(PGH_murphy_dif2_JDwNA$JostsD)

# Jost's D with mmod 
mmod_murphy_JD_wNA <- pairwise_D(Genind4)

mmod_murphy_JD_mat_wNA <- as.matrix(mmod_murphy_JD_wNA)

# Fst and Neis D  
Stmp_fst4 <- stamppFst(Glight4, nboots = 0)
Stmp_popND4 <- stamppNeisD(Glight4)
Stmp_indND4 <- stamppNeisD(Glight4, pop = FALSE)

### Let's see how different the estimates from the different packages are?

# Make the StAMPP and mmod results matrices with the upper triangle NA

# Farleigh et al. (2021) data
mmod_mat <- as.matrix(mmod)
mmod_mat[upper.tri(mmod_mat)] <- NA

Stmp_indND_mat <- Stmp_indND
Stmp_indND_mat[upper.tri(Stmp_indND_mat)] <- NA

Stmp_popND_mat <- Stmp_popND
Stmp_popND_mat[upper.tri(Stmp_popND_mat)] <- NA

# Wiens data 
mmod_wiens_JD_mat <- as.matrix(mmod_wiens_JD)
mmod_wiens_JD_mat[upper.tri(mmod_wiens_JD_mat)] <- NA

Stmp_indND2_mat <- Stmp_indND2
Stmp_indND2_mat[upper.tri(Stmp_indND2_mat)] <- NA

Stmp_popND2_mat <- Stmp_popND2
Stmp_popND2_mat[upper.tri(Stmp_popND2_mat)] <- NA


# Murphy data
mmod_murphy_JD_mat <- as.matrix(mmod_murphy_JD)
mmod_murphy_JD_mat[upper.tri(mmod_murphy_JD_mat)] <- NA
Stmp_indND4_mat <- Stmp_indND4
Stmp_indND4_mat[upper.tri(Stmp_indND4_mat)] <- NA

Stmp_popND4_mat <- Stmp_popND4
Stmp_popND4_mat[upper.tri(Stmp_popND4_mat)] <- NA


## Jost's D
# Farleigh et al. (2021) data
max(abs(PGH_dif_HL$JostsD-mmod_mat), na.rm = T)
mean(PGH_dif_HL$JostsD-mmod_mat, na.rm = T)
median(PGH_dif_HL$JostsD-mmod_mat, na.rm = T)


# Wiens data 
max(abs(PGH_wiens_dif_vcf$JostsD-mmod_wiens_JD_mat), na.rm = T)
mean(PGH_wiens_dif_vcf$JostsD-mmod_wiens_JD_mat, na.rm = T)
median(PGH_wiens_dif_vcf$JostsD-mmod_wiens_JD_mat, na.rm = T)

# Murphy data needs to be looped, mmod outputs populations in different order than PGH 
murph_pop <- unique(murphy_pops$population)

# get the comparisons 
murph_comps <- combn(murph_pop, m = 2)

JD_abs_dif <- c()
JD_dif <- c()

for(i in 1:ncol(murph_comps)){
  
  tmp_comp <- murph_comps[,i]
  
  tmp_pop1 <- tmp_comp[1]
  tmp_pop2 <- tmp_comp[2]
  
  idx_pop1 <- which(colnames(mmod_murphy_JD_mat) == tmp_pop1)
  idx_pop2 <- which(rownames(mmod_murphy_JD_mat) == tmp_pop2)
  
  tmp_mmod <- mmod_murphy_JD_mat[idx_pop2, idx_pop1]
  tmp_mmod <- c(tmp_mmod, mmod_murphy_JD_mat[idx_pop1, idx_pop2])
  
  tmp_mmod <- na.omit(tmp_mmod)
  
  idx_pop1_PGH <- which(colnames(PGH_murphy_JD) == tmp_pop1)
  idx_pop2_PGH <- which(colnames(PGH_murphy_JD) == tmp_pop2)
  
  tmp_PGH <- PGH_murphy_JD[idx_pop2_PGH, idx_pop1_PGH]
  tmp_PGH <- c(tmp_PGH,PGH_murphy_JD[idx_pop1_PGH, idx_pop2_PGH])
  
  tmp_PGH <- na.omit(tmp_PGH)
  
  JD_abs_dif <- c(JD_abs_dif, abs(tmp_PGH-tmp_mmod))
  JD_dif <- c(JD_dif, tmp_PGH-tmp_mmod)
  
  remove(idx_pop1,idx_pop1_PGH,idx_pop2,idx_pop2_PGH,tmp_PGH,tmp_mmod)
}

max(JD_abs_dif)
mean(JD_dif)
median(JD_dif)

## Fst
# Farleigh et al. (2021) data
max(abs(PGH_dif_HL$Fst-Stmp_fst), na.rm = T)
mean(PGH_dif_HL$Fst-Stmp_fst, na.rm = T)
median(PGH_dif_HL$Fst-Stmp_fst, na.rm = T)

# Wiens data 
max(abs(PGH_wiens_dif_vcf$Fst-Stmp_fst2), na.rm = T)
mean(PGH_wiens_dif_vcf$Fst-Stmp_fst2, na.rm = T)
median(PGH_wiens_dif_vcf$Fst-Stmp_fst2, na.rm = T)

# Murphy with NAs
max(abs(PGH_murphy_dif2_JDwNA$Fst-Stmp_fst4), na.rm = T)
mean(PGH_murphy_dif2_JDwNA$Fst-Stmp_fst4, na.rm = T)
median(PGH_murphy_dif2_JDwNA$Fst-Stmp_fst4, na.rm = T)

## Population Nei's D
# Farleigh et al. (2021) data
max(abs(PGH_dif_HL$NeisD_pop-Stmp_popND), na.rm = T)
mean(PGH_dif_HL$NeisD_pop-Stmp_popND, na.rm = T)
median(PGH_dif_HL$NeisD_pop-Stmp_popND, na.rm = T)

# Wiens data 
max(abs(PGH_wiens_dif_vcf$NeisD_pop-Stmp_popND2_mat), na.rm = T)
mean(PGH_wiens_dif_vcf$NeisD_pop-Stmp_popND2_mat, na.rm = T)
median(PGH_wiens_dif_vcf$NeisD_pop-Stmp_popND2_mat, na.rm = T)

# Murphy with NA
max(abs(PGH_murphy_dif2_JDwNA$NeisD_pop-Stmp_popND4_mat), na.rm = T)
mean(PGH_murphy_dif2_JDwNA$NeisD_pop-Stmp_popND4_mat, na.rm = T)
median(PGH_murphy_dif2_JDwNA$NeisD_pop-Stmp_popND4_mat, na.rm = T)

## Individual Nei's D
# Farleigh et al. (2021) data
max(abs(PGH_dif_HL$NeisD_ind-Stmp_indND), na.rm = T)
mean(PGH_dif_HL$NeisD_ind-Stmp_indND, na.rm = T)
median(PGH_dif_HL$NeisD_ind-Stmp_indND, na.rm = T)

# Wiens data 
max(abs(PGH_wiens_dif_vcf$NeisD_ind-Stmp_indND2_mat), na.rm = T)
mean(PGH_wiens_dif_vcf$NeisD_ind-Stmp_indND2_mat, na.rm = T)
median(PGH_wiens_dif_vcf$NeisD_ind-Stmp_indND2_mat, na.rm = T)

# Murphy with NA
max(abs(PGH_murphy_dif2_JDwNA$NeisD_ind-Stmp_indND4_mat), na.rm = T)
mean(PGH_murphy_dif2_JDwNA$NeisD_ind-Stmp_indND4_mat, na.rm = T)
median(PGH_murphy_dif2_JDwNA$NeisD_ind-Stmp_indND4_mat, na.rm = T)


##########################
##### Heterozygosity #####
##########################

PGH_het_HL <- Heterozygosity(data = HornedLizard_VCF, pops = HornedLizard_Pop)

PGH_het_HL_wNA <- Heterozygosity(data = vcf_wNA, pops = HornedLizard_Pop)

PGH_het_wiens_vcf <- Heterozygosity(data = wiens_files[2], pops = wiens_pops)

PGH_het_wiens_vcfR <- Heterozygosity(data = wiens_vcf, pops = wiens_pops)

PGH_het_wiens_geno <- Heterozygosity(data = wiens_files[3], pops = wiens_pops)

PGH_het_murphy <- Heterozygosity(data = murphy_vcf, pops = murphy_pops)

PGH_het_murphy_wNA <- Heterozygosity(data = murphy_vcf2, pops = murphy_pops)

# Check to make sure the heterozygosity estimates match regardless of file type

table(PGH_het_wiens_geno$Ho_perpop == PGH_het_wiens_vcf$Ho_perpop)
table(PGH_het_wiens_geno$Ho_perloc == PGH_het_wiens_vcf$Ho_perloc)
table(PGH_het_wiens_geno$He_perpop == PGH_het_wiens_vcf$He_perpop)
table(PGH_het_wiens_geno$He_perloc == PGH_het_wiens_vcf$He_perloc)
table(PGH_het_wiens_geno$PHt == PGH_het_wiens_vcf$PHt)
table(PGH_het_wiens_geno$Hs_exp  == PGH_het_wiens_vcf$Hs_exp)
table(PGH_het_wiens_geno$Hs_obs  == PGH_het_wiens_vcf$Hs_obs)
table(PGH_het_wiens_geno$IR  == PGH_het_wiens_vcf$IR)
table(PGH_het_wiens_geno$HL  == PGH_het_wiens_vcf$HL)

table(PGH_het_wiens_geno$Ho_perpop == PGH_het_wiens_vcfR$Ho_perpop)
table(PGH_het_wiens_geno$Ho_perloc == PGH_het_wiens_vcfR$Ho_perloc)
table(PGH_het_wiens_geno$He_perpop == PGH_het_wiens_vcfR$He_perpop)
table(PGH_het_wiens_geno$He_perloc == PGH_het_wiens_vcfR$He_perloc)
table(PGH_het_wiens_geno$PHt == PGH_het_wiens_vcfR$PHt)
table(PGH_het_wiens_geno$Hs_exp  == PGH_het_wiens_vcfR$Hs_exp)
table(PGH_het_wiens_geno$Hs_obs  == PGH_het_wiens_vcfR$Hs_obs)
table(PGH_het_wiens_geno$IR  == PGH_het_wiens_vcfR$IR)
table(PGH_het_wiens_geno$HL  == PGH_het_wiens_vcfR$HL)

# Check to make sure the heterozygosity estimates match regardless of missing data
table(PGH_het_HL$Ho_perpop == PGH_het_HL_wNA$Ho_perpop)
table(PGH_het_HL$Ho_perloc == PGH_het_HL_wNA$Ho_perloc)
table(PGH_het_HL$He_perpop == PGH_het_HL_wNA$He_perpop)
table(PGH_het_HL$He_perloc == PGH_het_HL_wNA$He_perloc)
table(PGH_het_HL$PHt == PGH_het_HL_wNA$PHt)
table(PGH_het_HL$Hs_exp  == PGH_het_HL_wNA$Hs_exp)
table(PGH_het_HL$Hs_obs  == PGH_het_HL_wNA$Hs_obs)
table(PGH_het_HL$IR  == PGH_het_HL_wNA$IR)
table(PGH_het_HL$HL  == PGH_het_HL_wNA$HL)

table(PGH_het_murphy$Ho_perpop == PGH_het_murphy_wNA$Ho_perpop)
table(PGH_het_murphy$Ho_perloc == PGH_het_murphy_wNA$Ho_perloc)
table(PGH_het_murphy$He_perpop == PGH_het_murphy_wNA$He_perpop)
table(PGH_het_murphy$He_perloc == PGH_het_murphy_wNA$He_perloc)
table(PGH_het_murphy$PHt == PGH_het_murphy_wNA$PHt)
table(PGH_het_murphy$Hs_exp  == PGH_het_murphy_wNA$Hs_exp)
table(PGH_het_murphy$Hs_obs  == PGH_het_murphy_wNA$Hs_obs)
table(PGH_het_murphy$IR  == PGH_het_murphy_wNA$IR)
table(PGH_het_murphy$HL  == PGH_het_murphy_wNA$HL)

### Check to make sure they match with other programs

## Observed heterozygosity
# Horned lizard
Hstat_HL <- genind2hierfstat(Genind)
Hstat_HL_wNA <- genind2hierfstat(Genind_wNA)

Hstat_hets_HL <- basic.stats(Hstat_HL)
Hstat_Ho_HL <- colMeans(Hstat_hets_HL$Ho)

Hstat_hets_HL_wNA <- basic.stats(Hstat_HL_wNA)
Hstat_Ho_HL_wNA <- colMeans(Hstat_hets_HL_wNA$Ho)

# Weins
Hstat_weins <- genind2hierfstat(Genind2)

Hstat_het_weins <- basic.stats(Hstat_weins)
Hstat_Ho_weins <- colMeans(Hstat_het_weins$Ho)

# Murphy 
Hstat_murphy <- genind2hierfstat(Genind3)
Hstat_murphy_wNA <- genind2hierfstat(Genind4)

Hstat_hets_murphy <- basic.stats(Hstat_murphy)
Hstat_Ho_murphy <- colMeans(Hstat_hets_murphy$Ho, na.rm = T)

Hstat_hets_murphy_wNA <- basic.stats(Hstat_murphy_wNA)
Hstat_Ho_murphy_wNA <- colMeans(Hstat_hets_murphy_wNA$Ho, na.rm = T)

# Check if the horned lizard ones match
PGH_het_HL$Ho_perpop[,1]-Hstat_Ho_HL
PGH_het_HL_wNA$Ho_perpop[,1]-Hstat_Ho_HL_wNA

# Check if weins matched
PGH_het_wiens_vcf$Ho_perpop[,1]-Hstat_Ho_weins

# Check if murphy data matched, we do this manually because they are ordered differently

PGH_het_murphy$Ho_perpop[order(match(PGH_het_murphy$Ho_perpop$Pop, names(Hstat_Ho_murphy))),1] - Hstat_Ho_murphy
PGH_het_murphy_wNA$Ho_perpop[order(match(PGH_het_murphy_wNA$Ho_perpop$Pop, names(Hstat_Ho_murphy_wNA))),1] - Hstat_Ho_murphy_wNA


## Expected heterozygosity (He)
# Horned lizard
He_HL <- Hs(Genind)
He_HL_wNA <- Hs(Genind_wNA)

# Weins
He_weins <- Hs(Genind2)

# Murphy 
He_murphy <- Hs(Genind3)
He_murphy_wNA <- Hs(Genind4)

# Check if the horned lizard ones match
PGH_het_HL$He_perpop$Expected.Heterozygosity-He_HL
PGH_het_HL_wNA$He_perpop$Expected.Heterozygosity-He_HL_wNA

# Check if the weins ones match
PGH_het_wiens_vcf$He_perpop$Expected.Heterozygosity-He_weins

# Check if the murphy ones match

PGH_het_murphy$He_perpop[order(match(PGH_het_murphy$He_perpop$Pop, names(He_murphy))),1]-He_murphy
PGH_het_murphy_wNA$He_perpop[order(match(PGH_het_murphy_wNA$He_perpop$Pop, names(He_murphy_wNA))),1]-He_murphy_wNA

#############################################
### Test statistics originally from Rhh #####
#############################################

source("./DC_hetero.R")

HL_rhh <- as.data.frame(t(extract.gt(HornedLizard_VCF, return.alleles = T)))
HL_rhh <- HL_rhh %>% 
  separate_wider_delim(everything(), delim = "/", names_sep = ".")

HL_rhh_final <- cbind(HornedLizard_Pop[1:2], HL_rhh)

HL_ir <- ir(HL_rhh_final)
HL_hl <- hl(HL_rhh_final)
HL_mlh <- mlh(HL_rhh_final)
HL_oh <- oh(HL_rhh_final)
HL_sh <- sh(HL_rhh_final)

# Use the functions from Daren Card's script to test IR and HL
max(PGH_het_HL$IR[,1]-HL_ir, na.rm = T)
mean(PGH_het_HL$IR[,1]-HL_ir, na.rm = T)

max(PGH_het_HL$HL[,1]-HL_hl, na.rm = T)
mean(PGH_het_HL$HL[,1]-HL_hl, na.rm = T)

# Try it with missing data 
HL_rhh_wNA <- as.data.frame(t(extract.gt(vcf_wNA, return.alleles = T)))
HL_rhh_wNA <- HL_rhh_wNA %>% 
  separate_wider_delim(everything(), delim = "/", names_sep = ".")

HL_rhh_wNA_final <- cbind(HornedLizard_Pop[1:2], HL_rhh_wNA)

HL_wNA_ir <- ir(HL_rhh_wNA_final)
HL_wNA_hl <- hl(HL_rhh_wNA_final)

max(PGH_het_HL$IR[,1]-HL_ir, na.rm = T)
mean(PGH_het_HL$IR[,1]-HL_ir, na.rm = T)

max(PGH_het_HL_wNA$HL[,1]-HL_wNA_hl, na.rm = T)
mean(PGH_het_HL_wNA$HL[,1]-HL_wNA_hl, na.rm = T)

max(PGH_het_HL_wNA$IR[,1]-HL_wNA_ir, na.rm = T)
mean(PGH_het_HL_wNA$IR[,1]-HL_wNA_ir, na.rm = T)
