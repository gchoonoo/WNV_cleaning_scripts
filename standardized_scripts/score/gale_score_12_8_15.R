##########################################################

# clinical score 12/8/15

##########################################################

# setwd("/Users/choonoo/WNV")

# save.image("~/gale_cs_dec8.RData")

# load("~/gale_cs_dec8.RData")

library(gdata)

##########################################################

# read in data
read.xls(xls="./Gale_12_8_15/Gale_Lab_scores_12_8_15.xlsx",sheet=1) -> gale_lund_weight_cs_pheno

# 1. Clean mating: "X" changed to "x"
gale_lund_weight_cs_pheno[,"Mating"] <- gsub("X","x",gale_lund_weight_cs_pheno[,"Mating"])

# check mating for line 18 and 39

# if true, then line 18 has correct mating
sum(gale_lund_weight_cs_pheno[which(gale_lund_weight_cs_pheno[,"UW_Line"] == 18),"Mating"] == "8042x16513") == length(gale_lund_weight_cs_pheno[which(gale_lund_weight_cs_pheno[,"UW_Line"] == 18),"Mating"])

# set warning if line 39 has wrong mating
if(sum(gale_lund_weight_cs_pheno[which(gale_lund_weight_cs_pheno[,"UW_Line"] == 39),"Mating"]  == "16211x16557") != length(gale_lund_weight_cs_pheno[which(gale_lund_weight_cs_pheno[,"UW_Line"] == 39),"Mating"])){
  print("Warning: Check the mating on UW Line 39")
}

# set warning if "8026x8050" mating exists
if(length(gale_lund_weight_cs_pheno[which(gale_lund_weight_cs_pheno[,"Mating"] == "8026x8050"),"Mating"]) != 0){
  print("Warning: Mating '8026x8050' is actually '8026x5080'")
}

# check if any duplicated rows
sum(duplicated(gale_lund_weight_cs_pheno))

# interger 0 if all matings have unique uw line
unique(gale_lund_weight_cs_pheno[,"UW_Line"])[sapply(1:length(unique(gale_lund_weight_cs_pheno[,"UW_Line"])), function(xx)length(unique(gale_lund_weight_cs_pheno[as.vector(unlist(sapply(unique(gale_lund_weight_cs_pheno[,"UW_Line"])[xx],function(x)which(x==gale_lund_weight_cs_pheno[,"UW_Line"])))),"Mating"]))) > 1]

# Specific data correction #2, clean line 82
gale_lund_weight_cs_pheno[which(gale_lund_weight_cs_pheno[,"UW_Line"] == 82),"Mating"] <- "16557x3154"

# 2. Clean ID
gale_lund_weight_cs_pheno[,"ID"] <- paste(gale_lund_weight_cs_pheno[,"Mating"], gale_lund_weight_cs_pheno[,"RIX_ID"], sep="_")

# observe duplicated rows, line 100
gale_lund_weight_cs_pheno[as.vector(unlist(sapply(unique(gale_lund_weight_cs_pheno[duplicated(gale_lund_weight_cs_pheno[,1]),1]),function(x)which(x==gale_lund_weight_cs_pheno[,1])))),]

# remove empty cols
empty_cols = c("CV1", "CV2", "Notes")

gale_lund_weight_cs_pheno[,-which(names(gale_lund_weight_cs_pheno) %in% empty_cols)] -> gale_lund_weight_cs_pheno_v2

# 3. Add Virus, Lab, Data_Altered, and Notes columns
gale_lund_weight_cs_pheno_v2$Virus = NA
gale_lund_weight_cs_pheno_v2[grep("M",gale_lund_weight_cs_pheno_v2[,"UWID"]),"Virus"] <- "Mock"
gale_lund_weight_cs_pheno_v2[-grep("M",gale_lund_weight_cs_pheno_v2[,"UWID"]),"Virus"] <- "WNV"

gale_lund_weight_cs_pheno_v2$Lab = "Gale"
gale_lund_weight_cs_pheno_v2$Data_Altered = NA
gale_lund_weight_cs_pheno_v2$Notes = NA

# check if any empty rows
days = c("D0","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15","D16","D17","D18","D19","D20","D21","D22","D23","D24","D25","D26","D27","D28")

all_na_rows = sapply(1:dim(gale_lund_weight_cs_pheno_v2)[1],function(x)sum(is.na(gale_lund_weight_cs_pheno_v2[x,days])) == length(gale_lund_weight_cs_pheno_v2[x,days]))
sum(all_na_rows)
gale_lund_weight_cs_pheno_v2[all_na_rows,] # 1 NA 44 jumpy mice, has already been flagged

# make version 3 copy
gale_lund_weight_cs_pheno_v3 = gale_lund_weight_cs_pheno_v2

# check if there are any question marks in score data
sum(sapply(days,function(x)sum(gale_lund_weight_cs_pheno_v2[,x] == "?",na.rm=T) > 0))

# 4. Add Death_Date, Death_Euthanized, Death_FoundInCage, Died_of_Virus, and Died_from_Anesthesia columns from weight data

# read in cleaned weight data
read.xls(xls="./Gale_12_8_15/processed_weight/Gale_Lab_weights_12_8_15_GC_flags.xlsx", sheet=1) -> gale_weight

# merge death data from weight with clinical score data
merge(gale_lund_weight_cs_pheno_v3, gale_weight[,c("ID","Death_Date","Death_Euthanized", "Death_FoundInCage", "Died_of_Virus", "Died_from_Anesthesia")], by="ID",all.x=T) -> gale_cs_death

# check if any death date's are empty
gale_cs_death[which(is.na(gale_cs_death[,"Death_Date"])),]

# add internal putative death day column
gale_cs_death$putative_death_day <- sapply(1:dim(gale_cs_death)[1],function(x)max(as.numeric(as.character(gsub("D","",names(gale_cs_death[x,days][which(!is.na(gale_cs_death[x,days]))]))))))

# check if any timepoints are not equal to putative death day and death date is not annotated
gale_cs_death[gale_cs_death[,"Timepoint"] - gale_cs_death[,"putative_death_day"] !=0 & is.na(gale_cs_death[,"Death_Date"]),]

# save copy
gale_cs_death_v3 = gale_cs_death

# 5. Add Flag_weight_day where Timepoint is >= 3 more than putative death day (day the weights are recorded up to), (internal check)

gale_cs_death_v3$Flag_Death_Day = NA

gale_cs_death_v3[which(gale_cs_death_v3[,"Timepoint"] - gale_cs_death_v3[,"putative_death_day"]>=3 & is.na(gale_cs_death_v3[,"Death_Euthanized"]) & is.na(gale_cs_death_v3[,"Death_FoundInCage"]) & is.na(gale_cs_death_v3[,"Died_of_Virus"]) & is.na(gale_cs_death_v3[,"Died_from_Anesthesia"])),"Flag_Death_Day"]

length(which(gale_cs_death_v3[,"Timepoint"] - gale_cs_death_v3[,"putative_death_day"]>=3 & is.na(gale_cs_death_v3[,"Death_Euthanized"]) & is.na(gale_cs_death_v3[,"Death_FoundInCage"]) & is.na(gale_cs_death_v3[,"Died_of_Virus"])&is.na(gale_cs_death_v3[,"Died_from_Anesthesia"]))) # there are no flags for death day

gale_cs_death_v3$Flag_Death_Day = FALSE

# 6. Add Flag_weight_date where death day > 28, (internal check)
gale_cs_death_v3$Flag_Death_Date = NA

# if no death day greater than 28, set flag to false
if(length(which(gale_cs_death_v3[,"Death_Euthanized"] > 28 | gale_cs_death_v3[,"Death_FoundInCage"] > 28 | gale_cs_death_v3[,"Died_of_Virus"] > 28 | gale_cs_death_v3[,"Died_from_Anesthesia"] > 28)) == 0){
  gale_cs_death_v3$Flag_Death_Date = FALSE
}


# 7. Alter data clinical score = 0 to NA past death day

length(which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Death_Euthanized"]))

length(which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Death_FoundInCage"]))

length(which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Died_from_Anesthesia"]))

length(which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Died_of_Virus"]))

for(i in 1:length(which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Died_from_Anesthesia"]))){
  print(i)
  gale_cs_death_v3[which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Died_from_Anesthesia"])[i],"Data_Altered"] <- "Yes"
  
  gale_cs_death_v3[which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Died_from_Anesthesia"])[i],"Notes"] <- "CS past death day set to NA"
  
  cols = gale_cs_death_v3[which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Died_from_Anesthesia"])[i],"Died_from_Anesthesia"] < as.numeric(as.character(gsub("D","",days)))
  
  gale_cs_death_v3[which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Died_from_Anesthesia"])[i],days][cols] <- NA
}

for(i in 1:length(which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Died_of_Virus"]))){
  print(i)
  gale_cs_death_v3[which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Died_of_Virus"])[i],"Data_Altered"] <- "Yes"
  
  gale_cs_death_v3[which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Died_of_Virus"])[i],"Notes"] <- "CS past death day set to NA"
  
  cols = gale_cs_death_v3[which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Died_of_Virus"])[i],"Died_of_Virus"] < as.numeric(as.character(gsub("D","",days)))
  
  gale_cs_death_v3[which(gale_cs_death_v3[,"putative_death_day"] != gale_cs_death_v3[,"Died_of_Virus"])[i],days][cols] <- NA
}

gale_cs_death_v3[which(gale_cs_death_v3[,"Data_Altered"] == "Yes"),]

# 8. Order names based on data dictionary
mycols <- rep("NULL", 4)
mycols[1] <- NA

read.xls(xls="/Users/choonoo/WNV/WNV_Data_Dictionary.xlsx", sheet="Clinical Score Data", colClasses=mycols) -> data_dict

# internal cols
internal_cols = c("putative_death_day","Flag_Death_Day", "Flag_Death_Date")
internal_cols_v2 = which(names(gale_cs_death_v3) %in% internal_cols)

# order columns by data dictionary and add internal column at the end
gale_cs_death_v3[,c(as.vector(unlist(sapply(data_dict[,1],function(x)which(x==names(gale_cs_death_v3))))),internal_cols_v2)] -> gale_cs_death_v4

# 9. Remove line 51
gale_cs_death_v4[-which(gale_cs_death_v4[,"UW_Line"] == 51),] -> gale_cs_death_v5

dim(gale_cs_death_v5)[1] == length(unique(gale_cs_death_v5[,"ID"]))

gale_cs_death_v5[which(is.na(gale_cs_death_v5[,"Death_Date"])),]

gale_cs_death_v5[,"Date_Infected"] <- as.Date(gale_cs_death_v5[,"Date_Infected"],format="%m/%d/%Y")

dup_ids = gale_cs_death_v5[which(duplicated(gale_cs_death_v5[,1])),1]

gale_cs_death_v5$Duplicated = NA
gale_cs_death_v5[as.vector(unlist(sapply(dup_ids,function(x)which(x==gale_cs_death_v5[,1])))),"Duplicated"] <- TRUE

gale_cs_death_v5[-as.vector(unlist(sapply(dup_ids,function(x)which(x==gale_cs_death_v5[,1])))),"Duplicated"] <- FALSE

gale_cs_death_v5[,"Death_Date"] <- as.Date(gale_cs_death_v5[,"Death_Date"],format="%Y-%m-%d")

# write processed data to file
write.table(x=gale_cs_death_v5, file="./Gale_12_8_15/processed_score/Gale_Lab_scores_12_8_15_GC.txt", sep="\t", quote=F, row.names=F)

