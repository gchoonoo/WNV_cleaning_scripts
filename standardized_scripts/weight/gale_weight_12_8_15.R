##########################################################

# Gale Weight 12/8/15

##########################################################

# setwd("/Users/choonoo/WNV")

# save.image("~/gale_weight_dec8.RData")

# load("~/gale_weight_dec8.RData")

library(gdata)

##########################################################

# Functions

##########################################################

# clean NA

clean_na = function(data_set){
  
  for (i in 1:dim(data_set)[2]){
    print(i)    
    if(sum(na.omit(data_set[,i] == "")) > 0){
      data_set[which(data_set[,i] == ""),i] <- NA
    }
    
  }
  return(data_set)
}

##########################################################

# weight cleaning steps

##########################################################

# read in data and flag NA, duplicated, no mock, high weight change, and putative death day
read.xls(xls="/Users/choonoo/WNV/Gale_12_8_15/Gale_Lab_weights_12_8_15.xlsx", sheet=1) -> new_gale_weight

# 1. Clean mating: "X" changed to "x"
dim(unique(new_gale_weight[,c("Mating","UW_Line")]))[1] == length(unique(new_gale_weight[,"UW_Line"]))

# clean mating
new_gale_weight[,"Mating"] <- gsub("X","x",new_gale_weight[,"Mating"])

# interger 0 if all matings have unique uw line
unique(new_gale_weight[,"UW_Line"])[sapply(1:length(unique(new_gale_weight[,"UW_Line"])), function(xx)length(unique(new_gale_weight[as.vector(unlist(sapply(unique(new_gale_weight[,"UW_Line"])[xx],function(x)which(x==new_gale_weight[,"UW_Line"])))),"Mating"]))) > 1]

new_gale_weight[which(new_gale_weight[,"UW_Line"] == 82),"Mating"]

# Specific data correction #2, clean line 82
new_gale_weight[which(new_gale_weight[,"UW_Line"] == 82),"Mating"] <- "16557x3154"

# 2. Recreate ID with new mating: Mating + RIX_ID
new_gale_weight[,"ID"] <- paste(new_gale_weight[,"Mating"], new_gale_weight[,"RIX_ID"], sep="_")

# check if any duplicated rows
sum(duplicated(new_gale_weight))

# check if any all NA weight starting at D0
days = c("D0","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15","D16","D17","D18","D19","D20","D21","D22","D23","D24","D25","D26","D27","D28","D29","D30","D31","D32","D45")
na_now = sapply(1:dim(new_gale_weight)[1],function(x)sum(is.na(new_gale_weight[x,days])) == length(new_gale_weight[x,days]))

sum(na_now)

new_gale_weight[na_now,] # died from anesthesia, no need to check

# 3. Add Lab column
new_gale_weight$Lab = "Gale"

# annotate putative death day and check early deaths with no notes, some return -Inf which can be interpreted as NA
putative_death_day = sapply(1:dim(new_gale_weight)[1],function(x){
  print(x)
  days[max(which(!is.na(new_gale_weight[x,c(days)])), na.rm=T)]
})

new_gale_weight$putative_death_day <- as.numeric(as.character(gsub("D","",putative_death_day)))

# check if any timepoints are greater than or equal to 3 days than putative death day
new_gale_weight[which(new_gale_weight[,"Timepoint"]-new_gale_weight$putative_death_day >= 3),]

# annotate actual death day
new_gale_weight$Actual_death_day = NA

# actual death day = timepoint if no died early notes or equals the died early day
new_gale_weight[which(!is.na(new_gale_weight[,"Died_from_anesthesia"])),"Actual_death_day"] <- new_gale_weight[which(!is.na(new_gale_weight[,"Died_from_anesthesia"])),"Died_from_anesthesia"]

new_gale_weight[which(!is.na(new_gale_weight[,"Animal_was_euthansized_at_day"])),"Actual_death_day"] <- new_gale_weight[which(!is.na(new_gale_weight[,"Animal_was_euthansized_at_day"])),"Animal_was_euthansized_at_day"]

new_gale_weight[which(!is.na(new_gale_weight[,"Died_of_virus_at_day"])),"Actual_death_day"] <- new_gale_weight[which(!is.na(new_gale_weight[,"Died_of_virus_at_day"])),"Died_of_virus_at_day"]

new_gale_weight[which(is.na(new_gale_weight[,"Actual_death_day"])),"Actual_death_day"] <- new_gale_weight[which(is.na(new_gale_weight[,"Actual_death_day"])),"Timepoint"]

# 4. Add Death date based on date of infection to timepoint or died early notes
new_gale_weight$Death_Date = as.Date(new_gale_weight[,"Date_Infected"],format="%m/%d/%Y")+new_gale_weight[,"Actual_death_day"]

# 5. Add Death_FoundInCage column
new_gale_weight$Death_FoundInCage = NA

# 6. Change column name "Died_of_virus_at_day" to "Died_of_Virus" 
names(new_gale_weight)[which(names(new_gale_weight) == "Died_of_virus_at_day")] <- "Died_of_Virus"

# 7. Change column name "Animal_was_euthansized_at_day" to "Death_Euthanized"
names(new_gale_weight)[which(names(new_gale_weight) == "Animal_was_euthansized_at_day")] <- "Death_Euthanized"

# 8. Change column name "Died_from_anesthesia" to "Died_from_Anesthesia"
names(new_gale_weight)[which(names(new_gale_weight) == "Died_from_anesthesia")] <- "Died_from_Anesthesia"

# 9. Remove columns: CV1, CV2, Interesting_Phenotype, Notes (empty)
remove_cols = c("CV1", "CV2", "Interesting_phenotype", "Notes", "Actual_death_day")

new_gale_weight[,-which(names(new_gale_weight) %in% remove_cols)] -> new_gale_weight_v2

# 10. Add columns Data_Altered and Notes
new_gale_weight_v2$Data_Altered = NA
new_gale_weight_v2$Notes = NA

# 11. Add D0_Percentage column, set to baseline 0
new_gale_weight_v2$D0_Percentage = 0

## 12. Recalculate weight change percentages

# make version 3 copy of data
new_gale_weight_v3 = new_gale_weight_v2

# percentages
percent_days = days[-1]

for(i in percent_days){
  print(i)
  d2_p = sapply(1:dim(new_gale_weight_v3)[1],function(x){
    print(x)
    (diff(c(new_gale_weight_v3[x,"D0"], new_gale_weight_v3[x,i]))/new_gale_weight_v3[x,"D0"])*100
  })
  
  new_gale_weight_v3[,which(names(new_gale_weight_v3) == paste(i,"Percentage",sep="_"))] <- d2_p
}

days_percent = c("D1_Percentage","D2_Percentage","D3_Percentage","D4_Percentage","D5_Percentage","D6_Percentage","D7_Percentage","D8_Percentage","D9_Percentage","D10_Percentage","D11_Percentage","D12_Percentage","D13_Percentage","D14_Percentage","D15_Percentage","D16_Percentage","D17_Percentage","D18_Percentage","D19_Percentage","D20_Percentage","D21_Percentage","D22_Percentage","D23_Percentage","D24_Percentage","D25_Percentage","D26_Percentage","D27_Percentage","D28_Percentage","D29_Percentage","D30_Percentage","D31_Percentage","D32_Percentage","D45_Percentage")

# 13. Add Flag_Weight_Drop column (True if animals lost 20% body weight)
flag_weight_drop = sapply(1:dim(new_gale_weight_v3)[1],function(x)sum(new_gale_weight_v3[x,days_percent] <= -20, na.rm=T) > 0)

new_gale_weight_v3$Flag_Weight_Drop = flag_weight_drop

# 14. Add Flag_Identical_Weights (True if identical weights on consecutive measurements)
identical_weights = sapply(1:dim(new_gale_weight_v3)[1],function(x)sum(abs(diff(na.omit(as.numeric(as.character(new_gale_weight_v3[x,days_percent]))))) == 0) > 0)

new_gale_weight_v3$Flag_Identical_Weights = identical_weights

# 15. Add Flag_Per_Day_Weight_Change (True if weight change > 10% on consecutive measurements)
weight_change_flag = sapply(1:dim(new_gale_weight_v3)[1],function(x)sum(diff(na.omit(as.numeric(as.character(new_gale_weight_v3[x,days_percent])))) >= 10) > 0)

new_gale_weight_v3$Flag_Per_Day_Weight_Change = weight_change_flag

# 16. Add Flag_weight_day where Timepoint is >= 3 more than putative death day (day the weights are recorded up to), (internal check)
new_gale_weight_v3$Flag_Death_Day = NA

new_gale_weight_v3[which(new_gale_weight_v3[,"Timepoint"] - new_gale_weight_v3[,"putative_death_day"]>=3 & is.na(new_gale_weight_v3[,"Death_Euthanized"]) & is.na(new_gale_weight_v3[,"Death_FoundInCage"]) & is.na(new_gale_weight_v3[,"Died_of_Virus"])&is.na(new_gale_weight_v3[,"Died_from_Anesthesia"])),"Flag_Death_Day"]

# if length is 0 (all early deaths have been documented), set flag death day to false
length(which(new_gale_weight_v3[,"Timepoint"] - new_gale_weight_v3[,"putative_death_day"]>=3 & is.na(new_gale_weight_v3[,"Death_Euthanized"]) & is.na(new_gale_weight_v3[,"Death_FoundInCage"]) & is.na(new_gale_weight_v3[,"Died_of_Virus"])&is.na(new_gale_weight_v3[,"Died_from_Anesthesia"]))) #none

new_gale_weight_v3$Flag_Death_Day = FALSE

# 17. Add Flag_weight_date where death day > 28, (internal check)
new_gale_weight_v3$Flag_Death_Date = NA

# if no death day extends past 28, set to false
if(length(which(new_gale_weight_v3[,"Death_Euthanized"] > 28 | new_gale_weight_v3[,"Death_FoundInCage"] > 28 | new_gale_weight_v3[,"Died_of_Virus"] > 28 | new_gale_weight_v3[,"Died_from_Anesthesia"] > 28)) == 0){
  new_gale_weight_v3$Flag_Death_Date = FALSE
}

# 18. Order columns by data dictionary

# read first column of weight sheet in data dictionary
mycols <- rep("NULL", 4)
mycols[1] <- NA

read.xls(xls="/Users/choonoo/WNV/WNV_Data_Dictionary.xlsx", sheet="Weight Data", colClasses=mycols) -> data_dict

# internal cols
internal_cols = c("putative_death_day","Flag_Death_Day", "Flag_Death_Date")
internal_cols_v2 = which(names(new_gale_weight_v3) %in% internal_cols)

# order columns by data dictionary and add internal column at the end
new_gale_weight_v3[,c(as.vector(unlist(sapply(data_dict[,1],function(x)which(x==names(new_gale_weight_v3))))),internal_cols_v2)] -> new_gale_weight_v4

# 19. Remove line 51
new_gale_weight_v4[-which(new_gale_weight_v4[,"UW_Line"] == 51),] -> new_gale_weight_v5

# check no duplicates
dim(new_gale_weight_v5)[1] == length(unique(new_gale_weight_v5[,"ID"]))

# observe duplicates
dup_ids = new_gale_weight_v5[which(duplicated(new_gale_weight_v5[,1])),1]

new_gale_weight_v5$Duplicated = NA
new_gale_weight_v5[as.vector(unlist(sapply(dup_ids,function(x)which(x==new_gale_weight_v5[,1])))),"Duplicated"] <- TRUE

new_gale_weight_v5[-as.vector(unlist(sapply(dup_ids,function(x)which(x==new_gale_weight_v5[,1])))),"Duplicated"] <- FALSE

# line 100?
new_gale_weight_v5[which(new_gale_weight_v5[,"UW_Line"] == 100),]

new_gale_weight_v5[which(new_gale_weight_v5$Flag_Death_Day),]
new_gale_weight_v5[which(new_gale_weight_v5$Flag_Death_Date),]
new_gale_weight_v5[which(new_gale_weight_v5$Flag_Weight_Drop),]
new_gale_weight_v5[which(new_gale_weight_v5$Flag_Identical_Weights),]
new_gale_weight_v5[which(new_gale_weight_v5$Flag_Per_Day_Weight_Change),]

new_gale_weight_v5[which(new_gale_weight_v5[,"Duplicated"]),]

# write processed version to file
write.table(file="./Gale_12_8_15/processed_weight/Gale_Lab_weights_12_8_15_GC_flags.txt", x=new_gale_weight_v5, sep="\t", quote=F, row.names=F)

# fix flag per day
read.xls(xls="./Gale_12_8_15/processed_weight/Gale_Lab_weights_12_8_15_GC_flags.xlsx", sheet=1) -> gale_data

gale_data$Flag_Per_Day_Weight_Change = NA

# check weights that changed more than 10% day to day
weight_change_flag = sapply(1:dim(gale_data)[1],function(x)sum(abs(diff(na.omit(as.numeric(as.character(gale_data[x,51:79]))))) >= 10) > 0)

gale_data$Flag_Per_Day_Weight_Change = weight_change_flag

# write to table
write.table(file="./Gale_12_8_15/processed_weight/Gale_Lab_weights_12_11_15_GC_flags.txt", x=gale_data, sep="\t", quote=F, row.names=F)
