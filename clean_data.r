###############
# Read in Data
###############

library(gdata)

# Gale: Read in weight sheet
# Lund: Read in all weight sheets and combine
# Lund Data Example:
read.xls(xls="102 weights 3_21_16.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_102
read.xls(xls="101 weights 3_21_16.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_101
read.xls(xls="096 weights 3_21_16.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_96
read.xls(xls="108 weights 3_21_16.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_108
read.xls(xls="107 weights 3_21_16.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_107
read.xls(xls="100 weights 3_21_16.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_100
read.xls(xls="099 weights 3_21_16.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_99
read.xls(xls="106 weights 3_15_16.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_106
read.xls(xls="105 weights 3_15_16.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_105
read.xls(xls="061 weights 3_15_16.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_61
read.xls(xls="090 weights 3_15_16.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_90
read.xls(xls="089 weights 3_15_16.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_89
read.xls(xls="045 weights 11_5_15.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_45

# Clean days
days = gsub("d","D",gsub("\\.","",names(line_102)[9:37]))

# Formatting function: removes percentages in second half of each file
format_data = function(data){
  
  names(data) <- c("UW_Line","ID","Mating","RIX_ID","Timepoint","Date_Infected","Death_Euthanized","Death_FoundInCage",days)
  
  data[1:which(data[,1] == "% of initial weight")-1,] -> data
  
  return(data)
}

# Format all data
format_data(line_102) -> lund1
format_data(line_101) -> lund2
format_data(line_96) -> lund3
format_data(line_108) -> lund4
format_data(line_107) -> lund5
format_data(line_100) -> lund6
format_data(line_99) -> lund7
format_data(line_106) -> lund8
format_data(line_105) -> lund9
format_data(line_61) -> lund10
format_data(line_90) -> lund11
format_data(line_89) -> lund12
format_data(line_45) -> lund13

# Combine Data
rbind(lund1, lund2, lund3, lund4, lund5, lund6, lund7, lund8, lund9, lund10, lund11, lund12, lund13) -> new_lund_line

# Cleaning steps as outlined by Readme

# Cleaning Blanks to NA function
clean_na = function(data_set){
  
  for (i in 1:dim(data_set)[2]){
    print(i)
    if(sum(na.omit(data_set[,i] == "") > 0)){
      
      data_set[which( data_set[,i] == ""),i] <- NA
    }
    
  }
  return(data_set)
}

clean_na(new_lund_line) -> new_lund_line

# 1. Clean mating: "X" changed to "x"
new_lund_line[,"Mating"] <- gsub("X","x",new_lund_line[,"Mating"])

# Check if all matings have unique uw line, if not returns lines that has duplicated mating
unique(new_lund_line[,"UW_Line"])[sapply(1:length(unique(new_lund_line[,"UW_Line"])), function(xx)length(unique(new_lund_line[as.vector(unlist(sapply(unique(new_lund_line[,"UW_Line"])[xx],function(x)which(x==new_lund_line[,"UW_Line"])))),"Mating"]))) > 1]

# Clean mating if necessary: see Gale and Lund specific cleaning steps in Read me

# 2. Clean ID with new mating: Mating + RIX_ID
new_lund_line[,"ID"] <- paste(new_lund_line[,"Mating"],new_lund_line[,"RIX_ID"],sep="_")

# 3. Add/remove columns and edit column names: Add Lab, Data_Altered
  # Gale: Add Death_FoundInCage column
  # Change column name "Died_of_virus_at_day" to "Died_of_Virus" 
  # Change column name "Animal_was_euthansized_at_day" to "Death_Euthanized"
  # Change column name "Died_from_anesthesia" to "Died_from_Anesthesia"
  # Remove columns: CV1, CV2, Interesting_Phenotype (empty)
# Lund: Add UWID, Cohort, Sex, Died_of_Virus, Died_from_Anesthesia, and Notes columns
new_lund_line$Lab = "Lund" 
new_lund_line$UWID = NA
new_lund_line$Cohort = NA
new_lund_line$Sex = NA
new_lund_line$Died_of_Virus = NA
new_lund_line$Died_from_Anesthesia = NA
new_lund_line$Data_Altered = NA
new_lund_line$Notes = NA

# 4. Check the Virus column is annotated and clean timepoint
unique(new_lund_line$Virus)

# add virus column
new_lund_line$Virus = NA
new_lund_line[grep("m",new_lund_line[,"Timepoint"]),"Virus"] <- "Mock"
new_lund_line[-grep("m",new_lund_line[,"Timepoint"]),"Virus"] <- "WNV"

# clean timepoint
new_lund_line[,"Timepoint"] <- as.numeric(as.character(gsub("d","",gsub("m","",new_lund_line[,"Timepoint"]))))

# 5. Check for any duplicated or NA rows
# check if any duplicated rows
sum(duplicated(new_lund_line[,"ID"]))

new_lund_line[,"ID"][which(duplicated(new_lund_line[,"ID"]))] # line 61, see specified notes, this line not recieved

# check if any all NA weight starting at D0
days = c("D0","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15","D16","D17","D18","D19","D20","D21","D22","D23","D24","D25","D26","D27","D28")
na_now = sapply(1:dim(new_lund_line)[1],function(x)sum(is.na(new_lund_line[x,days])) == length(new_lund_line[x,days]))

sum(na_now)

new_lund_line[na_now,] # Flag if no died early notes

# 6. Annotate Death_date based on date of infection to timepoint or died early notes, annotate putative death day

# Lund: Change Death_Euthanized and Death_FoundInCage to day instead of date
new_lund_line[,"Death_Euthanized"] <- as.Date(new_lund_line[,"Death_Euthanized"], format="%m/%d/%y")

new_lund_line[,"Death_FoundInCage"] <- as.Date(new_lund_line[,"Death_FoundInCage"], format="%m/%d/%y")

new_lund_line[,"Date_Infected"] <- as.Date(new_lund_line[,"Date_Infected"], format="%m/%d/%y")

new_lund_line[,"Death_FoundInCage"] <- as.vector(difftime(new_lund_line[,"Death_FoundInCage"],new_lund_line[,"Date_Infected"],units="days"))

new_lund_line[,"Death_Euthanized"] <- as.vector(difftime(new_lund_line[,"Death_Euthanized"],new_lund_line[,"Date_Infected"],units="days"))

# Annotate death date
new_lund_line$Death_Date = NA

new_lund_line[which(is.na(new_lund_line[,"Date_Infected"])),"UW_Line"]

new_lund_line[which(!is.na(new_lund_line[,"Death_Euthanized"])),"Death_Date"] <- as.character(new_lund_line[which(!is.na(new_lund_line[,"Death_Euthanized"])),"Date_Infected"]+new_lund_line[which(!is.na(new_lund_line[,"Death_Euthanized"])),"Death_Euthanized"])

new_lund_line[which(!is.na(new_lund_line[,"Death_FoundInCage"])),"Death_Date"] <- as.character(new_lund_line[which(!is.na(new_lund_line[,"Death_FoundInCage"])),"Date_Infected"]+new_lund_line[which(!is.na(new_lund_line[,"Death_FoundInCage"])),"Death_FoundInCage"])

new_lund_line[which(is.na(new_lund_line[,"Death_FoundInCage"]) & is.na(new_lund_line[,"Death_Euthanized"])),"Death_Date"] <- as.character(new_lund_line[which(is.na(new_lund_line[,"Death_FoundInCage"]) & is.na(new_lund_line[,"Death_Euthanized"])),"Date_Infected"]+new_lund_line[which(is.na(new_lund_line[,"Death_FoundInCage"]) & is.na(new_lund_line[,"Death_Euthanized"])),"Timepoint"])

# annotate putative death day
new_lund_line$putative_death_day <- sapply(1:dim(new_lund_line)[1],function(x)max(as.numeric(as.character(gsub("D","",names(new_lund_line[x,days][which(!is.na(new_lund_line[x,days]))]))))))

# Check if any still NA, annotate death date as date of infection to timepoint if the difference between the timepoint and putative death day is <= 3, otherwise this is flagged below
new_lund_line[which(is.na(new_lund_line[,"Death_Date"])),]

# 7. Calculate weight change percentages, add D0_Percentage column and set to baseline 0

# add weight percentage columns
ncol(new_lund_line) -> ncols
new_lund_line[,(ncols+1):(ncols+29)] <- NA
names(new_lund_line)[(ncols+1):(ncols+29)] <- paste(days, "Percentage",sep="_")

# Make sure D0_Percentage column exists and is set to baseline = 0
new_lund_line[,"D0_Percentage"] <- 0

# change weights to numeric
for(i in days){
  print(i)
  new_lund_line[,i] <- as.numeric(as.character(new_lund_line[,i]))
}

# calculate percentages
for(i in days[-1]){
  print(i)
  d2_p = sapply(1:dim(new_lund_line)[1],function(x){
    print(x)
    (diff(c(new_lund_line[x,"D0"], new_lund_line[x,i]))/new_lund_line[x,"D0"])*100
  })
  
  new_lund_line[,paste0(i,"_Percentage")] <- d2_p
}

# 8. Add Flag_Weight_Drop column (True if animals lost 20% body weight)
days_percent = paste0(days,"_Percentage")

flag_weight_drop = sapply(1:dim(new_lund_line)[1],function(x)sum(new_lund_line[x,days_percent] <= -20, na.rm=T) > 0)

new_lund_line$Flag_Weight_Drop = flag_weight_drop

# 9. Add Flag_Identical_Weights (True if identical weights on consecutive measurements)
identical_weights = sapply(1:dim(new_lund_line)[1],function(x)sum(diff(na.omit(as.numeric(as.character(new_lund_line[x,days_percent])))) == 0) > 0)

new_lund_line$Flag_Identical_Weights = identical_weights

# 10. Add Flag_Per_Day_Weight_Change (True if weight change > 10% on consecutive measurements)
weight_change_flag = sapply(1:dim(new_lund_line)[1],function(x)sum(abs(diff(na.omit(as.numeric(as.character(new_lund_line[x,days_percent]))))) >= 10) > 0)

new_lund_line$Flag_Per_Day_Weight_Change = weight_change_flag

# 11. Add Flag_weight_day where Timepoint is >= 3 more than putative death day (day the weights are recorded up to) (internal check)
new_lund_line$Flag_Death_Day = NA

if(length(which(new_lund_line[,"Timepoint"] - new_lund_line[,"putative_death_day"]>=3 & is.na(new_lund_line[,"Death_Euthanized"]) & is.na(new_lund_line[,"Death_FoundInCage"]))) == 0){
  new_lund_line$Flag_Death_Day = FALSE
}else{
  new_lund_line[which(new_lund_line[,"Timepoint"] - new_lund_line[,"putative_death_day"]>=3 & is.na(new_lund_line[,"Death_Euthanized"]) & is.na(new_lund_line[,"Death_FoundInCage"])),"Flag_Death_Day"] <- TRUE
  
  new_lund_line[-which(new_lund_line[,"Timepoint"] - new_lund_line[,"putative_death_day"]>=3 & is.na(new_lund_line[,"Death_Euthanized"]) & is.na(new_lund_line[,"Death_FoundInCage"])),"Flag_Death_Day"] <- FALSE
}

# 12. Add Flag_weight_date where death day > timepoint, (internal check)
new_lund_line$Death_Date_greater_timepoint = NA

if(length(which(new_lund_line[,"Death_Euthanized"] > new_lund_line[,"Timepoint"] | new_lund_line[,"Death_FoundInCage"] > new_lund_line[,"Timepoint"] | new_lund_line[,"Died_of_Virus"] > new_lund_line[,"Timepoint"] | new_lund_line[,"Died_from_Anesthesia"] > new_lund_line[,"Timepoint"])) >0){
  
  new_lund_line[which(new_lund_line[,"Death_Euthanized"] > new_lund_line[,"Timepoint"] | new_lund_line[,"Death_FoundInCage"] > new_lund_line[,"Timepoint"] | new_lund_line[,"Died_of_Virus"] > new_lund_line[,"Timepoint"] | new_lund_line[,"Died_from_Anesthesia"] > new_lund_line[,"Timepoint"]),"Death_Date_greater_timepoint"] <- TRUE
  
  new_lund_line[-which(new_lund_line[,"Death_Euthanized"] > new_lund_line[,"Timepoint"] | new_lund_line[,"Death_FoundInCage"] > new_lund_line[,"Timepoint"] | new_lund_line[,"Died_of_Virus"] > new_lund_line[,"Timepoint"] | new_lund_line[,"Died_from_Anesthesia"] > new_lund_line[,"Timepoint"]),"Death_Date_greater_timepoint"] <- FALSE
  
}else{
  new_lund_line$Death_Date_greater_timepoint = FALSE
}

# 13. Order columns according to the data dictionary
# Read in data dictionary
mycols <- rep("NULL", 4)
mycols[1] <- NA

read.xls(xls="WNV_Data_Dictionary.xlsx", sheet="Weight Data", colClasses=mycols, stringsAsFactors=F) -> data_dict

# internal cols
internal_cols = c("putative_death_day","Flag_Death_Day", "Death_Date_greater_timepoint")
internal_cols_v2 = which(names(new_lund_line) %in% internal_cols)

# order columns by data dictionary and add internal column at the end
new_lund_line[,c(as.vector(unlist(sapply(data_dict[,1],function(x)which(x==names(new_lund_line))))),internal_cols_v2)] -> new_lund_line_v2

# 14. Annotate array_exists column
# Read in array annotation
read.table(file="all_expression_ids.txt") -> exp_id

new_lund_line_v2$array_exists = NA

if(length(as.vector(unlist(sapply(as.character(exp_id[,1]),function(x)which(x==as.character(new_lund_line_v2[,1])))))) >0){
  new_lund_line_v2[as.vector(unlist(sapply(as.character(exp_id[,1]),function(x)which(x==as.character(new_lund_line_v2[,1]))))),"array_exists"] <- "Yes"
  
  new_lund_line_v2[-as.vector(unlist(sapply(as.character(exp_id[,1]),function(x)which(x==as.character(new_lund_line_v2[,1]))))),"array_exists"] <- "No"
}else{
  new_lund_line_v2$array_exists = "No"
}


# Lund: Remove beginning 0 in UW_Line
new_lund_line_v2[which(substring(new_lund_line_v2[,"UW_Line"],1,1) == "0"),"UW_Line"] <- substring(new_lund_line_v2[which(substring(new_lund_line_v2[,"UW_Line"],1,1) == "0"),"UW_Line"],2,nchar(new_lund_line_v2[which(substring(new_lund_line_v2[,"UW_Line"],1,1) == "0"),"UW_Line"]))

# 15. Make any specific alterations listed in the Readme. Note if you change weights, need to also re-calculate weight percentages. Record these changes in Data_Altered and Notes column

# Most of these are already saved in the old versions of the data, which is updated in the next step.

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 61 & new_lund_line_v2[,"Virus"] == "Mock" & new_lund_line_v2[,"Timepoint"] == 28),"Notes"] <- "Not receiving, line complete"

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 89 & new_lund_line_v2[,"Timepoint"] > 12 & new_lund_line_v2[,"Virus"] == "WNV"),"Notes"] <- "Jumpy, no weight"

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 90 & new_lund_line_v2[,"RIX_ID"] ==76),"Death_Euthanized"] <- "11"

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 90 & new_lund_line_v2[,"RIX_ID"] ==76),"Notes"] <- "Changed death euthanized to 11"

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 90 & new_lund_line_v2[,"RIX_ID"] ==76),"Data_Altered"] <- "Yes"

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 90 & new_lund_line_v2[,"RIX_ID"] ==76),"Death_Date"] <- new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 90 & new_lund_line_v2[,"RIX_ID"] ==76),"Date_Infected"]+11

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 100 & (new_lund_line_v2[,"RIX_ID"] ==82 | new_lund_line_v2[,"RIX_ID"] ==83)),"Death_Euthanized"] <- c(10,9)

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 100 & (new_lund_line_v2[,"RIX_ID"] ==82 | new_lund_line_v2[,"RIX_ID"] ==83)),"Data_Altered"] <- "Yes"

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 100 & (new_lund_line_v2[,"RIX_ID"] ==82 | new_lund_line_v2[,"RIX_ID"] ==83)),"Notes"] <- "Changed Death euthanized to one day less"

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 100 & (new_lund_line_v2[,"RIX_ID"] ==82 | new_lund_line_v2[,"RIX_ID"] ==83)),"Death_Date"] <- new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 100 & (new_lund_line_v2[,"RIX_ID"] ==82 | new_lund_line_v2[,"RIX_ID"] ==83)),"Date_Infected"]+c(10,9)

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 102 & (new_lund_line_v2[,"RIX_ID"] ==62 | new_lund_line_v2[,"RIX_ID"] ==64)),"D7"] <- c(35.24,35.30)

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 102 & (new_lund_line_v2[,"RIX_ID"] ==62 | new_lund_line_v2[,"RIX_ID"] ==64)),"Notes"] <- "changed 25.xx to 35.xx"

new_lund_line_v2[which(new_lund_line_v2[,"UW_Line"] == 102 & (new_lund_line_v2[,"RIX_ID"] ==62 | new_lund_line_v2[,"RIX_ID"] ==64)),"Data_Altered"] <- "Yes"

# calculate percentages
for(i in days[-1]){
  print(i)
  d2_p = sapply(1:dim(new_lund_line_v2)[1],function(x){
    print(x)
    (diff(c(new_lund_line_v2[x,"D0"], new_lund_line_v2[x,i]))/new_lund_line_v2[x,"D0"])*100
  })
  
  new_lund_line_v2[,paste0(i,"_Percentage")] <- d2_p
}

# 16. Update Data: Add flags_checked column
# Gale: Cumulative data, make sure old lines are the same. Set flag for old ID's that they have already been verified
# Lund: Need to read in old data and update any overlapping lines and add the new lines. Set flag for old ID's that have already been verified
# read in full data
read.xls(xls="Lund_Weight_22-Mar-2016_final.xlsx", sheet=1, stringsAsFactors=F) -> lund_weight_full

# remove any lines from full that are in updated data
remove_lines = intersect(new_lund_line_v2[,"UW_Line"], lund_weight_full[,"UW_Line"])

# update line 45
length(remove_lines) == 0

lund_weight_full[-which(lund_weight_full[,"UW_Line"] == 45),] -> lund_weight_full

lund_weight_full$flags_checked = TRUE

new_lund_line_v2$flags_checked = FALSE

# clean date format
lund_weight_full[,"Death_Date"] <- as.character(lund_weight_full[,"Death_Date"])
lund_weight_full[,"Date_Infected"] <- as.Date(lund_weight_full[,"Date_Infected"],format="%Y-%m-%d")

new_lund_line_v2[,"Death_Date"] <- as.character(new_lund_line_v2[,"Death_Date"])
new_lund_line_v2[,"Date_Infected"] <- as.Date(as.character(new_lund_line_v2[,"Date_Infected"]),format="%Y-%m-%d")

# Merge new data and old data
merge(lund_weight_full, new_lund_line_v2, by=intersect(names(new_lund_line_v2),names(lund_weight_full)), all=T) -> lund_weight_full_v2

lund_weight_full_v2[which(lund_weight_full_v2[,"Data_Altered"] == "Yes"),"flags_checked"] <- TRUE

lund_weight_full_v2[which(lund_weight_full_v2[,"Notes"] == "Not receiving, line complete"),"flags_checked"] <- TRUE

lund_weight_full_v2[which(lund_weight_full_v2[,"flags_checked"]),c("putative_death_day")] <- lund_weight_full_v2[which(lund_weight_full_v2[,"flags_checked"]),c("Timepoint")]

lund_weight_full_v2[which(lund_weight_full_v2[,"flags_checked"]),c("Flag_Death_Day")] <- FALSE

lund_weight_full_v2[which(lund_weight_full_v2[,"flags_checked"]),c("Death_Date_greater_timepoint")] <- FALSE

# Save cleaned file
write.table(file="./Lund_Weight_cleaned.txt", x=lund_weight_full_v2, sep="\t", quote=F, row.names=F, na="")
