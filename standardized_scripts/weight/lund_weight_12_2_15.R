##########################################################

# lund weight 12/2/15

##########################################################

# setwd("/Users/choonoo/WNV")

# save.image("~/lund_weight_dec2.RData")

# load("~/lund_weight_dec2.RData")

library(gdata)

##########################################################

clean_na = function(data_set){
  
  for (i in 1:dim(data_set)[2]){
    print(i)
    if(sum(na.omit(data_set[,i] == "") > 0)){
      
      data_set[which( data_set[,i] == ""),i] <- NA
    }
    
  }
  return(data_set)
}
##########################################################

# read in updated data 11/30

##########################################################

read.xls(xls="./Lund_Weight_12_2_15/012 weights 7_2_15.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)), stringsAsFactors=F) -> line_12

read.xls(xls="./Lund_Weight_12_2_15/022 weights 7_2_15.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)),stringsAsFactors=F) -> line_22

read.xls(xls="./Lund_Weight_12_2_15/024 weights 11_9_15.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)),stringsAsFactors=F) -> line_24

read.xls(xls="./Lund_Weight_12_2_15/039 weights 7_5_15.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)),stringsAsFactors=F) -> line_39

read.xls(xls="./Lund_Weight_12_2_15/043 weights 7_5_15.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)),stringsAsFactors=F) -> line_43

read.xls(xls="./Lund_Weight_12_2_15/056 weights 7_6_15.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)),stringsAsFactors=F) -> line_56

read.xls(xls="./Lund_Weight_12_2_15/067 weights 7_5_15.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)),stringsAsFactors=F) -> line_67

read.xls(xls="./Lund_Weight_12_2_15/082 weights 11_5_15.xlsx", sheet=1, skip=2,colClasses=c("NULL",rep(NA,37)),stringsAsFactors=F) -> line_82

# clean score days
days = gsub("d","D",gsub("\\.","",names(line_82)[9:37]))

# format files, remove percentages in second half of each file
format_data = function(data){
  
  names(data) <- c("UW_Line","ID","Mating","RIX_ID","Timepoint","Date_Infected","Death_Euthanized","Death_FoundInCage",days)
  
  data[1:which(data[,1] == "% of initial weight")-1,] -> data
  
  return(data)
}

format_data(line_12) -> lund1
format_data(line_22) -> lund2
format_data(line_24) -> lund3
format_data(line_39) -> lund4
format_data(line_43) -> lund5
format_data(line_56) -> lund6
format_data(line_67) -> lund7
format_data(line_82) -> lund8


# collapse
rbind(lund1, lund2, lund3, lund4, lund5, lund6, lund7,lund8) -> new_lund_line

# Cleaning steps as outlined by Readme

# 1. Clean mating: "X" changed to "x"
new_lund_line[,"Mating"] <- gsub("X","x",new_lund_line[,"Mating"])

# 2. Recreate ID with new mating: Mating + RIX_ID
new_lund_line[,"ID"] <- paste(new_lund_line[,"Mating"],new_lund_line[,"RIX_ID"],sep="_")

# 3. Add Lab column
new_lund_line$Lab = "Lund"

# 4. Change Death_Euthanized and Death_FoundInCage to day instead of date
new_lund_line[,"Death_Euthanized"] <- as.Date(new_lund_line[,"Death_Euthanized"], format="%m/%d/%y")

new_lund_line[,"Death_FoundInCage"] <- as.Date(new_lund_line[,"Death_FoundInCage"], format="%m/%d/%y")

new_lund_line[,"Date_Infected"] <- as.Date(new_lund_line[,"Date_Infected"], format="%m/%d/%y")

new_lund_line[,"Death_FoundInCage"] <- as.vector(difftime(new_lund_line[,"Death_FoundInCage"],new_lund_line[,"Date_Infected"],units="days"))

new_lund_line[,"Death_Euthanized"] <- as.vector(difftime(new_lund_line[,"Death_Euthanized"],new_lund_line[,"Date_Infected"],units="days"))

# add virus column
new_lund_line$Virus = NA
new_lund_line[grep("m",new_lund_line[,"Timepoint"]),"Virus"] <- "Mock"
new_lund_line[-grep("m",new_lund_line[,"Timepoint"]),"Virus"] <- "WNV"

# clean timepoint
new_lund_line[,"Timepoint"] <- as.numeric(as.character(gsub("d","",gsub("m","",new_lund_line[,"Timepoint"]))))

# 5. Add UWID, Cohort, Sex, Died_of_Virus, Died_from_Anesthesia, Data_Altered, and Notes columns
new_lund_line$UWID = NA
new_lund_line$Cohort = NA
new_lund_line$Sex = NA
new_lund_line$Died_of_Virus = NA
new_lund_line$Died_from_Anesthesia = NA
new_lund_line$Data_Altered = NA
new_lund_line$Notes = NA

# 6. Add Death date based on date of infection to timepoint or died early notes
new_lund_line$Death_Date = NA

new_lund_line[which(!is.na(new_lund_line[,"Death_Euthanized"])),"Death_Date"] <- as.character(new_lund_line[which(!is.na(new_lund_line[,"Death_Euthanized"])),"Date_Infected"]+new_lund_line[which(!is.na(new_lund_line[,"Death_Euthanized"])),"Death_Euthanized"])

new_lund_line[which(!is.na(new_lund_line[,"Death_FoundInCage"])),"Death_Date"] <- as.character(new_lund_line[which(!is.na(new_lund_line[,"Death_FoundInCage"])),"Date_Infected"]+new_lund_line[which(!is.na(new_lund_line[,"Death_FoundInCage"])),"Death_FoundInCage"])

new_lund_line[which(is.na(new_lund_line[,"Death_FoundInCage"]) & is.na(new_lund_line[,"Death_Euthanized"])),"Death_Date"] <- as.character(new_lund_line[which(is.na(new_lund_line[,"Death_FoundInCage"]) & is.na(new_lund_line[,"Death_Euthanized"])),"Date_Infected"]+new_lund_line[which(is.na(new_lund_line[,"Death_FoundInCage"]) & is.na(new_lund_line[,"Death_Euthanized"])),"Timepoint"])

# Specific data correction, #4 update line 67 animals 43,44,45
# 43
new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 43),"Death_Date"] <- "2015-02-06"

new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 43),"Death_FoundInCage"] <- round(difftime(new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 43),"Death_Date"],new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 43),"Date_Infected"],units="days"),0)

new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 43),"Data_Altered"] <- "Yes"

new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 43),"Notes"] <- "Changed Death Date 2025-02-06 to 2015-02-06"

#44
new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 44),"Death_Date"] <- "2015-02-01"

new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 44),"Death_FoundInCage"] <- round(difftime(new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 44),"Death_Date"],new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 44),"Date_Infected"],units="days"),0)

new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 44),"Data_Altered"] <- "Yes"

new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 44),"Notes"] <- "Changed Death Date 2015-03-01 to 2015-02-01"

#45
new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 45),"Death_Date"] <- "2015-02-01"

new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 45),"Death_Euthanized"] <- round(difftime(new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 45),"Death_Date"],new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 45),"Date_Infected"],units="days"),0)

new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 45),"Data_Altered"] <- "Yes"

new_lund_line[which(new_lund_line[,"UW_Line"] == "067" & new_lund_line[,"RIX_ID"] == 45),"Notes"] <- "Changed Death Date 2015-03-01 to 2015-02-01"

# order columns
new_lund_line[,c("ID","Mating","RIX_ID","UW_Line","UWID","Cohort","Sex","Virus","Date_Infected","Death_Date","Timepoint","Death_Euthanized","Death_FoundInCage","Died_of_Virus","Died_from_Anesthesia","Lab","Data_Altered","Notes",days)] -> new_lund_line_v2

# 7. Recalculate weight change percentages

# add weight percentage columns
new_lund_line_v2[,48:76] <- NA
names(new_lund_line_v2)[48:76] <- paste(names(new_lund_line_v2)[19:47], "Percentage",sep="_")

# set baseline
new_lund_line_v2[,"D0_Percentage"] <- 0

# change weights to numeric
for(i in 19:47){
  print(i)
  new_lund_line_v2[,i] <- as.numeric(as.character(new_lund_line_v2[,i]))
  
}

# calculate percentages
for(i in 20:47){
  print(i)
  d2_p = sapply(1:dim(new_lund_line_v2)[1],function(x){
    print(x)
    (diff(c(new_lund_line_v2[x,"D0"], new_lund_line_v2[x,i]))/new_lund_line_v2[x,"D0"])*100
  })
  
  new_lund_line_v2[,i+29] <- d2_p
}

# 8. Add Flag_Weight_Drop column (True if animals lost 20% body weight)

# check high weight drop
flag_weight_drop = sapply(1:dim(new_lund_line_v2)[1],function(x)sum(new_lund_line_v2[x,48:76] <= -20, na.rm=T) > 0)

new_lund_line_v2$Flag_Weight_Drop = flag_weight_drop

# 9. Add Flag_Identical_Weights (True if identical weights on consecutive measurements)

# check identical weights
identical_weights = sapply(1:dim(new_lund_line_v2)[1],function(x)sum(diff(na.omit(as.numeric(as.character(new_lund_line_v2[x,48:76])))) == 0) > 0)

new_lund_line_v2$Flag_Identical_Weights = identical_weights

# 10. Add Flag_Per_Day_Weight_Change (True if weight change > 10% on consecutive measurements)

# check weights that changed more than 10% day to day
weight_change_flag = sapply(1:dim(new_lund_line_v2)[1],function(x)sum(diff(na.omit(as.numeric(as.character(new_lund_line_v2[x,48:76])))) >= 10) > 0)

new_lund_line_v2$Flag_Per_Day_Weight_Change = weight_change_flag

# make version 3 copy
new_lund_line_v3 = new_lund_line_v2

# 11. Add Flag_weight_day where Timepoint is >= 3 more than putative death day (day the weights are recorded up to), (internal check)
new_lund_line_v3$putative_death_day <- sapply(1:dim(new_lund_line_v3)[1],function(x)max(as.numeric(as.character(gsub("D","",names(new_lund_line_v3[x,days][which(!is.na(new_lund_line_v3[x,days]))]))))))

new_lund_line_v3$Flag_Death_Day = NA

new_lund_line_v3[which(new_lund_line_v3[,"Timepoint"] - new_lund_line_v3[,"putative_death_day"]>=3 & is.na(new_lund_line_v3[,"Death_Euthanized"]) & is.na(new_lund_line_v3[,"Death_FoundInCage"])),"Flag_Death_Day"] <- TRUE

new_lund_line_v3[-which(new_lund_line_v3[,"Timepoint"] - new_lund_line_v3[,"putative_death_day"]>=3 & is.na(new_lund_line_v3[,"Death_Euthanized"]) & is.na(new_lund_line_v3[,"Death_FoundInCage"])),"Flag_Death_Day"] <- FALSE

# 12. Add Flag_weight_date where death day > 28, (internal check)
new_lund_line_v3$Flag_Death_Date = NA

new_lund_line_v3[which(new_lund_line_v3[,"Death_Euthanized"] > 28 | new_lund_line_v3[,"Death_FoundInCage"] > 28),"Flag_Death_Date"] <- TRUE

new_lund_line_v3[-which(new_lund_line_v3[,"Death_Euthanized"] > 28 | new_lund_line_v3[,"Death_FoundInCage"] > 28),"Flag_Death_Date"] <- FALSE

new_lund_line_v3[which(new_lund_line_v3$Flag_Death_Date),]

# 13. Reorder by data dictionary
mycols <- rep("NULL", 4)
mycols[1] <- NA

read.xls(xls="/Users/choonoo/WNV/WNV_Data_Dictionary.xlsx", sheet="Weight Data", colClasses=mycols, stringsAsFactors=F) -> data_dict

# internal cols
internal_cols = c("putative_death_day","Flag_Death_Day", "Flag_Death_Date")
internal_cols_v2 = which(names(new_lund_line_v3) %in% internal_cols)

# order columns by data dictionary and add internal column at the end
new_lund_line_v3[,c(as.vector(unlist(sapply(data_dict[,1],function(x)which(x==names(new_lund_line_v3))))),internal_cols_v2)] -> new_lund_line_v4

# 14. Remove beginning 0 in UW_Line
new_lund_line_v4[,"UW_Line"] <- gsub("0","",new_lund_line_v4[,"UW_Line"])

# write processed file
write.table(file="./Lund_Weight_12_2_15/Lund_Weight_11_30_15_upload_GC.txt", x=new_lund_line_v4, sep="\t", quote=F, row.names=F)

##########################################################

# 15. Add updated data to full cleaned data set

##########################################################

# read in full data
read.xls(xls="./Lund_Weight_12_2_15/Lund_Weight_11_18_2015_GC.xlsx",sheet=1, stringsAsFactors=F) -> lund_weight_full

# read in updated data from 11/9
read.xls(xls="./Lund_Weight_12_2_15/Lund_Weight_11_5_and_9_2015_upload_GC.xlsx",sheet=1,stringsAsFactors=F) -> lund_weight_11_9

# merge full and updated 11/9
rbind(lund_weight_full, lund_weight_11_9) -> lund_weight_merge

# remove any lines from full that are in updated data 11/30
remove_lines = intersect(new_lund_line_v4[,"UW_Line"], lund_weight_merge[,"UW_Line"])

lund_weight_merge[-as.vector(unlist(sapply(remove_lines,function(x)which(x==lund_weight_merge[,"UW_Line"])))),] -> lund_weight_merge_v2

# update date format
lund_weight_merge_v2[,"Date_Infected"] <- as.Date(lund_weight_merge_v2[,"Date_Infected"],format="%Y-%m-%d")

# clean 11/30, if no questions, add to full
rbind(lund_weight_merge_v2, new_lund_line_v4) -> lund_weight_full_12_3

# write processed file, has all corrections
write.table(file="./Lund_Weight_12_2_15/Lund_Weight_12_3_15_GC.txt", x=lund_weight_full_12_3, sep="\t", quote=F, row.names=F)

dim(lund_weight_full_12_3)[1] == length(unique(lund_weight_full_12_3[,"ID"]))
