##########################################################

# lund clinical score 12/3/15

##########################################################

# setwd("/Users/choonoo/WNV")

# save.image("~/lund_cs_dec3.RData")

# load("~/lund_cs_dec3.RData")

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

# process updated clinical score
read.xls(xls="./Lund_CS_12_2_15/CS_Lund_UW019.xlsx", sheet=1,stringsAsFactors=F) -> line_19

##########################################################
# Code to read in other file formats
##########################################################

# read.xls(xls="./Lund_CS_11_9_2015/072 clinical score.xlsx", sheet=1,stringsAsFactors=F) -> line_72
# 
# read.xls(xls="./Lund_CS_11_9_2015/082 clinical scores 11_5_15.xlsx", sheet=1, skip=1, colClasses=c("NULL",rep(NA,35),"NULL","NULL","NULL"),stringsAsFactors=F) -> line_82
# 
# line_72[1:which(line_72[,1] == "***CS only for this strain, no weights")-1,] -> line_72_v2
# 
# # clean columns
# 
# # remove sex column
# line_72_v2[,-5] -> line_72_v3
# 
# # add UWID, lab, data altered, notes
# line_72_v3[,37:40] <- NA
# names(line_72_v3)[37:40] <- c("UWID","Lab","Data_Altered","Notes")
# 
# names(line_72_v3)
# 
# names(line_82)[7:35] <- names(line_72_v3)[8:36]
# names(line_50)[7:35] <- names(line_72_v3)[8:36]
# 
# names(line_82)[1] <- "UW_Line"
# names(line_82)[3] <- "Mating"
# names(line_82)[4] <- "RIX_ID"
# names(line_82)[6] <- "Date_Infected"
# 
# names(line_50)[1] <- "UW_Line"
# names(line_50)[3] <- "Mating"
# names(line_50)[4] <- "RIX_ID"
# names(line_50)[6] <- "Date_Infected"
# 
# line_50[,36:40] <- NA
# names(line_50)[36:40]<- c("UWID","Virus","Lab","Data_Altered","Notes")
# 
# line_82[,36:40] <- NA
# names(line_82)[36:40]<- c("UWID","Virus","Lab","Data_Altered","Notes")
# 
# line_72_v3[,c(1,2,3,4,37,5,6,7,38,39,40,8:36)] -> line_72_v4
# 
# line_82[,as.vector(unlist(sapply(names(line_72_v4),function(x)which(x==names(line_82)))))] -> line_82_v2
# 
# line_50[,as.vector(unlist(sapply(names(line_72_v4),function(x)which(x==names(line_50)))))] -> line_50_v2
# 
# 
# # cbind
# rbind(line_50_v2, line_72_v4, line_82_v2) -> new_cs
# 
# clean_na(new_cs) -> new_cs_v2

##########################################################

# make copy
new_cs = line_19

# 1. Clean mating: "X" changed to "x"
new_cs[,"Mating"] <- gsub("X","x",new_cs[,"Mating"])

# 2. Clean ID
new_cs[,"ID"] <- paste(new_cs[,"Mating"], new_cs[,"RIX_ID"],sep="_")

# 3. Add UWID, Lab, Data_Altered, and Notes columns
new_cs$Lab <- "Lund"
new_cs$UWID <- NA
new_cs$Data_Altered <- NA
new_cs$Notes <- NA

# 4. Remove Sex column
new_cs[,!names(new_cs)%in%"Sex"] -> new_cs_v2

# 5. Add Death_Date, Death_Euthanized, Death_FoundInCage, Died_of_Virus, and Died_from_Anesthesia columns from weight data

# read in cleaned weight sheets
read.xls(xls="./Lund_Weight_12_2_15/Lund_Weight_12_3_15_GC.xlsx", sheet=1,stringsAsFactors=F) -> lund_weight

# merge clinical score with weight with deaths columns only
merge(new_cs_v2, lund_weight[,c("ID","Death_Date","Death_Euthanized", "Death_FoundInCage", "Died_of_Virus", "Died_from_Anesthesia")], by="ID",all.x=T) -> new_cs_v3

# clean timepoint
new_cs_v3[,"Timepoint"] <- as.numeric(as.character(gsub("d","",new_cs_v3[,"Timepoint"])))

# format date infected
new_cs_v3[,"Date_Infected"] <- as.character(as.Date(new_cs_v3[,"Date_Infected"], format="%m/%d/%y"))

##########################################################
# Code to extract notes from Day columns if necessary
##########################################################

# new_cs_v4[which(new_cs_v4[,"D8"] == "Found in cage 3/25/15"),"Data_Altered"] <- "Yes"
# 
# new_cs_v4[which(new_cs_v4[,"D8"] == "Found in cage 3/25/15"),"Notes"] <- "Found in cage 3/25/15"
# 
# new_cs_v4[which(new_cs_v4[,"D8"] == "Found in cage 3/25/15"),"Death_Date"] <- "2015-03-25"
# 
# new_cs_v4[which(new_cs_v4[,"D8"] == "Found in cage 3/25/15"),"Death_FoundInCage"] <- difftime(new_cs_v4[which(new_cs_v4[,"D8"] == "Found in cage 3/25/15"),"Death_Date"],new_cs_v4[which(new_cs_v4[,"D8"] == "Found in cage 3/25/15"),"Date_Infected"])
# 
# new_cs_v4[which(new_cs_v4[,"D8"] == "Found in cage 3/25/15"),"D8"] <- NA

# merge with older clinical score data if missing death date

# merge(new_cs_v4, lund_cs_final_11_16_15_order[,c("ID","Virus","Date_Infected","Death_Date","Death_Euthanized","Death_FoundInCage")], by="ID",all.x=T) -> new_cs_v5
# 
# new_cs_v5[which(is.na(new_cs_v5[,"Virus.x"])),"Virus.x"] <- new_cs_v5[which(is.na(new_cs_v5[,"Virus.x"])),"Virus.y"]
# 
# 
# new_cs_v5[which(is.na(new_cs_v5[,"Death_Date.x"])),"Death_Date.x"] <- new_cs_v5[which(is.na(new_cs_v5[,"Death_Date.x"])),"Death_Date.y"]
# 
# new_cs_v5[which(is.na(new_cs_v5[,"Date_Infected.x"])),"Date_Infected.x"] <- new_cs_v5[which(is.na(new_cs_v5[,"Date_Infected.x"])),"Date_Infected.y"]
# 
# new_cs_v5[which(is.na(new_cs_v5[,"Death_Euthanized.x"])),"Death_Euthanized.x"] <- new_cs_v5[which(is.na(new_cs_v5[,"Death_Euthanized.x"])),"Death_Euthanized.y"]
# 
# new_cs_v5[which(is.na(new_cs_v5[,"Death_FoundInCage.x"])),"Death_FoundInCage.x"] <- new_cs_v5[which(is.na(new_cs_v5[,"Death_FoundInCage.x"])),"Death_FoundInCage.y"]
# 
# new_cs_v5[which(new_cs_v5[,"UW_Line"] == 72),]
# 
# new_cs_v5[,-c(46:50)] -> new_cs_v6
# 
# names(new_cs_v6) <- gsub(".x","",names(new_cs_v6))
# 
# new_cs_v6[which(new_cs_v6[,"UW_Line"] == 72),]

# # D8, D10, D11
# new_cs_v6[which(new_cs_v6[,"D8"] == "Found in cage 3/25/15" | new_cs_v6[,"D8"] == "3, euthanized"),"Data_Altered"] <- "Yes"
# 
# new_cs_v6[which(new_cs_v6[,"D8"] == "Found in cage 3/25/15" | new_cs_v6[,"D8"] == "3, euthanized"),"Notes"] <- "D8 Euthanized"
# 
# new_cs_v6[which(new_cs_v6[,"D8"] == "Found in cage 3/25/15" | new_cs_v6[,"D8"] == "3, euthanized"),"D8"] <- c(NA, 3)
# 
# new_cs_v6[which(new_cs_v6[,"D10"] == "1, euthanized"),"Data_Altered"] <- "Yes"
# 
# new_cs_v6[which(new_cs_v6[,"D10"] == "1, euthanized"),"Notes"] <- "D10 Euthanized"
# 
# new_cs_v6[which(new_cs_v6[,"D10"] == "1, euthanized"),"D10"] <- 1
# 
# new_cs_v6[which(new_cs_v6[,"D11"] == "found in cage"),"Data_Altered"] <- "Yes"
# 
# new_cs_v6[which(new_cs_v6[,"D11"] == "found in cage"),"Notes"] <- "D11 Found in cage"
# 
# new_cs_v6[which(new_cs_v6[,"D11"] == "found in cage"),"D11"] <- NA
# 
# # change to days
# new_cs_v6[-which(is.na(new_cs_v6[,"Death_Euthanized"])),"Death_Euthanized"] <- as.vector(difftime(new_cs_v6[-which(is.na(new_cs_v6[,"Death_Euthanized"])),"Death_Date"],new_cs_v6[-which(is.na(new_cs_v6[,"Death_Euthanized"])),"Date_Infected"],units="days"))

# new_cs_v6[-which(is.na(new_cs_v6[,"Death_FoundInCage"])),"Death_FoundInCage"] <- as.vector(difftime(new_cs_v6[-which(is.na(new_cs_v6[,"Death_FoundInCage"])),"Death_Date"],new_cs_v6[-which(is.na(new_cs_v6[,"Death_FoundInCage"])),"Date_Infected"],units="days"))

##########################################################

# make copy if not using above code
new_cs_v7 = new_cs_v3

# Create early death flags

# add putative death day
days = c("D0","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12","D13","D14","D15","D16","D17","D18","D19","D20","D21","D22","D23","D24","D25","D26","D27","D28")

new_cs_v7$putative_death_day <- sapply(1:dim(new_cs_v7)[1],function(x)max(as.numeric(as.character(gsub("D","",names(new_cs_v7[x,days][which(!is.na(new_cs_v7[x,days]))]))))))

# 6. Add Flag_weight_day where Timepoint is >= 3 more than putative death day (day the weights are recorded up to), (internal check)
new_cs_v7$Flag_Death_Day = NA

new_cs_v7[which(new_cs_v7[,"Timepoint"] - new_cs_v7[,"putative_death_day"]>=3 & is.na(new_cs_v7[,"Death_Euthanized"]) & is.na(new_cs_v7[,"Death_FoundInCage"])),"Flag_Death_Day"] <- TRUE

length(which(new_cs_v7[,"Timepoint"] - new_cs_v7[,"putative_death_day"]>=3 & is.na(new_cs_v7[,"Death_Euthanized"]) & is.na(new_cs_v7[,"Death_FoundInCage"]))) # no flags, set to false

new_cs_v7$Flag_Death_Day = FALSE

# 7. Add Flag_weight_date where death day > 28, (internal check)
new_cs_v7$Flag_Death_Date = NA

if(length(which(new_cs_v7[,"Death_Euthanized"] > 28 | new_cs_v7[,"Death_FoundInCage"] > 28 | new_cs_v7[,"Died_of_Virus"] > 28 | new_cs_v7[,"Died_from_Anesthesia"] > 28)) == 0){
  new_cs_v7$Flag_Death_Date <- FALSE
}


# 8. Alter data clinical score = 0 to NA past death day

length(which(new_cs_v7[,"putative_death_day"] != new_cs_v7[,"Death_Euthanized"]))

length(which(new_cs_v7[,"putative_death_day"] != new_cs_v7[,"Death_FoundInCage"]))

length(which(new_cs_v7[,"putative_death_day"] != new_cs_v7[,"Died_from_Anesthesia"]))

length(which(new_cs_v7[,"putative_death_day"] != new_cs_v7[,"Died_of_Virus"]))

# clinical score past death day for death found in cage
for(i in 1:length(which(new_cs_v7[,"putative_death_day"] != new_cs_v7[,"Death_FoundInCage"]))){
  print(i)
  new_cs_v7[which(new_cs_v7[,"putative_death_day"] != new_cs_v7[,"Death_FoundInCage"])[i],"Data_Altered"] <- "Yes"
  
  new_cs_v7[which(new_cs_v7[,"putative_death_day"] != new_cs_v7[,"Death_FoundInCage"])[i],"Notes"] <- "CS past death day set to NA"
  
  cols = new_cs_v7[which(new_cs_v7[,"putative_death_day"] != new_cs_v7[,"Death_FoundInCage"])[i],"Death_FoundInCage"] < as.numeric(as.character(gsub("D","",days)))
  
  new_cs_v7[which(new_cs_v7[,"putative_death_day"] != new_cs_v7[,"Death_FoundInCage"])[i],days][cols] <- NA
}

new_cs_v7[which(new_cs_v7[,"Data_Altered"] == "Yes"),]

# 9. Order names based on data dictionary
mycols <- rep("NULL", 4)
mycols[1] <- NA

read.xls(xls="/Users/choonoo/WNV/WNV_Data_Dictionary.xlsx", sheet="Clinical Score Data", colClasses=mycols) -> data_dict

# internal cols
internal_cols = c("putative_death_day","Flag_Death_Day", "Flag_Death_Date")
internal_cols_v2 = which(names(new_cs_v7) %in% internal_cols)

# order columns by data dictionary and add internal column at the end
new_cs_v7[,c(as.vector(unlist(sapply(data_dict[,1],function(x)which(x==names(new_cs_v7))))),internal_cols_v2)] -> new_cs_v8


write.table(x=new_cs_v8, file="./Lund_CS_12_2_15/Lund_Score_11_30_2015_upload_GC.txt", sep="\t", quote=F, row.names=F)


##########################################################

# read in full version
read.xls(xls="./Lund_CS_12_2_15/Lund_Score_11_23_15.xlsx", sheet=1,stringsAsFactors=F) -> cs_full

# read in updated 11/9
read.xls(xls="./Lund_CS_12_2_15/Lund_Score_11_5_and_9_2015_upload_GC.xlsx", sheet=1,stringsAsFactors=F) -> cs_11_9

# merge
rbind(cs_full, cs_11_9) -> cs_merge

# remove any lines from full that are in updated data 11/30
remove_lines = intersect(cs_merge[,"UW_Line"], new_cs_v7[,"UW_Line"])

# remove new 11/30 upload
cs_merge[-as.vector(unlist(sapply(remove_lines,function(x)which(x==cs_merge[,"UW_Line"])))),] -> cs_merge_v2

# clean 11/30, if no questions, add to full
rbind(cs_merge_v2, new_cs_v8) -> lund_cs_full_12_3

# write processed file, has all corrections
write.table(file="./Lund_CS_12_2_15/Lund_CS_12_3_15_GC.txt", x=lund_cs_full_12_3, sep="\t", quote=F, row.names=F)

dim(lund_cs_full_12_3)[1] == length(unique(lund_cs_full_12_3[,"ID"]))




