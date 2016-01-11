##########################################################

# WNV: qPCR data nov 24

##########################################################

# setwd("/Users/choonoo/WNV")

# save.image("~/gale_qpcr_nov_24.RData")

# load("~/gale_qpcr_nov_24.RData")

library(gdata)

##########################################################

# Functions

##########################################################

# clean data

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

# new data qpcr

##########################################################

# read in qpcr byline
read.xls(xls="/Users/choonoo/WNV/Gale_11_23_15/Gale_qPCR_byLine_11_23_15.xlsx", sheet=1) -> qpcr_data

#check dimensions
dim(qpcr_data)

# no duplicates
sum(duplicated(qpcr_data))

# check if any lines have more than one mating

# line 50 has 2 different matings
unique(qpcr_data[,"UW_Line"])[sapply(1:length(unique(qpcr_data[,"UW_Line"])), function(xx)length(unique(qpcr_data[as.vector(unlist(sapply(unique(qpcr_data[,"UW_Line"])[xx],function(x)which(x==qpcr_data[,"UW_Line"])))),"Mating"]))) > 1]

qpcr_data[which(qpcr_data[,"UW_Line"] == 50),"Mating"]

# clean line 50 mating 3393x8052 to 16680x8016

# 1. Add Data_Altered and Notes columns
qpcr_data$Data_Altered = NA
qpcr_data$Notes = NA

# Specific qPCR Data Corrections #1
qpcr_data[which(qpcr_data[,"UW_Line"] == 50),"Data_Altered"] <- "Yes"

qpcr_data[which(qpcr_data[,"UW_Line"] == 50),"Notes"] <- "Mating 3393x8052 changed to 16680x8016"

qpcr_data[which(qpcr_data[,"UW_Line"] == 50),"Mating"] <- "16680x8016" 

# now that we cleaned line 50, remove duplicates
sum(duplicated(qpcr_data))

# save as version 2
qpcr_data[!duplicated(qpcr_data),] -> qpcr_data_v2

# check number of unique matings equals number of unique lines
length(unique(qpcr_data_v2[,"UW_Line"])) == length(unique(qpcr_data_v2[,"Mating"]))

# 2. Add Virus column from time points
qpcr_data_v2$Virus = NA
qpcr_data_v2[grep("M",qpcr_data_v2[,"Timepoint"]),"Virus"] <- "Mock"
qpcr_data_v2[-grep("M",qpcr_data_v2[,"Timepoint"]),"Virus"] <- "WNV"

# 3. Remove 'M' from Time points and convert to numeric values
qpcr_data_v2[,"Timepoint"] <- as.numeric(as.character(gsub("M","",qpcr_data_v2[,"Timepoint"])))

# check experiment name
sum(names(summary(qpcr_data_v2[,"Experiment"])) == c("IFIT1","IFITM1", "IFNb1", "IL12b", "WNV")) == 5

# 4. Add Group column: UW Line, Timepoint, Virus, Tissue, Experiment separated by "_"
qpcr_data_v2$Group <- paste(qpcr_data_v2[,"UW_Line"],qpcr_data_v2[,"Timepoint"],qpcr_data_v2[,"Virus"],qpcr_data_v2[,"Tissue"],qpcr_data_v2[,"Experiment"], sep="_")

# 5. Add Lab column
qpcr_data_v2$Lab = "Gale"

# check number of duplicates is 0
sum(duplicated(qpcr_data_v2)) == 0

# check number of empty fc.mean is 0
sum(is.na(qpcr_data_v2[,"fc.mean"])) == 0

# check by mouse columns: "N","dCt.mean","dCt.sd"
# check byline calculations: "baseline.dCt","ddCt.mean"    "ddCt.sd", "fc.mean", "fc.sd"

##########################################################

# read in by mouse file
read.xls(xls="/Users/choonoo/WNV/Gale_11_23_15/Gale_qPCR_byMouse_11_23_15.xlsx", sheet=1) -> qpcr_data_mouse

# check if any CT < 15, none in this case
length(which(qpcr_data_mouse[,"Ct"] < 15))

# save version 2 copy
qpcr_data_mouse_v2 = qpcr_data_mouse

# 1. Remove 'M' from Time points and convert to numeric values
qpcr_data_mouse_v2[,"Timepoint"] <- as.numeric(as.character(gsub("M","",qpcr_data_mouse_v2[,"Timepoint"])))

names(summary(qpcr_data_mouse_v2[,"Mating"]))

# 2. Change "Condition" column name to "Virus"
names(qpcr_data_mouse_v2)[which(names(qpcr_data_mouse_v2) == "Condition")] <- "Virus"

# check tissue names
names(summary(qpcr_data_mouse_v2[,"Tissue"]))

# check experiment names
sum(names(summary(qpcr_data_mouse_v2[,"Experiment"])) == c("IFIT1","IFITM1", "IFNb1", "IL12b", "WNV")) == 5

# check if any lines have 2 different matings, none in this case
unique(qpcr_data_mouse_v2[,"UW_Line"])[sapply(1:length(unique(qpcr_data_mouse_v2[,"UW_Line"])), function(xx)length(unique(qpcr_data_mouse_v2[as.vector(unlist(sapply(unique(qpcr_data_mouse_v2[,"UW_Line"])[xx],function(x)which(x==qpcr_data_mouse_v2[,"UW_Line"])))),"Mating"]))) > 1]

# check line 50 mating
unique(qpcr_data_mouse_v2[which(qpcr_data_mouse_v2[,"UW_Line"] == 50),"Mating"])

# 3. Add Group: UW Line, Timepoint, Virus, Tissue, Experiment separated by "_"
qpcr_data_mouse_v2$Group <- paste(qpcr_data_mouse_v2[,"UW_Line"],qpcr_data_mouse_v2[,"Timepoint"], qpcr_data_mouse_v2[,"Virus"],qpcr_data_mouse_v2[,"Tissue"],qpcr_data_mouse_v2[,"Experiment"], sep="_")

# 4. Add Data_Altered and Notes columns
qpcr_data_mouse_v2$Data_Altered = NA
qpcr_data_mouse_v2$Notes = NA
qpcr_data_mouse_v2$Lab = "Gale"

##########################################################

# check numbers that go in byline

# dct mean
dct_mean = aggregate(formula=qpcr_data_mouse_v2[,"dCt"]~qpcr_data_mouse_v2[,"Group"], data=qpcr_data_mouse_v2, FUN=mean)
names(dct_mean) <- c("Group","dCt.mean.V2")
dct_mean[order(dct_mean[,1]),] -> dct_mean_order

byline_dct_mean = unique(qpcr_data_v2[,c("Group","dCt.mean")])
byline_dct_mean[order(byline_dct_mean[,1]),] -> byline_dct_mean_order

check_dct_mean = sapply(1:dim(dct_mean_order)[1],function(x)isTRUE(all.equal(dct_mean_order[x,2],byline_dct_mean_order[x,2])))

# all dct mean correct
sum(check_dct_mean) == dim(dct_mean_order)[1]

# code to correct
# save correct means
#byline_dct_mean_order[!check_dct_mean,1] -> flag_dct_mean
#dct_mean_order[!check_dct_mean,2] -> correct_dct_mean

# document change
#qpcr_data_v2[as.vector(unlist(sapply(flag_dct_mean,function(x)which(x==qpcr_data_v2[,"Group"])))),"Data_Altered"] <- 'Yes'

# overwrite with correction
#qpcr_data_v2[as.vector(unlist(sapply(flag_dct_mean,function(x)which(x==qpcr_data_v2[,"Group"])))),"dCt.mean"] <- correct_dct_mean

# check n
data.frame(summary(as.factor(qpcr_data_mouse_v2[,"Group"]),maxsum=8000)) -> bymouse_n
names(bymouse_n) <- c("N")

bymouse_n[,2] <- row.names(bymouse_n)
names(bymouse_n)[2] <- "Group"

bymouse_n[order(bymouse_n[,"Group"]),] -> bymouse_n_v2

qpcr_data_v2[,c("Group","N")] -> byline_n
byline_n[order(byline_n[,"Group"]),] -> byline_n_v2

# all n correct
sum(bymouse_n_v2[,1] == byline_n_v2[,2]) == dim(bymouse_n_v2)[1]

# code to correct
#bymouse_n_v2[-which(bymouse_n_v2[,1] == byline_n_v2[,2]),] -> correct_n

# document
#qpcr_data_v2[as.vector(unlist(sapply(correct_n[,2],function(x)which(x==qpcr_data_v2[,"Group"])))),"Data_Altered"] <- 'Yes'

#qpcr_data_v2[as.vector(unlist(sapply(correct_n[,2],function(x)which(x==qpcr_data_v2[,"Group"])))),"Notes"] <- paste0(qpcr_data_v2[as.vector(unlist(sapply(correct_n[,2],function(x)which(x==qpcr_data_v2[,"Group"])))),"Notes"], ", corrected N")

# make correction
#qpcr_data_v2[as.vector(unlist(sapply(correct_n[,2],function(x)which(x==qpcr_data_v2[,"Group"])))),"N"] <- correct_n[,1]

# check dct sd
dct_sd = aggregate(formula=qpcr_data_mouse_v2[,"dCt"]~qpcr_data_mouse_v2[,"Group"], data=qpcr_data_mouse_v2, FUN=sd)
names(dct_sd) <- c("Group","dCt.sd.V2")
dct_sd[order(dct_sd[,1]),] -> dct_sd_order

byline_dct_sd = qpcr_data_v2[,c("Group","dCt.sd")]
byline_dct_sd[order(byline_dct_sd[,1]),] -> byline_dct_sd_order

check_dct_sd = sapply(1:dim(dct_sd_order)[1],function(x)isTRUE(all.equal(dct_sd_order[x,2],byline_dct_sd_order[x,2])))

# all dct sd correct
sum(check_dct_sd) == dim(dct_sd_order)[1]

# code for correction
#byline_dct_sd_order[!check_dct_sd,1] -> flag_dct_sd
#dct_sd_order[!check_dct_sd,2] -> correct_dct_sd

# document
#qpcr_data_v2[as.vector(unlist(sapply(flag_dct_sd,function(x)which(x==qpcr_data_v2[,"Group"])))),"Data_Altered"] <- 'Yes'

#qpcr_data_v2[as.vector(unlist(sapply(flag_dct_sd,function(x)which(x==qpcr_data_v2[,"Group"])))),"Notes"] <- paste0(qpcr_data_v2[as.vector(unlist(sapply(flag_dct_sd,function(x)which(x==qpcr_data_v2[,"Group"])))),"Notes"], ", corrected dCt.sd")

# make correction
#qpcr_data_v2[as.vector(unlist(sapply(flag_dct_sd,function(x)which(x==qpcr_data_v2[,"Group"])))),"dCt.sd"] <- correct_dct_sd

# check baseline. dct

# add group g to annotate baseline
qpcr_data_v2$Group_g <- paste(qpcr_data_v2[,"UW_Line"],qpcr_data_v2[,"Tissue"], qpcr_data_v2[,"Experiment"],sep="_")

qpcr_data_v2$baseline.dCt.V2 = NA

# calculate baseline, baseline is 12 for this data
baseline = 12

for(i in unique(qpcr_data_v2[,"Group_g"])){
  print(i)
  qpcr_data_v2[which(qpcr_data_v2[,"Group_g"] == i),"baseline.dCt.V2"] <- qpcr_data_v2[which(qpcr_data_v2[,"Group_g"] == i & qpcr_data_v2[,"Virus"] == "Mock" & qpcr_data_v2[,"Timepoint"] == baseline),"dCt.mean"]
}

# check all baseline correct
sum(qpcr_data_v2[,"baseline.dCt"] == qpcr_data_v2[,"baseline.dCt.V2"]) == dim(qpcr_data_v2)[1]

# code to correct

# document 
#qpcr_data_v2[which(qpcr_data_v2[,"baseline.dCt"] != qpcr_data_v2[,"baseline.dCt.V2"]),c("Data_Altered")] <- 'Yes'

#qpcr_data_v2[which(qpcr_data_v2[,"baseline.dCt"] != qpcr_data_v2[,"baseline.dCt.V2"]),c("Notes")] <- paste0(qpcr_data_v2[which(qpcr_data_v2[,"baseline.dCt"] != qpcr_data_v2[,"baseline.dCt.V2"]),c("Notes")], ", corrected baseline.dCt")

# make correction
#qpcr_data_v2[which(qpcr_data_v2[,"baseline.dCt"] != qpcr_data_v2[,"baseline.dCt.V2"]),"baseline.dCt"] <- qpcr_data_v2[which(qpcr_data_v2[,"baseline.dCt"] != qpcr_data_v2[,"baseline.dCt.V2"]),"baseline.dCt.V2"]


# calculate ddct mean to check 
qpcr_data_v2$ddCt.mean.V2 <- as.numeric(as.character(qpcr_data_v2[,"dCt.mean"])) - as.numeric(as.character(qpcr_data_v2[,"baseline.dCt"]))

check_ddct_mean = sapply(1:dim(qpcr_data_v2)[1],function(x)isTRUE(all.equal(qpcr_data_v2[x,"ddCt.mean"], qpcr_data_v2[x,"ddCt.mean.V2"], tolerance=5.5e-8)))

# all ddct mean correct
sum(check_ddct_mean) == dim(qpcr_data_v2)[1]

all.equal(as.numeric(as.character(qpcr_data_v2[!check_ddct_mean,"ddCt.mean"] )),as.numeric(as.character(qpcr_data_v2[!check_ddct_mean,"ddCt.mean.V2"])))

# code to correct
# document
#qpcr_data_v2[!check_ddct_mean,"Data_Altered"] <- 'Yes'

#qpcr_data_v2[!check_ddct_mean,"Notes"] <- "corrected ddCt.mean"

# make correction
#qpcr_data_v2[!check_ddct_mean,"ddCt.mean"] <- qpcr_data_v2[!check_ddct_mean,"ddCt.mean.V2"]

#8. Remove ddCt.sd and fc.sd

# take out V2's if calculated correctly, and ddct.sd, fc.sd since not using
remove_cols = c("ddCt.sd", "fc.sd", "baseline.dCt.V2", "ddCt.mean.V2")

qpcr_data_v2[,-as.vector(unlist(sapply(remove_cols,function(x)which(x==names(qpcr_data_v2)))))] -> qpcr_data_v3

#6. Add baseline.dCt.sd column (this is the dCt.sd value for the baseline animal)
qpcr_data_v3$baseline.dCt.sd = NA

# compute baseline sd, use baseline = 12 saved above
for(i in unique(qpcr_data_v3[,"Group_g"])){
  print(i)
  qpcr_data_v3[which(qpcr_data_v3[,"Group_g"] == i),"baseline.dCt.sd"] <- qpcr_data_v3[which(qpcr_data_v3[,"Group_g"] == i & qpcr_data_v3[,"Virus"] == "Mock" & qpcr_data_v3[,"Timepoint"] == baseline),"dCt.sd"]
}

# 7. Add ddCt.se = sqrt(dCt.sd^2/N1 + baseline.dCt.sd^2/N2)
ddCt.se = sapply(1:dim(qpcr_data_v3)[1],function(x){
  sqrt(
    (qpcr_data_v3[x,"dCt.sd"]^2/qpcr_data_v3[x,"N"])+                                                                  (qpcr_data_v3[x,"baseline.dCt.sd"]^2/qpcr_data_v3[which(qpcr_data_v3[,"Group_g"] == qpcr_data_v3[x,"Group_g"] & qpcr_data_v3[,"Virus"] == "Mock" & qpcr_data_v3[,"Timepoint"] == baseline),"N"])
  )
})

qpcr_data_v3$ddCt.se <- ddCt.se

# Compute fc
qpcr_data_v3$fc.mean.V2 <- 2^-qpcr_data_v3[,"ddCt.mean"]

check_fc_mean = sapply(1:dim(qpcr_data_v3)[1],function(x)isTRUE(all.equal(qpcr_data_v3[x,"fc.mean"], qpcr_data_v3[x,"fc.mean.V2"])))

# check if all fc mean correct
sum(check_fc_mean) == dim(qpcr_data_v3)[1]

# code to correct

# document
#qpcr_data_v3[!check_fc_mean,"Data_Altered"] <- 'Yes'

#qpcr_data_v3[!check_fc_mean,"Notes"] <- paste0(qpcr_data_v3[!check_fc_mean,"Notes"], ", corrected fc.mean")

# make correction
#qpcr_data_v3[!check_fc_mean,"fc.mean"] <- qpcr_data_v3[!check_fc_mean,"fc.mean.V2"]

# remove fc.meanV2 if all correct
qpcr_data_v3[,!(names(qpcr_data_v3) %in% "fc.mean.V2")] -> qpcr_data_v4

# order by group
qpcr_data_v4[order(qpcr_data_v4[,"Group"]),] -> qpcr_data_v4_final

# double check all calculations are correct

# dct mean
isTRUE(all.equal(dct_mean[,2],qpcr_data_v4_final[,"dCt.mean"]))

# N
isTRUE(all.equal(bymouse_n[,"N"],qpcr_data_v4_final[,"N"]))

# dct sd
isTRUE(all.equal(dct_sd[,2],qpcr_data_v4_final[,"dCt.sd"]))

# baseline
isTRUE(all.equal(qpcr_data_v2[order(qpcr_data_v2[,"Group"]),"baseline.dCt.V2"],qpcr_data_v4_final[,"baseline.dCt"]))

# ddctmean
isTRUE(all.equal(qpcr_data_v2[order(qpcr_data_v2[,"Group"]),"ddCt.mean.V2"],qpcr_data_v4_final[,"ddCt.mean"]))

# fc mean
isTRUE(all.equal(qpcr_data_v3[order(qpcr_data_v3[,"Group"]),"fc.mean.V2"],qpcr_data_v4_final[,"fc.mean"]))

# baseline sd, computed using corrected dct sd's at timepoint 12

# ddct se, computed using correct dct sd and baseline sd

# remove group g
qpcr_data_v4_final[,!(names(qpcr_data_v4_final) %in% "Group_g")] -> qpcr_data_final_format



# order by data dictionary

# read in first column of qpcr by line
mycols <- rep("NULL", 4)
mycols[1] <- NA

read.xls(xls="/Users/choonoo/WNV/WNV_Data_Dictionary.xlsx", sheet="qPCR Data - By Line", colClasses=mycols) -> data_dict

qpcr_data_final_format[,as.vector(unlist(sapply(data_dict[,1],function(x)which(x==names(qpcr_data_final_format)))))] -> qpcr_data_final_format_v2

# read in first column of qpcr by mouse
mycols <- rep("NULL", 4)
mycols[1] <- NA

read.xls(xls="/Users/choonoo/WNV/WNV_Data_Dictionary.xlsx", sheet="qPCR Data - By Mouse", colClasses=mycols) -> data_dict_mouse

qpcr_data_mouse_v2[,as.vector(unlist(sapply(data_dict_mouse[,1],function(x)which(x==names(qpcr_data_mouse_v2)))))] -> qpcr_data_mouse_v3

# check both byline and bymouse have the same UW lines
names(summary(as.factor(qpcr_data_mouse_v3[,"UW_Line"]))) == names(summary(as.factor(qpcr_data_final_format_v2[,"UW_Line"])))

# write processed data to files
write.table(file="./Gale_11_23_15/processed_qpcr/Gale_qPCR_byLine_GC_11_23_15.txt", x=qpcr_data_final_format_v2, sep="\t", quote=F, row.names=F)

write.table(file="./Gale_11_23_15/processed_qpcr/Gale_qPCR_byMouse_GC_11_23_15.txt", x=qpcr_data_mouse_v3, sep="\t", quote=F, row.names=F)

##########################################################

# add to full qpcr

##########################################################

# full byline
read.xls(xls="/Users/choonoo/U19 WNV Data Clean Archive Shared/9-Nov-2015/Gale_qPCR_byLine_9-Nov-2015_final.xlsx", sheet=1) -> qpcr_line

names(qpcr_line)
names(qpcr_data_final_format_v2)
rbind(qpcr_line, qpcr_data_final_format_v2) -> qpcr_11_25_15

write.table(file="./Gale_11_23_15/processed_qpcr/Gale_qPCR_byLine_GC_full_11_25_15.txt", x=qpcr_11_25_15, sep="\t", quote=F, row.names=F)

# full bymouse
read.xls(xls="/Users/choonoo/U19 WNV Data Clean Archive Shared/9-Nov-2015/Gale_qPCR_byMouse_9-Nov-2015_final.xlsx", sheet=1) -> qpcr_bymouse

rbind(qpcr_bymouse, qpcr_data_mouse_v3) -> qpcr_bymouse_11_25_15

write.table(file="./Gale_11_23_15/processed_qpcr/Gale_qPCR_byMouse_GC_full_11_25_15.txt", x=qpcr_bymouse_11_25_15, sep="\t", quote=F, row.names=F)


