##########################################################

# WNV: qPCR data nov 24: 3032, 8002

##########################################################

# setwd("/Users/choonoo/WNV")

# save.image("~/gale_qpcr_3032_RI.RData")

# load("~/gale_qpcr_3032_RI.RData")

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

# new data qpcr, 3032

##########################################################

# read in qpcr
read.xls(xls="/Users/choonoo/WNV/Gale_11_23_15/Gale_3032_qPCR_byLine.xlsx", sheet=1) -> qpcr_data

#check dimensions
dim(qpcr_data)

# check if any duplicates
sum(duplicated(qpcr_data))

# no mating in this data
#table(qpcr_data[,"Mating"], qpcr_data[,"UW_Line"])

#names(summary(qpcr_data[,"Mating"]))

#unique(qpcr_data[,"UW_Line"])[sapply(1:length(unique(qpcr_data[,"UW_Line"])), function(xx)length(unique(qpcr_data[as.vector(unlist(sapply(unique(qpcr_data[,"UW_Line"])[xx],function(x)which(x==qpcr_data[,"UW_Line"])))),"Mating"]))) > 1]

# add data altered and notes column
qpcr_data$Data_Altered = NA
qpcr_data$Notes = NA

sum(duplicated(qpcr_data))

# remove duplicated rows and save as version 2
qpcr_data[!duplicated(qpcr_data),] -> qpcr_data_v2

# clean timepoint and add virus column
qpcr_data_v2$Virus = NA
qpcr_data_v2[grep("M",qpcr_data_v2[,"Timepoint"]),"Virus"] <- "Mock"
qpcr_data_v2[-grep("M",qpcr_data_v2[,"Timepoint"]),"Virus"] <- "WNV"

qpcr_data_v2[,"Timepoint"] <- as.numeric(as.character(gsub("M","",qpcr_data_v2[,"Timepoint"])))

# check experiment, all viral titer in this data set
names(summary(qpcr_data_v2[,"Experiment"]))
      
# annotate grouping uw line, timepoint, virus, tissue, experiment
qpcr_data_v2$Group <- paste(qpcr_data_v2[,"UW_Line"],qpcr_data_v2[,"Timepoint"],qpcr_data_v2[,"Virus"],qpcr_data_v2[,"Tissue"],qpcr_data_v2[,"Experiment"], sep="_")

# check there are no duplicated rows
sum(duplicated(qpcr_data_v2)) == 0

# check there are no empy FC means
sum(is.na(qpcr_data_v2[,"fc.mean"])) == 0

# check columns: "N","dCt.mean","dCt.sd" from Mouse file

# calculate and check in By line: "baseline.dCt","ddCt.mean"    "ddCt.sd", "fc.mean", "fc.sd"

##########################################################

# read in by mouse file
read.xls(xls="/Users/choonoo/WNV/Gale_11_23_15/Gale_3032_qPCR_byMouse.xlsx", sheet=1) -> qpcr_data_mouse

names(summary(as.factor(qpcr_data_mouse[,"UW_Line"])))

names(summary(qpcr_data_mouse[,"Timepoint"]))

# ct.mean appears to be CT value
# check no CT.mean < 15
length(which(qpcr_data_mouse[,"Ct.mean"] < 15))

# make version 2 copy
qpcr_data_mouse_v2 = qpcr_data_mouse

# no mating
#names(summary(qpcr_data_mouse_v2[,"Mating"]))

# adding virus column and clean timepoint
qpcr_data_mouse_v2$Virus = NA

qpcr_data_mouse_v2[grep("M",qpcr_data_mouse_v2[,"Timepoint"]),"Virus"] <- "Mock"
qpcr_data_mouse_v2[-grep("M",qpcr_data_mouse_v2[,"Timepoint"]),"Virus"] <- "WNV"

qpcr_data_mouse_v2[,"Timepoint"] <- as.numeric(as.character(gsub("M","",qpcr_data_mouse_v2[,"Timepoint"])))

# more tissue types (cecum, large/small intestine, stomach)
names(summary(qpcr_data_mouse_v2[,"Tissue"]))

# viral titer only
names(summary(qpcr_data_mouse_v2[,"Experiment"]))

# annotate group
qpcr_data_mouse_v2$Group <- paste(qpcr_data_mouse_v2[,"UW_Line"],qpcr_data_mouse_v2[,"Timepoint"], qpcr_data_mouse_v2[,"Virus"],qpcr_data_mouse_v2[,"Tissue"],qpcr_data_mouse_v2[,"Experiment"], sep="_")

# add data alteration, notes and lab
qpcr_data_mouse_v2$Data_Altered = NA
qpcr_data_mouse_v2$Notes = NA
qpcr_data_mouse_v2$Lab = "Gale"

##########################################################

# check numbers in byline

# dct mean
dct_mean = aggregate(formula=qpcr_data_mouse_v2[,"dCt"]~qpcr_data_mouse_v2[,"Group"], data=qpcr_data_mouse_v2, FUN=mean)
names(dct_mean) <- c("Group","dCt.mean.V2")
dct_mean[order(dct_mean[,1]),] -> dct_mean_order

byline_dct_mean = unique(qpcr_data_v2[,c("Group","dCt.mean")])
byline_dct_mean[order(byline_dct_mean[,1]),] -> byline_dct_mean_order

check_dct_mean = sapply(1:dim(dct_mean_order)[1],function(x)isTRUE(all.equal(dct_mean_order[x,2],byline_dct_mean_order[x,2])))

# check if all dct mean correct
sum(check_dct_mean) == dim(dct_mean_order)[1]

# code to update if not all correct
#byline_dct_mean_order[!check_dct_mean,1] -> flag_dct_mean
#dct_mean_order[!check_dct_mean,2] -> correct_dct_mean

# update documentation
#qpcr_data_v2[as.vector(unlist(sapply(flag_dct_mean,function(x)which(x==qpcr_data_v2[,"Group"])))),"Data_Altered"] <- 'Yes'

#qpcr_data_v2[as.vector(unlist(sapply(flag_dct_mean,function(x)which(x==qpcr_data_v2[,"Group"])))),"Notes"] <- paste0(qpcr_data_v2[as.vector(unlist(sapply(flag_dct_mean,function(x)which(x==qpcr_data_v2[,"Group"])))),"Notes"], ", corrected dCt.mean")

# overwrite with correct values
#qpcr_data_v2[as.vector(unlist(sapply(flag_dct_mean,function(x)which(x==qpcr_data_v2[,"Group"])))),"dCt.mean"] <- correct_dct_mean

# check n
data.frame(summary(as.factor(qpcr_data_mouse_v2[,"Group"]),maxsum=8000)) -> bymouse_n
names(bymouse_n) <- c("N")

bymouse_n[,2] <- row.names(bymouse_n)
names(bymouse_n)[2] <- "Group"

bymouse_n[order(bymouse_n[,"Group"]),] -> bymouse_n_v2

qpcr_data_v2[,c("Group","N")] -> byline_n
byline_n[order(byline_n[,"Group"]),] -> byline_n_v2

# check if all n correct
sum(bymouse_n_v2[,1] == byline_n_v2[,2]) == dim(bymouse_n_v2)[1]

# code to correct
#bymouse_n_v2[-which(bymouse_n_v2[,1] == byline_n_v2[,2]),] -> correct_n

# document
#qpcr_data_v2[as.vector(unlist(sapply(correct_n[,2],function(x)which(x==qpcr_data_v2[,"Group"])))),"Data_Altered"] <- 'Yes'

#qpcr_data_v2[as.vector(unlist(sapply(correct_n[,2],function(x)which(x==qpcr_data_v2[,"Group"])))),"Notes"] <- paste0(qpcr_data_v2[as.vector(unlist(sapply(correct_n[,2],function(x)which(x==qpcr_data_v2[,"Group"])))),"Notes"], ", corrected N")

# correct values
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

# code to update if not all correct
# document
#qpcr_data_v2[as.vector(unlist(sapply(flag_dct_sd,function(x)which(x==qpcr_data_v2[,"Group"])))),"Data_Altered"] <- 'Yes'

#qpcr_data_v2[as.vector(unlist(sapply(flag_dct_sd,function(x)which(x==qpcr_data_v2[,"Group"])))),"Notes"] <- paste0(qpcr_data_v2[as.vector(unlist(sapply(flag_dct_sd,function(x)which(x==qpcr_data_v2[,"Group"])))),"Notes"], ", corrected dCt.sd")

#byline_dct_sd_order[!check_dct_sd,1] -> flag_dct_sd
#dct_sd_order[!check_dct_sd,2] -> correct_dct_sd

#qpcr_data_v2[as.vector(unlist(sapply(flag_dct_sd,function(x)which(x==qpcr_data_v2[,"Group"])))),"dCt.sd"] <- correct_dct_sd

# check baseline. dct
qpcr_data_v2$Group_g <- paste(qpcr_data_v2[,"UW_Line"],qpcr_data_v2[,"Tissue"], qpcr_data_v2[,"Experiment"],sep="_")

qpcr_data_v2$baseline.dCt.V2 = NA

# calculate baseline, 3032 uses baseline D4
baseline = 4
for(i in unique(qpcr_data_v2[,"Group_g"])){
  print(i)
  qpcr_data_v2[which(qpcr_data_v2[,"Group_g"] == i),"baseline.dCt.V2"] <- qpcr_data_v2[which(qpcr_data_v2[,"Group_g"] == i & qpcr_data_v2[,"Virus"] == "Mock" & qpcr_data_v2[,"Timepoint"] == baseline),"dCt.mean"]
}

# check if all baseline correct
sum(qpcr_data_v2[,"baseline.dCt"] == qpcr_data_v2[,"baseline.dCt.V2"]) == dim(qpcr_data_v2)[1]

# code to correct
#qpcr_data_v2[which(qpcr_data_v2[,"baseline.dCt"] != qpcr_data_v2[,"baseline.dCt.V2"]),c("Data_Altered")] <- 'Yes'

#qpcr_data_v2[which(qpcr_data_v2[,"baseline.dCt"] != qpcr_data_v2[,"baseline.dCt.V2"]),c("Notes")] <- paste0(qpcr_data_v2[which(qpcr_data_v2[,"baseline.dCt"] != qpcr_data_v2[,"baseline.dCt.V2"]),c("Notes")], ", corrected baseline.dCt")

# correct data
#qpcr_data_v2[which(qpcr_data_v2[,"baseline.dCt"] != qpcr_data_v2[,"baseline.dCt.V2"]),"baseline.dCt"] <- qpcr_data_v2[which(qpcr_data_v2[,"baseline.dCt"] != qpcr_data_v2[,"baseline.dCt.V2"]),"baseline.dCt.V2"]

# check ddct mean
qpcr_data_v2$ddCt.mean.V2 <- as.numeric(as.character(qpcr_data_v2[,"dCt.mean"])) - as.numeric(as.character(qpcr_data_v2[,"baseline.dCt"]))

check_ddct_mean = sapply(1:dim(qpcr_data_v2)[1],function(x)isTRUE(all.equal(qpcr_data_v2[x,"ddCt.mean"], qpcr_data_v2[x,"ddCt.mean.V2"])))

# all ddct mean correct
sum(check_ddct_mean) == dim(qpcr_data_v2)[1]

# code to correct
#qpcr_data_v2[!check_ddct_mean,"Data_Altered"] <- 'Yes'

#qpcr_data_v2[!check_ddct_mean,"Notes"] <- "corrected ddCt.mean"

# overwrite correct values
#qpcr_data_v2[!check_ddct_mean,"ddCt.mean"] <- qpcr_data_v2[!check_ddct_mean,"ddCt.mean.V2"]

# take out V2s if they were calculated correct, ddct.sd, fc.sd
names(qpcr_data_v2)

remove_cols = c("ddCt.sd","fc.sd","baseline.dCt.V2","ddCt.mean.V2")

qpcr_data_v2[,!names(qpcr_data_v2) %in% remove_cols] -> qpcr_data_v3

# compute baseline sd, use baseline timepoint 4 saved above
qpcr_data_v3$baseline.dCt.sd = NA

for(i in unique(qpcr_data_v3[,"Group_g"])){
  print(i)
  qpcr_data_v3[which(qpcr_data_v3[,"Group_g"] == i),"baseline.dCt.sd"] <- qpcr_data_v3[which(qpcr_data_v3[,"Group_g"] == i & qpcr_data_v3[,"Virus"] == "Mock" & qpcr_data_v3[,"Timepoint"] == baseline),"dCt.sd"]
}

# compute ddCt.se, baseline 4 saved above
ddCt.se = sapply(1:dim(qpcr_data_v3)[1],function(x){
  sqrt(
    (qpcr_data_v3[x,"dCt.sd"]^2/qpcr_data_v3[x,"N"])+                                                                  (qpcr_data_v3[x,"baseline.dCt.sd"]^2/qpcr_data_v3[which(qpcr_data_v3[,"Group_g"] == qpcr_data_v3[x,"Group_g"] & qpcr_data_v3[,"Virus"] == "Mock" & qpcr_data_v3[,"Timepoint"] == baseline),"N"])
  )
})

qpcr_data_v3$ddCt.se <- ddCt.se

# Compute fc
qpcr_data_v3$fc.mean.V2 <- 2^-qpcr_data_v3[,"ddCt.mean"]

# check if fc calculated correct
check_fc_mean = sapply(1:dim(qpcr_data_v3)[1],function(x)isTRUE(all.equal(qpcr_data_v3[x,"fc.mean"], qpcr_data_v3[x,"fc.mean.V2"])))

# all fc mean correct
sum(check_fc_mean) == dim(qpcr_data_v3)[1]

# code to correct
#qpcr_data_v3[!check_fc_mean,"Data_Altered"] <- 'Yes'

#qpcr_data_v3[!check_fc_mean,"Notes"] <- paste0(qpcr_data_v3[!check_fc_mean,"Notes"], ", corrected fc.mean")

#qpcr_data_v3[!check_fc_mean,"fc.mean"] <- qpcr_data_v3[!check_fc_mean,"fc.mean.V2"]

# remove fc v2 if calculated correct
qpcr_data_v3[,!names(qpcr_data_v3) %in% 'fc.mean.V2'] -> qpcr_data_v4

# order by group
qpcr_data_v4[order(qpcr_data_v4[,"Group"]),] -> qpcr_data_v4_final

# double check calculations

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

# baseline sd, computed using corrected dct sd's at timepoint 4

# ddct se, computed using correct dct sd and baseline sd

# remove group g
qpcr_data_v4_final[,!names(qpcr_data_v4_final) %in% 'Group_g'] -> qpcr_data_final_format

# add lab
qpcr_data_final_format$Lab = "Gale"

# read in data dictionary

# read first column in of by line data dictionary
mycols <- rep("NULL", 4)
mycols[1] <- NA

read.xls(xls="/Users/choonoo/WNV/WNV_Data_Dictionary.xlsx", sheet="qPCR Data - By Line", colClasses=mycols) -> data_dict

# remove mating
order_cols = data_dict[,1]
order_cols[-which(order_cols=="Mating")] -> order_cols_v2

qpcr_data_final_format[,as.vector(unlist(sapply(order_cols_v2,function(x)which(x==names(qpcr_data_final_format)))))] -> qpcr_data_final_format_v2

# order bymouse

# read first column of by mouse in data dictionary
mycols <- rep("NULL", 4)
mycols[1] <- NA

read.xls(xls="/Users/choonoo/WNV/WNV_Data_Dictionary.xlsx", sheet="qPCR Data - By Mouse", colClasses=mycols) -> data_dict_mouse

# remove mating
order_cols_m = data_dict_mouse[,1]
order_cols_m[-which(order_cols_m=="Mating")] -> order_cols_m_v2

# change ct.mean name and dct linear
names(qpcr_data_mouse_v2)[which(names(qpcr_data_mouse_v2) == "Ct.mean")] <- "Ct"

names(qpcr_data_mouse_v2)[which(names(qpcr_data_mouse_v2) == "dCt_linear")] <- "dCt.linear"

qpcr_data_mouse_v2[,as.vector(unlist(sapply(order_cols_m_v2,function(x)which(x==names(qpcr_data_mouse_v2)))))] -> qpcr_data_mouse_v3

# check both byline and bymouse have the same UW lines
names(summary(as.factor(qpcr_data_mouse_v3[,"UW_Line"]))) == names(summary(as.factor(qpcr_data_final_format_v2[,"UW_Line"])))

# save processed data to files
write.table(file="./Gale_11_23_15/processed_qpcr/Gale_3032_qPCR_byLine_GC.txt", x=qpcr_data_final_format_v2, sep="\t", quote=F, row.names=F)

write.table(file="./Gale_11_23_15/processed_qpcr/Gale_3032_qPCR_byMouse_GC.txt", x=qpcr_data_mouse_v3, sep="\t", quote=F, row.names=F)

