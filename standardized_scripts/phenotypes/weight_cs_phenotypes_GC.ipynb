##########################################################

# Phenotypes (Version 12/7/15)

##########################################################

# setwd("/Users/choonoo/WNV")
# save.image("./pheno_12_7_15.RData")
# load("./pheno_12_7_15.RData")

library(gdata)

library(R.utils)

##########################################################

#setwd('Documents/MyDocuments/SystemsImmunogenetics/WNV')

# read in dec 7 release data
read.xls(xls="./8-Dec-2015/Lund_Weight_7-Dec-2015_final.xlsx", sheet=1) -> lund_weight

read.xls(xls="./8-Dec-2015/Lund_Scores_7-Dec-2015_final.xlsx", sheet=1) -> lund_score

read.xls(xls="./8-Dec-2015/Gale_Weight_7-Dec-2015_final.xlsx", sheet=1) -> gale_weight

read.xls(xls="./8-Dec-2015/Gale_Scores_7-Dec-2015_final.xlsx", sheet=1) -> gale_score

## Merge Lund Weight and CS Data

## Fix CS column names
days = c("D0", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12", "D13", "D14", "D15", "D16", "D17", "D18", "D19", "D20", "D21", "D22", "D23", "D24", "D25", "D26", "D27", "D28")
colnames(lund_score)[colnames(lund_score) %in% days] = paste(colnames(lund_score)[colnames(lund_score) %in% days], "cs", sep="_")

cs_days = paste(days,"cs",sep="_")

lund_score[,c("ID","Mating","RIX_ID","UW_Line","Virus","Date_Infected","Death_Date","Timepoint","Lab",cs_days)] -> lund_score_subset

lund_weight[,c("ID","Mating","RIX_ID","UW_Line","Virus","Date_Infected","Death_Date","Timepoint","Lab",days,paste0(days,"_Percentage"))] -> lund_weight_subset

merge_cols = intersect(names(lund_weight_subset), names(lund_score_subset))

all_lund_weight_cs = merge(lund_weight_subset, lund_score_subset, by=merge_cols, all=T)

## Merge Gale Weight and CS Data

## Fix CS column names
colnames(gale_score)[colnames(gale_score) %in% days] = paste(colnames(gale_score)[colnames(gale_score) %in% days], "cs", sep="_")

gale_score[,c("ID","Mating","RIX_ID","UW_Line","Virus","Date_Infected","Death_Date","Timepoint","Lab",cs_days)] -> gale_score_subset

gale_weight[,c("ID","Mating","RIX_ID","UW_Line","Virus","Date_Infected","Death_Date","Timepoint","Lab",days,paste0(days,"_Percentage"))] -> gale_weight_subset

all_gale_weight_cs = merge(gale_weight_subset, gale_score_subset, by=merge_cols, all=T)

### Merge both labs
all_weight_cs = rbind(all_lund_weight_cs, all_gale_weight_cs)

#### Calculate Phenotypes

## Calculate Max Clinical Scores
cs_days7 = c("D0_cs","D1_cs","D2_cs", "D3_cs", "D4_cs", "D5_cs", "D6_cs","D7_cs")
cs_days8_14 = c("D8_cs", "D9_cs", "D10_cs", "D11_cs","D12_cs","D13_cs","D14_cs")
cs_days15_28 = c("D15_cs", "D16_cs", "D17_cs", "D18_cs", "D19_cs", "D20_cs","D21_cs", "D22_cs", "D23_cs", "D24_cs", "D25_cs", "D26_cs", "D27_cs", "D28_cs")

all_weight_cs$max_cs = apply(all_weight_cs[,cs_days], 1, max, na.rm=T)
all_weight_cs$max_cs_d7 = apply(all_weight_cs[,cs_days7], 1, max, na.rm=T)
all_weight_cs$max_cs_d8_d14 = apply(all_weight_cs[,cs_days8_14], 1, max, na.rm=T)
all_weight_cs$max_cs_d15_d28 = apply(all_weight_cs[,cs_days15_28], 1, max, na.rm=T)

all_weight_cs$max_cs[is.infinite(all_weight_cs$max_cs)] = NA
all_weight_cs$max_cs_d7[is.infinite(all_weight_cs$max_cs_d7)] = NA
all_weight_cs$max_cs_d8_d14[is.infinite(all_weight_cs$max_cs_d8_d14)] = NA
all_weight_cs$max_cs_d15_d28[is.infinite(all_weight_cs$max_cs_d15_d28)] = NA

per_day_weight_change = function(w, start_day=1) {
	#print(w)
	max_day = max(which(!is.na(w)))
	day_span = max_day+start_day-1
	if (!is.infinite(max_day) & day_span > 1) {
		return(w[max_day]/(day_span-1))
	} else {
		return(NA)
	}
}

days7 = c("D0", "D1", "D2", "D3", "D4", "D5", "D6", "D7")
days8_14 = c("D8", "D9", "D10", "D11", "D12", "D13", "D14")
days15_28 = c("D15", "D16", "D17", "D18", "D19", "D20", "D21", "D22", "D23", "D24", "D25", "D26", "D27", "D28")

days_percent = paste0(days[-1],"_Percentage")
days7_percent = paste0(days7[-1],"_Percentage")
days8_14_percent = paste0(days8_14,"_Percentage")
days15_28_percent = paste0(days15_28,"_Percentage")

all_weight_cs$avg_w_change = apply(all_weight_cs[,days_percent], 1, per_day_weight_change)
all_weight_cs$avg_w_change_d7 = apply(all_weight_cs[,days7_percent], 1, per_day_weight_change)
all_weight_cs$avg_w_change_d8_d14 = apply(all_weight_cs[,days8_14_percent], 1, per_day_weight_change, start_day=8)
all_weight_cs$avg_w_change_d15_d28 = apply(all_weight_cs[,days15_28_percent], 1, per_day_weight_change, start_day=15)

all_weight_cs$max_w_change = apply(all_weight_cs[,days_percent], 1, min, na.rm=T)
all_weight_cs$max_w_change_d7 = apply(all_weight_cs[,days7_percent], 1, min, na.rm=T)
all_weight_cs$max_w_change_d8_d14 = apply(all_weight_cs[,days8_14_percent], 1, min, na.rm=T)
all_weight_cs$max_w_change_d15_d28 = apply(all_weight_cs[,days15_28_percent], 1, min, na.rm=T)

all_weight_cs$max_w_change[is.infinite(all_weight_cs$max_w_change)] = NA
all_weight_cs$max_w_change_d7[is.infinite(all_weight_cs$max_w_change_d7)] = NA
all_weight_cs$max_w_change_d8_d14[is.infinite(all_weight_cs$max_w_change_d8_d14)] = NA
all_weight_cs$max_w_change_d15_d28[is.infinite(all_weight_cs$max_w_change_d15_d28)] = NA

## Aggregate average weight loss by Line and Virus
## Before D8
all_weight_cs$d7_w_change_gt5percent_mock = NA
agg_w_change = aggregate(all_weight_cs$avg_w_change_d7, by=list(all_weight_cs$UW_Line, all_weight_cs$Virus), mean, na.rm=T)
colnames(agg_w_change) = c('line', 'virus', 'w_change')

lines = unique(agg_w_change$line)
for (l in lines) {
	mock_change = agg_w_change[agg_w_change$line==l & agg_w_change$virus=='Mock', 'w_change']
	wnv_change = agg_w_change[agg_w_change$line==l & agg_w_change$virus=='WNV', 'w_change']
	if (length(mock_change) > 0 & length(wnv_change) > 0) { 
		if (!is.na(mock_change) & !is.na(wnv_change)) {
			delta = wnv_change-mock_change
			print(l)
			print(delta)
			if (delta < -5) {
				all_weight_cs$d7_w_change_gt5percent_mock[all_weight_cs$UW_Line==l] = TRUE
			} else {
				all_weight_cs$d7_w_change_gt5percent_mock[all_weight_cs$UW_Line==l] = FALSE
			}
		}
	}
}

## D8 - D14
all_weight_cs$d8_d14_w_change_gt5percent_mock = NA
agg_w_change = aggregate(all_weight_cs$avg_w_change_d8_d14, by=list(all_weight_cs$UW_Line, all_weight_cs$Virus), mean, na.rm=T)
colnames(agg_w_change) = c('line', 'virus', 'w_change')

lines = unique(agg_w_change$line)
for (l in lines) {
	mock_change = agg_w_change[agg_w_change$line==l & agg_w_change$virus=='Mock', 'w_change']
	wnv_change = agg_w_change[agg_w_change$line==l & agg_w_change$virus=='WNV', 'w_change']
	if (length(mock_change) > 0 & length(wnv_change) > 0) { 
		if (!is.na(mock_change) & !is.na(wnv_change)) {
			delta = wnv_change-mock_change
			print(l)
			print(delta)
			if (delta < -5) {
				all_weight_cs$d8_d14_w_change_gt5percent_mock[all_weight_cs$UW_Line==l] = TRUE
			} else {
				all_weight_cs$d8_d14_w_change_gt5percent_mock[all_weight_cs$UW_Line==l] = FALSE
			}
		}
	}
}

## After D14
all_weight_cs$d15_d28_w_change_gt5percent_mock = NA
agg_w_change = aggregate(all_weight_cs$avg_w_change_d15_d28, by=list(all_weight_cs$UW_Line, all_weight_cs$Virus), mean, na.rm=T)
colnames(agg_w_change) = c('line', 'virus', 'w_change')

lines = unique(agg_w_change$line)
for (l in lines) {
	mock_change = agg_w_change[agg_w_change$line==l & agg_w_change$virus=='Mock', 'w_change']
	wnv_change = agg_w_change[agg_w_change$line==l & agg_w_change$virus=='WNV', 'w_change']
	if (length(mock_change) > 0 & length(wnv_change) > 0) { 
		if (!is.na(mock_change) & !is.na(wnv_change)) {
			delta = wnv_change-mock_change
			print(l)
			print(delta)
			if (delta < -5) {
				all_weight_cs$d15_d28_w_change_gt5percent_mock[all_weight_cs$UW_Line==l] = TRUE
			} else {
				all_weight_cs$d15_d28_w_change_gt5percent_mock[all_weight_cs$UW_Line==l] = FALSE
			}
		}
	}
}


## Aggregate max weight loss by Line and Virus
## Before D8
all_weight_cs$d7_max_w_change_gt5percent_mock = NA
agg_w_change = aggregate(all_weight_cs$max_w_change_d7, by=list(all_weight_cs$UW_Line, all_weight_cs$Virus), mean, na.rm=T)
colnames(agg_w_change) = c('line', 'virus', 'w_change')

lines = unique(agg_w_change$line)
for (l in lines) {
	mock_change = agg_w_change[agg_w_change$line==l & agg_w_change$virus=='Mock', 'w_change']
	wnv_change = agg_w_change[agg_w_change$line==l & agg_w_change$virus=='WNV', 'w_change']
	if (length(mock_change) > 0 & length(wnv_change) > 0) { 
		if (!is.na(mock_change) & !is.na(wnv_change)) {
			delta = wnv_change-mock_change
			print(l)
			print(delta)
			if (delta < -5) {
				all_weight_cs$d7_max_w_change_gt5percent_mock[all_weight_cs$UW_Line==l] = TRUE
			} else {
				all_weight_cs$d7_max_w_change_gt5percent_mock[all_weight_cs$UW_Line==l] = FALSE
			}
		}
	}
}

## D8 - D14
all_weight_cs$d8_d14_max_w_change_gt5percent_mock = NA
agg_w_change = aggregate(all_weight_cs$max_w_change_d8_d14, by=list(all_weight_cs$UW_Line, all_weight_cs$Virus), mean, na.rm=T)
colnames(agg_w_change) = c('line', 'virus', 'w_change')

lines = unique(agg_w_change$line)
for (l in lines) {
	mock_change = agg_w_change[agg_w_change$line==l & agg_w_change$virus=='Mock', 'w_change']
	wnv_change = agg_w_change[agg_w_change$line==l & agg_w_change$virus=='WNV', 'w_change']
	if (length(mock_change) > 0 & length(wnv_change) > 0) { 
		if (!is.na(mock_change) & !is.na(wnv_change)) {
			delta = wnv_change-mock_change
			print(l)
			print(delta)
			if (delta < -5) {
				all_weight_cs$d8_d14_max_w_change_gt5percent_mock[all_weight_cs$UW_Line==l] = TRUE
			} else {
				all_weight_cs$d8_d14_max_w_change_gt5percent_mock[all_weight_cs$UW_Line==l] = FALSE
			}
		}
	}
}

## After D14
all_weight_cs$d15_d28_max_w_change_gt5percent_mock = NA
agg_w_change = aggregate(all_weight_cs$max_w_change_d15_d28, by=list(all_weight_cs$UW_Line, all_weight_cs$Virus), mean, na.rm=T)
colnames(agg_w_change) = c('line', 'virus', 'w_change')

lines = unique(agg_w_change$line)
for (l in lines) {
	mock_change = agg_w_change[agg_w_change$line==l & agg_w_change$virus=='Mock', 'w_change']
	wnv_change = agg_w_change[agg_w_change$line==l & agg_w_change$virus=='WNV', 'w_change']
	if (length(mock_change) > 0 & length(wnv_change) > 0) { 
		if (!is.na(mock_change) & !is.na(wnv_change)) {
			delta = wnv_change-mock_change
			print(l)
			print(delta)
			if (delta < -5) {
				all_weight_cs$d15_d28_max_w_change_gt5percent_mock[all_weight_cs$UW_Line==l] = TRUE
			} else {
				all_weight_cs$d15_d28_max_w_change_gt5percent_mock[all_weight_cs$UW_Line==l] = FALSE
			}
		}
	}
}


## Aggregate max CS by Line
## Before D8
all_weight_cs$d7_max_cs_gt1 = NA
agg_max_cs = aggregate(all_weight_cs$max_cs_d7[all_weight_cs$Virus=='WNV'], by=list(all_weight_cs$UW_Line[all_weight_cs$Virus=='WNV']), max, na.rm=T)
colnames(agg_max_cs) = c('line', 'max_cs')

lines = unique(agg_max_cs$line)
for (l in lines) {
	max_cs = agg_max_cs[agg_max_cs$line==l, 'max_cs']
	if (length(max_cs) > 0) { 
		if (!is.na(max_cs)) {
			if (max_cs > 1) {
				all_weight_cs$d7_max_cs_gt1[all_weight_cs$UW_Line==l] = TRUE
			} else {
				all_weight_cs$d7_max_cs_gt1[all_weight_cs$UW_Line==l] = FALSE
			}
		}
	}
}

## Between D8 and D14
all_weight_cs$d8_d14_max_cs_gt1 = NA
agg_max_cs = aggregate(all_weight_cs$max_cs_d8_d14[all_weight_cs$Virus=='WNV'], by=list(all_weight_cs$UW_Line[all_weight_cs$Virus=='WNV']), max, na.rm=T)
colnames(agg_max_cs) = c('line', 'max_cs')

lines = unique(agg_max_cs$line)
for (l in lines) {
	max_cs = agg_max_cs[agg_max_cs$line==l, 'max_cs']
	if (length(max_cs) > 0) { 
		if (!is.na(max_cs)) {
			if (max_cs > 1) {
				all_weight_cs$d8_d14_max_cs_gt1[all_weight_cs$UW_Line==l] = TRUE
			} else {
				all_weight_cs$d8_d14_max_cs_gt1[all_weight_cs$UW_Line==l] = FALSE
			}
		}
	}
}

## After D14
all_weight_cs$d15_d28_max_cs_gt1 = NA
agg_max_cs = aggregate(all_weight_cs$max_cs_d15_d28[all_weight_cs$Virus=='WNV'], by=list(all_weight_cs$UW_Line[all_weight_cs$Virus=='WNV']), max, na.rm=T)
colnames(agg_max_cs) = c('line', 'max_cs')

lines = unique(agg_max_cs$line)
for (l in lines) {
	max_cs = agg_max_cs[agg_max_cs$line==l, 'max_cs']
	if (length(max_cs) > 0) { 
		if (!is.na(max_cs)) {
			if (max_cs > 1) {
				all_weight_cs$d15_d28_max_cs_gt1[all_weight_cs$UW_Line==l] = TRUE
			} else {
				all_weight_cs$d15_d28_max_cs_gt1[all_weight_cs$UW_Line==l] = FALSE
			}
		}
	}
}


all_weight_cs$phenotype = NA
all_weight_cs$phenotype[is.na(all_weight_cs$d7_max_cs_gt1) & is.na(all_weight_cs$d8_d14_max_cs_gt1) & is.na(all_weight_cs$d15_d28_max_cs_gt1) & is.na(all_weight_cs$d7_max_w_change_gt5percent_mock) & is.na(all_weight_cs$d8_d14_max_w_change_gt5percent_mock) & is.na(all_weight_cs$d15_d28_max_w_change_gt5percent_mock)] = 'unknown'
all_weight_cs$phenotype[any(!all_weight_cs$d7_max_cs_gt1 & !all_weight_cs$d8_d14_max_cs_gt1 & !all_weight_cs$d15_d28_max_cs_gt1 & !all_weight_cs$d7_max_w_change_gt5percent_mock & !all_weight_cs$d8_d14_max_w_change_gt5percent_mock & !all_weight_cs$d15_d28_max_w_change_gt5percent_mock) & !any(all_weight_cs$d7_max_cs_gt1 & all_weight_cs$d8_d14_max_cs_gt1 & all_weight_cs$d15_d28_max_cs_gt1 & all_weight_cs$d7_max_w_change_gt5percent_mock & all_weight_cs$d8_d14_max_w_change_gt5percent_mock & all_weight_cs$d15_d28_max_w_change_gt5percent_mock)] = 'no disease'
all_weight_cs$phenotype[all_weight_cs$d7_max_cs_gt1 | all_weight_cs$d7_max_w_change_gt5percent_mock] = 'early susceptibility'
all_weight_cs$phenotype[(!all_weight_cs$d7_max_cs_gt1 | is.na(all_weight_cs$d7_max_cs_gt1)) & (!all_weight_cs$d7_max_w_change_gt5percent_mock | is.na(all_weight_cs$d7_max_w_change_gt5percent_mock)) & (all_weight_cs$d8_d14_max_cs_gt1 | all_weight_cs$d8_d14_max_w_change_gt5percent_mock) & (!all_weight_cs$d15_d28_max_cs_gt1 | is.na(all_weight_cs$d15_d28_max_cs_gt1)) & (!all_weight_cs$d15_d28_max_w_change_gt5percent_mock | is.na(all_weight_cs$d15_d28_max_w_change_gt5percent_mock))] = 'normal disease'
all_weight_cs$phenotype[(all_weight_cs$d15_d28_max_cs_gt1 | all_weight_cs$d15_d28_max_w_change_gt5percent_mock) & (!all_weight_cs$d7_max_cs_gt1 | is.na(all_weight_cs$d7_max_cs_gt1)) & (!all_weight_cs$d7_max_w_change_gt5percent_mock | is.na(all_weight_cs$d7_max_w_change_gt5percent_mock))] = 'extended disease'

table(all_weight_cs$phenotype)

# Change line 4 back to no disease, score is outlier
all_weight_cs[which(all_weight_cs[,"UW_Line"] == '4'),"phenotype"] <- "no disease"

pheno_list = list("no disease" = sort(unique(all_weight_cs$UW_Line[all_weight_cs$phenotype=='no disease'])),'early susceptibility'=sort(unique(all_weight_cs$UW_Line[all_weight_cs$phenotype=='early susceptibility'])),'normal disease'=sort(unique(all_weight_cs$UW_Line[all_weight_cs$phenotype=='normal disease'])),'extended disease'=sort(unique(all_weight_cs$UW_Line[all_weight_cs$phenotype=='extended disease'])))

save(pheno_list, file="./Phenotypes/phenotype_list.RData")

write.table(file="./Phenotypes/phenotype_data.txt",x=all_weight_cs,quote=F,row.names=F,sep="\t")


# update phenotypes in data matrix
read.xls(xls="./Data_Matrices/Lund_Gale_Data_Matrix_12_7_15.xlsx", sheet=1,stringsAsFactors=F) -> data_matrix

pheno_data = unique(all_weight_cs[,c("UW_Line","phenotype")])

pheno_data[order(pheno_data[,1]),] -> pheno_data_v2

merge(data_matrix, pheno_data_v2,by="UW_Line",all.x=T) -> data_matrix_v2

# Change line 51 to unknown, dropped this line
data_matrix_v2[which(is.na(data_matrix_v2[,"phenotype"])),"phenotype"] <- "Unknown"

# change line 51 to dropped
data_matrix_v2[which((data_matrix_v2[,"UW_Line"] == 51)),"Overall.Progress...Weight.CS"] <- "Dropped"


data_matrix_v2[which(toupper(data_matrix_v2[,"Phenotype"]) != toupper(data_matrix_v2[,"phenotype"])),"Phenotype_Note"] <- c("Changed from normal to extended","Changed from normal to extended", "Dropped this line","Changed from no disease to extended","Changed from normal to early susceptibility","Changed from unknown to extended")

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

update = data_matrix_v2[which(toupper(data_matrix_v2[,"Phenotype"]) != toupper(data_matrix_v2[,"phenotype"])),"phenotype"]

data_matrix_v2[which(toupper(data_matrix_v2[,"Phenotype"]) != toupper(data_matrix_v2[,"phenotype"])),"Phenotype"] <- as.vector(sapply(update, simpleCap))

data_matrix_v2[,!names(data_matrix_v2)%in%"phenotype"] -> data_matrix_v3

write.table(file="./Data_Matrices/Lund_Gale_Data_Matrix_12_9_15.txt",x=data_matrix_v3,quote=F,row.names=F,sep="\t")
