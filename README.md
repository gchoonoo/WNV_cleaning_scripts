# WNV SIG Data Release: March 23, 2016

### Overall Data Status:

http://church.ohsu.edu:3838/choonoo/WNV_Data_Matrix_App/

--------------------------------------------------


## Histology Data Cleaning Steps: 

1. slide_label column used to create Virus, Timepoint, UW_Line columns
2. Mating and RIX_ID columns created and populated based on matching UWID in weight data file
3. Column names are standardized (see WNV_Data_Dictionary.xlsx spreadsheet)
4. Unique IDs are created by appending the Mating and RIX_ID
5. Duplicates identified
6. ND values set to missing values; Notes updated to reflect ND values
7. Missing values set to 0 if subtotals are present
8. Subtotals and Total score recalculated
9. Tissue column added
10. GI_Lesions column added for non-brain samples (0/1 values)

### Specific Histology Data Corrections / Modifications

1. Duplicates removed (line 55, D4 animals)
2. For line 36, slide_labels were corrected for D12 mocks.
3. For line 30, slide_labels were corrected (2.3M changed to 3.2M; and 3.12 changed to 2.12 for RIX_ID 74)
4. For line 58, slide label 3.12 changed to 3.12M for RIX_ID 156
5. For line 39, remove duplicate animals 143 (1.12M) and 145 (2.12M).
6. For line 39, animal 142 (3.12) scores set to NA. Needs to be re-scored, since scores from duplicate samples do not agree.
7. For line 38, scores for animals 101 (2.7) and 102 (3.7) set to NA. Needs to be re-scored, since scores from duplicate samples do not agree.
8. Removed line 51
9. For additional line 9 animals, IDs changed from 1,2,3 to 4,5,6.

### Cleaned Data Files

1. Gale_Histology_21-Mar-2016_final.xlsx

### Histology Plots

http://church.ohsu.edu:3838/mooneymi/wnv_histology_lineplots/

--------------------------------------------------

## Weight Data Cleaning Steps:

1. Clean mating: "X" changed to "x"
2. Clean ID with new mating: Mating + RIX_ID
3. Add/remove columns and edit column names: Add Lab, Data_Altered
	Gale: Add Death_FoundInCage column
			Change column name "Died_of_virus_at_day" to "Died_of_Virus" 
			Change column name "Animal_was_euthansized_at_day" to "Death_Euthanized"
			Change column name "Died_from_anesthesia" to "Died_from_Anesthesia"
			Remove columns: CV1, CV2, Interesting_Phenotype (empty)
	Lund: Add UWID, Cohort, Sex, Died_of_Virus, Died_from_Anesthesia, and Notes columns
4. Check the Virus column is annotated and clean timepoint
5. Check for any duplicated or NA rows
6. Annotate Death_date based on date of infection to timepoint or died early notes, annotate putative death day
	Lund: Change Death_Euthanized and Death_FoundInCage to day instead of date
7. Calculate weight change percentages, add D0_Percentage column and set to baseline 0
8. Add Flag_Weight_Drop column (True if animals lost 20% body weight)
9. Add Flag_Identical_Weights (True if identical weights on consecutive measurements)
10. Add Flag_Per_Day_Weight_Change (True if weight change > 10% on consecutive measurements)
11. Add Flag_weight_day where Timepoint is >= 3 more than putative death day (day the weights are recorded up to) (internal check)
12. Add Flag_weight_date where death day > timepoint, (internal check)
13. Order columns according to the data dictionary
14. Annotate array_exists column
	Lund: Remove beginning 0 in UW_Line
15. Update Data: Add flags_checked column
16. Make any specific alterations listed below

### Specific Weight Data Corrections / Modifications

Gale Lab

1. For line 18, UWID 3.13M changed to 3.12M (RIX_ID 62)
2. For line 95, animal 23, D8 weight changed from 27.74 to 24.74
3. For line 100, changed D2 RIX_IDs from 1 to 7, 2 to 15, and 3 to 20
4. For line 24 D2 animals, weights were switched between mocks and infected animals; notes added
5. For line 102, D12 animals swapped infection status
6. For line 82, clean mating "" to "16557x3154"
7. Remove line 51

Lund Lab

1. For line 63, typos fixed for animal 73
2. For line 55, typos fixed for animals 122 and 133
3. For line 8, d7 animal RIX_ID changed from 174 to 171
4. For line 67, death dates changed for animals 43, 44, and 45
5. For line 38, typo fixed for animal 7
6. For line 54, typos fixed for animals 54 and 58
7. For line 46, typo fixed for animal 43
8. For line 29, death dates changed for animals 163, 164, and 165
9. For line 57, mating changed from 3415x6012 to 3415x16012 to be consistent with flow data
10. For line 18, mating changed for D21 and D28 animals from 8042x16513 to 18042x16513 to be consistent with flow data
11. For line 67, death dates updated for animals 43, 44, and 45
12. For line 16, animal 291 d8 weight changed to 27.46
13. Line 61, D28 mock, set Notes to "Not receiving, line complete
14. Line 89 D12, set Notes to "Jumpy, no weight"
15. Line 90, animal 76 set Death_Euthanized to 11 and updated Death_Date
16. Line 100, animal 82 and 83 Death_Euthanized to 10 and 9 and updated Death_Date
17. Line 102, animal 62 and 64 changed D7 weight to 35.24 and 35.30 and updated weight percentages

### Cleaned Data Files

1. Gale_Weight_22-Mar-2016_final.xlsx
2. Lund_Weight_22-Mar-2016_final.xlsx

### Weight Plots

http://church.ohsu.edu:3838/mooneymi/wnv_pheno_time_series/

--------------------------------------------------


## Clinical Scores Data Cleaning Steps:

1. Clean mating: "X" changed to "x"
2. Clean ID
3. Add/remove columns and edit column names: Add Lab, Data_Altered, Notes
	Gale: Add Virus
	Lund: Add UWID, remove Sex 
4. Annotate died early columns: Add Death_Date, Death_Euthanized, Death_FoundInCage, Died_of_Virus, and Died_from_Anesthesia columns from weight data
5. Add Flag_weight_day where Timepoint is >= 3 more than putative death day (day the weights are recorded up to), (internal check)
6. Add Flag_weight_date where death day > 28, (internal check)
7. Alter data clinical score = 0 to NA past death day
8. Order columns according to the data dictionary
9. "array_exists" column added to indicate availability of gene expression data
10. Update Data: Add flags_checked column
11. Make any specific alterations listed below

### Specific Clinical Score Data Corrections / Modifications

Gale Lab

1. Remove line 51
2. For line 100, changed D2 RIX_IDs from 1 to 7, 2 to 15, and 3 to 20
3. For line 102, D12 animals swapped infection status
4. For line 82, clean mating "" to "16557x3154"

Lund Lab

1. Remove line 45, timepoint 21
2. For line 57, mating changed from 3415x6012 to 3415x16012 to be consistent with flow data
3. For line 18, mating changed for D21 and D28 animals from 8042x16513 to 18042x16513 to be consistent with flow data
4. For line 72, death annotations added for animals 54 and 57, since there is no weight data for this line yet


### Cleaned Data Files

1. Gale_Scores_22-Mar-2016_final.xlsx
2. Lund_Scores_22_Mar_2016_final.xlsx

### Clinical Scores Plots

http://church.ohsu.edu:3838/mooneymi/wnv_pheno_time_series/

--------------------------------------------------


## qPCR Data Cleaning Steps:

ByMouse Data File

1. Remove 'M' from Time points and convert to numeric values
2. Change "Condition" column name to "Virus"
3. Add Group: UW Line, Timepoint, Virus, Tissue, Experiment separated by "_"
4. Add Data_Altered, Notes, and Lab columns

ByLine Data File

1. Add Data_Altered and Notes columns
2. Add Virus column from time points
3. Remove 'M' from Time points and convert to numeric values.
4. Add Group column: UW Line, Timepoint, Virus, Tissue, Experiment separated by "_"
5. Add Lab column
6. Add baseline.dCt.sd column (this is the dCt.sd value for the baseline animal)
7. Add ddCt.se = sqrt(dCt.sd^2/N1 + baseline.dCt.sd^2/N2)
8. Remove ddCt.sd and fc.sd

### Specific qPCR Data Corrections / Modifications

ByMouse Data File

1. Remove samples with Ct < 15 (1 row, line 11)

ByLine Data File

1. UW Line 18 mating 18042x16513 changed to 8042x16513 to stay consistent with ByMouse file
2. Remove duplicated rows (these were line 18 after correction in #1)
3. Corrected N, dCt.mean, dCt.sd, baseline.dCt, ddCt.mean, and fc.mean for line 18

ByLine Data File (11/23/15)

1. UW Line 50 mating 3393x8052 to 16680x8016 to stay consistent with weight file
2. Remove duplicated rows (these were line 50 after correction above)

3032 ByMouse Data File

1. Note: Ct.mean column name changed to Ct, no ID, mating or RIX ID

8002 ByMouse Data File

1. Note: Ct.mean column name changed to Ct, no ID, mating or RIX ID

3032 ByLine Data File

1. Note: no mating, baseline is D4 mock

8002 ByLine Data File

1. Note: no mating, baseline is D7 mock

### Cleaned Data Files

1. Gale_qPCR_byMouse_23-Mar-2016_final.xlsx
2. Gale_qPCR_byLine_23-Mar-2016_final.xlsx
3. Gale_qPCR_dCt_Heritability_23-Mar-2016_final.xlsx

### qPCR Plots

http://church.ohsu.edu:3838/mooneymi/wnv_qpcr_barplots/

http://church.ohsu.edu:3838/mooneymi/wnv_qpcr_heatmaps/

--------------------------------------------------


## Flow Cytometry Data Cleaning Steps: 

1. Non-numeric characters (special characters; e.g. ¥) where numeric values are expected are set to NA (missing)
2. Column names are standardized (see WNV_Data_Dictionary.xlsx spreadsheet)
3. Matings, RIX_IDs, and time points are compared to a master list (from weight data). Mating formatting is fixed (no spaces; lowercase x)
4. Unique IDs are recreated by appending the Mating and RIX_ID
5. A Virus column is created from time points, and time points are then converted to numeric values (removing 'd' and 'm')
6. Add Data_Altered and Notes columns

### Specific Flow Data Corrections / Modifications (other than formatting mentioned above):

1. For files uploaded on June 17th, 2015, the Time variable on the Treg and T-cell panels was assumed to be 100. ***Note: as "Expt" files have been uploaded, I've noticed that the Time=100 assumption is not always correct. It may be best to upload all the older data as "Expt" files.
2. Data in files uploaded on June 17th, 2015 is overwritten by data in Experiment files
3. In Expt 7 and 8, animals labelled "mock" were removed
4. For animals 7, 8, and 32 from UW line 18, the time points were changed to d28 mock
5. For animals 4 and 6 from UW line 18 (18042x16513), the time points were changed to d28
6. For Expt 60, the mating for line 42 was corrected (8008x8016)
7. For Expt 64, the mating for line 68 was corrected (4410x3460) for d21 animals
8. For Expt 12, 13, and 14, the mating for line 12 was corrected (3252x8042)
9. For Expt 70, 71, and 72, the mating for line 71 was corrected (3564x8027)
10. For Expt 37, matings for line 18 were corrected (18042x16513)
11. For Expt 86, matings for line 50 were corrected (16680x8016)
12. For Expt 88, timepoints for line 86 were corrected (animals 49, 54, 58 changed to d12)

### Cleaned Data Files

1. Lund_Flow_21-Mar-2016_final.xlsx
    - This file contains only the cell population percentages in the original files
2. Lund_Flow_Full_21-Mar-2016_final.xlsx
    - This file contains the percentages, the calculated cell counts ('count' appended to variable name), and ratios for the ICS panel ('ratio' for ratio of percentages and 'ratio_count' for ratio of counts)
3. Lund_Flow_Heritability_21-Mar-2016_final.xlsx 
    - Heritability estimates for all flow variables, calculated using all animals and mocks only

### Flow Cytometry Plots

http://church.ohsu.edu:3838/mooneymi/wnv_flow_heatmaps/

http://church.ohsu.edu:3838/mooneymi/wnv_flow_boxplots/

http://church.ohsu.edu:3838/mooneymi/wnv_flow_lineplots/

--------------------------------------------------