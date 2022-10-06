#1)Generate a R function to calculate FDR by using a data frame of all p-values of probe sets as parameter, and return a new data frame for all FDR of corresponding probeset. You CANNOT use any existed FDR-calculation functions to do the analysis below.
#2)Write R code to read in the p-value list file, process the p-values with your own FDR calculation function, and output the FDR for all probe sets.
#3)Paste your R codes and output FDRs into the report file to submit to assignment.

# p-values for FDR in R
# using Benjamini-Hochberg procedure (BH)
# BH is one of the most popular correction of multiple testing for biomedical data

# 11098 Sample size (N) with p-value

#FORMULA

#FDR = p-value * N/Rank@i


#STEPS
# 1. Sort p-value in ascending order
# 2. Assign rank 1 to N(11098)
# 3. apply formular

##############################################################################
#1)Generate a R function to calculate FDR by using a data frame of all p-values of probe sets as parameter,
#and return a new data frame for all FDR of corresponding probeset. You CANNOT use any existed FDR-calculation 
#functions to do the analysis below.
##################################################################################
df1 <- read.table("Ibrahim_Odunayo_Salaudeen.pValue.tsv", header=T, sep="\t")
as.data.frame.matrix(df1)
head(df1)
tail(df1)
df2 <- df1[ , c('ProbeSet', 'p.value')]
df2[0]
#data frame with 0 columns and 11098 rows

attach(df2)

#order or sort dataframe by p.value

df3 <- df2[order(p.value), ]
df3
df3["Rank"] <- seq(1:11098)
df3


#####################################################################################
#2)Write R code to read in the p-value list file, process the p-values with your own 
#FDR calculation function, and output the FDR for all probe sets.

#for loop to calculate the False Discovery rate using the FDR Formula
########################################################################################

library(data.table)
setDT(df3)[, FDR := p.value* 11098/Rank, by = Rank]
df3


#################################################################################
#3)Paste your R codes and output FDRs into the report file to submit to assignment.
#################################################################################
write.table(df3, file="FDR-result.tsv", sep="\t", quote=F)
