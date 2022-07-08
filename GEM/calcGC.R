#Calculate GC lambda (output in print) and save new file with added adjusted P values for marginal, 1df, 2df.


library(tidyverse)

dir1<-"/scratch/mf91122/LipidsxVeg/GEM"

phe<-c("LDL", "HDL", "Tot_Chol", "TAGs")
exp<-c("Veg1", "Veg2")
mod<-c("M1", "M2")

for (p in 1:length(phe)){
for (e in 1:length(exp)){
for (m in 1:length(mod)){

#p=1;e=1;m=1
chr<-list()
for (i in 1:22){
chr[[i]]<-read.table(paste(dir1, "/", mod[m], "/", phe[p], "/", 
        exp[e], "/chr", i, sep=""),
                header=T)
}

stat<-do.call(rbind,chr)
stat<-as_tibble(stat)

#1df
chisq1df <- qchisq(1-stat$P_Value_Interaction,1)
lambda1df<-median(chisq1df)/qchisq(0.5,1)
newchisq1df<-chisq1df/lambda1df
stat$adjP_Value_Interaction<-pchisq(newchisq1df, df=1, lower.tail=FALSE)

#2df
chisq2df <- qchisq(1-stat$P_Value_Joint,1)
lambda2df<-median(chisq2df)/qchisq(0.5,1) 
newchisq2df<-chisq2df/lambda2df
stat$adjP_Value_Joint<-pchisq(newchisq2df, df=1, lower.tail=FALSE)

#Marginal
chisqMar <- qchisq(1-stat$P_Value_Marginal,1)
lambdaMar<-median(chisqMar)/qchisq(0.5,1)
newchisqMar<-chisqMar/lambdaMar
stat$adjP_Value_Marginal<-pchisq(newchisqMar, df=1, lower.tail=FALSE)

print(paste(phe[p], exp[e], mod[m], 
	lambda1df, lambda2df, lambdaMar))

write_tsv(stat,paste(dir1, "/adj/",mod[m], "-", phe[p],
                        "-", exp[e],"-adj.tsv", sep=""),
                col_names=TRUE)

}}}

#Process for FUMA
x<-read.delim("/scratch/mf91122/LipidsxVeg/GEM/adj/chr2-M2-TAGs-Veg1-adj.tsv", sep="\t", header=F)
colnames(x)<-c("SNPID", "RSID", "CHR", "POS", "Non_Effect_Allele", "Effect_Allele", "N_Samples", "AF", "N_Veg1_1", "AF_Veg1_1", "N_Veg1_0", "AF_Veg1_0", "Beta_Marginal", "SE_Beta_Marginal", "Beta_G.Veg1", "SE_Beta_G.Veg1", "P_Value_Marginal", "P_Value_Interaction", "P_Value_Joint", "adjP_Value_Interaction", "adjP_Value_Joint", "adjP_Value_Marginal")
write.table(x, "/scratch/mf91122/LipidsxVeg/GEM/adj/chr2-M2-TAGs-Veg1-adj.txt", fileEncoding="ascii", quote=F, row.names=F)


