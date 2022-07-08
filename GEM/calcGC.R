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
