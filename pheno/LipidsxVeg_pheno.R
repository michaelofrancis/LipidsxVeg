suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(RNOmni))

setwd("/work/kylab/mike/LipidsxVeg/pheno")

source("manyColsToDummy.R")

withdrawn<-read.csv("w48818_20210809.csv", header = FALSE)

QCids<-read.table("/scratch/mf91122/LipidsxVeg/pheno/bd_QC-keep.txt",header=TRUE)

###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Load data=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#Load UK Biobank datasets-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
source('/work/kylab/mike/PUFA-GWAS/pheno/load_UKBphenotables.R') #20 min

#Phenotypes  ------------------------------------------------------------------------------
#Covariates 
#Model 1: sex, age, age squared, genotyping array, and assessment center indicators (sites of recruitment); 
#lipid medication, socioeconomic status measured by Townsend deprivation index;  


pheno<-bd%>%select(f.eid, f.21003.0.0, f.31.0.0, 
			f.189.0.0, f.30800.0.0,
                   f.30890.0.0, f.30600.0.0, 
                   f.54.0.0, f.22000.0.0,
		f.30690.0.0, f.30780.0.0, f.30760.0.0, f.30870.0.0,
                f.30850.0.0, f.30830.0.0    
		)

colnames(pheno)<-c("IID", "Age", "Sex",  
			"Townsend", "Oestradiol",
                   "VitaminD_Blood", "Albumin", 
                   "Assessment_center", "Geno_batch",
			"Tot_Chol", "LDL", "HDL", "TAGs",
			"Test", "SHBG"
                    )

pheno2<-bd_join4%>%select(f.eid, f.74.0.0, f.30050.0.0,
                          f.21001.0.0, f.20116.0.0, f.20160.0.0,
                          f.6177.0.0,f.6153.0.0,
				f.3166.0.0,
                          f.2724.0.0
                          )
colnames(pheno2)<-c("IID","Fasting_time", "Mean_corpuscular_haemoglobin",
                    "BMI", "SmokeStatus","Ever_smoked",                    
                    "lipid_med", "lipid_med_plushormones",
			"blood_draw_time",
                    "Menopause"
                    )

new<-left_join(pheno, pheno2, by="IID")
new<-as_tibble(new)

#Remove withdrawn participants------------------------------------
new<-new[!(new$IID %in% withdrawn$V1), ]

#QC participants via output of UKB_participantQC.R----------------

new<-new[(new$IID %in% QCids$IID),]

#Age squared----------------------------
new$Age2<-new$Age^2

#Make dummy 0/1 cols for each assessment center----------------------
#table(pheno$Assessment_center)
centers<-unique(new$Assessment_center)
centercols<-paste("center", 1:22, sep="")
new[centercols]<-0

for (i in 1:length(centers)){
    new[new$Assessment_center==centers[i],][centercols[i]]<-1
}

new<-new%>%select(-Assessment_center)
new

#Genotype batch
new$Geno_batch1<-0
new$Geno_batch1[new$Geno_batch>0]<-1
#sum(pheno$Geno_batch1) #[1] 438313
new$Geno_batch<-new$Geno_batch1
new<-new%>%select(-Geno_batch1)
#table(new$Geno_batch) #it worked


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#Switch sex values to numeric
new$Sex<-mapvalues(as.character(new$Sex), 
                     c("Male", "Female"), c(0,1))
new$Sex<-as.numeric(new$Sex)



#Statins
statincols<-c(sprintf("f.20003.0.%s", 0:47))
statincodes<-c(1141146234,1141192414,1140910632,1140888594,1140864592,
	1141146138,1140861970,1140888648,1141192410,
	1141188146,1140861958,1140881748,1141200040)

manyColsToDummy(statincodes, bd_join4[,statincols], "statinoutput")
statinoutput$statins<-rowSums(statinoutput) 
statinoutput$statins[statinoutput$statins>1]<-1

statinoutput$IID<-bd_join4$f.eid

statinoutput<-statinoutput%>%select(IID, statins)

new<-left_join(new, statinoutput, by="IID")


#ADD VEGETARIAN COLUMNS

VEG<-read.table("/scratch/mf91122/LipidsxVeg/pheno/Vegetarian/vegQC2_04032021.txt", 
		header=T, stringsAsFactors=F)

new<-inner_join(new, VEG, by="IID")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#EXTRA PROCESSING FOR T COLUMNS---------------
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# SHBG Transformations---------------------------------------------
new$SHBG2<-new$SHBG^2




###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###WRITE OUTPUT=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

outdir="/scratch/mf91122/LipidsxVeg/pheno"
new$FID<-new$IID

new<-new%>%select(FID, everything())


participants<-new%>%select(FID,IID)

#Participant list
write.table(participants, 
	paste(outdir, "/LipidsxVeg_phenoQC_IDS.txt",sep=""), 
	row.names=FALSE, quote=FALSE)

write.table(new, 
	paste(outdir, "/LipidsxVeg_pheno.txt", sep=""),
	row.names=FALSE, quote=FALSE)

write.csv(new,
        paste(outdir, "/LipidsxVeg_pheno.csv", sep=""),
        row.names=FALSE, quote=FALSE)


paste(colnames(new[3:50]), " ", 
	sapply(new[3:50], mean, na.rm=T), "(", 
	sapply(new[3:50], sd, na.rm=T), ")", sep="")
