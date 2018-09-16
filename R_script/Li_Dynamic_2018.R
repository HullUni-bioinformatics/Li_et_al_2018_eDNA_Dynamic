####START####


####Set working directory under your home directory ('R'is my home directory)
setwd("~/R")

###Dynamic experiment 2017###

###packages require in this script###
##The package(reshape) should after package(dplyr), otherwise the rename funtion in package(reshape) will not work


#plyr
library(plyr)
#libraries pakage,including ggplot2
library(tidyverse)
#cowplot
library(cowplot)
#scales
library(scales)
#reshape
library(reshape)
library(reshape2)
#gtable
library(gtable)
#data.table
library(dtplyr)
#gridExtra
library(gridExtra)
#geom_text_repel
library(ggrepel)
#Chi square indepandent
library(rcompanion)
#Dunn test for multiple comparisons after Kruskal-Wallis
library(FSA)
#PCA plotting
library(factoextra)
#Species Prediction And Diversity Estimation
library(SpadeR)

library(car)

library(vegan)

library(ade4)
#Mantel test

###metaBEAT command###

#metaBEAT.py \
#-Q Querymap.txt -R REFmap.txt --cluster --clust_match 1 --clust_cov 3 --blast --min_ident 1 \
#-m 12S -n 5 -E -v -o Dyn_12S-trim_30-merged-nonchimeras-cluster_1c3-blast-min_ident_1.0 &> log

#Avoid false positive approaches#


#1.Low-frequency noise threshold (0.002)



Dynamic <- read.csv(file = "Dynamic_2018/12S_Dyn_May_2018.csv")


for (i in 1:nrow(Dynamic)){Dynamic$Leuciscus_leuciscus[i]=Dynamic$Leuciscus_leuciscus[i]+Dynamic$Leuciscus_idus[i]; 
Dynamic$Leuciscus_idus[i]<-0}
Dynamic$Leuciscus_idus <- NULL

for (i in 1:nrow(Dynamic)){Dynamic$Carassius_carassius[i]=Dynamic$Carassius_auratus[i]+Dynamic$Carassius_carassius[i]; 
Dynamic$Carassius_auratus[i]<-0}
Dynamic$Carassius_auratus <- NULL


#Summary the the total read counts for each sample
Dynamic$SUM <- rowSums(Dynamic[7:ncol(Dynamic)])

#Reshape the dataset


#Using the 1:6, and 20 variables as ID

Dynamic_rs <- melt(Dynamic,id=c(1:6,ncol(Dynamic)))

#Rename the variables

Dynamic_rs <- rename (Dynamic_rs,c("variable"="Species","value"="Reads"))



#Calculate each species reads percentage for each sample

Dynamic_rs$Ratio <- Dynamic_rs$Reads/Dynamic_rs$SUM

Threshold <- 0.002 #based on the positive reads in samples percentage

for (a in 1:nrow(Dynamic_rs)) {
  if (!is.na(Dynamic_rs$Ratio[a])) {
    if(Dynamic_rs$Ratio[a] < Threshold){Dynamic_rs$Reads[a] <- 0} 
  }
}



###Control samples###
#Check the Control samples, then remove all of control samples, positive species and unassigned species

Dynamic_rs_control <- Dynamic_rs[which(Dynamic_rs$Position == 'Ctrl'
                                      &Dynamic_rs$Replicate != 'POS1'
                                      &Dynamic_rs$Replicate != 'POS2'),]


Dynamic_rs_control_filter <- Dynamic_rs_control [which(Dynamic_rs_control $Species != 'Astatotilapia_calliptera' 
                                                     & Dynamic_rs_control$Species != 'unassigned'),]


ggplot(Dynamic_rs_control_filter, aes(x = Species, y= Reads, fill = Replicate))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Pond,ncol =3, scales = "free")+
  labs(x="Species", y="Read counts", title= "Control samples_12S_T0.002")+theme_bw()+
  theme(text=element_text(size=14),axis.text.x = element_text(angle =45, vjust = 1, hjust = 0.95),
        plot.title = element_text(size = rel(0.8),face = "bold", hjust = 0.5))

ggsave("Control samples_12S_T0.002.jpeg",path = "Dynamic_2018/Figures/", width = 10, height = 8, units = "in", dpi=500)






###Samples###

levels(Dynamic_rs$Species)

Dynamic_rs_filter <- Dynamic_rs [which(Dynamic_rs$Reads >0),]



Dynamic_rs_sample <- Dynamic_rs[which(Dynamic_rs$Position != 'Ctrl' 
                                      & Dynamic_rs$Species != 'Astatotilapia_calliptera' 
                                      & Dynamic_rs$Species != 'Alburnus_alburnus'
                                      & Dynamic_rs$Species != 'Blicca_bjoerkna'
                                      & Dynamic_rs$Species != 'Hypophthalmichthys_molitrix'
                                      & Dynamic_rs$Species != 'unassigned'),]


Dynamic_rs_sample <- ddply(Dynamic_rs_sample, 'SampleID', mutate, Percent_reads = Reads/sum(Reads))


#rename the species
Dynamic_rs_sample$Species <-mapvalues(Dynamic_rs_sample$Species, c("Abramis_brama","Barbus_barbus","Carassius_carassius","Squalius_cephalus","Leuciscus_leuciscus",
                                                                 "Rutilus_rutilus","Scardinius_erythrophthalmus","Tinca_tinca"), 
                                                                 c("BRE","BAR","CAR","CHU","DAC","ROA","RUD","TEN"))



####E1####


Dynamic_rs_sample_E1 <- Dynamic_rs_sample[which(Dynamic_rs_sample$Pond == 'E1'),]


Dynamic_rs_sample_E1$Day <- factor(Dynamic_rs_sample_E1$Day, levels = (c("D0","D2","D4","D6","D8","D10","D12","D14")))




Dynamic_rs_sample_E1$Species <- factor(Dynamic_rs_sample_E1$Species, levels = (c("RUD","CHU","BAR","BRE","CAR","ROA","TEN","DAC")))


 # ggplot(Dynamic_rs_sample_E1,aes(x=Replicate,y=Percent_reads,fill=Species))+
 #  geom_bar(stat="identity",position="stack",width = 0.8)+facet_wrap(~Position+Day,ncol=8,labeller=labeller(.multi_line = FALSE))+
 #  scale_fill_manual(values=c("firebrick3","green4","blue","orange2","tan4","hotpink","gray10","purple3"))+
 #  scale_y_continuous(labels=percent)+
 #  labs(x="Replicate", y= "Average read counts percentage")+theme_bw()+
 #  theme(text=element_text(size=13))

 
 


names(Dynamic_rs_sample_E1)

Dynamic_E1_mean <- aggregate(Percent_reads~Pond+DP+Day+Position+Species, data = Dynamic_rs_sample_E1, FUN=mean)


Dynamic_E1_SD <- aggregate(Percent_reads~Pond+DP+Day+Position+Species, data = Dynamic_rs_sample_E1, FUN=sd)


Dynamic_E1_SD <- rename(Dynamic_E1_SD,c("Percent_reads"="SD"))

Dynamic_E1 <- merge(Dynamic_E1_mean, Dynamic_E1_SD, by=c("Pond","DP","Day","Position","Species"))


Dynamic_E1$Species <- factor(Dynamic_E1$Species, levels = (c("RUD","CHU","BAR","BRE","CAR","ROA","TEN","DAC")))



# ggplot(Dynamic_E1,aes(x=Day,y=Percent_reads,ymin=Percent_reads-SD,ymax=Percent_reads+SD,fill=Species))+
#  geom_bar(stat="identity",position="dodge")+facet_wrap(~Position,ncol=1)+
#  geom_errorbar(position = position_dodge(0.9), width = 0.5)+
#  scale_fill_manual(values=c("firebrick3","green4","blue","orange2","tan4","hotpink","gray10","purple3"))+
#  scale_y_continuous(labels=percent)+
#  labs(x="Sampling day", y= "Average read counts percentage")+theme_bw()+
#  theme(text=element_text(size=13))



Dynamic_E1$Day<- factor(Dynamic_E1$Day, levels = (c("D14","D12","D10","D8","D6","D4","D2","D0")))



Density <- read.csv(file = 'Dynamic_2018/Fish density.csv')


Density <- ddply(Density, .(Pond,Position), mutate, Percent_reads = Value/sum(Value))

Density$Value <- NULL

Density_E1 <- Density[which(Density$Pond == "E1"),]

Dynamic_Den_E1<-rbind(Dynamic_E1,Density_E1)


Dynamic_Den_E1$Species <- factor(Dynamic_Den_E1$Species, levels = (c("RUD","CHU","BAR","BRE","CAR","ROA","TEN","DAC")))


####E1_COM####

# ggplot(Dynamic_Den_E1,aes(x=Day,y=Percent_reads,fill=Species))+
#           geom_bar(stat="identity",position="stack",width = 0.8)+facet_grid(Position ~ ., switch = "both",scales="free_y",space = "free")+
#           coord_flip()+
#           scale_y_continuous(labels=percent,expand = c(0,0),limits = c(0, 1.04), breaks=seq(0, 1, 0.25))+
#           scale_fill_manual(values=c("firebrick3","green4","blue","orange2","tan4","purple3","gray10","hotpink"))+
#           labs(x="Sampling position", y= "Species composition")+
#           guides(fill = guide_legend(ncol = 8,reverse=TRUE))+
#           theme_bw()+theme(text=element_text(size=15),legend.position="bottom",legend.text = element_text(size=12),
#                            strip.background = element_blank(),panel.border = element_blank(),strip.text.y=element_text(angle=180),
#                            strip.placement= "outside",axis.ticks.length = unit(0.1, "cm"),
#                            legend.title = element_text(size=12),axis.text.y =element_text(size=10),
#                            panel.grid.minor=element_blank(),panel.grid.major =element_blank())






FIG2_E1<- ggplot(Dynamic_Den_E1,aes(x=Position,y=Percent_reads,fill=Species))+
          geom_bar(stat="identity",position="stack",width = 0.8)+facet_grid(Day ~ ., switch = "both",scales="free_y",space = "free")+
          coord_flip()+
          scale_y_continuous(labels=percent,expand = c(0,0),limits = c(0, 1.04), breaks=seq(0, 1, 0.25))+
          scale_fill_manual(values=c("firebrick3","green4","blue","orange2","tan4","hotpink","gray10","purple3"))+
          labs(x="Sampling day", y= "Species composition")+
          guides(fill = guide_legend(ncol = 8,reverse=TRUE))+
          theme_bw()+theme(text=element_text(size=15),legend.position="bottom",legend.text = element_text(size=9),
                           strip.background = element_blank(),panel.border = element_blank(),strip.text.y=element_text(angle=180),
                           strip.placement= "outside",axis.ticks.length = unit(0.1, "cm"),
                           legend.title = element_text(size=10),axis.text.y =element_text(size=8.5),axis.text.x =element_text(size=8.5),
                           panel.grid.minor=element_blank(),panel.grid.major =element_blank(),legend.key.size =  unit(0.2, "in"))



legend_FIG2_E1 <- get_legend(FIG2_E1)
 
 
FIG2_E1 <- FIG2_E1 + theme(legend.position="none")




Dyn_E1_read_mean <- aggregate(Reads~Pond+DP+Day+Position+Species, data = Dynamic_rs_sample_E1, FUN=mean)



BCD_E1 <- Dyn_E1_read_mean [c(2,5,6)]



BCD_E1_wide <- reshape(BCD_E1, idvar= c("Species"), timevar= "DP", direction = "wide")


BCD_E1_wide.rownames <- data.frame(BCD_E1_wide[,-1], row.names=BCD_E1_wide[,1])



##standardization to propostion
BCD_E1_wide.std <-decostand(BCD_E1_wide.rownames, "total",MARGIN=2)

 

BCD_E1_P1 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("P1", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_P1_data <- data.frame(Position=c("P1"),Bray=BCD_E1_P1)


BCD_E1_P2 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("P2", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_P2_data <- data.frame(Position=c("P2"),Bray=BCD_E1_P2)


BCD_E1_P3 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("P3", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_P3_data <- data.frame(Position=c("P3"),Bray=BCD_E1_P3)


BCD_E1_P4 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("P4", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_P4_data <- data.frame(Position=c("P4"),Bray=BCD_E1_P4)


BCD_E1_P5 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("P5", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_P5_data <- data.frame(Position=c("P5"),Bray=BCD_E1_P5)


BCD_E1_data <- rbind(BCD_E1_P1_data,BCD_E1_P2_data,BCD_E1_P3_data,BCD_E1_P4_data,BCD_E1_P5_data)

FIG2_BCD_pos_E1<- ggplot(BCD_E1_data,aes(x=Position,y=Bray))+
                  geom_boxplot(outlier.shape = NA)+
                  geom_jitter(size=2,width =0.2)+
                  scale_y_continuous(limits=c(0,1))+
                  labs(x="Sampling position",y= "")+
                  geom_smooth(method = "loess", se=TRUE, color="black", aes(group=1), size=1.0,lty="F1")+theme_bw()+
                  theme(text=element_text(size=15),panel.grid.minor=element_blank(),panel.grid.major =element_blank(),
                        legend.title = element_text(size=8))



# Shapiro-Wilk test of normality
shapiro.test(BCD_E1_data$Bray)
qqnorm(BCD_E1_data$Bray)


#Not normal distribution, use Kruskal-Wallis one-way analysis of variance

#kruskal.test for non-parametric test https://rcompanion.org/rcompanion/d_06.html
kruskal.test(Bray ~ Position, data = BCD_E1_data)



DT_BCD_E1 <- dunnTest(Bray ~ Position, data = BCD_E1_data, method="bh") 


DT_BCD_E1 <- DT_BCD_E1$res

DT_BCD_E1 


ST_BCD_E1 <-cldList(comparison = DT_BCD_E1$Comparison, p.value = DT_BCD_E1$P.adj, threshold= 0.05)

ST_BCD_E1$Pond <- "E1"


ST_BCD_E1$y <- 0.9




FIG2_BCD_pos_E1 <- FIG2_BCD_pos_E1+geom_text(data=ST_BCD_E1,aes(Group,y, label=Letter),hjust = 0.5, size=5, parse = TRUE)






BCD_E1_D0 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D0", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D0_data <- data.frame(Day=c("D0"),Bray=BCD_E1_D0)


BCD_E1_D2 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D2", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D2_data <- data.frame(Day=c("D2"),Bray=BCD_E1_D2)


BCD_E1_D4 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D4", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D4_data <- data.frame(Day=c("D4"),Bray=BCD_E1_D4)


BCD_E1_D6 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D6", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D6_data <- data.frame(Day=c("D6"),Bray=BCD_E1_D6)


BCD_E1_D8 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D8", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D8_data <- data.frame(Day=c("D8"),Bray=BCD_E1_D8)


BCD_E1_D10 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D10", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D10_data <- data.frame(Day=c("D10"),Bray=BCD_E1_D10)


BCD_E1_D12 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D12", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D12_data <- data.frame(Day=c("D12"),Bray=BCD_E1_D12)


BCD_E1_D14 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D14", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D14_data <- data.frame(Day=c("D14"),Bray=BCD_E1_D14)


BCD_E1_Day <- rbind(BCD_E1_D0_data,BCD_E1_D2_data,BCD_E1_D4_data,BCD_E1_D6_data,
                    BCD_E1_D8_data,BCD_E1_D10_data,BCD_E1_D12_data,BCD_E1_D14_data)


####E1_BCD####



FIG2_BCD_E1<-ggplot(BCD_E1_Day,aes(x=Day,y=Bray))+
            geom_boxplot(outlier.shape = NA)+
            geom_jitter(size=2,width =0.2)+
            scale_y_continuous(limits=c(0,1))+
            labs(x="Sampling day",y= "Bray-Curtis dissimilarity")+
            geom_smooth(method = "loess", se=TRUE, color="black", aes(group=1), size=1.0,lty="F1")+theme_bw()+
            theme(text=element_text(size=15),panel.grid.minor=element_blank(),panel.grid.major =element_blank(),
                  legend.title = element_text(size=8))




# Shapiro-Wilk test of normality
shapiro.test(BCD_E1_Day$Bray)
qqnorm(BCD_E1_Day$Bray)


#Not normal distribution, use Kruskal-Wallis one-way analysis of variance

#kruskal.test for non-parametric test https://rcompanion.org/rcompanion/d_06.html
kruskal.test(Bray ~ Day, data = BCD_E1_Day)



DT_BCD_E1 <- dunnTest(Bray ~ Day, data = BCD_E1_Day, method="bh") 


DT_BCD_E1 <- DT_BCD_E1$res

DT_BCD_E1 


ST_BCD_E1 <-cldList(comparison = DT_BCD_E1$Comparison, p.value = DT_BCD_E1$P.adj, threshold= 0.05)

ST_BCD_E1$Pond <- "E1"


ST_BCD_E1$y <- 0.9



ST_BCD_E1$Group <-mapvalues(ST_BCD_E1$Group, c("D","D1"), c("D0","D10"))

FIG2_BCD_E1 <- FIG2_BCD_E1+geom_text(data=ST_BCD_E1,aes(Group,y, label=Letter),hjust = 0.5, size=5, parse = TRUE)




######stage_E1#########
BCD_E1_D0 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D0", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D0_st <- data.frame(Day=c("Before introduction"),Bray=BCD_E1_D0)


BCD_E1_D2 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D2", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D2_st <- data.frame(Day=c("Introduction"),Bray=BCD_E1_D2)


BCD_E1_D4 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D4", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D4_st <- data.frame(Day=c("Introduction"),Bray=BCD_E1_D4)


BCD_E1_D6 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D6", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D6_st <- data.frame(Day=c("Introduction"),Bray=BCD_E1_D6)


BCD_E1_D8 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D8", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D8_st <- data.frame(Day=c("Introduction"),Bray=BCD_E1_D8)


BCD_E1_D10 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D10", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D10_st <- data.frame(Day=c("Removal"),Bray=BCD_E1_D10)


BCD_E1_D12 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D12", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D12_st <- data.frame(Day=c("Removal"),Bray=BCD_E1_D12)


BCD_E1_D14 <- as.vector(vegdist(t(BCD_E1_wide.std[,grepl("D14", colnames(BCD_E1_wide.std))]), method="bray"))
BCD_E1_D14_st <- data.frame(Day=c("Removal"),Bray=BCD_E1_D14)


BCD_E1_stage <- rbind(BCD_E1_D0_st,BCD_E1_D2_st,BCD_E1_D4_st,BCD_E1_D6_st,
                    BCD_E1_D8_st,BCD_E1_D10_st,BCD_E1_D12_st,BCD_E1_D14_st)




FIGS1_BCD_E1<-ggplot(BCD_E1_stage,aes(x=Day,y=Bray))+
              geom_boxplot(outlier.shape = NA)+
              geom_jitter(size=2,width =0.2)+
              scale_y_continuous(limits=c(0,1))+
              labs(x="Sampling stage",y= "Bray-Curtis dissimilarity")+theme_bw()+
              theme(text=element_text(size=15),panel.grid.minor=element_blank(),panel.grid.major =element_blank(),
                   axis.text.x =element_text(size=10),axis.text.y =element_text(size=10))

# Shapiro-Wilk test of normality
shapiro.test(BCD_E1_stage$Bray)
qqnorm(BCD_E1_stage$Bray)


#Not normal distribution, use Kruskal-Wallis one-way analysis of variance

#kruskal.test for non-parametric test https://rcompanion.org/rcompanion/d_06.html
kruskal.test(Bray ~ Day, data = BCD_E1_stage)



DT_BCD_E1 <- dunnTest(Bray ~ Day, data = BCD_E1_stage, method="bh") 


DT_BCD_E1 <- DT_BCD_E1$res

DT_BCD_E1 


ST_BCD_E1 <-cldList(comparison = DT_BCD_E1$Comparison, p.value = DT_BCD_E1$P.adj, threshold= 0.05)

ST_BCD_E1$Pond <- "E1"


ST_BCD_E1$y <- 0.9

ST_BCD_E1$Group <-mapvalues(ST_BCD_E1$Group, c("Beforeintroduction"), c("Before introduction"))




FIGS1_BCD_E1 <- FIGS1_BCD_E1+geom_text(data=ST_BCD_E1,aes(Group,y, label=Letter),hjust = 0.5, size=5, parse = TRUE)



###E4####



Dynamic_rs_sample_E4 <- Dynamic_rs_sample[which(Dynamic_rs_sample$Pond == 'E4'),]


Dynamic_rs_sample_E4$Day <- factor(Dynamic_rs_sample_E4$Day, levels = (c("D0","D2","D4","D6","D8","D10","D12","D14")))




Dynamic_rs_sample_E4$Species <- factor(Dynamic_rs_sample_E4$Species, levels = (c("RUD","DAC","CHU","BAR","BRE","CAR","ROA","TEN")))




# ggplot(Dynamic_rs_sample_E4,aes(x=Replicate,y=Percent_reads,fill=Species))+
#   geom_bar(stat="identity",position="stack",width = 0.8)+facet_wrap(~Position+Day,ncol=8,labeller=labeller(.multi_line = FALSE))+
#   scale_fill_manual(values=c("firebrick3","purple3","green4","blue","orange2","tan4","hotpink","gray10"))+
#   scale_y_continuous(labels=percent)+
#   labs(x="Replicate", y= "Average read counts percentage")+theme_bw()+
#   theme(text=element_text(size=13))
# 




names(Dynamic_rs_sample_E4)

Dynamic_E4_mean <- aggregate(Percent_reads~Pond+DP+Day+Position+Species, data = Dynamic_rs_sample_E4, FUN=mean)


Dynamic_E4_SD <- aggregate(Percent_reads~Pond+DP+Day+Position+Species, data = Dynamic_rs_sample_E4, FUN=sd)


Dynamic_E4_SD <- rename(Dynamic_E4_SD,c("Percent_reads"="SD"))

Dynamic_E4 <- merge(Dynamic_E4_mean, Dynamic_E4_SD, by=c("Pond","DP","Day","Position","Species"))


Dynamic_E4$Species <- factor(Dynamic_E4$Species, levels = (c("RUD","DAC","CHU","BAR","BRE","CAR","ROA","TEN")))






# 
# ggplot(Dynamic_E4,aes(x=Day,y=Percent_reads,ymin=Percent_reads-SD,ymax=Percent_reads+SD,fill=Species))+
#   geom_bar(stat="identity",position="dodge")+facet_wrap(~Position,ncol=1)+
#   geom_errorbar(position = position_dodge(0.9), width = 0.5)+
#   scale_fill_manual(values=c("firebrick3","purple3","green4","blue","orange2","tan4","hotpink","gray10"))+
#   scale_y_continuous(labels=percent)+
#   labs(x="Sampling day", y= "Average read counts percentage")+theme_bw()+
#   theme(text=element_text(size=13))




Dynamic_E4$Day<- factor(Dynamic_E4$Day, levels = (c("D14","D12","D10","D8","D6","D4","D2","D0")))


Density_E4 <- Density[which(Density$Pond == "E4"),]

Dynamic_Den_E4<-rbind(Dynamic_E4,Density_E4)

Dynamic_Den_E4$Species <- factor(Dynamic_Den_E4$Species, levels = (c("RUD","DAC","CHU","BAR","BRE","CAR","ROA","TEN")))

# ####E4_COM####
# 
# FIG2_E4<- ggplot(Dynamic_Den_E4,aes(x=Day,y=Percent_reads,fill=Species))+
#           geom_bar(stat="identity",position="stack",width = 0.8)+facet_grid(Position ~ ., switch = "both",scales="free_y",space = "free")+
#           coord_flip()+
#           scale_y_continuous(labels=percent,expand = c(0,0),limits = c(0, 1.04), breaks=seq(0, 1, 0.25))+
#           scale_fill_manual(values=c("firebrick3","hotpink","green4","blue","orange2","tan4","purple3","gray10"))+
#           labs(x="", y= "Species composition")+
#           theme_bw()+theme(text=element_text(size=15),legend.position="none",legend.text = element_text(size=12),
#                            strip.background = element_blank(),panel.border = element_blank(),strip.text.y=element_text(angle=180),
#                            strip.placement= "outside",axis.ticks.length = unit(0.1, "cm"),
#                            legend.title = element_text(size=12),axis.title.y=element_blank(),axis.text.y =element_text(size=10),
#                            panel.grid.minor=element_blank(),panel.grid.major =element_blank())



FIG2_E4<- ggplot(Dynamic_Den_E4,aes(x=Position,y=Percent_reads,fill=Species))+
          geom_bar(stat="identity",position="stack",width = 0.8)+facet_grid(Day~ ., switch = "both",scales="free_y",space = "free")+
          coord_flip()+
          scale_y_continuous(labels=percent,expand = c(0,0),limits = c(0, 1.04), breaks=seq(0, 1, 0.25))+
          scale_fill_manual(values=c("firebrick3","purple3","green4","blue","orange2","tan4","hotpink","gray10"))+
          labs(x="", y= "Species composition")+
          theme_bw()+theme(text=element_text(size=15),legend.position="none",legend.text = element_text(size=9),
                           strip.background = element_blank(),panel.border = element_blank(),strip.text.y=element_text(angle=180),
                           strip.placement= "outside",axis.ticks.length = unit(0.1, "cm"),
                           legend.title = element_text(size=10),axis.title.y=element_blank(),axis.text.y =element_text(size=8.5),axis.text.x =element_text(size=8.5),
                           panel.grid.minor=element_blank(),panel.grid.major =element_blank())





####Fig2_com####

Com <- grid.arrange(FIG2_E1,FIG2_E4, ncol = 2, nrow = 1,
                    widths = c(2.7,2.6), heights = c(2.5))



Com2 <- grid.arrange(Com, legend_FIG2_E1,nrow = 2,
                     heights = c(2.5,0.2))

Com2_final <- ggdraw(Com2)+draw_plot_label(c("(a)","(b)"), c(0, 0.51), c(1.0, 1.0), size = 10,fontface = "plain")




ggsave("Fig2_Com.jpeg",Com2_final,path = "Dynamic_2018/Figures/",  width = 9, height = 8, units = "in", dpi=500)

ggsave("Fig2_Com.pdf",Com2_final,path = "Dynamic_2018/Figures/",  width = 9, height = 8, units = "in", dpi=800)





Dyn_E4_read_mean <- aggregate(Reads~Pond+DP+Day+Position+Species, data = Dynamic_rs_sample_E4, FUN=mean)



BCD_E4 <- Dyn_E4_read_mean [c(2,5,6)]



BCD_E4_wide <- reshape(BCD_E4, idvar= c("Species"), timevar= "DP", direction = "wide")


BCD_E4_wide.rownames <- data.frame(BCD_E4_wide[,-1], row.names=BCD_E4_wide[,1])



##standardization to propostion
BCD_E4_wide.std <-decostand(BCD_E4_wide.rownames, "total",MARGIN=2)



BCD_E4_P1 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("P1", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_P1_data <- data.frame(Position=c("P1"),Bray=BCD_E4_P1)


BCD_E4_P2 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("P2", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_P2_data <- data.frame(Position=c("P2"),Bray=BCD_E4_P2)


BCD_E4_P3 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("P3", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_P3_data <- data.frame(Position=c("P3"),Bray=BCD_E4_P3)


BCD_E4_P4 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("P4", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_P4_data <- data.frame(Position=c("P4"),Bray=BCD_E4_P4)


BCD_E4_P5 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("P5", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_P5_data <- data.frame(Position=c("P5"),Bray=BCD_E4_P5)


BCD_E4_data <- rbind(BCD_E4_P1_data,BCD_E4_P2_data,BCD_E4_P3_data,BCD_E4_P4_data,BCD_E4_P5_data)

FIG2_BCD_pos_E4<- ggplot(BCD_E4_data,aes(x=Position,y=Bray))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=2,width =0.2)+
  scale_y_continuous(limits=c(0,1))+
  labs(x="Sampling position",y= "")+
  geom_smooth(method = "loess", se=TRUE, color="black", aes(group=1), size=1.0,lty="F1")+theme_bw()+
  theme(text=element_text(size=15),panel.grid.minor=element_blank(),panel.grid.major =element_blank(),
        legend.title = element_text(size=8))



# Shapiro-Wilk test of normality
shapiro.test(BCD_E4_data$Bray)
qqnorm(BCD_E4_data$Bray)


#Not normal distribution, use Kruskal-Wallis one-way analysis of variance

#kruskal.test for non-parametric test https://rcompanion.org/rcompanion/d_06.html
kruskal.test(Bray ~ Position, data = BCD_E4_data)



DT_BCD_E4 <- dunnTest(Bray ~ Position, data = BCD_E4_data, method="bh") 


DT_BCD_E4 <- DT_BCD_E4$res

DT_BCD_E4 


ST_BCD_E4 <-cldList(comparison = DT_BCD_E4$Comparison, p.value = DT_BCD_E4$P.adj, threshold= 0.05)

ST_BCD_E4$Pond <- "E4"


ST_BCD_E4$y <- 0.9




FIG2_BCD_pos_E4 <- FIG2_BCD_pos_E4+geom_text(data=ST_BCD_E4,aes(Group,y, label=Letter),hjust = 0.5, size=5, parse = TRUE)






BCD_E4_D0 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D0", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D0_data <- data.frame(Day=c("D0"),Bray=BCD_E4_D0)


BCD_E4_D2 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D2", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D2_data <- data.frame(Day=c("D2"),Bray=BCD_E4_D2)


BCD_E4_D4 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D4", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D4_data <- data.frame(Day=c("D4"),Bray=BCD_E4_D4)


BCD_E4_D6 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D6", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D6_data <- data.frame(Day=c("D6"),Bray=BCD_E4_D6)


BCD_E4_D8 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D8", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D8_data <- data.frame(Day=c("D8"),Bray=BCD_E4_D8)


BCD_E4_D10 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D10", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D10_data <- data.frame(Day=c("D10"),Bray=BCD_E4_D10)


BCD_E4_D12 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D12", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D12_data <- data.frame(Day=c("D12"),Bray=BCD_E4_D12)


BCD_E4_D14 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D14", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D14_data <- data.frame(Day=c("D14"),Bray=BCD_E4_D14)


BCD_E4_Day <- rbind(BCD_E4_D0_data,BCD_E4_D2_data,BCD_E4_D4_data,BCD_E4_D6_data,
                    BCD_E4_D8_data,BCD_E4_D10_data,BCD_E4_D12_data,BCD_E4_D14_data)


####E4_BCD####

# FIG2_BCD_E4<- ggplot(BCD_E4_data,aes(x=Position,y=Bray))+
#               geom_boxplot(outlier.shape = NA)+
#               geom_jitter(size=2,width =0.2)+
#               scale_y_continuous(limits=c(0,1))+
#               labs(x="Sampling position",y= "Bray-Curtis dissimilarity")+theme_bw()+
#               theme(text=element_text(size=15),panel.grid.minor=element_blank(),panel.grid.major =element_blank(),
#                     legend.title = element_text(size=12))

FIG2_BCD_E4<-ggplot(BCD_E4_Day,aes(x=Day,y=Bray))+
              geom_boxplot(outlier.shape = NA)+
              geom_jitter(size=2,width =0.2)+
              scale_y_continuous(limits=c(0,1))+
              labs(x="Sampling day",y= "Bray-Curtis dissimilarity")+
              geom_smooth(method = "loess", se=TRUE, color="black", aes(group=1), size=1.0,lty="F1")+theme_bw()+
              theme(text=element_text(size=15),panel.grid.minor=element_blank(),panel.grid.major =element_blank(),
                    legend.title = element_text(size=8))



# Shapiro-Wilk test of normality
shapiro.test(BCD_E4_Day$Bray)
qqnorm(BCD_E4_Day$Bray)


#Not normal distribution, use Kruskal-Wallis one-way analysis of variance

#kruskal.test for non-parametric test https://rcompanion.org/rcompanion/d_06.html
kruskal.test(Bray ~ Day, data = BCD_E4_Day)



DT_BCD_E4 <- dunnTest(Bray ~ Day, data = BCD_E4_Day, method="bh") 


DT_BCD_E4 <- DT_BCD_E4$res

DT_BCD_E4 


ST_BCD_E4 <-cldList(comparison = DT_BCD_E4$Comparison, p.value = DT_BCD_E4$P.adj, threshold= 0.05)

ST_BCD_E4$Pond <- "E4"


ST_BCD_E4$y <- 0.9


ST_BCD_E4$Group <-mapvalues(ST_BCD_E4$Group, c("D","D1"), c("D0","D10"))


FIG2_BCD_E4 <- FIG2_BCD_E4+geom_text(data=ST_BCD_E4,aes(Group,y, label=Letter),hjust = 0.5, size=5, parse = TRUE)




######stage_E4#########
BCD_E4_D0 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D0", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D0_st <- data.frame(Day=c("Before introduction"),Bray=BCD_E4_D0)


BCD_E4_D2 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D2", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D2_st <- data.frame(Day=c("Introduction"),Bray=BCD_E4_D2)


BCD_E4_D4 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D4", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D4_st <- data.frame(Day=c("Introduction"),Bray=BCD_E4_D4)


BCD_E4_D6 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D6", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D6_st <- data.frame(Day=c("Introduction"),Bray=BCD_E4_D6)


BCD_E4_D8 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D8", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D8_st <- data.frame(Day=c("Introduction"),Bray=BCD_E4_D8)


BCD_E4_D10 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D10", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D10_st <- data.frame(Day=c("Removal"),Bray=BCD_E4_D10)


BCD_E4_D12 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D12", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D12_st <- data.frame(Day=c("Removal"),Bray=BCD_E4_D12)


BCD_E4_D14 <- as.vector(vegdist(t(BCD_E4_wide.std[,grepl("D14", colnames(BCD_E4_wide.std))]), method="bray"))
BCD_E4_D14_st <- data.frame(Day=c("Removal"),Bray=BCD_E4_D14)


BCD_E4_stage <- rbind(BCD_E4_D0_st,BCD_E4_D2_st,BCD_E4_D4_st,BCD_E4_D6_st,
                      BCD_E4_D8_st,BCD_E4_D10_st,BCD_E4_D12_st,BCD_E4_D14_st)




FIGS1_BCD_E4<-ggplot(BCD_E4_stage,aes(x=Day,y=Bray))+
              geom_boxplot(outlier.shape = NA)+
              geom_jitter(size=2,width =0.2)+
              scale_y_continuous(limits=c(0,1))+
              labs(x="Sampling stage",y= "")+theme_bw()+
              theme(text=element_text(size=15),panel.grid.minor=element_blank(),panel.grid.major =element_blank(),
                    legend.title = element_text(size=8),axis.text.x =element_text(size=10),axis.text.y =element_text(size=10))

# Shapiro-Wilk test of normality
shapiro.test(BCD_E4_stage$Bray)
qqnorm(BCD_E4_stage$Bray)


#Not normal distribution, use Kruskal-Wallis one-way analysis of variance

#kruskal.test for non-parametric test https://rcompanion.org/rcompanion/d_06.html
kruskal.test(Bray ~ Day, data = BCD_E4_stage)



DT_BCD_E4 <- dunnTest(Bray ~ Day, data = BCD_E4_stage, method="bh") 


DT_BCD_E4 <- DT_BCD_E4$res

DT_BCD_E4 


ST_BCD_E4 <-cldList(comparison = DT_BCD_E4$Comparison, p.value = DT_BCD_E4$P.adj, threshold= 0.05)

ST_BCD_E4$Pond <- "E4"


ST_BCD_E4$y <- 0.9

ST_BCD_E4$Group <-mapvalues(ST_BCD_E4$Group, c("Beforeintroduction"), c("Before introduction"))




FIGS1_BCD_E4 <- FIGS1_BCD_E4+geom_text(data=ST_BCD_E4,aes(Group,y, label=Letter),hjust = 0.5, size=5, parse = TRUE)



# 
# Com <- grid.arrange(FIG2_E1,FIG2_E4, ncol = 2, nrow = 1,
#                           widths = c(2.7,2.6), heights = c(2.5))
# 
# 
# 
# Com2 <- grid.arrange(legend_FIG2_E1,Com, nrow = 2,
#                      heights = c(0.3,2.5))
# 
# 
# 
# 
# Bray <- grid.arrange(FIG2_BCD_E1,FIG2_BCD_E4, ncol = 2, nrow = 1,
#                      widths = c(2.7,2.7), heights = c(2.5))
# 
# 
# FIG2 <- grid.arrange(Com2,Bray, ncol = 1, nrow = 2,
#              widths = c(2.7), heights = c(3.0,2.3))
# 
# FIG2_final <- ggdraw(FIG2)+draw_plot_label(c("A", "C","B","D"), c(0, 0, 0.51,0.51), c(0.95, 0.46, 0.95, 0.46), size = 15)
# 
# 
# ggsave("Fig2.jpeg",FIG2_final,path = "Dynamic_2018/Figures/",  width = 11, height = 12, units = "in", dpi=800)







# 
# Bray <- grid.arrange(FIG2_BCD_E1,FIG2_BCD_E4, ncol = 1, nrow = 2,
#                      widths = c(2.7), heights = c(2.5,2.5))
# 
# Bray_final <- ggdraw(Bray)+draw_plot_label(c("A","B"), c(0, 0), c(1.0, 0.5), size = 12)
# 
# 
# ggsave("Bray.jpeg",Bray_final,path = "Dynamic_2018/Figures/",  width = 5, height = 6, units = "in", dpi=800)

Bray <- grid.arrange(FIG2_BCD_E1, FIG2_BCD_pos_E1, FIG2_BCD_E4, FIG2_BCD_pos_E4, ncol = 2, nrow = 2,
                     widths = c(2.7,2.5), heights = c(2.5,2.5))

Bray_final <- ggdraw(Bray)+draw_plot_label(c("(a1)",'(b1)',"(a2)",'(b2)'), c(0, 0,0.51,0.51), c(1.0, 0.52), size = 12,fontface = "plain")


###Fig5_Bray#####

ggsave("Fig5_Bray.jpeg",Bray_final,path = "Dynamic_2018/Figures/",  width = 10, height = 8, units = "in", dpi=500)


ggsave("Fig5_Bray.pdf",Bray_final,path = "Dynamic_2018/Figures/",  width = 10, height = 8, units = "in", dpi=800)




####FigS3_Bray_stage####

Bray_stage <- grid.arrange(FIGS1_BCD_E1,FIGS1_BCD_E4, ncol = 2, nrow = 1,
                     widths = c(2.7,2.7), heights = c(2.5))

Bray_stage_final <- ggdraw(Bray_stage)+draw_plot_label(c("(a)","(b)"), c(0, 0.5), c(1.0, 1.0), size = 12,fontface = "plain")


ggsave("FigS3_Bray_stage.jpeg",Bray_stage_final,path = "Dynamic_2018/Figures/",  width = 8, height = 5, units = "in", dpi=500)

ggsave("FigS3_Bray_stage.pdf",Bray_stage_final,path = "Dynamic_2018/Figures/",  width = 8, height = 5, units = "in", dpi=800)



###mantel test####



manfun <-function(x, y) {
  man=(mantel.rtest(x, y,nrepet = 9999))}


Pos <- as.matrix(read.table("Dynamic_2018/Distance_between_positions.txt", head=T))
Pos_dis <- as.dist(Pos)




Dyn_E1_read_mean <- aggregate(Reads~Pond+DP+Day+Position+Species, data = Dynamic_rs_sample_E1, FUN=mean)



BCD_E1 <- Dyn_E1_read_mean [c(2,5,6)]



BCD_E1_wide <- reshape(BCD_E1, idvar= c("Species"), timevar= "DP", direction = "wide")


BCD_E1_wide.rownames <- data.frame(BCD_E1_wide[,-1], row.names=BCD_E1_wide[,1])



##standardization to propostion
BCD_E1_wide.std <-decostand(BCD_E1_wide.rownames, "total",MARGIN=2)




BCD_E1_D0 <- as.dist(vegdist(t(BCD_E1_wide.std[,grepl("D0", colnames(BCD_E1_wide.std))]), method="bray"))

BCD_E1_D2 <- as.dist(vegdist(t(BCD_E1_wide.std[,grepl("D2", colnames(BCD_E1_wide.std))]), method="bray"))

BCD_E1_D4 <- as.dist(vegdist(t(BCD_E1_wide.std[,grepl("D4", colnames(BCD_E1_wide.std))]), method="bray"))

BCD_E1_D6 <- as.dist(vegdist(t(BCD_E1_wide.std[,grepl("D6", colnames(BCD_E1_wide.std))]), method="bray"))

BCD_E1_D8 <- as.dist(vegdist(t(BCD_E1_wide.std[,grepl("D8", colnames(BCD_E1_wide.std))]), method="bray"))

BCD_E1_D10 <- as.dist(vegdist(t(BCD_E1_wide.std[,grepl("D10", colnames(BCD_E1_wide.std))]), method="bray"))

BCD_E1_D12 <- as.dist(vegdist(t(BCD_E1_wide.std[,grepl("D12", colnames(BCD_E1_wide.std))]), method="bray"))

BCD_E1_D14 <- as.dist(vegdist(t(BCD_E1_wide.std[,grepl("D14", colnames(BCD_E1_wide.std))]), method="bray"))


Mantel_E1_list <- list(BCD_E1_D0,BCD_E1_D2,BCD_E1_D4,BCD_E1_D6,BCD_E1_D8,BCD_E1_D10,BCD_E1_D12,BCD_E1_D14,Pos_dis) 





Dayclass <- c("D0","D2","D4","D6","D8","D10","D12","D14","Distance")



M_E1_result <- matrix(nrow = 9, ncol = 9) 




for (i in 1:8){
  for (j in (i+1):9){
    M_E1_result[j,i] <- round(manfun(Mantel_E1_list[[i]],Mantel_E1_list[[j]])$obs,3)
    M_E1_result[i,j] <- round(manfun(Mantel_E1_list[[i]], Mantel_E1_list[[j]])$pvalue,3)
  }
  
}


# for (i in 1:8){
#   for (j in (i+1):9){
#     M_E1_result[j,i] <- paste("r=",round(manfun(Mantel_E1_list[[i]],Mantel_E1_list[[j]])$obs,3))
#     M_E1_result[i,j] <- paste("P=",round(manfun(Mantel_E1_list[[i]], Mantel_E1_list[[j]])$pvalue,3))
#   }
# 
# }



# for (i in 1:8){
#   for (j in (i+1):9){
#     M_E1_result[j,i] <- round(manfun(Mantel_E1_list[[i]],Mantel_E1_list[[j]])$obs,3)
#     M_E1_result[i,j] <- round(manfun(Mantel_E1_list[[i]],Mantel_E1_list[[j]])$obs,3)
#   }
#   
# }

rownames(M_E1_result) <-Dayclass

colnames(M_E1_result) <-Dayclass


M_E1_result [is.na(M_E1_result)] <- 1


heat_E1 <- melt(M_E1_result, na.rm = TRUE)   # transform the matrix to columns




# Fix the order
heat_E1$Var1 <- factor(heat_E1$Var1, levels=c("D0","D2","D4","D6","D8","D10","D12","D14","Distance"), ordered = TRUE)
heat_E1$Var2 <- factor(heat_E1$Var2, levels=c("D0","D2","D4","D6","D8","D10","D12","D14","Distance"), ordered = TRUE)



HM_E1 <-ggplot(heat_E1, aes(Var2,Var1, fill=value)) + geom_tile(colour = "white")+
        scale_fill_gradient2(name="") +
        coord_fixed()+
        theme_minimal()+ # minimal theme
        theme(text=element_text(size=9.5),axis.text.y = element_text(hjust = 0.5),
              axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.minor=element_blank(),
              panel.grid.major =element_blank(),legend.title = element_text(size=8),legend.position="bottom")+
        geom_text(aes(label = value), color = "black", size = 2.5)



legend_HM_E1 <- get_legend(HM_E1)




HM_E1<- HM_E1 + theme(legend.position="none")




Dyn_E4_read_mean <- aggregate(Reads~Pond+DP+Day+Position+Species, data = Dynamic_rs_sample_E4, FUN=mean)



BCD_E4 <- Dyn_E4_read_mean [c(2,5,6)]



BCD_E4_wide <- reshape(BCD_E4, idvar= c("Species"), timevar= "DP", direction = "wide")


BCD_E4_wide.rownames <- data.frame(BCD_E4_wide[,-1], row.names=BCD_E4_wide[,1])



##standardization to propostion
BCD_E4_wide.std <-decostand(BCD_E4_wide.rownames, "total",MARGIN=2)




BCD_E4_D0 <- as.dist(vegdist(t(BCD_E4_wide.std[,grepl("D0", colnames(BCD_E4_wide.std))]), method="bray"))

BCD_E4_D2 <- as.dist(vegdist(t(BCD_E4_wide.std[,grepl("D2", colnames(BCD_E4_wide.std))]), method="bray"))

BCD_E4_D4 <- as.dist(vegdist(t(BCD_E4_wide.std[,grepl("D4", colnames(BCD_E4_wide.std))]), method="bray"))

BCD_E4_D6 <- as.dist(vegdist(t(BCD_E4_wide.std[,grepl("D6", colnames(BCD_E4_wide.std))]), method="bray"))

BCD_E4_D8 <- as.dist(vegdist(t(BCD_E4_wide.std[,grepl("D8", colnames(BCD_E4_wide.std))]), method="bray"))

BCD_E4_D10 <- as.dist(vegdist(t(BCD_E4_wide.std[,grepl("D10", colnames(BCD_E4_wide.std))]), method="bray"))

BCD_E4_D12 <- as.dist(vegdist(t(BCD_E4_wide.std[,grepl("D12", colnames(BCD_E4_wide.std))]), method="bray"))

BCD_E4_D14 <- as.dist(vegdist(t(BCD_E4_wide.std[,grepl("D14", colnames(BCD_E4_wide.std))]), method="bray"))


Mantel_E4_list <- list(BCD_E4_D0,BCD_E4_D2,BCD_E4_D4,BCD_E4_D6,BCD_E4_D8,BCD_E4_D10,BCD_E4_D12,BCD_E4_D14,Pos_dis) 





Dayclass <- c("D0","D2","D4","D6","D8","D10","D12","D14","Distance")



M_E4_result <- matrix(nrow = 9, ncol = 9) 


for (i in 1:8){
  for (j in (i+1):9){
    M_E4_result[j,i] <- round(manfun(Mantel_E4_list[[i]],Mantel_E4_list[[j]])$obs,3)
    M_E4_result[i,j] <- round(manfun(Mantel_E4_list[[i]], Mantel_E4_list[[j]])$pvalue,3)
  }

}




# for (i in 1:8){
#   for (j in (i+1):9){
#     M_E4_result[j,i] <- paste("r=",round(manfun(Mantel_E4_list[[i]],Mantel_E4_list[[j]])$obs,3))
#     M_E4_result[i,j] <- paste("P=",round(manfun(Mantel_E4_list[[i]], Mantel_E4_list[[j]])$pvalue,3))
#   }
#   
# }



# for (i in 1:8){
#   for (j in (i+1):9){
#     M_E4_result[j,i] <- round(manfun(Mantel_E4_list[[i]],Mantel_E4_list[[j]])$obs,3)
#     M_E4_result[i,j] <- round(manfun(Mantel_E4_list[[i]],Mantel_E4_list[[j]])$obs,3)
#   }
#   
# }

rownames(M_E4_result) <-Dayclass

colnames(M_E4_result) <-Dayclass


M_E4_result [is.na(M_E4_result)] <- 1


heat_E4 <- melt(M_E4_result, na.rm = TRUE)   # transform the matrix to columns


names(heat_E4)

# Fix the order
heat_E4$Var1 <- factor(heat_E4$Var1, levels=(c("D0","D2","D4","D6","D8","D10","D12","D14","Distance")))
heat_E4$Var2 <- factor(heat_E4$Var2, levels=(c("D0","D2","D4","D6","D8","D10","D12","D14","Distance")))



HM_E4 <-ggplot(heat_E4, aes(Var2, Var1, fill=value)) + geom_tile(colour = "white")+
        scale_fill_gradient2(name="") +
        coord_fixed()+
        theme_minimal()+ # minimal theme
        theme(text=element_text(size=9.5),axis.text.y = element_text(hjust = 0.5),
              axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.minor=element_blank(),
              panel.grid.major =element_blank(),legend.title = element_text(size=8))+
        geom_text(aes(label = value), color = "black", size = 2.5)











###Fig4_HM#####
HM <- grid.arrange(HM_E1,HM_E4, ncol = 2, nrow = 1,
                     widths = c(2.7,3.38), heights = c(2.5))


HM2_final <- ggdraw(HM)+draw_plot_label(c("(a)","(b)"), c(0, 0.45), c(1.0, 1.0), size = 10,fontface = "plain")

ggsave("Fig4_HM.jpeg",HM2_final,path = "Dynamic_2018/Figures/",  width = 8, height = 3.5, units = "in", dpi=500)

ggsave("Fig4_HM.pdf",HM2_final,path = "Dynamic_2018/Figures/",  width = 8, height = 3.5, units = "in", dpi=800)

# 
# FIG2 <- grid.arrange(Com2,HM2, ncol = 1, nrow = 2,
#                      widths = c(2.7), heights = c(3.0,2.3))
# 
# FIG2_final <- ggdraw(FIG2)+draw_plot_label(c("A", "C","B","D"), c(0, 0, 0.51,0.51), c(0.95, 0.46, 0.95, 0.46), size = 15)
# 
# 
# ggsave("Fig3.jpeg",FIG2_final,path = "Dynamic_2018/Figures/",  width = 11, height = 12, units = "in", dpi=800)


#NMDS####
# E1PA.dist=vegdist(t(BCD_E1_wide.std),method="bray")  # estimate the beta diversity
# E1PA.pcoa=cmdscale(E1PA.dist)           # PCoA analysis
# E1PA.pcoa2 <- data.frame(E1PA.pcoa)
# 
# E1PA.pcoa2$Rep <- rownames(E1PA.pcoa2)
# 
# ggplot(E1PA.pcoa2, aes(x=X1, y=X2))+
#   geom_point(shape = 21,size=3, stroke=1)+
#   theme_bw()+ xlab("PCoA1") + ylab("PCoA2") + 
#   geom_text_repel(aes(label = E1PA.pcoa2$Rep), show.legend=F)





BCD_E1_wide.NMDS <- t(BCD_E1_wide.std)


NMDS_E1_data <- metaMDS(BCD_E1_wide.NMDS, distance="bray", k=2, trymax=2000)


NMDS_E1_data.scores <- as.data.frame(scores(NMDS_E1_data))

NMDS_E1_data.scores$Rep <- rownames(NMDS_E1_data.scores)

NMDS_E1_data.scores$Rep <- as.factor(NMDS_E1_data.scores$Rep)


levels(NMDS_E1_data.scores$Rep) <- gsub("Reads.", "", levels(NMDS_E1_data.scores$Rep))



NMDS_E1<- ggplot(NMDS_E1_data.scores, aes(x=NMDS1, y=NMDS2))+
  geom_point(size=3)+
  theme_bw()+
  geom_text_repel(aes(label = NMDS_E1_data.scores$Rep), show.legend=F)





BCD_E4_wide.NMDS <- t(BCD_E4_wide.std)


NMDS_E4_data <- metaMDS(BCD_E4_wide.NMDS, distance="bray", k=2, trymax=2000)


NMDS_E4_data.scores <- as.data.frame(scores(NMDS_E4_data))

NMDS_E4_data.scores$Rep <- rownames(NMDS_E4_data.scores)

NMDS_E4_data.scores$Rep <- as.factor(NMDS_E4_data.scores$Rep)


levels(NMDS_E4_data.scores$Rep) <- gsub("Reads.", "", levels(NMDS_E4_data.scores$Rep))




NMDS_E4<-ggplot(NMDS_E4_data.scores, aes(x=NMDS1, y=NMDS2))+
  geom_point(size=3)+
  theme_bw()+
  geom_text_repel(aes(label = NMDS_E4_data.scores$Rep), show.legend=F)



NMDS <- grid.arrange(NMDS_E1,NMDS_E4, ncol = 2, nrow = 1,
                     widths = c(2.7,2.7), heights = c(2.5))

NMDS_final <- ggdraw(NMDS)+draw_plot_label(c("A","B"), c(0, 0.5), c(1.0, 1.0), size = 12)


ggsave("NMDS.jpeg",NMDS_final,path = "Dynamic_2018/Figures/",  width = 8, height = 5, units = "in", dpi=800)




#SO####
Dyn_E1_SO <- Dyn_E1_read_mean[which(Dyn_E1_read_mean$Reads >0),] 



Dyn_E1_SO2 = data.frame(Pond=factor(Dyn_E1_SO$Pond),Day=Dyn_E1_SO$Day,
                             Position=factor(Dyn_E1_SO$Position),Species=factor(Dyn_E1_SO$Species))


Dyn_E1_SO_data <-dcast(Dyn_E1_SO2, Pond+Species~Day)


Dyn_E4_SO <- Dyn_E4_read_mean[which(Dyn_E4_read_mean$Reads >0),] 



Dyn_E4_SO2 = data.frame(Pond=factor(Dyn_E4_SO$Pond),Day=Dyn_E4_SO$Day,
                        Position=factor(Dyn_E4_SO$Position),Species=factor(Dyn_E4_SO$Species))


Dyn_E4_SO_data <-dcast(Dyn_E4_SO2, Pond+Species~Day)


Dyn_SO_data<- rbind(Dyn_E1_SO_data,Dyn_E4_SO_data)

write.csv(Dyn_SO_data,file = "Dynamic_2018/Dyn_SO_data.csv",row.names = FALSE)


Dyn_SO_rs <- melt(Dyn_SO_data,id=c(1:2,ncol(Dyn_SO_data)))

#Rename the variables

Dyn_SO_rs <- rename (Dyn_SO_rs,c("variable"="Day","value"="SN"))



ggplot(Dyn_SO_rs, aes(x = Day, y=SN/5, fill=Species))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~Pond)+
  scale_y_continuous(limits = c(0, 1), breaks=seq(0, 1.0, 0.1))+
  labs(x="Sampling day", y="Site occupancy")+
  scale_fill_manual(values=c("firebrick3","green4","blue","orange2","tan4","hotpink","gray10","purple3"))+
  guides(fill = guide_legend(ncol = 8, byrow = TRUE))+
  theme_bw()+
  theme(text=element_text(size=15),legend.position="bottom",legend.title = element_text(size=12),legend.text = element_text(size=12),
        panel.grid.minor=element_blank(),panel.grid.major =element_blank(),strip.text=element_blank())



#PERMANOVA####



E1.dis <- as.matrix(vegdist(t(BCD_E1_wide.std), method="bray"))

BCD_E1_meta <- Dyn_E1_read_mean [c(1,3,4,5,6)]

BCD_E1_meta <- aggregate(Reads~Pond+Day+Position, data = BCD_E1_meta, FUN=sum)


E1.adois <- adonis(E1.dis~BCD_E1_meta$Day+BCD_E1_meta$Position)$aov.tab


E4.dis <- as.matrix(vegdist(t(BCD_E4_wide.std), method="bray"))

BCD_E4_meta <- Dyn_E4_read_mean [c(1,3,4,5,6)]

BCD_E4_meta <- aggregate(Reads~Pond+Day+Position, data = BCD_E4_meta, FUN=sum)


E4.adois <- adonis(E4.dis~BCD_E4_meta$Day+BCD_E4_meta$Position)$aov.tab




#https://github.com/pmartinezarbizu/pairwiseAdonis
library(pairwiseAdonis)




E1_merge <- cbind(t(BCD_E1_wide.std),BCD_E1_meta)


pairwise.adonis(E1_merge[,1:8],E1_merge$Day)

pairwise.adonis(E1_merge[,1:8],E1_merge$Position)


E4_merge <- cbind(t(BCD_E4_wide.std),BCD_E4_meta)


pairwise.adonis(E4_merge[,1:8],E4_merge$Day)

pairwise.adonis(E4_merge[,1:8],E4_merge$Position)


###addtional species####

Dynamic_E1_add <- Dynamic_E1[which(Dynamic_E1$Species == 'RUD' |
                                     Dynamic_E1$Species == 'CHU'),]

Dynamic_E4_add <- Dynamic_E4[which(Dynamic_E4$Species == 'RUD' |
                                     Dynamic_E4$Species == 'DAC'),]

Dynamic_add <- rbind (Dynamic_E1_add,Dynamic_E4_add)


Dynamic_add$Day<- factor(Dynamic_add$Day, levels = (c("D0","D2","D4","D6","D8","D10","D12","D14")))




# ggplot(Dynamic_add,aes(x=Day,y=Percent_reads,shape=Species))+
#   geom_point(size=3)+facet_wrap(~Pond+Position,nrow=2,labeller=labeller(.multi_line = FALSE))+
#   labs(x="Sampling posiiton", y= "Proportional abundance")+
#   scale_shape_manual(values=c(16,0,2))+
#   theme_bw()+theme(text=element_text(size=12),strip.background = element_blank(),
#                    panel.grid.minor=element_blank(),panel.grid.major =element_blank())


    
add <- ggplot(Dynamic_add,aes(x=Day,y=Percent_reads,ymin=Percent_reads-SD,ymax=Percent_reads+SD,shape=Species,group=Species))+
      geom_point(size=2)+geom_errorbar(width=0.2)+
      facet_wrap(~Pond+Position,nrow=2,labeller=labeller(.multi_line = FALSE))+
      scale_shape_manual(values=c(16,0,2))+geom_line(size=0.5,linetype = 3)+
      labs(x="Sampling day", y= "Proportional abundance")+
      theme_bw()+theme(text=element_text(size=12),strip.background = element_blank(),
                       panel.grid.minor=element_blank(),panel.grid.major =element_blank(),
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                       legend.title = element_text(size=10),legend.text = element_text(size=9))



add_final <- ggdraw(add)+draw_plot_label(c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)"), c(0.08,0.24,0.4,0.56,0.72), c(0.983,0.983,0.983,0.983,0.983,0.536,0.536,0.536,0.536,0.536), size = 12,fontface = "plain")

#Fig3_add_pos####
ggsave("Fig3_add_pos.jpeg",add_final,path = "Dynamic_2018/Figures/",  width = 8, height = 6, units = "in", dpi=500)

ggsave("Fig3_add_pos.pdf",add_final,path = "Dynamic_2018/Figures/",  width = 8, height = 6, units = "in", dpi=800)





####Pond fish summary####


Pond_summary <- read.csv(file = "Dynamic_2018/Fish information summary 2016.csv")

Pond_summary$Species <-mapvalues(Pond_summary$Species, c("Abramis brama","Barbus barbus","Carassius carassius","Squalius cephalus","Leuciscus leuciscus",
                                                         "Rutilus rutilus","Tinca tinca"), 
                                                      c("BRE","BAR","CAR","CHU","DAC","ROA","TEN"))


Pond_summary$Pond <-mapvalues(Pond_summary$Pond, c("E1","E4"), 
                              c("(a)","(b)"))

# FigS1 AND S2 ------------------------------------------------------------



Pond_summary$Species <-factor(Pond_summary$Species,levels = (c("CHU","BAR","BRE","CAR","ROA","TEN","DAC")))



ggplot(Pond_summary, aes(x= Date, y= Number,colour=Species,group=Species))+
geom_point(stat="identity",size=3.5) + geom_line(size=0.8,linetype = 3)+
facet_wrap(~Pond,ncol=2)+
scale_colour_manual(values=c("green4","blue","orange2","tan4","hotpink","gray10","purple3"))+
scale_x_discrete(limits=c("Jun-16","Jul-16","Aug-16",'Sep-16','Oct-16','Nov-16'))+
labs(x="Month", y="Fish abundance")+ theme_bw()+
theme(text=element_text(size=15),strip.background = element_blank(),strip.text=element_text(hjust=0),
      panel.grid.minor=element_blank(),panel.grid.major =element_blank(),
      legend.text = element_text(size=9),legend.title = element_text(size=10),
      axis.text.x =element_text(size=10),axis.text.y =element_text(size=10))



ggsave("FigS1_number.jpeg",path = "Dynamic_2018/Figures/",  width = 8, height = 4.5, units = "in", dpi=800)







ggplot(Pond_summary, aes(x= Date, y= Biomass,colour=Species,group=Species))+
geom_point(stat="identity",size=3.5) + geom_line(size=0.8,linetype = 3)+
facet_wrap(~Pond,ncol=2)+
scale_colour_manual(values=c("green4","blue","orange2","tan4","hotpink","gray10","purple3"))+
scale_x_discrete(limits=c("Jun-16","Jul-16","Aug-16",'Sep-16','Oct-16','Nov-16'))+
labs(x="Month", y="Fish biomass (kg)")+ theme_bw()+
theme(text=element_text(size=15),strip.background = element_blank(),strip.text=element_text(hjust=0),
      panel.grid.minor=element_blank(),panel.grid.major =element_blank(),
      legend.text = element_text(size=9),legend.title = element_text(size=10),
      axis.text.x =element_text(size=10),axis.text.y =element_text(size=10))



ggsave("FigS2_biomass.jpeg",path = "Dynamic_2018/Figures/",  width = 8, height = 4.5, units = "in", dpi=800)




