library(tidyverse)
library(reshape)
library(ggpubr)

#Allelic Diversity
pi_est<-read.table("...\\allsamples_merged_sorted2.pi.windowed.pi", header = T)

gghistogram(pi_est$PI, bins=20, fill ="grey70")+
  labs(x="Allelic diversity")

#Tajima's D
taj_est <- read.table("...\\tajimasd.Tajima.D",header=T)

hist(taj_est$TajimaD, br=20)

ggarrange(labels="AUTO",
          ggboxplot(pi_est$PI, fill="orange")+
            labs(y="Diversity", x="")+
            theme(axis.text.x = element_blank(), 
                  axis.ticks.x =element_blank()), 
          ggboxplot(taj_est$TajimaD, fill="purple")+
            labs(y="Tajima's D", x="")+
            theme(axis.text.x = element_blank(), 
                  axis.ticks.x =element_blank()), 
          nrow=1)

taj_est2<-subset(taj_est, taj_est$TajimaD>min(boxplot.stats(taj_est$TajimaD)$out))

#Fst
fst_est<-read.table("...\\sample2.tsv", header=T)
length(unique(subset(fst_est, WEIR_AND_COCKERHAM_FST<0)$CHROM))
length(unique(subset(fst_est, WEIR_AND_COCKERHAM_FST==1)$CHROM))

fst_est[fst_est$WEIR_AND_COCKERHAM_FST < 0, "WEIR_AND_COCKERHAM_FST"] <- 0

ggboxplot(fst_est$WEIR_AND_COCKERHAM_FST, fill="lightblue")+
  labs(y="Fst values", x="")+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x =element_blank())

length(unique(fst_est$CHROM))

ggplot(fst_est, aes(POS, WEIR_AND_COCKERHAM_FST))+
  geom_point()+
  theme_minimal()

ggplot(subset(fst_est, CHROM=="VZTU01030148.1"), aes(POS, WEIR_AND_COCKERHAM_FST))+
  geom_point()+
  theme_minimal()

#count of genes with fixed snps
length(unique(subset(fst_est, WEIR_AND_COCKERHAM_FST==1)$CHROM))

fst_means<-fst_est %>%
  group_by(CHROM) %>%
  summarise(MEAN = mean(WEIR_AND_COCKERHAM_FST), 
            COUNT=n(), 
            COUNT1=sum(WEIR_AND_COCKERHAM_FST==1))

ggplot(fst_means, aes(CHROM, MEAN))+
  geom_point()+
  labs(y="Mean Fst per Gene", x="Gene Order")+
  theme_minimal()

