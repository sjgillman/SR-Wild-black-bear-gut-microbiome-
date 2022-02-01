####### Supplemental Text 2 Script for Statistical Analysis in R #######
### Statistical Analysis for Gillman et al 2020



library(microbiome) ## data analysis
library(qiime2R) # import data
library(phyloseq) # also the basis of data object. Data analysis and visualization
library(vegan) # some utility tools
library(data.table) # alternative to data.frame
library(dplyr) # data handling
library(tidyverse)
library(ggpubr) ## plotting 
library(ggplot2)
library(mctoolsr)
library(picante) ## faith's PD
library(Rmisc)

setwd("~/Desktop/Projects/Bear/Bear-R/CLEAN/FINAL")

#### Import & create phyloseq dataframe with qiime2R and QIIME2 artifacts #####
## Following Tutorial: Integrating QIIME2 and R for data visualization and analysis using qiime2R by J. Bisanz
## you will need
# 1.) Metafile.tsv (alpha_tableR.tsv) -the alpha_table file will need to have the second row removed and the # infront of SampleID removed for it to read okay
# 2.) taxonomy.qza
# 3.) table.qza
# 4.) rooted.qza

## import artifacts & metadata file
metadata<-read_tsv("Metafile.tsv")
SVs<-read_qza("table.qza")
taxonomy<-read_qza("taxonomy.qza")
taxtable<-taxonomy$data %>% as_tibble() %>% separate(Taxon, sep=";", c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) #convert the table into a tabular split version
tree<-read_qza("rooted-tree.qza")


## Create the phyloseq object
phy_obj<-phyloseq(
  otu_table(SVs$data, taxa_are_rows = T), 
  phy_tree(tree$data), 
  tax_table(as.data.frame(taxtable) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), #moving the taxonomy to the way phyloseq wants it
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("SampleID")))



## view data table
datatable(tax_table(phy_obj))

##### Clean Taxonomy table #####
## Rename NAs to last known group
tax.clean <- data.frame(tax_table(phy_obj))
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}
## import new taxonomy table
tax_table(phy_obj) <- as.matrix(tax.clean)

## view
datatable(tax_table(phy_obj))

###### Rename uncultured
tax.clean2 <- data.frame(tax_table(phy_obj))

for (i in 1:7){ tax.clean2[,i] <- as.character(tax.clean2[,i])}
for (i in 1:nrow(tax.clean2)){
  if (tax.clean2[i,2] == "uncultured"){
    kingdom <- paste("Kingdom_", tax.clean2[i,1], sep = "")
    tax.clean2[i, 2:7] <- kingdom
  } else if (tax.clean2[i,3] == "uncultured"){
    phylum <- paste("Phylum_", tax.clean2[i,2], sep = "")
    tax.clean2[i, 3:7] <- phylum
  } else if (tax.clean2[i,4] == "uncultured"){
    class <- paste("Class_", tax.clean2[i,3], sep = "")
    tax.clean2[i, 4:7] <- class
  } else if (tax.clean2[i,5] == "uncultured"){
    order <- paste("Order_", tax.clean2[i,4], sep = "")
    tax.clean2[i, 5:7] <- order
  } else if (tax.clean2[i,6] == "uncultured"){
    family <- paste("Family_", tax.clean2[i,5], sep = "")
    tax.clean2[i, 6:7] <- family
  } else if (tax.clean2[i,7] == ""){
    tax.clean2$Species[i] <- paste("Genus",tax.clean.bear$Genus[i], sep = "_")
  }
}

## import new taxonomy table
tax_table(phy_obj) <- as.matrix(tax.clean2)

## view new table
datatable(tax_table(phy_obj))


## save phyloseq object
saveRDS(phy_obj, "~/Desktop/Projects/Bear/Bear-R/CLEAN/FINAL/physeq.rds")

## if you ever want to pull back in
phy_obj<- readRDS("physeq.rds")

####  Alpha Diversity ####
##  Equal sample sums
set.seed(9242) ## ensures rarifies the same each time script is run

summary(sample_sums(phy_obj)) ## helps determine depth for rarifying 

## rarefying: we already know our depth: 1050 so rarefy to that
phyb.rar <- rarefy_even_depth(phy_obj, sample.size = 1050)
## lost one sample: S100J
summary(sample_sums(phyb.rar)) ##checking to see they all have the same sequence depth
any(taxa_sums(phy_obj)== 0) # making sure we dont have any sequences not in at least one samples
###  run this is you do have 0's:ps1a <- prune_taxa(taxa_sums(phyb.rar) > 0, phyb.rar)


##  pull metadata from physeq object
sam.meta <- meta(phyb.rar)
sam.meta

## put variables is particular order
sam.meta$GIT<-factor(sam.meta$GIT, levels=c("Jejunum", "Colon"))
sam.meta$AgeClass<-factor(sam.meta$AgeClass, levels=c("Yearling", "Subadult", "Adult", "Unknown"))

## Add the rownames as a new colum for easy integration later.
sam.meta$sam_name <- rownames(sam.meta)

#### Non-phylogenetic diversities: Shannon ####
## calculated with microbiome package
## 
div_shan<- microbiome::alpha(phyb.rar, index = "shannon")
## can run index= "all" if you desire all alpha indices

## Add the rownames to diversity table
div_shan$sam_name <- rownames(div_shan)


#### Non-phylogenetic diversities: Simpson ####
## calculated with microbiome package
div_sim<- microbiome::alpha(phyb.rar, index="diversity_inverse_simpson") 

## Add the rownames to diversity table
div_sim$sam_name <- rownames(div_sim)


#### Phylogenetic diversity: Faith's PD #####
#Phylogenetic diversity is calculated using the picante package.

## pull ASV table
phyb.rar.asvtab <- as.data.frame(phyb.rar@otu_table)

## pull tree
phyb.rar.tree <- phyb.rar@phy_tree

## We first need to check if the tree is rooted or not 

phyb.rar@phy_tree
###rooted so we are good to go

## Getting the data ready
div_pd <- pd(t(phyb.rar.asvtab), phyb.rar.tree,include.root=T) 
# t(ou_table) transposes the table for use in picante and the
#tree file comes from the first code  we used to read tree
#file (see making a phyloseq object section)

## Add the rownames to diversity table
div_pd$sam_name <- rownames(div_pd)

## STEP 4p. merge all of the alphas into one file
merged_table<-merge(div_pd,div_shan, by = "sam_name", all=T)
merged_table2<-merge(merged_table,sam.meta, by = "sam_name", all=T)
alpha_table <- merge(merged_table2,div_sim, by = "sam_name", all=T)

datatable(alpha_table) ## this now has all alpha measures in one datatable!

##### produce summary tables for diversity indices for  age-class ####
## note: only one yearling female will not analyze 
## females



##### Community composition  ######

## relative abundance
pseq.rel <- microbiome::transform(phyb.rar, "compositional")

## merge to phylum rank
phlyum <- tax_glom(pseq.rel, taxrank = "Phylum")
ntaxa(phlyum)
#23
## melt
phylum_melt<- psmelt(phlyum)


unique(phylum_melt$Phylum)
#23
## get summary statistics phyla GIT

p_abund<-summarySE(phylum_melt, measurevar = "Abundance", groupvars =c("Phylum", "GIT"))

##remove 0 abundance
p_abund$Abundance[p_abund$Abundance==0] <- NA
p_abund<-p_abund[complete.cases(p_abund$Abundance),]
p_abund<- p_abund %>% 
  mutate_if(is.numeric, round, digits = 5)
unique(p_abund$Phylum)

## ageclass

## remove unknown ages
phylum_melt2<-phylum_melt
phylum_melt2$AgeClass[phylum_melt2$Ages==0] <- NA
phylum_melt2<-phylum_melt2[complete.cases(phylum_melt2$Age),]
age_p_abund<-summarySE(phylum_melt2, measurevar = "Abundance", groupvars =c("Phylum", "AgeClass"))
age_p_abund<- age_p_abund %>% 
  mutate_if(is.numeric, round, digits = 5)

##remove 0 abundance
age_p_abund$Abundance[age_p_abund$Abundance==0] <- NA
age_p_abund<-age_p_abund[complete.cases(age_p_abund$Abundance),]
unique(age_p_abund$Phylum)
## genus 
## merge to phylum rank
genus<- tax_glom(pseq.rel, taxrank = "Genus")
ntaxa(genus)
#
## melt
genus_melt<- psmelt(genus)

## get summary statistics genus GIT

g_abund<-summarySE(genus_melt, measurevar = "Abundance", groupvars =c("Genus", "GIT"))

##remove 0 abundance
g_abund$Abundance[g_abund$Abundance==0] <- NA
g_abund<-g_abund[complete.cases(g_abund$Abundance),]
g_abund<- g_abund %>% 
  mutate_if(is.numeric, round, digits = 5)
unique(g_abund$Genus)

## ageclass
genus_melt2<-genus_melt
genus_melt2$AgeClass[genus_melt2$Ages==0] <- NA
genus_melt2<-genus_melt2[complete.cases(genus_melt2$Age),]
age_g_abund<-summarySE(genus_melt2, measurevar = "Abundance", groupvars =c("Genus", "AgeClass"))

##remove 0 abundance
age_g_abund$Abundance[age_g_abund$Abundance==0] <- NA
age_g_abund<-age_g_abund[complete.cases(age_g_abund$Abundance),]
age_g_abund<- age_g_abund %>% 
  mutate_if(is.numeric, round, digits = 5)
unique(age_g_abund$Genus)


#### LMER Analysis for Alpha diversity #####

library(lmerTest)
library(lme4)
library(car)
library(emmeans)
library(sjmisc)
library(sjPlot)
library(tidyverse)



### bear samples without age-class (n=4) have been removed from the analysis
## rename
alpha_table0<-alpha_table

alpha_table0$GIT<-factor(alpha_table0$GIT, levels=c("Jejunum", "Colon"))
alpha_table0$AgeClass<-factor(alpha_table0$AgeClass, levels=c("Yearling", "Subadult", "Adult"))



## remove unknown ages
alpha_table0$AgeClass[alpha_table0$Ages==0] <- NA
alpha_table0<-alpha_table0[complete.cases(alpha_table0$Age),]



#### LMER PD ####
## histogram
ggplot(alpha_table0,aes(x=PD))+geom_histogram()
ggplot(alpha_table0,aes(x=log(PD)))+geom_histogram()+ ggtitle("log transformed PD values") 

## skewed right so log transform

## create intial models with maximum likelihood
pd_lme1<-lmer(log(PD)~GIT+Sex+AgeClass+(1|Subject),data=alpha_table0, REML=F)
pd_lme2<-lmer(log(PD)~GIT+GIT*Sex+Sex+AgeClass+(1|Subject),data=alpha_table0, REML=F)
pd_lme3<-lmer(log(PD)~GIT*AgeClass+Sex+GIT+AgeClass+(1|Subject), data=alpha_table0, REML=F)
pd_lme4<-lmer(log(PD)~GIT+Sex+AgeClass+GIT*AgeClass+GIT*Sex+(1|Subject), data=alpha_table0, REML=F)

## compare to determine best model
anova(pd_lme1, pd_lme2, pd_lme3, pd_lme4)

## best model run witm REML=T
pd_lme1<-lmer(log(PD)~GIT+Sex+AgeClass+(1|Subject),data=alpha_table0, REML=T)
summary(pd_lme1)
## run final model to get 
Anova(pd_lme1)
## sig difference between GIT site and for age-class

## get R^2
library(performance)
performance::r2(pd_lme1)

## Estimated Marginal Means
library(emmeans)

emmeans(pd_lme1, pairwise~AgeClass,lmer.df="satterthwaite", adjust="tukey")
## sig difference is between adults and subadults

##Assumption 1 - Linearity

## Graphically, plotting the model residuals (the difference 
#between the observed value and the model-estimated value) vs 
#the predictor 
alpha_table1<-alpha_table0
alpha_table1$log<-log(alpha_table1$PD)
Plot.Model.F.Linearity<-plot(resid(pd_lme1),alpha_table1$log)

## Assumption 2 Homogeneity of Variance
#Regression models assume that variance of the residuals
#is equal across groups. 

#extracts the residuals and places them in a new column in our original data table
alpha_table1$lme10<- residuals(pd_lme1) 
alpha_table1$baslme10 <-abs(alpha_table1$lme10) #creates a new column with the absolute value of the residuals
alpha_table1$lme102 <- alpha_table1$baslme10^2 #squares the absolute values of the residuals to provide the more robust estimate
pd_leven <- lm(lme102 ~ Subject, data=alpha_table1) #ANOVA of the squared residuals
anova(pd_leven) #displays the results

##visually
plot(pd_lme1) #creates a fitted vs residual plot

##Assumption 3: The residuals of the model are normally distributed.

#QQ plots 
library(lattice)

qqmath(pd_lme1, id=0.05)
## overall looks good!!

#### LMER Shannon ####
## histogram
ggplot(alpha_table0,aes(x=diversity_shannon))+geom_histogram() 


## create intial models with maximum likelihood
shan_lme1<-lmer(diversity_shannon~GIT+Sex+AgeClass+(1|Subject),data=alpha_table0, REML=F)
shan_lme2<-lmer(diversity_shannon~GIT+GIT*Sex+Sex+AgeClass+(1|Subject),data=alpha_table0, REML=F)
shan_lme3<-lmer(diversity_shannon~GIT*AgeClass+Sex+GIT+AgeClass+(1|Subject), data=alpha_table0, REML=F)
shan_lme4<-lmer(diversity_shannon~GIT+Sex+AgeClass+GIT*AgeClass+GIT*Sex+(1|Subject), data=alpha_table0, REML=F)

## compare to determine best model
anova(shan_lme1, shan_lme2, shan_lme3, shan_lme4)

## best model run witm REML=T
shan_lme1<-lmer(diversity_shannon~GIT+Sex+AgeClass+(1|Subject),data=alpha_table0, REML=T)

## run final model to get 
Anova(shan_lme1)
## no sig
## get R^2
performance::r2(shan_lme1)


##Assumption 1 - Linearity

## Graphically, plotting the model residuals (the difference 
#between the observed value and the model-estimated value) vs 
#the predictor 
Plot.Model.F.Linearity<-plot(resid(shan_lme1),alpha_table0$diversity_shannon)

## Assumption 2 Homogeneity of Variance
#Regression models assume that variance of the residuals
#is equal across groups. 

#extracts the residuals and places them in a new column in our original data table
alpha_table2<-alpha_table0
alpha_table2$lme10<- residuals(shan_lme1) 
alpha_table2$baslme10 <-abs(alpha_table2$lme10) #creates a new column with the absolute value of the residuals
alpha_table2$lme102 <- alpha_table2$baslme10^2 #squares the absolute values of the residuals to provide the more robust estimate
shan_leven <- lm(lme102 ~ Subject, data=alpha_table2) #ANOVA of the squared residuals
anova(shan_leven) #displays the results

##visually
plot(shan_lme1) #creates a fitted vs residual plot

##Assumption 3: The residuals of the model are normally distributed.

#QQ plots 
qqmath(shan_lme1, id=0.05)
## overall looks good!!

#####  LMER Simpson ####
## histogram
ggplot(alpha_table0,aes(x=diversity_inverse_simpson))+geom_histogram()
ggplot(alpha_table0,aes(x=log(diversity_inverse_simpson)))+geom_histogram()
#### transform with log as it is skewed

## create intial models with maximum likelihood
sim_lme1<-lmer(log(diversity_inverse_simpson)~GIT+Sex+AgeClass+(1|Subject),data=alpha_table0, REML=F)
sim_lme2<-lmer(log(diversity_inverse_simpson)~GIT+GIT*Sex+Sex+AgeClass+(1|Subject),data=alpha_table0, REML=F)
sim_lme3<-lmer(log(diversity_inverse_simpson)~GIT*AgeClass+Sex+GIT+AgeClass+(1|Subject), data=alpha_table0, REML=F)
sim_lme4<-lmer(log(diversity_inverse_simpson)~GIT+Sex+AgeClass+GIT*AgeClass+GIT*Sex+(1|Subject), data=alpha_table0, REML=F)

## determine best fit model
anova(sim_lme1, sim_lme2, sim_lme3, sim_lme4)

## lme1 best run with REML=T
sim_lme1<-lmer(log(diversity_inverse_simpson)~GIT+Sex+AgeClass+(1|Subject),data=alpha_table0, REML=T)
summary(sim_lme1)

## run final model to get 
Anova(sim_lme1)
## no sig

## get R^2
performance::r2(sim_lme1)

##Assumption 1 - Linearity

## Graphically, plotting the model residuals (the difference 
#between the observed value and the model-estimated value) vs 
#the predictor 
alpha_table3<-alpha_table0
alpha_table3$log<-log(alpha_table3$diversity_inverse_simpson)
Plot.Model.F.Linearity<-plot(resid(sim_lme1),alpha_table3$log)

## Assumption 2 Homogeneity of Variance
#Regression models assume that variance of the residuals
#is equal across groups. 

#extracts the residuals and places them in a new column in our original data table
alpha_table3$lme10<- residuals(sim_lme1) 
alpha_table3$baslme10 <-abs(alpha_table3$lme10) #creates a new column with the absolute value of the residuals
alpha_table3$lme102 <- alpha_table3$baslme10^2 #squares the absolute values of the residuals to provide the more robust estimate
sim_leven <- lm(lme102 ~ Subject, data=alpha_table3) #ANOVA of the squared residuals
anova(sim_leven) #displays the results

##visually
plot(sim_lme1) #creates a fitted vs residual plot



##Assumption 3: The residuals of the model are normally distributed.

#QQ plots 

qqmath(sim_lme1, id=0.05)
## some deviation but overall looks good!!


##### Beta Diversity ####

#### UniFrac ####
## remove unknowns
new_obj.rel= subset_samples(pseq.rel, AgeClass != "Unknown")
## weighted
unifrac.dist <- UniFrac(new_obj.rel, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

WU_permanova <- adonis(unifrac.dist ~ GIT+Sex+AgeClass, strata=alpha_table1$Subject,data = alpha_table1)
WU_permanova 
## GIT sig
##### Checking the homogeneity condition

#GIT
permutest(betadisper(unifrac.dist, alpha_table1$GIT), strata=Subject)
##  not significant!
permutest(betadisper(unifrac.dist, alpha_table1$Sex), strata=Subject)
## not significant
permutest(betadisper(unifrac.dist, alpha_table1$AgeClass), strata=Subject)
## not significant

## unweighted
ununifrac.dist <- UniFrac(new_obj, 
                          weighted = FALSE, 
                          parallel = FALSE, 
                          fast = TRUE)

unWU_permanova <- adonis(ununifrac.dist ~ GIT+Sex+AgeClass, strata=alpha_table1$Subject,data = alpha_table1)
unWU_permanova 
## GIT significant
#GIT
permutest(betadisper(ununifrac.dist, alpha_table1$GIT), strata=Subject)
## significance for GIT so permanova result may be potentially explained by that.

permutest(betadisper(ununifrac.dist, alpha_table1$AgeClass), strata=Subject)
## not significant


#### Community composition  Visualization GIT ####
## bar plot to show mean relative abundance of each community at phyla and genus level.
## phyla

p_abund$Phylum <- as.character(p_abund$Phylum)
unique(p_abund$Phylum)


#simple way to rename phyla with < 1% abundance
p_abund$Phylum[p_abund$Abundance <= 0.01] <- "Other"
unique(p_abund$Phylum)

## put phyla in order fso it plots most abundant on bottom
p_abund$Phylum <- factor(p_abund$Phylum, levels = c("Other","Actinobacteria", "Epsilonbacteraeota", "Proteobacteria", "Firmicutes"))


## plot
spatial_plot <- ggplot(data=p_abund, aes(x=GIT, y=Abundance, fill=Phylum, width=.8))+
  coord_flip()
p1<-spatial_plot + geom_bar(aes(),stat="identity", position="stack", width =1) +
  scale_fill_manual(values = c("black","gray","deeppink4", "cyan", "forestgreen"))+
  theme_bw()+
  theme(legend.position="bottom",axis.title=element_text(size=9, family="Arial"),legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size=5, family="Arial"),legend.title.align=0.5,
        legend.key.width=unit(0.2,'cm'),legend.spacing.x = unit(.2, 'cm'),
        legend.spacing.y=unit(.2, 'cm'),legend.title = element_text(size=9, family = "Arial"),
        axis.text.x =element_blank(),axis.ticks.x =element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  guides(fill=guide_legend(nrow=1, byrow=TRUE, color=guide_legend(nrow=1),reverse = TRUE, title.position = "top"))+
  xlab("")+ ylab("")
p1

## genus
unique(g_abund$Genus)
g_abund$Genus <- as.character(g_abund$Genus)

#simple way to rename phyla with < 1% abundance
g_abund$Genus[g_abund$Abundance <= 0.01] <- "Other"
unique(g_abund$Genus)
write_csv(g_abund, "g_abund.csv")
### put in order you want 
g_abund$Genus <- factor(g_abund$Genus, levels = c("Other","Helicobacter","Moraxella","Weissella", "Staphylococcus","Sporosarcina","Family_Enterobacteriaceae","Family_Pasteurellaceae","Escherichia-Shigella","Terrisporobacter",
                                                  "Paeniclostridium","Romboutsia" ,"Lactococcus", "Turicibacter",  
                                                  "Enterococcus", "Cellulosilyticum","Bacillus","Family_Peptostreptococcaceae",
                                                  "Lactobacillus", "Streptococcus", 
                                                  "Clostridium sensu stricto 1","Sarcina")) 


spatial_plot2 <- ggplot(data=g_abund, aes(x=GIT, y=Abundance, fill=Genus, width=.8))+
  coord_flip()
p2<-spatial_plot2 + geom_bar(aes(),stat="identity", position="stack", width =1) +
  scale_fill_manual(values = c("black","lightpink","turquoise4","cyan2","steelblue4","lightblue","lightslateblue","blue1","cadetblue3","gold1","khaki1","yellow4","lightgoldenrod1","limegreen", "olivedrab2", "springgreen4","lightgreen","palegreen4", "darkseagreen2","green","forestgreen","darkgreen"))+
  theme_bw()+
  theme(legend.position="bottom",axis.title=element_text(size=9, family="Arial"),legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size=5, family="Arial"),
        legend.key.width=unit(0.2,'cm'),legend.spacing.x = unit(.2, 'cm'),legend.title.align=0.5,
        legend.spacing.y=unit(.2, 'cm'),legend.title = element_text(size=9, family="Arial"),
        axis.text =element_text(color="black", family="Arial"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  guides(fill=guide_legend(nrow=3, byrow=TRUE, color=guide_legend(nrow=3),reverse = TRUE, title.position="top"))+
  xlab("")+ ylab("Mean Relative Abundance %")
p2
## plot together & save
tiff('Community.tiff', units="in", width=7, height=4, res=300, compression = 'lzw')
see::plots(p1,p2, n_columns = 1, tags=c("A", "B"))
dev.off()                                     

#### Community composition  Visualization Age class ####
# put phyla in order fso it plots most abundant on bottom
age_p_abund$Phylum <- as.character(age_p_abund$Phylum)
#simple way to rename phyla with < 1% abundance
age_p_abund$Phylum[age_p_abund$Abundance <= 0.01] <- "Other"
unique(age_p_abund$Phylum)

age_p_abund$Phylum <- factor(age_p_abund$Phylum, levels = c("Other", "Actinobacteria","Tenericutes" ,"Epsilonbacteraeota", "Proteobacteria", "Firmicutes"))


## plot
spatial_plot3 <- ggplot(data=age_p_abund, aes(x=AgeClass, y=Abundance, fill=Phylum, width=.8))+
  coord_flip()
p3<-spatial_plot3 + geom_bar(aes(),stat="identity", position="stack", width =1) +
  scale_fill_manual(values = c("black","gray", "yellow","deeppink4", "cyan", "forestgreen"))+
  theme_bw()+
  theme(legend.position="bottom",axis.title=element_text(size=9, family="Arial"),legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size=5, family="Arial"),legend.title.align=0.5,
        legend.key.width=unit(0.2,'cm'),legend.spacing.x = unit(.2, 'cm'),
        legend.spacing.y=unit(.2, 'cm'),legend.title = element_text(size=9, family = "Arial"),
        axis.text.x =element_blank(),axis.ticks.x =element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  guides(fill=guide_legend(nrow=1, byrow=TRUE, color=guide_legend(nrow=1),reverse = TRUE, title.position = "top"))+
  xlab("")+ ylab("")
p3

## genus
age_g_abund$Genus <- as.character(age_g_abund$Genus)

#simple way to rename phyla with < 1% abundance
age_g_abund$Genus[age_g_abund$Abundance <= 0.01] <- "Other"
unique(age_g_abund$Genus)
### put in order you want 
age_g_abund$Genus <- factor(age_g_abund$Genus, levels = c("Other","Helicobacter","Family_Neisseriaceae","Bibersteinia", "Moraxella","Weissella", "Staphylococcus","Sporosarcina","Family_Pasteurellaceae","Escherichia-Shigella","Terrisporobacter",
                                                          "Mycoplasma", "Leuconostoc", "Paeniclostridium","Romboutsia" ,"Lactococcus", "Turicibacter",  
                                                          "Enterococcus", "Cellulosilyticum","Bacillus","Family_Peptostreptococcaceae",
                                                          "Lactobacillus", "Streptococcus", 
                                                          "Clostridium sensu stricto 1","Sarcina")) 


spatial_plot4 <- ggplot(data=age_g_abund, aes(x=AgeClass, y=Abundance, fill=Genus, width=.8))+
  coord_flip()
p4<-spatial_plot4 + geom_bar(aes(),stat="identity", position="stack", width =1) +
  scale_fill_manual(values = c("black","lightpink","purple1","mediumorchid","turquoise4","cyan2","steelblue4","lightblue","lightslateblue","blue1","cadetblue3","gold1","orange3","khaki1","yellow4","lightgoldenrod1","limegreen", "olivedrab2", "springgreen4","lightgreen","palegreen4", "darkseagreen2","green","forestgreen","darkgreen"))+
  theme_bw()+
  theme(legend.position="bottom",axis.title=element_text(size=9, family="Arial"),legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size=5, family="Arial"),
        legend.key.width=unit(0.2,'cm'),legend.spacing.x = unit(.2, 'cm'),legend.title.align=0.5,
        legend.spacing.y=unit(.2, 'cm'),legend.title = element_text(size=9, family="Arial"),
        axis.text =element_text(color="black", family="Arial"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  guides(fill=guide_legend(nrow=3, byrow=TRUE, color=guide_legend(nrow=3),reverse = TRUE, title.position="top"))+
  xlab("")+ ylab("Mean Relative Abundance %")
p4
## plot together & save
tiff('CommunityAge.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
see::plots(p3,p4, n_columns = 1, tags=c("A", "B"))
dev.off()                                


#### Phylogenetic differences via heat trees of GIT sites ####
## following metacoder example analysis

library(metacoder)
library(taxa)

## be sure to set.seed to ensure the plots are the same if you go back to recreate
set.seed(9242)

## Convert rarified phyloseq object to taxmap 
## subset between sites
Col<-subset_samples(phyb.rar, GIT=="Colon")
Je<-subset_samples(phyb.rar, GIT=="Jejunum")

## convert
ctm_obj <- parse_phyloseq(Col)
jtm_obj <- parse_phyloseq(Je)

# get rid of low counts
ctm_obj$data$tax_data <- zero_low_counts(ctm_obj, data = "otu_table", min_count = 5)
jtm_obj$data$tax_data <- zero_low_counts(jtm_obj, data = "otu_table", min_count = 5)

##Check observations
cno_reads <- rowSums(ctm_obj$data$tax_data[, ctm_obj$data$sample_data$sample_id]) == 0
sum(cno_reads)
jno_reads <- rowSums(jtm_obj$data$tax_data[, jtm_obj$data$sample_data$sample_id]) == 0
sum(jno_reads)


##remove
ctm_obj <- filter_obs(ctm_obj, data = "tax_data", ! cno_reads, drop_taxa = TRUE)
print(ctm_obj)


jtm_obj <- filter_obs(jtm_obj, data = "tax_data", ! jno_reads, drop_taxa = TRUE)
print(jtm_obj)

## calculate abundance
ctm_obj$data$tax_abund <- calc_taxon_abund(ctm_obj, "tax_data",
                                           cols = ctm_obj$data$sample_data$sample_id)

jtm_obj$data$tax_abund <- calc_taxon_abund(jtm_obj, "tax_data",
                                           cols = jtm_obj$data$sample_data$sample_id)


## counts per sample type
ctm_obj$data$tax_occ <- calc_n_samples(ctm_obj, "tax_abund", groups = ctm_obj$data$sample_data$GIT, cols = ctm_obj$data$sample_data$sample_id)

jtm_obj$data$tax_occ <- calc_n_samples(jtm_obj, "tax_abund", groups = jtm_obj$data$sample_data$GIT, cols = jtm_obj$data$sample_data$sample_id)


##plot
heat_tree(jtm_obj, 
          node_label = ifelse(n_obs == 0,"", taxon_names),
          node_size = n_obs,
          node_color = Jejunum, 
          node_size_axis_label = "ASV count",
          node_color_axis_label = "Samples with reads",
          node_color_range=c("grey74","khaki1","green", "deepskyblue" ),
          node_size_range = c(0.005, 0.03),
          edge_size_range=c(0.0005, 0.013),
          edge_label_size_range = c(10, 14),
          node_label_max = 200, edge_label_max =200,
          initial_layout = "re", layout = "da",
          output_file = "jejunum_heat_tree.pdf")


#### Colon
heat_tree(ctm_obj, 
          node_label = ifelse(n_obs == 0,"", taxon_names),
          node_size = n_obs,
          node_color = Colon, 
          node_size_axis_label = "ASV count",
          node_color_axis_label = "Samples with reads",
          node_color_range=c("grey74","khaki1","green", "deepskyblue" ),
          node_size_range = c(0.005, 0.03),
          edge_size_range=c(0.0005, 0.013),
          edge_label_size_range = c(10, 14),
          initial_layout = "re", layout = "da",
          output_file = "colon_heat_tree.pdf")

####  weighted Unifrac Visualization ####
weighted<-phyloseq::ordinate(new_obj.rel, "PCoA", "unifrac", weighted=TRUE)


tiff('weighted.tiff', units="in", width=7, height=4, res=300, compression = 'lzw')
phyloseq::plot_ordination(new_obj.rel, weighted, color="GIT", shape="AgeClass")+
  geom_point(size=2)+ scale_color_manual(values=c("cadetblue", "darkgreen"))+scale_shape_manual(values=c(4,2,1))+
  theme( 
    legend.position="bottom",
    axis.ticks = element_blank(),
    axis.text.x = element_text(family="Arial",size=8, color="black"),axis.text.y = element_text(family="Arial",size=8, color="black"),
    legend.text = element_text(size=8, family="Arial"),
    axis.title.x = element_text(family="Arial",size=8), 
    axis.title.y = element_text(size=8, family="Arial"), 
    panel.background = element_blank(), legend.title = element_blank(), 
    plot.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=.5))+ggtitle("Weighted Unifrac")
dev.off()




### BETA PLOTS ####



phyb.rar2<-subset_samples(phyb.rar, AgeClass != "Unknown")

sam.meta2 <- meta(phyb.rar2)
sam.meta2$SampleID <- rownames(sam.meta2)
sam.meta2$AgeClass<-factor(sam.meta2$AgeClass, levels=c("Yearling", "Subadult", "Adult"))

sample_data(phyb.rar2)<- sam.meta2

library(plyr)  ## or dplyr (transform -> mutate)
library(reshape2)
install.packages("remotes")
remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)



year<-subset_samples(phyb.rar2, AgeClass=="Yearling")

suba<-subset_samples(phyb.rar2, AgeClass=="Subadult")
adu<-subset_samples(phyb.rar2, AgeClass=="Adult")

merged<-merge_phyloseq(year,suba, adu)


###GIT
p = merged
m = "wunifrac"
s = "SampleID"
d = "GIT"

# calc distances
wu = phyloseq::distance(p, m)
wu.m = melt(as.matrix(wu))

# remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
sd = sample_data(p) %>%
  select(s, d) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")
wu.sd$Type2<-factor(wu.sd$Type2,levels=c("Jejunum", "Colon"))
wu.sd$Type1<-factor(wu.sd$Type1,levels=c("Jejunum", "Colon"))

dfgt<-wu.sd

dfgt$Compared <- apply(dfgt, 1, function(x) grepl(x['Type1'], x['Type2'],fixed=TRUE))
dfgt$Compared <- as.character(dfgt$Compared)
dfgt$Compared[dfgt$Compared == "TRUE"] <- "Within"
dfgt$Compared[dfgt$Compared == "FALSE"] <- "Between"
dfgt$Try[dfgt$Type1 == "Jejunum"& dfgt$Type2=="Jejunum"] <- "Jejunum"
dfgt$Try[dfgt$Type1 == "Jejunum"& dfgt$Type2=="Colon"] <- "Jejunum vs Colon"
dfgt$Try[dfgt$Type1 == "Colon"& dfgt$Type2=="Colon"] <- "Colon"
dfgt$Try[dfgt$Type1 == "Colon"& dfgt$Type2=="Jejunum"] <- "Jejunum vs Colon"



tiff('GITWeightedFinal.tiff', units="in", width=6, height=5, res=1200, compression = 'lzw')
ggplot(dfgt, aes(x = Try, y = value, fill=Try)) +
  theme_bw() +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  ggpattern::geom_boxplot_pattern(pattern=c("none","none", "stripe")) +
  scale_color_identity() + scale_fill_manual(values=c("gray", "white", "gray28"))+
  theme(plot.title=element_text(family= "Arial", hjust=.5),legend.position="none",panel.grid=element_blank(),legend.title=element_blank(), legend.text=element_blank(),strip.text.x=element_text(color="black", family="Arial", size=12),strip.background=element_rect(color="black", fill="white"),
        axis.text.y=element_text(family="Arial", color="black"),axis.title.y=element_text(family="Arial"),axis.text.x=element_text( hjust = .5, vjust=.8,family="Arial", color="black")) + 
  ggtitle(paste0(""))+ylab("Weighted UniFrac Distance")+xlab("")+
  scale_y_continuous(sec.axis = dup_axis(name=NULL))
dev.off()



###AGE

d = "AgeClass"

# calc distances
wu = phyloseq::distance(p, m)
wu.m = melt(as.matrix(wu))

# remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
sd = sample_data(p) %>%
  select(s, d) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")
#wu.sd$Type2<-factor(wu.sd$Type2,levels=c("Jejunum", "Colon"))
#wu.sd$Type1<-factor(wu.sd$Type1,levels=c("Jejunum", "Colon"))

df<-wu.sd

df$Compared <- apply(df, 1, function(x) grepl(x['Type1'], x['Type2'],fixed=TRUE))
df$Compared <- as.character(df$Compared)
df$Compared[df$Compared == "TRUE"] <- "Within"
df$Compared[df$Compared == "FALSE"] <- "Between"
df$Try[df$Type1 == "Yearling"& df$Type2=="Yearling"] <- "Yearling"
df$Try[df$Type1 == "Yearling"& df$Type2=="Subadult"] <- "Yearling vs Subadult"
df$Try[df$Type1 == "Yearling"& df$Type2=="Adult"] <- "Yearling vs Adult"
df$Try[df$Type1 == "Subadult"& df$Type2=="Adult"] <- "Subadult vs Adult"
df$Try[df$Type1 == "Subadult"& df$Type2=="Yearling"] <- "Yearling vs Subadult"
df$Try[df$Type1 == "Subadult"& df$Type2=="Subadult"] <- "Subadult"
df$Try[df$Type1 == "Adult"& df$Type2=="Subadult"] <- "Subadult vs Adult"
df$Try[df$Type1 == "Adult"& df$Type2=="Adult"] <- "Adult"
df$Try[df$Type1 == "Adult"& df$Type2=="Yearling"] <- "Yearling vs Adult"

# plot

df$Try<-factor(df$Try,levels=c("Yearling", "Subadult", "Adult", "Yearling vs Subadult", "Yearling vs Adult", "Subadult vs Adult"))
df$Wrap<-"Age"

tiff('AgeWeightedFinal.eps', units="in", width=6, height=5, res=1200, compression = 'lzw')
ggplot(df, aes(x = Try, y = value, fill=Try)) +
  theme_bw() +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  ggpattern::geom_boxplot_pattern(pattern=c("none","none","none", "stripe","stripe","stripe")) +
  scale_color_identity() + scale_fill_manual(values=c("gray", "white", "gray28","gray84", "gray66", "gray45"))+
  theme(plot.title =element_text(family="Arial", hjust = 0.5),legend.position="none",panel.grid=element_blank(),legend.title=element_blank(), legend.text=element_blank(),strip.text.x=element_text(color="black", family="Arial", size=12),strip.background=element_rect(color="black", fill="white"),
        axis.text.y=element_text(family="Arial", color="black"),axis.title.y=element_text(family="Arial"),axis.text.x=element_text(angle = 30, hjust = .8, vjust=.9,family="Arial", color="black")) + 
  ggtitle(paste0())+ylab("Weighted UniFrac Distance")+xlab("")+scale_y_continuous(sec.axis = dup_axis(name=NULL))

dev.off()



###SEX
d = "Sex"

# calc distances
wu = phyloseq::distance(p, m)
wu.m = melt(as.matrix(wu))

# remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
sd = sample_data(p) %>%
  select(s, d) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")
wu.sd$Type2<-factor(wu.sd$Type2,levels=c("Female", "Male"))
wu.sd$Type1<-factor(wu.sd$Type1,levels=c("Female", "Male"))



# plot
dfs<-wu.sd

dfs$Compared <- apply(dfs, 1, function(x) grepl(x['Type1'], x['Type2'],fixed=TRUE))
dfs$Compared <- as.character(dfs$Compared)
dfs$Compared[dfs$Compared == "TRUE"] <- "Within"
dfs$Compared[dfs$Compared == "FALSE"] <- "Between"
dfs$Try[dfs$Type1 == "Female"& dfs$Type2=="Female"] <- "Female"
dfs$Try[dfs$Type1 == "Female"& dfs$Type2=="Male"] <- "Female vs Male"
dfs$Try[dfs$Type1 == "Male"& dfs$Type2=="Male"] <- "Male"
dfs$Try[dfs$Type1 == "Male"& dfs$Type2=="Female"] <- "Female vs Male"


# plot

dfs$Try<-factor(dfs$Try,levels=c("Female", "Male", "Female vs Male"))
dfs$Wrap<-"Sex"

tiff('SexWeightedFinal.tiff', units="in", width=6, height=5, res=1200, compression = 'lzw')
ggplot(dfs, aes(x = Try, y = value, fill=Try)) +
  theme_bw() +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  ggpattern::geom_boxplot_pattern(pattern=c("none","none", "stripe")) +
  scale_color_identity() + scale_fill_manual(values=c("gray", "white", "gray28"))+
  theme(plot.title =element_text(family="Arial", hjust = 0.5),legend.position="none",panel.grid=element_blank(),legend.title=element_blank(), legend.text=element_blank(),strip.text.x=element_text(color="black", family="Arial", size=12),strip.background=element_rect(color="black", fill="white"),
        axis.text.y=element_text(family="Arial", color="black"),axis.title.y=element_text(family="Arial"),axis.text.x=element_text(hjust = .5, vjust=.9,family="Arial", color="black")) + 
  ggtitle(paste0())+ylab("Weighted UniFrac Distance")+xlab("")+scale_y_continuous(sec.axis = dup_axis(name=NULL))

dev.off()


