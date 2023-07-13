############################################################
###Procedure
#0. Package 
#1. Prepare the example data set
#2. Select and process the toxicity data
#3. Select chemicals to be analyzed
#4. Estimate acute SSDs for 4 distributions
#5. Estimate chronic SSDs for 4 distributions
#6. Calculate AICc differences and HC5 ratios
#7. prepare for visualization
#8. Visualization
############################################################



#### 0. Package ----
library(openxlsx)
library(tidyverse)
library(ssdtools)
library(ggplot2)
library(ggExtra)
library(dplyr)
library(EnvStats)
library(ggrepel)
library(cowplot)
library(mousetrap)



#### 1. Prepare the example data set ----
# This dataset "example.xlsx" includes 20,000 test records randomly selected from the "EnviroTox" database only for demonstration.
# All the data used in this study was collected from the "EnviroTox" database and please contact the authors if you like to exactly reproduce our results.

#import data
EnviroTox_test <- read.xlsx("example.xlsx", sheet="test")
EnviroTox_chem <- read.xlsx("example.xlsx", sheet="substance")
EnviroTox_taxo <- read.xlsx("example.xlsx", sheet="taxonomy")



#### 2. Select and process the toxicity data ----
EnviroTox_test_selected <- EnviroTox_test %>%
  filter (Test.statistic=="EC50" & Test.type=="A" | Test.statistic=="LC50" & Test.type=="A"  |
            Test.statistic=="NOEC" & Test.type=="C" | Test.statistic=="NOEL" & Test.type=="C") %>% 
  filter (Effect.is.5X.above.water.solubility =="0") %>%
  mutate (original.CAS = EnviroTox_chem[match(.$CAS, EnviroTox_chem$CAS),"original.CAS"] ) %>%
  mutate_at(vars(Effect.value), as.numeric) %>%
  mutate (Effect.value = replace(.$Effect.value, !is.na(.$Effect.value), .$Effect.value * 10^3) ) %>%  # transform unit (mg/L to ug/L)
  mutate (Unit = replace(Unit, Unit=="mg/L","µg/L")) %>%
  mutate (Substance=EnviroTox_chem[match (.$original.CAS, EnviroTox_chem$original.CAS) ,"Chemical.name"]) %>%
  separate (Substance, into=c("Short_name"),sep=";",extra="drop" ) 

## calculate geometric mean and select chemicals analyzed based on the number of species
EnviroTox_test_selected2 <- aggregate(EnviroTox_test_selected$Effect.value,
                                      by=list(original.CAS = EnviroTox_test_selected$original.CAS,
                                              Test.type=EnviroTox_test_selected$Test.type,
                                              Latin.name=EnviroTox_test_selected$Latin.name),
                                      function(x) geoMean(x) ) %>%
    dplyr::rename(Effect.value=x) %>%
    mutate (Trophic.Level = EnviroTox_taxo[match (.$Latin.name, EnviroTox_taxo$Latin.name) ,"Trophic.Level"] ) %>%
    mutate (Substance=EnviroTox_chem[match (.$original.CAS, EnviroTox_chem$original.CAS) ,"Chemical.name"]) %>%
    separate (Substance, into=c("Short_name"),sep=";",extra="drop" )  %>%
   group_by(original.CAS,Test.type) %>% 
   filter(n()>=10) 


## Organize information of chemicals and the toxicity
EnviroTox_ssd <- aggregate(x=as.numeric(EnviroTox_test_selected2$Effect.value),
                           by=list(original.CAS=EnviroTox_test_selected2$original.CAS, Test.type=EnviroTox_test_selected2$Test.type),
                           FUN=function(x) mean(log10( x ) ) )  %>%
  mutate(sd=aggregate(EnviroTox_test_selected2$Effect.value,
                      by=list(EnviroTox_test_selected2$original.CAS, EnviroTox_test_selected2$Test.type), function(x) sd(log10( x ) ) )[,3]   )   %>%
  dplyr::rename(mean=x) %>%
  mutate(HC5 = qlnorm (0.05, meanlog=log(10^mean), sdlog=log(10^sd) ) ) %>%
  mutate (Substance=EnviroTox_chem[match (.$original.CAS, EnviroTox_chem$original.CAS) ,"Chemical.name"]) %>%
  mutate (No_species = aggregate(EnviroTox_test_selected2$Latin.name,
                                 by=list(EnviroTox_test_selected2$original.CAS,EnviroTox_test_selected2$Test.type), function(x) length(unique(x)) )[,3]) %>%
  mutate (No_trophic=aggregate(EnviroTox_test_selected2$Trophic.Level,
                               by=list(EnviroTox_test_selected2$original.CAS,EnviroTox_test_selected2$Test.type), function(x) length(unique (x)) )[,3]) %>%
  filter(!is.na(sd)) %>%
  mutate(Test.type = replace(Test.type, Test.type=="A", "Acute")) %>%
  mutate(Test.type = replace(Test.type, Test.type=="C", "Chronic"))%>%
  pivot_wider(names_from=Test.type, values_from=c("mean","sd","HC5","No_species","No_trophic")) %>%
  mutate (ConsensusMoA = EnviroTox_chem[match (.$original.CAS, EnviroTox_chem$original.CAS), "Consensus.MOA"] ) %>%
  mutate (ASTER = EnviroTox_chem[match (.$original.CAS, EnviroTox_chem$original.CAS) ,"ASTER"] )

EnviroTox_ssd$ConsensusMoA <- replace (EnviroTox_ssd$ConsensusMoA, which(EnviroTox_ssd$ConsensusMoA=="N"),"Narcotic")
EnviroTox_ssd$ConsensusMoA <- replace (EnviroTox_ssd$ConsensusMoA, which(EnviroTox_ssd$ConsensusMoA=="U"),"Unclassified")
EnviroTox_ssd$ConsensusMoA <- replace (EnviroTox_ssd$ConsensusMoA, which(EnviroTox_ssd$ConsensusMoA=="S"),"Specifically acting")
EnviroTox_ssd$ConsensusMoA <- as.factor(EnviroTox_ssd$ConsensusMoA)


## Calculate bimodality coefficient (BC)

# acute data
BC_A <- EnviroTox_test_selected2 %>%
  filter(Test.type =="A") %>%
  group_by(original.CAS) %>%
  dplyr::summarize(BC = mousetrap::bimodality_coefficient(log10(Effect.value))) %>%
  mutate(Bimodal = ifelse(BC >0.555, "Bimodal","Not bimodal") )

# chronic data
BC_C <- EnviroTox_test_selected2 %>%
  filter(Test.type =="C") %>%
  group_by(original.CAS) %>%
  dplyr::summarize(BC = mousetrap::bimodality_coefficient(log10(Effect.value))) %>%
  mutate(Bimodal = ifelse(BC >0.555, "Bimodal","Not bimodal") )

## BC's criteria: 0.555  (Freeman et al., 2013; Pfister et al. 2013)
BC_CAS_A <- BC_A %>%
  filter(Bimodal == "Not bimodal")
BC_CAS_C <- BC_C %>%
  filter(Bimodal == "Not bimodal")


#### 3. Select chemicals to be analyzed ----
## Get the lists of chemicals to be used SSD estimation

## No ofspecies >= 10 and No of trophic groups >= 3 and "Not bimodal"
EnviroTox_ssd_HH_A <- EnviroTox_ssd %>%
  filter (No_trophic_Acute >= 3  ) %>%
  filter (No_species_Acute >= 10 ) %>%
  filter (original.CAS %in% BC_CAS_A$original.CAS) %>%
  separate (Substance, into=c("Short_name"), sep=";", extra="drop")

EnviroTox_ssd_HH_C <- EnviroTox_ssd %>%
  filter (No_trophic_Chronic >= 3  ) %>%
  filter (No_species_Chronic >= 10 ) %>%
  filter (original.CAS %in% BC_CAS_C$original.CAS) %>%
  separate (Substance, into=c("Short_name"), sep=";", extra="drop")

## Lists of chemicals (CAS) to be examined
StudyChemicals_A <- EnviroTox_ssd_HH_A$original.CAS
StudyChemicals_C <- EnviroTox_ssd_HH_C$original.CAS


## Prepare information of chemicals
EnviroTox_chem_rev <- EnviroTox_chem[,1:4]
d02_A <- left_join(EnviroTox_ssd_HH_A, EnviroTox_chem_rev, by="original.CAS")
d02_C <- left_join(EnviroTox_ssd_HH_C, EnviroTox_chem_rev, by="original.CAS")


#### 4. Estimate acute SSDs for 4 distributions ----

### Estimate SSDs using ssdtools package
temp.res <- data.frame(matrix(-9999, ncol= 1 + 8*4 + 6*4, nrow = length(StudyChemicals_A)))
head(temp.res)

for (i in 1:length(StudyChemicals_A)){
  Temp.data <- EnviroTox_test_selected2 %>% filter(Test.type == "A" & original.CAS==StudyChemicals_A[i]) 
  # fit distributions
  fits <- ssdtools::ssd_fit_dists(Temp.data, left = 'Effect.value',
                                  dists = c('lnorm', 'llogis', 'burrIII3', 'weibull'),
                                  at_boundary_ok=FALSE,computable=TRUE) 
  # plot distributions
  ssdtools::ssd_plot_cdf(fits)
  # goodness of fit table
  ssd_gof(fits) 
  
  set.seed(99)
  hc5 <- ssd_hc(fits, ci = TRUE, nboot = 1000,  average = FALSE, delta = 100) # 1000 bootstrap samples
  print(hc5) 
  
  temp.res[i,1] <- StudyChemicals_A[i]
  temp.res[i,2:33] <- c(ssd_gof(fits) [1,1:8], ssd_gof(fits) [2,1:8], ssd_gof(fits) [3,1:8],
                        ssd_gof(fits) [4,1:8])
  temp.res[i,34:57] <- c(hc5[1,1:6],hc5[2,1:6], hc5[3,1:6], hc5[4,1:6])
}

head(temp.res)
nrow(temp.res)
temp.res_acute_example <- temp.res

### Reshape the result data 

res_aic1 <- dplyr::select(temp.res_acute_example, c(1:9))
res_aic2 <- dplyr::select(temp.res_acute_example, c(1, 10:17))
res_aic3 <- dplyr::select(temp.res_acute_example, c(1, 18:25))
res_aic4 <- dplyr::select(temp.res_acute_example, c(1, 26:33))
colnames(res_aic2) <- colnames(res_aic1)
colnames(res_aic3) <- colnames(res_aic1)
colnames(res_aic4) <- colnames(res_aic1)

res_aic <- bind_rows(res_aic1, res_aic2, res_aic3, res_aic4)

res_hc1 <- dplyr::select(temp.res_acute_example, c(1, 34:39))
res_hc2 <- dplyr::select(temp.res_acute_example, c(1, 40:45))
res_hc3 <- dplyr::select(temp.res_acute_example, c(1, 46:51))
res_hc4 <- dplyr::select(temp.res_acute_example, c(1, 52:57))
colnames(res_hc2) <- colnames(res_hc1)
colnames(res_hc3) <- colnames(res_hc1)
colnames(res_hc4) <- colnames(res_hc1)

res_hc <- bind_rows(res_hc1, res_hc2, res_hc3, res_hc4)

colnames(res_aic) <- c("original.CAS", "dist", "ad", "ks", "cvm", "aic", "aicc", "bic", "delta")
colnames(res_hc) <- c("original.CAS","dist", "percent", "est", "se", "lcl", "ucl")

res_aic_wider <- res_aic%>%
  filter(!is.na(dist)) %>%
  group_by(original.CAS) %>%
  pivot_wider(id_cols = original.CAS, names_from=dist, values_from = c( "ad", "ks", "cvm", "aic", "aicc", "bic", "delta"))

res_hc_wider <- res_hc%>% 
  filter(!is.na(dist)) %>%
  pivot_wider(id_cols = original.CAS, names_from=dist, values_from = c( "percent", "est", "se", "lcl", "ucl"), values_fill = NA)

temp.res_acute_clean <- left_join(res_aic_wider, res_hc_wider, by = "original.CAS")

# Add the information of each compound such as MoAs
temp.res_acute_clean_2 <- left_join(temp.res_acute_clean,d02_A, by = "original.CAS")

# Remove rows that have NA
temp.res_acute_clean_3 <- temp.res_acute_clean_2 %>%
  filter(!is.na(aicc_lnorm) & !is.na(aicc_llogis) & !is.na(aicc_burrIII3) & !is.na(aicc_weibull)) %>%
  filter(!is.na(est_lnorm) & !is.na(est_llogis) & !is.na(est_burrIII3) & !is.na(est_weibull)) 


#### 5. Estimate chronic SSDs for 4 distributions ----

### Estimate SSDs using ssdtools package
temp.res <- data.frame(matrix(-9999, ncol= 1 + 8*4 + 6*4, nrow = length(StudyChemicals_C)))
head(temp.res)

# for logn, logl. burr, weibull 
for (i in 1:length(StudyChemicals_C)){
  Temp.data <-  EnviroTox_test_selected2 %>% filter(Test.type == "A" & original.CAS==StudyChemicals_C[i]) 
  fits <- ssdtools::ssd_fit_dists(Temp.data, left = 'Effect.value',
                                  dists = c('lnorm', 'llogis', 'burrIII3', 'weibull'),
                                  at_boundary_ok=FALSE,computable=TRUE) 
  # plot distributions
  ssdtools::ssd_plot_cdf(fits)
  # goodness of fit table
  ssd_gof(fits)
  
  set.seed(99)
  hc5 <- ssd_hc(fits, ci = TRUE, nboot = 1000,  average = FALSE, delta = 100) # 1000 bootstrap samples
  print(hc5) 

    temp.res[i,1] <- StudyChemicals_C[i]
  temp.res[i,2:33] <- c(ssd_gof(fits) [1,1:8], ssd_gof(fits) [2,1:8], ssd_gof(fits) [3,1:8],
                        ssd_gof(fits) [4,1:8])
  temp.res[i,34:57] <- c(hc5[1,1:6],hc5[2,1:6], hc5[3,1:6], hc5[4,1:6])
}

head(temp.res)
nrow(temp.res)
temp.res_chronic_example <- temp.res

### Reshape the result data 

res_aic1 <- dplyr::select(temp.res_chronic_example, c(1:9))
res_aic2 <- dplyr::select(temp.res_chronic_example, c(1, 10:17))
res_aic3 <- dplyr::select(temp.res_chronic_example, c(1, 18:25))
res_aic4 <- dplyr::select(temp.res_chronic_example, c(1, 26:33))
colnames(res_aic2) <- colnames(res_aic1)
colnames(res_aic3) <- colnames(res_aic1)
colnames(res_aic4) <- colnames(res_aic1)

res_aic <- bind_rows(res_aic1, res_aic2, res_aic3, res_aic4)

res_hc1 <- dplyr::select(temp.res_chronic_example, c(1, 34:39))
res_hc2 <- dplyr::select(temp.res_chronic_example, c(1, 40:45))
res_hc3 <- dplyr::select(temp.res_chronic_example, c(1, 46:51))
res_hc4 <- dplyr::select(temp.res_chronic_example, c(1, 52:57))

colnames(res_hc2) <- colnames(res_hc1)
colnames(res_hc3) <- colnames(res_hc1)
colnames(res_hc4) <- colnames(res_hc1)

res_hc <- bind_rows(res_hc1, res_hc2, res_hc3, res_hc4)

colnames(res_aic) <- c("original.CAS", "dist", "ad", "ks", "cvm", "aic", "aicc", "bic", "delta")
colnames(res_hc) <- c("original.CAS","dist", "percent", "est", "se", "lcl", "ucl")

res_aic_wider <- res_aic%>%
  filter(!is.na(dist)) %>%
  group_by(original.CAS) %>%
  pivot_wider(id_cols = original.CAS, names_from=dist, values_from = c( "ad", "ks", "cvm", "aic", "aicc", "bic", "delta"))

res_hc_wider <- res_hc%>% 
  filter(!is.na(dist)) %>%
  pivot_wider(id_cols = original.CAS, names_from=dist, values_from = c( "percent", "est", "se", "lcl", "ucl"), values_fill = NA)

temp.res_chronic_clean <- left_join(res_aic_wider, res_hc_wider, by = "original.CAS")

# Add the information of each compound such as MoAs
temp.res_chronic_clean_2 <- left_join(temp.res_chronic_clean, d02_C, by = "original.CAS")

# Remove rows that have NA
temp.res_chronic_clean_3 <- temp.res_chronic_clean_2 %>%
  filter(!is.na(aicc_lnorm) & !is.na(aicc_llogis) & !is.na(aicc_burrIII3) & !is.na(aicc_weibull)) %>%
  filter(!is.na(est_lnorm) & !is.na(est_llogis) & !is.na(est_burrIII3) & !is.na(est_weibull)) 


#### 6. Calculate AICc differences and HC5 ratios----

### Acute SSDs

# Process the data (acute SSD)
res_acute <- temp.res_acute_clean_3 %>%
  mutate(BestModel = case_when(delta_lnorm == 0 ~ "lnorm",
                               delta_llogis == 0 ~ "llogis",
                               delta_burrIII3 == 0 ~ "burrIII3",
                               delta_weibull == 0 ~ "weibull")) %>%
  mutate(AICcdif_llogis = aicc_llogis - aicc_lnorm,
         AICcdif_burrIII3 = aicc_burrIII3 - aicc_lnorm,
         AICcdif_weibull = aicc_weibull - aicc_lnorm) %>%
  mutate(HC5rat_llogis = ifelse(est_llogis == "NULL", NA,  
                                   ifelse(est_lnorm == "NULL", NA, as.numeric(unlist(est_llogis))/as.numeric(unlist(est_lnorm)))),
         HC5rat_burrIII3 = ifelse(est_burrIII3 == "NULL", NA,  
                                   ifelse(est_lnorm == "NULL", NA, as.numeric(unlist(est_burrIII3))/as.numeric(unlist(est_lnorm)))),
         HC5rat_weibull = ifelse(est_weibull == "NULL", NA,  
                                  ifelse(est_lnorm == "NULL", NA, as.numeric(unlist(est_weibull))/as.numeric(unlist(est_lnorm)))) ) 



### Chronic SSDs

# Process the data (chronic SSD)
res_chronic <- temp.res_chronic_clean_3 %>%
  mutate(BestModel = case_when(delta_lnorm == 0 ~ "lnorm",
                               delta_llogis == 0 ~ "llogis",
                               delta_burrIII3 == 0 ~ "burrIII3",
                               delta_weibull == 0 ~ "weibull")) %>%
  mutate(AICcdif_llogis = aicc_llogis - aicc_lnorm,
         AICcdif_burrIII3 = aicc_burrIII3 - aicc_lnorm,
         AICcdif_weibull = aicc_weibull - aicc_lnorm) %>%
  mutate(HC5rat_llogis = ifelse(est_llogis == "NULL", NA,  
                                ifelse(est_lnorm == "NULL", NA, as.numeric(unlist(est_llogis))/as.numeric(unlist(est_lnorm)))),
         HC5rat_burrIII3 = ifelse(est_burrIII3 == "NULL", NA,  
                                  ifelse(est_lnorm == "NULL", NA, as.numeric(unlist(est_burrIII3))/as.numeric(unlist(est_lnorm)))),
         HC5rat_weibull = ifelse(est_weibull == "NULL", NA,  
                                 ifelse(est_lnorm == "NULL", NA, as.numeric(unlist(est_weibull))/as.numeric(unlist(est_lnorm)))) ) 


#### 7. prepare for visualization ----

### Acute SSDs
res_acute_selected1 <- res_acute %>%
  select(original.CAS, AICcdif_llogis, AICcdif_burrIII3, AICcdif_weibull) %>%
  pivot_longer(cols = c(AICcdif_llogis, AICcdif_burrIII3, AICcdif_weibull), 
               names_to = "dist", names_prefix = "AICcdif_", values_to = "AICcdiff")

res_acute_selected2 <- res_acute %>%
  select(original.CAS, HC5rat_llogis, HC5rat_burrIII3, HC5rat_weibull) %>%
  pivot_longer(cols = c(HC5rat_llogis, HC5rat_burrIII3, HC5rat_weibull), 
               names_to = "dist", names_prefix = "HC5rat_", values_to = "HC5_ratio")

res_acute_visualization <- left_join(res_acute_selected1, res_acute_selected2, 
                                     by = c("original.CAS" = "original.CAS", "dist" = "dist"))  

res_acute_visualization2 <- left_join(res_acute_visualization, d02_A, by = "original.CAS")

res_acute_visualization3 <- res_acute_visualization2 %>%
  mutate(dist2 = case_when(dist == "llogis" ~ "Log-logistic",
                           dist == "burrIII3" ~ "Burr type III",
                           dist == "weibull" ~ "Weibull")) %>%
  mutate(Label1 = ifelse( abs(AICcdiff) > 20, Short_name ,"") ) %>%
  mutate(Label2 = ifelse( abs(log10(HC5_ratio)) > 1, Short_name ,""))



### Chronic SSDs
res_chronic_selected1 <- res_chronic %>%
  select(original.CAS, AICcdif_llogis, AICcdif_burrIII3, AICcdif_weibull) %>%
  pivot_longer(cols = c(AICcdif_llogis, AICcdif_burrIII3, AICcdif_weibull), 
               names_to = "dist", names_prefix = "AICcdif_", values_to = "AICcdiff")

res_chronic_selected2 <- res_chronic %>%
  select(original.CAS, HC5rat_llogis, HC5rat_burrIII3, HC5rat_weibull) %>%
  pivot_longer(cols = c(HC5rat_llogis, HC5rat_burrIII3, HC5rat_weibull), 
               names_to = "dist", names_prefix = "HC5rat_", values_to = "HC5_ratio")

res_chronic_visualization <- left_join(res_chronic_selected1, res_chronic_selected2, 
                                     by = c("original.CAS" = "original.CAS", "dist" = "dist"))  

res_chronic_visualization2 <- left_join(res_chronic_visualization, d02_C, by = "original.CAS")

res_chronic_visualization3 <- res_chronic_visualization2 %>%
  mutate(dist2 = case_when(dist == "llogis" ~ "Log-logistic",
                           dist == "burrIII3" ~ "Burr type III",
                           dist == "weibull" ~ "Weibull")) %>%
  mutate(Label1 = ifelse( abs(AICcdiff) > 20, Short_name ,"") ) %>%
  mutate(Label2 = ifelse( abs(log10(HC5_ratio)) > 1, Short_name ,""))



#### 8. Visualization ----

### Figure 1. Comparison of AICc differences

## Acute SSDs
Fig_AICdiff_Acute <- res_acute_visualization3 %>%
  mutate(across(dist2, factor, levels=c("Log-logistic",  "Burr type III", "Weibull"))) %>%
  ggplot(aes(dist2, AICcdiff, alpha =0.7)) +
  geom_boxplot()+
  geom_line(aes(group = original.CAS), size=0.5, color='gray', alpha=0.6)+
  geom_point(aes(color=ConsensusMoA,shape=ConsensusMoA, group=original.CAS),size=2.5, alpha = 0.7)+
  scale_fill_manual(values = c("Narcotic" = "olivedrab1", "Specifically acting" = "lightskyblue",
                               "Unclassified" = "palevioletred1"  ))+
  scale_shape_manual(values = c(8 , 15, 16)) +
  geom_text_repel(aes(label=Label1), size = 5, color = "darkblue",
                  min.segment.length = unit(0.01, "lines")) + #paper5, poster7
  geom_abline(slope=0, intercept=10, size=0.5, lty="dotted")+
  geom_abline(slope=0, intercept=-10, size=0.5, lty="dotted")+
  geom_abline(slope=0, intercept=20, size=0.5, lty="dashed")+
  geom_abline(slope=0, intercept=-20, size=0.5, lty="dashed")+
  theme_bw(base_size=20) + 
  theme(axis.text = element_text(color="black"), panel.grid=element_blank(),
        legend.position = 'bottom')+
  guides(alpha = "none") +
  ylim(-30,30)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  labs(x="Distribution", y="AICc difference")

## Chronic SSDs
Fig_AICdiff_Chronic <- res_chronic_visualization3 %>%
  mutate(across(dist2, factor, levels=c("Log-logistic",  "Burr type III", "Weibull"))) %>%
  ggplot(aes(dist2, AICcdiff, alpha =0.7)) +
  geom_boxplot()+
  geom_line(aes(group = original.CAS), size=0.5, color='gray', alpha=0.6)+
  geom_point(aes(color=ConsensusMoA,shape=ConsensusMoA, group=original.CAS),size=2.5, alpha = 0.7)+
  scale_fill_manual(values = c("Narcotic" = "olivedrab1", "Specifically acting" = "lightskyblue",
                               "Unclassified" = "palevioletred1"  ))+
  scale_shape_manual(values = c(8 , 15, 16)) +
  geom_text_repel(aes(label=Label1)) +
  geom_abline(slope=0, intercept=10, size=0.5, lty="dotted")+
  geom_abline(slope=0, intercept=-10, size=0.5, lty="dotted")+
  geom_abline(slope=0, intercept=20, size=0.5, lty="dashed")+
  geom_abline(slope=0, intercept=-20, size=0.5, lty="dashed")+
  theme_bw(base_size=20) +
  theme(axis.text = element_text(color="black"), panel.grid=element_blank(),
        legend.position = 'bottom')+
  guides(alpha = "none") +
  ylim(-30,30)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  labs(x="Distribution", y="AICc difference")

## Combine the two figures
legend_fig1 <- cowplot::get_legend(Fig_AICdiff_Acute)

fig1a <- Fig_AICdiff_Acute + theme(legend.position = "none")
fig1b <- Fig_AICdiff_Chronic + theme(legend.position = "none")
fig1ab <- cowplot::plot_grid(fig1a, fig1b, nrow=1, labels="AUTO")

# Figure 1
Figure1 <- cowplot::plot_grid(fig1ab,legend_fig1,nrow = 2,  rel_heights = c(6, 1))


### Figure 2. Comparison of HC5 ratios

## Acute SSDs
Fig_HC5diff_Acute <- res_acute_visualization3 %>%
  mutate(across(dist2, factor, levels=c("Log-logistic",  "Burr type III", "Weibull"))) %>%
  ggplot(aes(dist2, log10(HC5_ratio), alpha=0.7)) +
  geom_boxplot()+
  geom_line(aes(group = original.CAS), size=0.5, color='gray', alpha=0.6)+
  geom_point(aes(color=ConsensusMoA,shape=ConsensusMoA, group=original.CAS),size=3, alpha = 0.6)+
  scale_fill_manual(values = c("Narcotic" = "olivedrab1", "Specifically acting" = "lightskyblue",
                               "Unclassified" = "palevioletred1"  ))+
  scale_shape_manual(values = c(8 , 15, 16)) +
  geom_text_repel(aes(label=Label2), size = 5, color = "darkblue",
                  min.segment.length = unit(0.01, "lines")) + # direction = "x"
  geom_abline(slope=0, intercept=log10(2), size=0.5, lty="dotted")+
  geom_abline(slope=0, intercept=-log10(2), size=0.5, lty="dotted")+
  geom_abline(slope=0, intercept=1, size=0.5, lty="dashed")+
  geom_abline(slope=0, intercept=-1, size=0.5, lty="dashed")+
  theme_bw(base_size=20) +
  theme(axis.text = element_text(color="black"), panel.grid=element_blank(),
        legend.position = 'none')+
  guides(alpha = "none") +
  ylim(-2,0.75)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  labs(x="Distribution")+
  ylab(expression(paste(Log[10]," HC5 ratio")) )



Fig_HC5diff_Chronic <- res_chronic_visualization3 %>%
  mutate(across(dist2, factor, levels=c("Log-logistic",  "Burr type III", "Weibull"))) %>%
  ggplot(aes(dist2, log10(HC5_ratio), alpha=0.7)) +
  geom_boxplot()+
  geom_line(aes(group = original.CAS), size=0.5, color='gray', alpha=0.6)+
  geom_point(aes(color=ConsensusMoA,shape=ConsensusMoA, group=original.CAS),size=3, alpha = 0.6)+
  scale_fill_manual(values = c("Narcotic" = "olivedrab1", "Specifically acting" = "lightskyblue",
                               "Unclassified" = "palevioletred1"  ))+
  scale_shape_manual(values = c(8 , 15, 16)) +
  geom_text_repel(aes(label=Label2), size = 5, color = "darkblue",
                  min.segment.length = unit(0.01, "lines") ) + #direction = "y"
  geom_abline(slope=0, intercept=log10(2), size=0.5, lty="dotted")+
  geom_abline(slope=0, intercept=-log10(2), size=0.5, lty="dotted")+
  geom_abline(slope=0, intercept=1, size=0.5, lty="dashed")+
  geom_abline(slope=0, intercept=-1, size=0.5, lty="dashed")+
  theme_bw(base_size=20) +
  theme(axis.text = element_text(color="black"), panel.grid=element_blank(),
        legend.position = 'bottom')+
  guides(alpha = "none") +
  ylim(-2,0.75)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  labs(x="Distribution")+
  ylab(expression(paste(Log[10]," HC5 ratio")) )


## Combine the two figures
fig2a <- Fig_HC5diff_Acute + theme(legend.position = "none")
fig2b <- Fig_HC5diff_Chronic + theme(legend.position = "none")
fig2ab <- cowplot::plot_grid(fig2a, fig2b, nrow=1, labels="AUTO", label_size = 20)
legend_fig2 <- cowplot::get_legend(Fig_HC5diff_Chronic)

# Figure 2
Figure2 <- cowplot::plot_grid(fig2ab,legend_fig2,nrow = 2,  rel_heights = c(7, 1))



