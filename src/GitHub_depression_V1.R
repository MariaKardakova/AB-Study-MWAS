library(dplyr)
library(tidyverse)
library(qwraps2)
library(ANCOMBC)
library(DT)
library(magrittr)
library(microbiome)
library(vegan)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(Maaslin2)
library(phyloseqGraphTest)
library("phyloseq")
library("ggplot2")
library(meta)
library(jsonlite)
library(tidyr)
library(visdat)
library(janitor)
library(ggrepel)


#INPUTE_DATA

metadata_full <- read_csv("/Users/metadata.csv") 
bacterialist <- read.csv("/Users/microbiome_analysis.csv")

# matching number of individuals in microbiome and metadata data sets 
metadata_full <- metadata_full %>% filter(barcode%in% bacterialist$barcode)

# PEOPLE WITH Depression
Depression_data <- metadata_full %>%  filter(qMedicalHistoryDepression == 2, BMI<70, BMI>10, qWaistCircumferenceCm >30, age < 95) 

## PEOPLE WITHOUT Depression  (control)
noDepression_data <- metadata_full %>%  filter(qMedicalHistoryDepression == 1, BMI<70, BMI>10, qWaistCircumferenceCm >30, age <95)

# _______________----
## MERGE DATA ---- 
metadata_full <- bind_rows(Depression_data, noDepression_data) #6,911

# Bacteria list for MAASLIN2
bacterialist.maaslin <- bacterialist 

# Bacteria list for ANCOM-BC2 
bacterialist.ancom <- (bacterialist*100) # absolute bacterial abundance 10000 per sample for ANCOMBC2 analysis

# DATA SET 1 - DISCOVERY SET and DATA SET 2 - REPLICATION SET
#Split into two data sets 157 weeks - is a date before 31.12.2020 (157 weeks or less since the date of calculation of the whole dataset which is 2024-01-07) and after 31.12.2020 - two data loads from data provider

# _______________----
# DS and RS ----
metadata_RS <- metadata_full %>% select(cohort, weeks.since.testing.microbiome, everything()) %>% filter(weeks.since.testing.microbiome <157) %>%  arrange(barcode)  #Replication depression dataset 2,452
metadata_DS <- metadata_full %>% select(cohort, weeks.since.testing.microbiome, everything()) %>% filter(weeks.since.testing.microbiome >157) %>%  arrange(barcode)  #Discovery depression dataset 4,456

# _______________----
# BACTERIA NAMES ----
bacteriadata_names <- as.data.frame(colnames(bacterialist.ancom)) %>% 
  rename(bacterianame = 'colnames(bacterialist.ancom)') %>% 
  separate(bacterianame, sep="\\.", into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), remove = FALSE) %>% 
  mutate(Genus = paste(Genus, Species)) %>% 
  mutate(FamilyGenus = paste(Family, Genus, Species)) %>% 
  select(bacterianame, Genus, FamilyGenus) #975 bacteria associated names on the genus level

# _______________----
# Metadata For PCA ----

metadata_PCA <- metadata_full %>% mutate(qMedicalHistoryDepression = ifelse(qMedicalHistoryDepression == '1', "NoDepression", "Depression"))   #ifelse works good with 2 answer options "if gender = 1, else 2"
metadata_PCA <- metadata_PCA %>% mutate(qMedicalHistoryDepression = factor(qMedicalHistoryDepression, levels = c("NoDepression", "Depression")))

metadata_PCA <- metadata_PCA %>% mutate(Obesity.status = ifelse(Obesity.status == '1', "NotObese", "Obese"))   
metadata_PCA <- metadata_PCA %>% mutate(Obesity.status = factor(Obesity.status, levels = c("NotObese", "Obese")))

metadata_PCA <- metadata_PCA %>% mutate(qMedicalHistoryMigraine = ifelse(qMedicalHistoryMigraine == '1', "NoMigraine", "Migraine"))  
metadata_PCA <- metadata_PCA %>% mutate(qMedicalHistoryMigraine = factor(qMedicalHistoryMigraine, levels = c("NoMigraine", "Migraine"))) # NoMigraine 6391, Migraine 520

metadata_PCA <- metadata_PCA %>% mutate(qPArecreationalVigorousPresent = ifelse(qPArecreationalVigorousPresent == '1', "NoPhysicalActivity", "Active_RecreationalVigorous"))   # NoPhysicalActivity 2650, Active_RecreationalVigorous 2081
metadata_PCA <- metadata_PCA %>% mutate(qPArecreationalVigorousPresent = factor(qPArecreationalVigorousPresent, levels = c("NoPhysicalActivity", "Active_RecreationalVigorous"))) 

metadata_PCA <- metadata_PCA %>% mutate(qGender = ifelse(qGender == '1', "FEMALE", "MALE"))   
metadata_PCA <- metadata_PCA %>% arrange(barcode)
metadata_PCA <- metadata_PCA %>% mutate(qDrugHistoryProbioticsType_lactobacillus = ifelse(qDrugHistoryProbioticsType_lactobacillus == '1', "No probiotics.lacto", "probiotics.lacto"))  
metadata_PCA <- metadata_PCA %>% mutate(cohort = ifelse(cohort =="0", "DiscoverySet", "ReplicationSet")) %>% filter(!is.na(cohort)) # DiscoverySet 4459 ReplicationSet 2452
metadata_PCA$enterotype[metadata_PCA$enterotype == 0] <- "Ent2" #1522
metadata_PCA$enterotype[metadata_PCA$enterotype == 1] <- "Ent1" #3650
metadata_PCA$enterotype[metadata_PCA$enterotype == 2] <- "Ent3" #1739
metadata_PCA %>% count(enterotype)
metadata_PCA %>% count(reg)

metadata_PCA$weight.status[metadata_PCA$weight.status == 0] <- "Underweight" #300
metadata_PCA$weight.status[metadata_PCA$weight.status == 1] <- "HealthyWeight" #4029
metadata_PCA$weight.status[metadata_PCA$weight.status == 2] <- "Overweight" #1775
metadata_PCA$weight.status[metadata_PCA$weight.status == 3] <- "Obese" #807

metadata_PCA$qSmokingStatus[metadata_PCA$qSmokingStatus == 0] <- "NotSmoking" #3902
metadata_PCA$qSmokingStatus[metadata_PCA$qSmokingStatus == 1] <- "QuitSmoking" #2495
metadata_PCA$qSmokingStatus[metadata_PCA$qSmokingStatus == 2] <- "Smoking" #514


metadata_PCA <- metadata_PCA
bacterialist_PCA <- bacterialist %>% 
  rownames_to_column("barcode") %>% 
  filter(barcode%in%metadata_PCA$barcode) %>% 
  column_to_rownames("barcode")
metadata_PCA %>% count(qPArecreationalVigorousPresent)

# _______________----
## Stats for desriptive ----

# DS, population
metadata_DS %>% count(qMedicalHistoryDepression) %>% mutate(freq = (n / sum(n)) * 100) # DS, 3912 controls, 544 cases
metadata_DS %>% filter(qMedicalHistoryDepression ==1) %>% count(qGender) %>% mutate(freq = (n / sum(n)) * 100) 
metadata_DS %>% filter(qMedicalHistoryDepression ==2) %>% count(qGender) %>% mutate(freq = (n / sum(n)) * 100) 
metadata_DS %>% filter(qMedicalHistoryDepression ==1) %>% count(reg) %>% mutate(freq = (n / sum(n)) * 100) 
metadata_DS %>% filter(qMedicalHistoryDepression ==2) %>% count(reg) %>% mutate(freq = (n / sum(n)) * 100) 
# RS, population
metadata_RS %>% count(qMedicalHistoryDepression) %>% mutate(freq = (n / sum(n)) * 100) # RS, 2167 controls, 285 cases
metadata_RS %>% filter(qMedicalHistoryDepression ==1) %>% count(qGender) %>% mutate(freq = (n / sum(n)) * 100)
metadata_RS %>% filter(qMedicalHistoryDepression ==2) %>% count(qGender) %>% mutate(freq = (n / sum(n)) * 100) 
metadata_RS %>% filter(qMedicalHistoryDepression ==1) %>% count(reg) %>% mutate(freq = (n / sum(n)) * 100) 
metadata_RS %>% filter(qMedicalHistoryDepression ==2) %>% count(reg) %>% mutate(freq = (n / sum(n)) * 100)

# DS, age, cases
metadata_DS %>% filter(age>18) %>%filter(qMedicalHistoryDepression == 2) %>% filter(qGender==1) %>%  summarise(mean_age = mean(age), sd_age = sd(age))  #DS, male, mean age 
metadata_DS %>% filter(age>18) %>%filter(qMedicalHistoryDepression == 2) %>% filter(qGender==0) %>%  summarise(mean_age = mean(age), sd_age = sd(age))  #DS, female , mean age 

# DS, age, controls
metadata_DS %>% filter(age>18) %>%filter(qMedicalHistoryDepression == 1) %>% filter(qGender==1) %>%  summarise(mean_age = mean(age), sd_age = sd(age))  #DS, male, mean age 
metadata_DS %>% filter(age>18) %>%filter(qMedicalHistoryDepression == 1) %>% filter(qGender==0) %>%  summarise(mean_age = mean(age), sd_age = sd(age))  #DS, female , mean age 

# RS, age, cases
metadata_RS %>% filter(age>18) %>%filter(qMedicalHistoryDepression == 2) %>% filter(qGender==1) %>%  summarise(mean_age = mean(age), sd_age = sd(age))  #RS, male, mean age 
metadata_RS %>% filter(age>18) %>%filter(qMedicalHistoryDepression == 2) %>% filter(qGender==0) %>%  summarise(mean_age = mean(age), sd_age = sd(age))  #RS, female , mean age 

# RS, age, controls
metadata_RS %>% filter(age>18) %>%filter(qMedicalHistoryDepression == 1) %>% filter(qGender==1) %>%  summarise(mean_age = mean(age), sd_age = sd(age))  #RS, male, mean age 
metadata_RS %>% filter(age>18) %>%filter(qMedicalHistoryDepression == 1) %>% filter(qGender==0) %>%  summarise(mean_age = mean(age), sd_age = sd(age))  #RS, female , mean age 

metadata_DS %>% filter(age>18) %>% filter(qGender==0) %>% pull(age) %>% mean() #female, mean age 41.96233
metadata_RS %>% filter(age>18) %>% filter(qGender==0) %>% pull(age) %>% mean() #female, mean age 43.1786

metadata_DS %>% filter(age>18) %>% filter(qGender==1) %>% pull(age) %>% mean() #male, mean age 41.70152
metadata_RS %>% filter(age>18) %>% filter(qGender==1) %>% pull(age) %>% mean() #male, mean age 41.91492

metadata_DS %>% filter(qGender==0) %>% pull(BMI) %>% mean() #female, mean BMI 23.85919
metadata_RS %>% filter(qGender==0) %>% pull(BMI) %>% mean() #female, mean BMI 24.08314

metadata_DS %>% filter(qGender==1) %>% pull(BMI) %>% mean() #male, mean BMI 25.47665
metadata_RS %>% filter(qGender==1) %>% pull(BMI) %>% mean() #male, mean BMI 25.28928

metadata_DS %>% filter(qGender==0) %>% pull(age) %>% mean() #female, mean BMI 41.88494
metadata_RS %>% filter(qGender==0) %>% pull(age) %>% mean() #female, mean BMI 43.12817

metadata_DS %>% filter(qGender==1) %>% pull(age) %>% mean() #male, mean BMI 41.63517
metadata_RS %>% filter(qGender==1) %>% pull(age) %>% mean() #male, mean BMI 41.86478

#gender
metadata_DS %>% count(qGender) %>% mutate(freq = (n / sum(n)) * 100) # 2477 female (55.6%), 1979 male (44.4%)
metadata_RS %>% count(qGender) %>% mutate(freq = (n / sum(n)) * 100) # 1498 female (61.1%), 954 male (38.9%)

#depression prevalence
metadata_DS %>% count(qMedicalHistoryDepression) %>% mutate(freq = (n / sum(n)) * 100) # 544 Yes (12.2%), 3912 No (87.8%)
metadata_RS %>% count(qMedicalHistoryDepression) %>% mutate(freq = (n / sum(n)) * 100) # 285 Yes (11.6%), 2167 No (88.4%)

### gender % ----
metadata_DS %>% count(qGender) %>% mutate(freq = (n / sum(n)) * 100) #  female   2477  (55.6%), male   1986  (44.4%) total DS
metadata_RS %>% count(qGender) %>% mutate(freq = (n / sum(n)) * 100)  #  female   1498  (61.1%), male   954  (38.9%) total RS

### depression prevalence ----
metadata_DS %>% count(qMedicalHistoryDepression) %>%  mutate(freq = (n / sum(n)) * 100) # n=544, 12% 
metadata_RS %>% count(qMedicalHistoryDepression) %>%  mutate(freq = (n / sum(n)) * 100) # n=285, 12% 

### depression prevalence in female ----
metadata_DS %>% filter(qGender == 0) %>% count(qMedicalHistoryDepression) %>%  mutate(freq = (n / sum(n)) * 100) # n=348, 14% 
metadata_RS %>% filter(qGender == 0) %>% count(qMedicalHistoryDepression) %>%  mutate(freq = (n / sum(n)) * 100) # n=182, 12% 

### depression prevalence in male ----
metadata_DS %>% filter(qGender == 1) %>% count(qMedicalHistoryDepression) %>%  mutate(freq = (n / sum(n)) * 100) # n=196, 10% 
metadata_RS %>% filter(qGender == 1) %>% count(qMedicalHistoryDepression) %>%  mutate(freq = (n / sum(n)) * 100) # n=103, 11% 

metadata_full %>% filter(qGender == 0) %>% count(qMedicalHistoryDepression) %>%  mutate(freq = (n / sum(n)) * 100) # n=530, 13.3% 
metadata_full %>% filter(qGender == 1) %>% count(qMedicalHistoryDepression) %>%  mutate(freq = (n / sum(n)) * 100) # n=299, 10.2% 

### weight groups prevalence ----
metadata_DS %>% filter(qMedicalHistoryDepression == 1) %>% count(weight.status) %>%  mutate(freq = (n / sum(n)) * 100) 
metadata_DS %>% filter(qMedicalHistoryDepression == 2) %>% count(weight.status) %>%  mutate(freq = (n / sum(n)) * 100) 
metadata_RS %>% filter(qMedicalHistoryDepression == 1) %>% count(weight.status) %>%  mutate(freq = (n / sum(n)) * 100) 
metadata_RS %>% filter(qMedicalHistoryDepression == 2) %>% count(weight.status) %>%  mutate(freq = (n / sum(n)) * 100) 

### weight groups prevalence ----
metadata_full %>% count(qMedicalHistoryMigraine) %>%  mutate(freq = (n / sum(n)) * 100) #  n=520  7.52%
metadata_full %>% count(qMedicalHistoryBloating) %>%  mutate(freq = (n / sum(n)) * 100) #  n=1021  14.8%
metadata_full %>% count(qSmokingStatus) %>%  mutate(freq = (n / sum(n)) * 100) #  1 - Quit smoking n=2495, 36.1%
metadata_full %>% count(qSmokingStatus) %>%  mutate(freq = (n / sum(n)) * 100) #  2- Smoking now n=514  7.53%
metadata_full %>% count(qMedicalHistoryAsthma) %>%  mutate(freq = (n / sum(n)) * 100) #  n=381  5.51%
metadata_full %>% count(qMedicalHistorySmallIntestinalBacterialOvergrowth) %>%  mutate(freq = (n / sum(n)) * 100) #  n=133  1.92%
metadata_full %>% count(qPArecreationalVigorousPresent) %>%  mutate(freq = (n / sum(n)) * 100) #  n=2081  30.1%
metadata_full %>% count(qPArecreationalModeratePresent) %>%  mutate(freq = (n / sum(n)) * 100) #  n=2650  38.3%
metadata_full %>% count(qMedicalHistoryPolycysticOvarySyndrome) %>%  mutate(freq = (n / sum(n)) * 100) #  n=246  3.56%
metadata_full %>% count(qDrinkSweetBeverages) %>%  mutate(freq = (n / sum(n)) * 100) #  n=393  5.69%
metadata_full %>% count(qEatFoodRichInVitC) %>%  mutate(freq = (n / sum(n)) * 100) #  n=3956  57.3%
metadata_full %>% count(qMedicalHistoryLactoseIntolerance) %>%  mutate(freq = (n / sum(n)) * 100) #  n=484  7.00%
metadata_full %>% count(qMoreSalt) %>%  mutate(freq = (n / sum(n)) * 100) # 4 - always adding n=462  6.68%
metadata_full %>% count(qAlcoholUnitsPerWeek) %>%  mutate(freq = (n / sum(n)) * 100) # qAlcoholUnitsPerWeek3 Drinking more than 20 units of alcohol per week  n=307  4.43%
metadata_full %>% count(qPAtravelBicycleUse) %>%  mutate(freq = (n / sum(n)) * 100) #  n=3704  53.7%
metadata_full %>% count(qEatOliveOil) %>%  mutate(freq = (n / sum(n)) * 100) #  n=4272  61.8%
metadata_full %>% count(qEatFiber) %>%  mutate(freq = (n / sum(n)) * 100) #  n=5075  73.4%
metadata_full %>% count(qAlcoholStatus) %>%  mutate(freq = (n / sum(n)) * 100) 
# qAlcoholStatus4 Alcohol consumption from 1 to 2 units a week,  n=1656  24.0%
# qAlcoholStatus2 Alcohol consumption from 1 to 3 units a month and less,  n=1504  21.9%
metadata_full %>% count(qMedicalHistoryDiabetesType2) %>%  mutate(freq = (n / sum(n)) * 100) #  n=96  1.4%
metadata_full %>% count(qMedicalHistoryPsoriasis) %>%  mutate(freq = (n / sum(n)) * 100) #  n=269  7.53%
metadata_full %>% count(qEatNuts) %>%  mutate(freq = (n / sum(n)) * 100) #  n=3825  55%
metadata_full %>% count(qDrugHistoryProbioticsType_bifidobacterium) %>%  mutate(freq = (n / sum(n)) * 100) #  n=805  12%
metadata_full %>% count(qDrugHistoryProbioticsType_lactobacillus) %>%  mutate(freq = (n / sum(n)) * 100) #  n=950  14%
metadata_full %>% count(qAlcoholUnitsPerWeek) %>%  mutate(freq = (n / sum(n)) * 100) # qAlcoholUnitsPerWeek2 - from9to20units , n=981  14.2%
metadata_full %>% count(qMedicalHistoryHypothyroidism) %>%  mutate(freq = (n / sum(n)) * 100) #  n=449  6%
metadata_full %>% count(qMedicalHistoryHighBloodPressure) %>%  mutate(freq = (n / sum(n)) * 100) #  n=490  7.1%
metadata_full %>% count(short_sleep) %>%  mutate(freq = (n / sum(n)) * 100) #  n=1259  18%
metadata_full %>% count(long_sleep) %>%  mutate(freq = (n / sum(n)) * 100) #  n=67  1%
metadata_full %>% count(normal_sleep) %>%  mutate(freq = (n / sum(n)) * 100) #  n=5585  81%

### Manifestation date
metadata_full %>% count(qMedicalHistoryDepressionManifestationDate) %>%  mutate(freq = (n / sum(n)) * 100) # 400 individuals reported depression disease manifestation date

##ADDING DEPRESSION + Bloating group ----

metadata_full_Beta <- metadata_full %>%
  mutate(Depression_Bloating = case_when(
    qMedicalHistoryDepression == 1 ~ "1",
    qMedicalHistoryDepression == 2 & qMedicalHistoryBloating == 1 ~ "2",
    qMedicalHistoryDepression == 2 & qMedicalHistoryBloating == 2 ~ "3"
  ))

metadata_full_Beta <- metadata_full_Beta %>% filter(!is.na(Depression_Bloating))
metadata_full_Beta %>% count(Depression_Bloating) %>% mutate(freq = (n / sum(n)) * 100) 


# _______________----

## VIOLIN PLOT ----

library(ggpubr)
library(patchwork)

metadata_violin <- metadata_full %>% mutate(qMedicalHistoryDepression = ifelse(qMedicalHistoryDepression == '1', "controls", "DD cases")) 

# Fit a generalized linear model adjusting for age, gender, and region

glm_model <- glm(alpha_diversity ~ qMedicalHistoryDepression + age +  reg  + qGender + weeks.since.testing.microbiome, data = metadata_violin)
p_value <- summary(glm_model)$coefficients["qMedicalHistoryDepressionDD cases", "Pr(>|t|)"]

p_value <- summary(glm_model)$coefficients["qMedicalHistoryDepressionDD cases", "Pr(>|t|)"]
n_cases <- sum(metadata_violin$qMedicalHistoryDepression == "DD cases")
n_controls <- sum(metadata_violin$qMedicalHistoryDepression == "controls")

# Assuming your data is in a dataframe called 'data' and alpha-diversity is a numeric column

p2 <- ggplot(metadata_violin, aes(x = qMedicalHistoryDepression, y = alpha_diversity, fill = qMedicalHistoryDepression)) +
  geom_violin(trim = FALSE, alpha = 0.7) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Boxplot inside violin plot
  annotate("text", x = 1.2, y = min(metadata_violin$alpha_diversity), label = paste("n =", n_controls), size = 5) +  # Add number of controls
  annotate("text", x = 2.2, y = min(metadata_violin$alpha_diversity), label = paste("n =", n_cases), size = 5) +  # Add number of cases
  annotate("text", x = 1.5, y = max(metadata_violin$alpha_diversity), label = paste("Adjusted p =", signif(p_value, 3)), size = 5) +  # Add p-value
  scale_fill_manual(values = c("lightblue", "lightpink")) +  # Customize colors
  labs(x = "Group", y = "Alpha-Diversity", title = "Violin plot of alpha-diversity for DD cases and controls adjusted for age, gender, region, batch effect" ) + 
  theme_minimal()+
  theme(legend.position = "none")
p2

p1 <- ggplot(metadata_violin, aes(x = qMedicalHistoryDepression, y = alpha_diversity, fill = qMedicalHistoryDepression)) +
  geom_violin(trim = FALSE, alpha = 0.7) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Boxplot inside violin plot
  stat_compare_means(method = "wilcox.test", label = "p.format", hjust = -0.8, vjust = 3.1) +  # Add p-value using Wilcoxon test - no adjustment
  annotate("text", x = 1.2, y = min(metadata_violin$alpha_diversity), label = paste("n =", n_controls), size = 5) +  # Add number of controls
  annotate("text", x = 2.2, y = min(metadata_violin$alpha_diversity), label = paste("n =", n_cases), size = 5) +  # Add number of cases
  scale_fill_manual(values = c("lightblue", "lightpink")) +  # Customize fill colors
  labs(x = "Group", y = "Alpha-Diversity", title = "Violin plot of alpha-diversity for DD cases and controls" ) + 
  theme_minimal()+
  theme(legend.position = "none")
p1

combined_plot <- p1 + p2
combined_plot

# _______________----
# REGRESSION ----

mergedTotal <- read.csv("/Users/metadata.csv") #7,506

## MRS and MDD regression output ----
#MDD for 9 taxa was calculated in the MRS_MR.R script and added to the metadata file for regression analysis

regress_output <- list()

for ( i in 1:length(variables)){
  
  print(paste("Fitting ", i, variables[i])) 
  
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  #mod <- as.formula(sprintf("Microbiome_Risk_Score_MDD ~ %s", variables[i])) 
  mod <- as.formula(sprintf("Microbiome_Risk_Score_MDD ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort+BMI", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = mergedTotal ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL
rownames(allregression)

#### multiple testing -----

#remove Intersepr rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_MRS_MDD <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "BMI", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)
#allregression_MRS_DD_non_adj <- allregression %>% filter(!term %in% c("(Intercept)")) %>% arrange(p.value)
#p-value correction 
allregression_MRS_DD$Bonferroni =
  p.adjust(allregression_MRS_MDD$p.value,
           method = "bonferroni")
allregression_MRS_DD$BH_MRS_MDD =
  p.adjust(allregression_MRS_MDD$p.value,
           method = "BH")
#Add description
Description <- read.csv("/User/Metadata_names.csv") %>% select(term = Name, DisplayName , cases, percent)
allregression_MRS_MDD <- left_join(allregression_MRS_MDD, Description, by = "term") %>% select(DisplayName, cases, percent, everything())
write.csv(allregression_MRS_MDD, "/Users/allregression_multipletesting_MRS_MDD.csv")

## Alpha-diversity and MDD regression output ----

regress_output <- list()

for ( i in 1:length(variables)){
  
  print(paste("Fitting ", i, variables[i])) 
  
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  mod <- as.formula(sprintf("alpha_diversity ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = mergedTotal ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL

#### multiple testing -----

#Let's prepare spreadsheet for multiple testing Shannon diversity index
#remove Intersept rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_alpha <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)

#p-value correction 

allregression_alpha$Bonferroni =
  p.adjust(allregression_alpha$p.value,
           method = "bonferroni")

allregression_alpha$BH_Alpha =
  p.adjust(allregression_alpha$p.value,
           method = "BH")

#Add description
Description <- read.csv("/Users/Metadata_names.csv") %>% select(term = Name, DisplayName)

allregression_alpha <- left_join(allregression_alpha, Description, by = "term") %>% select(DisplayName, everything())

write.csv(allregression_alpha,"/Users/allregression_multipletesting_Alpha_diversity.csv")

# _______________----
# PERMANOVA ----

## Bray-Curtis ----

distance_matrix <- vegdist(bacteriadata, method = "bray") # bacteriadata is the abundance data with samples as rows and taxa as columns

# testing for differences in microbial composition based on depression status, adjusted for age and gender
permanova_result <- adonis2(distance_matrix ~ qMedicalHistoryDepression + age + reg + qGender + weeks.since.testing.microbiome, data = metadata_full, permutations = 999)
print(permanova_result)
ordination <- cmdscale(distance_matrix, k = 2, eig = TRUE) # k = 2 for 2D plotting

# Convert PCoA result into a data frame
ordination_df <- as.data.frame(ordination$points)
ordination_df$barcode <- rownames(ordination$points)

# Add metadata to the ordination data frame
ordination_df <- ordination_df %>%
  left_join(metadata_full, by = "barcode")

# Plot the PCoA with ggplot
pcoa_plot <- ggplot(ordination_df, aes(x = V1, y = V2, color = qMedicalHistoryDepression)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(title = "PCoA of Bray-Curtis Dissimilarity by Depression and Bloating Status",
       x = "PCoA1", y = "PCoA2") +
  scale_color_manual(values = c("1" = "#00BFC4", 
                                "2" = "#D55E00", 
                                "3" = "pink"))


## Aitchison Distance ----
library(compositions)
clr_data <- clr(bacteriadata + 1e-6) # Add a small pseudocount to avoid log(0)
aitchison_distance <- dist(clr_data, method = "euclidean") # Euclidean distance on the CLR-transformed data to obtain the Aitchison distance
permanova_aitchison <- adonis2(aitchison_distance ~ qMedicalHistoryDepression, data = metadata_full, permutations = 999) # Perform PERMANOVA using the Aitchison distance matrix 
write.csv(permanova_aitchison, "/Users/permanova_aitchison.csv")
print(permanova_aitchison)

# _______________----

# ___DISCOVERY SET___ ----

# MAASLIN 2 ----

bacterialist.maaslin.DS <- bacterialist.maaslin %>% rownames_to_column("barcode") #6,911, 976
bacterialist.maaslin.DS <- bacterialist.maaslin.DS %>% filter(barcode%in%metadata_DS$barcode) %>% column_to_rownames("barcode") #4,456
metadata_maaslin_DS <-  metadata_DS %>% column_to_rownames("barcode")

#for Depression + age + reg + gender + weeks.since.testing as a covariate

maaslin.data = Maaslin2(
  input_data = bacterialist.maaslin.DS, 
  input_metadata = metadata_maaslin_DS, 
  output = "/Users/Maaslin.Depression_DS_for_ancom", 
  #min_prevalence = 0.0,
  min_prevalence = 0.1,
  fixed_effects = c("qMedicalHistoryDepression", "age","reg", "qGender","weeks.since.testing.microbiome"),
  reference = c("qMedicalHistoryDepression,1","reg,en", "qGender,1","age","weeks.since.testing.microbiome"))

all_results <-read_tsv("/Users/Maaslin.Depression_DS_for_ancom/all_results.tsv")
all_results2 <- all_results %>% 
  filter(metadata == "qMedicalHistoryDepression") %>% 
  select(bacterianame = feature, N_DS = N , N.not.0_DS = N.not.0, estimate.maaslin_DS = coef, stderr.maaslin_DS = stderr, pval.maaslin_DS = pval, qval.maaslin_DS = qval ) %>% arrange(pval.maaslin_DS)
write_csv(all_results2, "/Users/Maaslin_Depression_DS.csv") #111

# PHYSEQ file for ANCOM-BC2 ----

# OTUmatrix

bacterialist_DS_0to1 = bacterialist.ancom_DS

bacterialist_DS_0to1 <- bacterialist_DS_0to1 %>% rownames_to_column("barcode")
bacterialist_DS_0to1 <- bacterialist_DS_0to1 %>% filter(barcode%in% metadata_DS$barcode) %>% column_to_rownames("barcode") #4,456 , 380
bacterialist_DS_0to1[bacterialist_DS_0to1 == 0] <- 1 #This is a step recommended by Author of ANCOM-BC2 method from the email

OTUmatrix__DS_0to1 <- bacterialist_DS_0to1 %>% t() #
OTU_0to1 = otu_table(OTUmatrix__DS_0to1, taxa_are_rows = TRUE)
view(OTUmatrix__DS_0to1)

# Metadata table 

SAMPLEdata <- metadata_DS %>% select(barcode, 
                                     reg, 
                                     weeks.since.testing.microbiome, 
                                     butirate.percentage, 
                                     age, qGender, BMI, 
                                     enterotype, 
                                     cohort, 
                                     weight.status, 
                                     qPAyoga, 
                                     alpha_diversity, 
                                     qMedicalHistoryDepression,
                                     qMoreSalt, 
                                     qEatFiber, 
                                     qDrinkSweetBeverages, 
                                     qEatBackedSweets, 
                                     qEatNonOilyFishRegularly, 
                                     qEatFoodRichInOmega3, 
                                     qEatFoodRichInVitC, 
                                     qEatFoodRichInVitD, 
                                     qAlcoholUnitsPerWeek, 
                                     qAlcoholStatus, 
                                     qDrugHistoryProbioticsType_lactobacillus,
                                     qDrugHistoryProbioticsType_bifidobacterium,
                                     qSmokingStatus, 
                                     qPArecreationalModeratePresent, 
                                     qPAworkVigorousPresent, 
                                     qPAtravelBicycleUse) %>% 
  tibble::column_to_rownames('barcode')
sampledata = sample_data(SAMPLEdata)
all.equal(colnames(OTUmatrix__DS_0to1), rownames(SAMPLEdata)) # TRUE, samples check
dim(bacterialist) #file size = 6911, 975  - all unique USER_ID's and removed barcode columns, renaimed user_id column (compare with row 9)

#TAXA DATA 

TAXmatrix_0to1 <- as.data.frame(colnames(bacterialist_DS_0to1)) %>% 
  rename(bacterianame = 'colnames(bacterialist_DS_0to1)') %>% 
  separate(bacterianame, sep="\\.", into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), remove = FALSE) %>% 
  mutate(Genus = paste(Genus, Species)) %>%  #we bind genus and species columns so we don't lose information about some more specific bacteria and our results and bacteria names match Maaslin2 analysis.
  select(-Species) %>% 
  tibble::column_to_rownames ('bacterianame') %>% as.matrix()
view(TAXmatrix_0to1)
#TAXmatrix to TAXdata format
TAX_0to1 = tax_table(TAXmatrix_0to1)
physeq_0to1 = phyloseq(OTU_0to1, TAX_0to1, sampledata)
save(physeq_0to1, file="/Users/phyloseq_Depression_for_ancom2_DS.RData")

#ANCOM-BC2 ----

#to read this file
load("/Users/phyloseq_Depression_for_ancom2_DS.RData") 

tse = mia::makeTreeSummarizedExperimentFromPhyloseq(physeq_0to1)

output = ancombc2(data = tse, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "age + reg + weeks.since.testing.microbiome + qGender + qMedicalHistoryDepression", rand_formula = NULL,
                  p_adj_method = "BH", pseudo_sens = TRUE,
                  prv_cut = 0.1, lib_cut = 1000, s0_perc = 0.05,
                  group = "qMedicalHistoryDepression", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(-1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2, 2, 1),
                                       solver = "ECOS",
                                       B = 10))
res_prim = output$res
tmp = output$res
res_prim <- res_prim %>% select(taxon, estimate.ancom = "lfc_qMedicalHistoryDepression2", 
                                stderr.ancom = "se_qMedicalHistoryDepression2", 
                                W = "W_qMedicalHistoryDepression2", 
                                pval.ancom = "p_qMedicalHistoryDepression2",
                                qval.ancom = "q_qMedicalHistoryDepression2",
                                diff = "diff_qMedicalHistoryDepression2") 

length(unique(res_prim$taxon)) # taxon items = 375
write.csv(res_prim,"/Users/ANCOMBC2.Depression_DS.csv") #all 375 taxa, 18 taxa significant q<0.05
# _______________----

# ___REPLICATION SET___ ----

# MAASLIN 2 ----

# bacterialist.maaslin.RS_full - 375 full list, sign = 32, q<0.05 

bacterialist.maaslin.RS_full <- bacterialist.maaslin %>% rownames_to_column("barcode")
bacterialist.maaslin.RS_full <- bacterialist.maaslin.RS_full %>% filter(barcode%in%metadata_RS$barcode) %>% column_to_rownames("barcode") #2,452, 975 taxa
metadata_maaslin_RS <- metadata_RS %>% column_to_rownames("barcode") #2,452, 237

#for Depression + age + reg + gender + weeks.since.testing.microbiome as a covariate

maaslin.data = Maaslin2(
  input_data = bacterialist.maaslin.RS_full, 
  input_metadata = metadata_maaslin_RS, 
  output = "/Users/Maaslin.Depression_RS_meta", 
  #min_prevalence = 0.0, #create a list of taxa for ANCOM-BC2 analysis
  min_prevalence = 0.1, #for main Maaslin analysis
  fixed_effects = c("qMedicalHistoryDepression", "age","reg", "qGender","weeks.since.testing.microbiome"),
  reference = c("qMedicalHistoryDepression,1","reg,en", "qGender,1","age","weeks.since.testing.microbiome"))

all_results <-read_tsv("/Users/Maaslin.Depression_RS_meta/all_results.tsv")
all_results2 <- all_results %>% 
  filter(metadata == "qMedicalHistoryDepression") %>% 
  select(bacterianame = feature, N_RS = N, N.not.0_RS = N.not.0, estimate.maaslin_RS = coef, stderr.maaslin_RS = stderr, pval.maaslin_RS = pval, qval.maaslin_RS = qval ) %>% arrange(pval.maaslin_RS)
write_csv(all_results2, "/Users/Depression/Maaslin2Depression_RS.csv") #116

#ANCOM-BC2 ----

# OTUmatrixvz

#ancom_RS - 321 full list

bacterialist_RS_0to1=ancom_RS

bacterialist_RS_0to1 <- bacterialist_RS_0to1 %>% rownames_to_column("barcode")
bacterialist_RS_0to1 <- bacterialist_RS_0to1 %>% filter(barcode%in% metadata_RS$barcode) %>% column_to_rownames("barcode") #2,452
bacterialist_RS_0to1[bacterialist_RS_0to1 == 0] <- 1 #This is a step recommended by Author of ANCOM-BC2 method from the email

OTUmatrix__RS_0to1 <- bacterialist_RS_0to1 %>% t() #
OTU_0to1 = otu_table(OTUmatrix__RS_0to1, taxa_are_rows = TRUE)
view(OTU_0to1) #2452

# Metadata table 

SAMPLEdata <- metadata_RS %>% select(barcode, reg, weeks.since.testing.microbiome, butirate.percentage, age, qGender, BMI, enterotype, cohort, weight.status, qPAyoga, alpha_diversity, qMedicalHistoryDepression , qMoreSalt, qEatFiber, qDrinkSweetBeverages, qEatBackedSweets, qEatNonOilyFishRegularly, qEatFoodRichInOmega3, qEatFoodRichInVitC, qEatFoodRichInVitD, qAlcoholUnitsPerWeek, qAlcoholStatus, BMI, qDrugHistoryProbioticsType_lactobacillus, qSmokingStatus, qPArecreationalModeratePresent, qPAworkVigorousPresent, qPAtravelBicycleUse) %>% 
  tibble::column_to_rownames('barcode')
sampledata = sample_data(SAMPLEdata)
all.equal(colnames(OTUmatrix__RS_0to1), rownames(SAMPLEdata)) # TRUE, samples check
dim(bacterialist) #file size = 6911  975 - all unique USER_ID's and removed barcode columns, renaimed user_id column (compare with row 9)

#TAXA DATA 

TAXmatrix_0to1 <- as.data.frame(colnames(bacterialist_RS_0to1)) %>% 
  rename(bacterianame = 'colnames(bacterialist_RS_0to1)') %>% 
  #filter(!bacterianame %in% c("user_id", "barcode_dots")) %>% #remove 2 columns
  separate(bacterianame, sep="\\.", into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), remove = FALSE) %>% 
  mutate(Genus = paste(Genus, Species)) %>%  #we bind genus and spisies colunms so we don't lose information about some more specific bacteria and our results and bacteria names match Maaslin analysis.
  select(-Species) %>% 
  tibble::column_to_rownames ('bacterianame') %>% as.matrix()

#TAXmatrix to TAXdata format

TAX_0to1 = tax_table(TAXmatrix_0to1)
physeq_0to1 = phyloseq(OTU_0to1, TAX_0to1, sampledata)
save(physeq_0to1, file="/Users/phyloseq_Depression_for_ancom2_RS.RData")

#to read this file
load("/Users/phyloseq_Depression_for_ancom2_RS.RData")

tse = mia::makeTreeSummarizedExperimentFromPhyloseq(physeq_0to1)
output = ancombc2(data = tse, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "age + reg + weeks.since.testing.microbiome + qGender + qMedicalHistoryDepression", rand_formula = NULL,
                  p_adj_method = "BH", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "qMedicalHistoryDepression", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(-1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2, 2, 1),
                                       solver = "BH",
                                       B = 10))
res_prim = output$res
tmp = output$res
res_prim <- res_prim %>% select(taxon, estimate.ancom_RS = "lfc_qMedicalHistoryDepression2", 
                                stderr.ancom_RS = "se_qMedicalHistoryDepression2", 
                                W_RS = "W_qMedicalHistoryDepression2", 
                                pval.ancom_RS = "p_qMedicalHistoryDepression2",
                                qval.ancom_RS = "q_qMedicalHistoryDepression2",
                                diff_RS = "diff_qMedicalHistoryDepression2") 
length(unique(res_prim$taxon)) # taxon items = 317

write.csv(res_prim,"/Users/ANCOMBC2.Depression_RS.csv") #317 taxa (full list analysis)

# __META-ANALYSIS DATA ----

# MASSLIN ----

maaslin_DS_meta <- read.csv("/Users/Maaslin_Depression_DS.csv")#111
maaslin_RS_meta <- read.csv("/Users/Maaslin.Depression_RS.csv")#116

maaslin_meta <- maaslin_DS_meta %>% left_join(maaslin_RS_meta, by = c("bacterianame"))

# ANCOM -----

ancombc_DS_meta <- read.csv("/Users/ANCOMBC2.Depression_DS.csv") %>% rename(Genus = taxon) %>% select(-1) #375
ancombc_RS_meta <- read.csv("/Users/ANCOMBC2.Depression_RS.csv") %>% select(-1)%>% rename(Genus = taxon)  #317
ancombc_DS_meta$Genus <- trimws(ancombc_DS_meta$Genus)
ancombc_RS_meta$Genus <- trimws(ancombc_RS_meta$Genus)

ancom_list <- left_join(maaslin_DS_meta, bacteriadata_names, by = c("bacterianame")) %>% select(bacterianame, FamilyGenus, Genus, N_DS, N.not.0_DS) #111
ancom_meta <- left_join(ancom_list, ancombc_DS_meta, by = c("Genus")) # 111 - we need only taxa with taxa level > 445 (10%) for Discovery Set
ancom_meta <- left_join(ancom_meta, ancombc_RS_meta, by = c("Genus")) # 111 - we need only taxa with taxa level > 240 (10%) for Replication Set

# TOTAL ----

total_meta <- left_join(maaslin_meta, Displayname.depression, by = c("bacterianame"))
total_meta <- left_join(total_meta, ancom_meta , by = c("bacterianame")) %>% select(-FamilyGenus) %>% select(Displayname, bacterianame, Genus, everything())

write.csv(total_meta, "/Users/METAANALYSIS_MODEL.csv")
# _______________----

# ____META-ANALYSIS____ ----

# MAASLIN METAGEN ----

maaslin_meta <- read_csv("/Users/METAANALYSIS_MODEL.csv") %>% select(-1)

maaslin_list <- maaslin_meta %>% select(
  bacterianame = c(Displayname),
  estimate_DS = c(estimate.maaslin_DS),  # Effect sizes DS
  estimate_RS = c(estimate.maaslin_RS),  # Effect sizes RS
  stderr_DS = c(stderr.maaslin_DS),
  stderr_RS = c(stderr.maaslin_RS),        # Standard errors
  pval_DS = c(pval.maaslin_DS),
  pval_RS = c(pval.maaslin_RS),
)
# Combine effect sizes and standard errors using inverse variance weighting
maaslin_list$combined_estimate <- (maaslin_list$estimate_DS / maaslin_list$stderr_DS^2 + 
                                     maaslin_list$estimate_RS / maaslin_list$stderr_RS^2) / 
  (1 / maaslin_list$stderr_DS^2 + 1 / maaslin_list$stderr_RS^2)
maaslin_list$combined_stderr <- sqrt(1 / (1 / maaslin_list$stderr_DS^2 + 1 / maaslin_list$stderr_RS^2))

# View the combined effect sizes and standard errors
print(maaslin_list[, c("bacterianame", "combined_estimate", "combined_stderr")])

# Combine p-values using Fisher's method 
maaslin_list$combined_pval <- 1 - pchisq(-2 * (log(maaslin_list$pval_DS) + log(maaslin_list$pval_RS)), df = 4)
# View the combined p-values
print(maaslin_list[, c("bacterianame", "combined_pval")])

# Perform a meta-analysis on the combined results
meta_combined_maaslin <- metagen(TE = maaslin_list$combined_estimate, 
                                 seTE = maaslin_list$combined_stderr, 
                                 studlab = maaslin_list$bacterianame, 
                                 comb.random = TRUE)  # Random-effects model

# View the meta-analysis summary
summary(meta_combined_maaslin)

# Choose significant data for the plot 
MaaslinMataA <- meta_combined_maaslin$data
MaaslinMataA$p_val <- meta_combined_maaslin$pval
view(meta_combined_maaslin$pval.random) # p = 1.452067e-05
view(meta_combined_maaslin$pval.fixed) # p = 0.0001423726
MaaslinMataA <- MaaslinMataA %>% select(Displayname = .studlab, Estimate.TE_MaAsLin2 = .TE, StdErr_MaAsLin2 = .seTE, p_val_MaAsLin2 = p_val )

significant_indices <- which(meta_combined_maaslin$pval < 0.05)
significant_results <- maaslin_list[significant_indices, ]
meta_result_significant_maaslin <- metagen(TE = significant_results$combined_estimate, 
                                           seTE = significant_results$combined_stderr, 
                                           studlab = significant_results$bacterianame, 
                                           comb.random = TRUE,
                                           data = significant_results)
summary(significant_results)
funnel(meta_result_significant_maaslin)

# Customize the forest plot
colors <- ifelse(meta_result_significant_maaslin$TE > 0, "red", "green")

forest(meta_result_significant_maaslin, 
       xlab = "Combined Effect Size (TE)",   
       main = "Forest plot for MaAsLin2 method with colour-coded genera and DD associatiation",  
       refline = 0,                 # Reference line at 0
       col.study = "black",         # Color for genus names
       col.square = colors,
       xlim = c(-0.5, 0.5),         # Adjust x-axis limits
       digits = 2,                  # Number of decimal places for effect sizes
       cex = 0.9)                   # Text size adjustment

# Create and save the forest plot
png("forest_plot_Maaslin2.png", width = 5100, height = 7700, res = 500)

forest(meta_result_significant_maaslin, 
       xlab = "Combined Effect Size (TE)",   
       main = "Forest plot for MaAsLin2 method with colour-coded genera and DD associatiation",  
       refline = 0,                 # Reference line at 0
       col.study = "black",         # Color for genus names
       col.square = colors,
       xlim = c(-0.5, 0.5),         # Adjust x-axis limits
       digits = 2,                  # Number of decimal places for effect sizes
       cex = 0.9)                   # Text size adjustment
# Use the color vector for the squares
dev.off()

# ANCOM METAGEN ----

ancom_meta <- read_csv("/Users/METAANALYSIS_MODEL.csv") %>% select(-1)
colnames(Full_list)

ancom_list <- ancom_meta %>% select(
  bacterianame = c(Displayname),
  estimate_DS = c(estimate.ancom),  # Effect sizes DS
  estimate_RS = c(estimate.ancom_RS),  # Effect sizes RS
  stderr_DS = c(stderr.ancom),
  stderr_RS = c(stderr.ancom_RS),        # Standard errors
  pval_DS = c(pval.ancom),
  pval_RS = c(pval.ancom_RS),
)
# Combine effect sizes and standard errors using inverse variance weighting
ancom_list$combined_estimate <- (ancom_list$estimate_DS / ancom_list$stderr_DS^2 + 
                                   ancom_list$estimate_RS / ancom_list$stderr_RS^2) / 
  (1 / ancom_list$stderr_DS^2 + 1 / ancom_list$stderr_RS^2)
ancom_list$combined_stderr <- sqrt(1 / (1 / ancom_list$stderr_DS^2 + 1 / ancom_list$stderr_RS^2))

# View the combined effect sizes and standard errors
print(ancom_list[, c("bacterianame", "combined_estimate", "combined_stderr")])

# Combine p-values using Fisher's method 
ancom_list$combined_pval <- 1 - pchisq(-2 * (log(ancom_list$pval_DS) + log(ancom_list$pval_RS)), df = 4)
# View the combined p-values
print(ancom_list[, c("bacterianame", "combined_pval")])

# Perform a meta-analysis on the combined results
meta_combined <- metagen(TE = ancom_list$combined_estimate, 
                         seTE = ancom_list$combined_stderr, 
                         studlab = ancom_list$bacterianame, 
                         comb.random = TRUE)  # Random-effects model

# View the meta-analysis summary
summary(meta_combined)

# Choose significant data for the plot 
AncomMataA <- meta_combined$data
AncomMataA$p_val <- meta_combined$pval
view(meta_combined$pval.random) # p = 0.00434101
view(meta_combined$pval.fixed) # p = 8.001241e-08
AncomMataA <- AncomMataA %>% select(Displayname = .studlab, Estimate.TE_ANCOMBC2 = .TE, StdErr_ANCOMBC2 = .seTE, p_val_ANCOMBC2 = p_val)

significant_indices <- which(meta_combined$pval < 0.05)
significant_results <- ancom_list[significant_indices, ]
meta_result_significant <- metagen(TE = significant_results$combined_estimate, 
                                   seTE = significant_results$combined_stderr, 
                                   studlab = significant_results$bacterianame, 
                                   comb.random = TRUE,
                                   data = significant_results)

summary(significant_results)
funnel(meta_result_significant)

# Customize the forest plot
colors <- ifelse(meta_result_significant$TE > 0, "red", "green")

forest(meta_result_significant, 
       xlab = "Combined Effect Size (TE)",   
       main = "Forest plot for ANCOM-BC2 method with colour-coded genera and DD associatiation",  
       refline = 0,                 # Reference line at 0
       col.study = "black",         # Color for genus names
       col.square = colors,
       xlim = c(-0.5, 0.5),         # Adjust x-axis limits
       digits = 2,                  # Number of decimal places for effect sizes
       cex = 0.9)                   # Text size adjustment
# Use the color vector for the squares
png("forest_plot_ANCOM-BC2.png", width = 5100, height = 7700, res = 500)
dev.off()

# _______________----

# Meta- analysis FDR ----

MetaAnalysis <- full_join(AncomMataA, MaaslinMataA, by = "Displayname") 

#Full list n=111

res_all <- MetaAnalysis 

#FDR for all metaanalysed taxa (n=111) 

res_all$BH_ANCOMBC2 =
  p.adjust(res_all$p_val_ANCOMBC2,
           method = "BH")
res_all$BH_MaAsLin2 =
  p.adjust(res_all$p_val_MaAsLin2,
           method = "BH")

res_all <- res_all %>% select(Displayname, Estimate.TE_ANCOMBC2, StdErr_ANCOMBC2, p_val_ANCOMBC2, BH_ANCOMBC2,everything())

write_csv(res_all, "/Users/MetaAnalysis_Maaslin2_ANCOMBC2_FDR.csv")

#MetaAnalysis full list ----

Meta_full <- res_all %>% left_join(Full_list, by = c("Displayname")) %>% 
  mutate(direction_category = case_when(
    MetaAnalysis$Estimate.TE_MaAsLin2 > 0 & MetaAnalysis$Estimate.TE_ANCOMBC2 > 0 ~ "Both Positive",
    MetaAnalysis$Estimate.TE_MaAsLin2 < 0 & MetaAnalysis$Estimate.TE_ANCOMBC2 < 0 ~ "Both Negative",
    MetaAnalysis$Estimate.TE_MaAsLin2 > 0 & MetaAnalysis$Estimate.TE_ANCOMBC2 < 0 ~ "Positive MaAsLin2, Negative ANCOMBC2",
    MetaAnalysis$Estimate.TE_MaAsLin2 < 0 & MetaAnalysis$Estimate.TE_ANCOMBC2 > 0 ~ "Negative MaAsLin2, Positive ANCOMBC2")) %>% 
  mutate(estimate_direction = case_when(direction_category == "Both Positive" ~ 1,
                                        direction_category == "Both Negative" ~ 2,
                                        direction_category == "Positive MaAsLin2, Negative ANCOMBC2" ~ 0,
                                        direction_category == "Negative MaAsLin2, Positive ANCOMBC2" ~ 0))


write_csv(Meta_full, "/Users/MetaAnalysis_DS_RS_Maaslin2_ANCOMBC2_Full_list.csv")
Meta_sign_names <- Meta_full%>%  filter (Meta_full$p_val_ANCOMBC2 < 0.05, Meta_full$p_val_MaAsLin2 <0.05 ) %>% 
  select(Displayname, bacterianame, estimate_direction)
write_csv(Meta_sign_names, "/Users/associated_bacteria.csv")
colnames(Meta_full)

# _______________----

# EULER PLOT ----

library(eulerr)

meta_analys <- read.csv("/Users/MetaAnalysis_Maaslin2_ANCOMBC2_FDR.csv")
replication_model <- read_csv("/Users/REPLICATION_MODEL.csv")

colnames(replication_model)

#FDR<0.05, from full list 
s3 <- list("ANCOM-BC2, MetaA" = meta_analys %>% filter(BH_ANCOMBC2<0.05) %>% pull("Displayname")%>% unique(),
           "MaAsLin2, MetaA" = meta_analys %>% filter(BH_MaAsLin2<0.05) %>% pull("Displayname") %>% unique(),
           "ANCOM-BC2, Repl" = replication_model %>% filter(pval.ancom_RS<0.05) %>% pull("Displayname")%>% unique(),
           "MaAsLin2, Repl" = replication_model %>% filter(pval.maaslin_RS<0.05) %>% pull("Displayname")%>% unique())

color_palette <- c("ANCOM-BC2, MetaA" = "#66c2a5", 
                   "MaAsLin2, MetaA" = "#b2e2cd", 
                   "ANCOM-BC2, Repl" = "#8da0cb", 
                   "MaAsLin2, Repl" = "#e78ac3")

plot(euler(s3, shape = "ellipse"), 
     quantities = TRUE, 
     fills = color_palette,  
     labels = list(font = 9,fontsize = 10),  
     main = expression("Meta-analysis and Replication for ANCOM-BC2 and MaAsLin2"))
# _______________----

## -log10 PLOT ----

#let's visualise only taxa < 0.05 in both ana
colnames(MAASLIN_ANCOM)
MAASLIN_ANCOM <- res_all %>% mutate(estimate_MAASLIN_direction = ifelse(Estimate.TE_MaAsLin2 >0, "positive","negative")) %>% 
  mutate(estimate_ANCOM_direction = ifelse(Estimate.TE_ANCOMBC2 >0, "positive","negative")) %>% 
  mutate(estimate_direction = ifelse(estimate_MAASLIN_direction == estimate_ANCOM_direction, "same", "different")) %>% 
  #filter(p_val_ANCOMBC2 <0.05) %>% 
  mutate(q_val_ANCOM_log10 = -log10(BH_ANCOMBC2)) %>%
  mutate(q_val_MAASLIN_log10 = -log10(BH_MaAsLin2)) 

# create column to label features reaching FDR 1E-4 in either method
MAASLIN_ANCOM <- MAASLIN_ANCOM %>% 
  filter(!is.na(q_val_ANCOM_log10)&!is.na(q_val_MAASLIN_log10)) %>% 
  mutate(labels = ifelse(BH_ANCOMBC2< 0.05 & BH_MaAsLin2 < 0.05, Displayname, "")) %>% 
  mutate(estimate_same_direction = ifelse(BH_ANCOMBC2 < 0.05 & BH_MaAsLin2 < 0.05, estimate_ANCOM_direction, 'not associated')) %>% 
  mutate(estimate_same_direction = factor(estimate_same_direction, levels = c('positive', 'negative', 'not associated')))

g3 <- ggplot(data=MAASLIN_ANCOM, aes(y=q_val_ANCOM_log10 , x=q_val_MAASLIN_log10,
                                     label=labels)) +
  geom_point(aes(fill = estimate_same_direction),size=1.5, shape=21) +
  geom_text_repel(min.segment.length=0, box.padding=0.05, size=3, color='grey25') +
  geom_vline(xintercept=-log10(0.05), color='grey50', linetype='dashed') +
  geom_hline(yintercept=-log10(0.05), color='grey50', linetype='dashed') +
  geom_vline(xintercept=-log10(0.1), color='grey90', linetype='dashed') +
  geom_hline(yintercept=-log10(0.1), color='grey90', linetype='dashed') +
  labs(x='-log10(MaAsLin2 meta-analysis FDR)', y='-log10(ANCOM-BC2 meta-analysis FDR)', fill = 'Effect direction') +
  scale_x_continuous(
    breaks=c(0, -log10(0.1), -log10(0.05), -log10(0.01), -log10(1E-4), -log10(1E-6), -log10(1E-8)),
    labels=c('0', paste(round(-log10(0.1),1), ' \n', '(0.1)', sep=''), 
             paste(round(-log10(0.05),1), ' \n', '(0.05)', sep=''), 
             paste(-log10(0.01), '\n', '(0.01)', sep=''), 
             paste(-log10(1E-4), '\n', '(1E-4)', sep=''), 
             paste(-log10(1E-6), '\n', '(1E-6)', sep=''), 
             paste(-log10(1E-8), '\n', '(1E-8)', sep='')),
    limits=c(0, -log10(min(MAASLIN_ANCOM[,c('BH_ANCOMBC2','BH_MaAsLin2')]))),
    minor_breaks=NULL  # This line closes the scale_x_continuous function
  ) +
  scale_fill_manual(values=c('red','blue','grey90')) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  theme_bw() +
  theme(legend.title=element_text(size=10), legend.text=element_text(size=8), legend.text.align=0,
        legend.position=c(0.85,0.1), legend.background=element_rect(fill='white', color='grey50'),
        axis.text.x=element_text(size=10), axis.text.y=element_text(size=10, vjust=0.8),
        axis.title=element_text(size=10),
        axis.title.x=element_text(vjust=-0.75), axis.title.y=element_text(vjust=3),
        plot.margin=margin(t=10, r=30, b=10, l=10, unit = "pt"))
ggsave(g3, "/Users/OUTPUT/MAASLIN_ANCOM_Q_depression.pdf", device='pdf', width=25, height=15)
g3

g3 <- ggplot(data=MAASLIN_ANCOM, aes(y=q_val_ANCOM_log10 , x=q_val_MAASLIN_log10, label=labels)) +
  geom_point(aes(fill = estimate_same_direction), size=3, shape=21) +
  geom_text_repel(
    min.segment.length=0, 
    box.padding=0.2, 
    size=2.5,  
    color='grey25', 
    max.overlaps = Inf,  
    force = 5,  
    nudge_y = 0.5,  
    nudge_x = 0.5  
  ) +
  geom_vline(xintercept=-log10(0.05), color='grey50', linetype='dashed') +
  geom_hline(yintercept=-log10(0.05), color='grey50', linetype='dashed') +
  geom_vline(xintercept=-log10(0.1), color='grey90', linetype='dashed') +
  geom_hline(yintercept=-log10(0.1), color='grey90', linetype='dashed') +
  labs(x='-log10(MaAsLin2 meta-analysis FDR)', y='-log10(ANCOM-BC2 meta-analysis FDR)', fill = 'Effect direction') +
  scale_x_continuous(
    breaks=c(0, -log10(0.1), -log10(0.05), -log10(0.01), 2.1),  # Set limit to 2.1
    labels=c('0', paste(round(-log10(0.1),1), ' \n', '(0.1)', sep=''), 
             paste(round(-log10(0.05),1), ' \n', '(0.05)', sep=''), 
             paste(-log10(0.01), '\n', '(0.01)', sep=''), 
             "2.1"),
    limits=c(0, 2.1),  # Limit the x-axis to 2.1
    minor_breaks=NULL  
  ) +
  scale_fill_manual(values=c('red','blue','grey90')) +
  coord_cartesian(clip="off") +  # Prevent clipping of points that are beyond the limit
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  theme_bw() +
  theme(legend.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.text.align=0,
        legend.position=c(0.85,0.1), 
        legend.background=element_rect(fill='white', color='grey50'),
        axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10, vjust=0.8),
        axis.title=element_text(size=10),
        axis.title.x=element_text(vjust=-0.75), 
        axis.title.y=element_text(vjust=3),
        plot.margin=margin(t=10, r=30, b=10, l=10, unit = "pt"))

# Save the plot
ggsave(g3, "/Users/MAASLIN_ANCOM_Q_depression.pdf", device='pdf', width=25, height=15)
g3


# FOREST PLOTS ---- 

combined_meta <- metagen(TE = c(meta_result_significant$TE, meta_result_significant_maaslin$TE),
                         seTE = c(meta_result_significant$seTE, meta_result_significant_maaslin$seTE),
                         studlab = c(meta_result_significant$studlab, meta_result_significant_maaslin$studlab),
                         comb.random = TRUE)

# Customize the forest plot
colors <- ifelse(combined_meta$TE > 0, "red", "green")

forest(combined_meta, 
       xlab = "Combined Effect Size (TE)",   
       main = "Forest plot for ANCOM-BC2 and MaAsLin2 method with colour-coded genera and DD associatiation",  
       refline = 0,                 # Reference line at 0
       col.study = "black",         # Color for genus names
       col.square = colors,
       xlim = c(-0.5, 0.5),         # Adjust x-axis limits
       digits = 2,                  # Number of decimal places for effect sizes
       cex = 0.9)                   # Text size adjustment
# Use the color vector for the squares


# Create and save the forest plot
png("forest_plot_full.png", width = 6000, height = 12000, res = 500)

forest(combined_meta, 
       xlab = "Combined Effect Size (TE)",   
       main = "Forest plot for ANCOM-BC2 and MaAsLin2 method with colour-coded genera and DD associatiation",  
       refline = 0,                 # Reference line at 0
       col.study = "black",         # Color for genus names
       col.square = colors,
       xlim = c(-0.5, 0.5),         # Adjust x-axis limits
       digits = 2,                  # Number of decimal places for effect sizes
       cex = 0.9)                   # Text size adjustment
# Use the color vector for the squares
dev.off()

# Separate plots for Maaslin and ANCOM

par(mfrow = c(1, 2))  # Creates a 1x2 grid of plots

colors <- ifelse(meta_result_significant_maaslin$TE > 0, "red", "green")

# Left plot: MaAsLin2 forest plot
forest(meta_result_significant_maaslin$TE,
       sei = meta_result_significant_maaslin$seTE, 
       slab = meta_result_significant_maaslin$studlab, 
       xlab = "Effect Size (TE) - MaAsLin2", 
       main = "MaAsLin2 Results",  
       col.square = colors, 
       refline = 0)

# Right plot: ANCOM-BC2 forest plot
forest(meta_result_significant$TE, 
       sei = meta_result_significant$seTE, 
       slab = meta_result_significant$studlab, 
       xlab = "Effect Size (TE) - ANCOM-BC2", 
       main = "ANCOM-BC2 Results", 
       col.square = "red", 
       refline = 0)

