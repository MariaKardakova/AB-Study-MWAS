library(tidyverse)
library(jsonlite)
library(tidyjson)
library(tidyr)
library(visdat)
library(semPlot)
#library(mediation)
library(stats)
library(ggplot2)
library(boot)
library(pROC)
library(gridExtra)
library(ggpubr)


#DEPRESSION ----

## Microbiome LIST REGRESSION ----

bacteriadata <- read_csv("/Users/MetaAnalysis_Maaslin2_ANCOMBC2_FDR.csv") 
Description <- read_csv("/Users/Displayname.depression.csv" )  
metadata <- read_csv("/Users/mergedDepression.csv" ) %>% dplyr::select(-1)   #7,506 individuals,  238 questions
bacterialist <- read_csv("/Users/differential_analysis_counts.csv") %>% mutate(barcode = paste0("X",barcode))
Description <- Description %>% dplyr::select(bacterianame, Displayname) 

#bacteriadata <- bacteriadata %>% 
#  filter(BH_ANCOMBC2 <= 0.05 | BH_MaAsLin2 <= 0.04)
bacteriadata <- bacteriadata %>% left_join(Description, by = "Displayname")
bacterialist <- bacterialist  %>% 
  filter(barcode %in% metadata$barcode) %>%  #6,942
  column_to_rownames("barcode") %>% 
  dplyr::select(all_of(bacteriadata$bacterianame)) %>% 
  rownames_to_column("barcode")

view(colnames(metadata))

metadata <- metadata %>% left_join(bacterialist, by = "barcode") #6,942 entries, 247 total columns

metadata$Romboutsia = as.numeric(metadata$Bacteria.Firmicutes.Clostridia.Clostridiales.Peptostreptococcaceae.Romboutsia)              
metadata$Clostridium.sensu_stricto = as.numeric(metadata$Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae_1.Clostridium_sensu_stricto_1)   
metadata$Agathobacter = as.numeric(metadata$Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Agathobacter)         
metadata$Victivallis = as.numeric(metadata$Bacteria.Lentisphaerae.Lentisphaeria.Victivallales.Victivallaceae.Victivallis)              
metadata$Flavonifractor = as.numeric(metadata$Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Flavonifractor)                
metadata$Bacteroides = as.numeric(metadata$Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Bacteroidaceae.Bacteroides)                
metadata$Coprococcus = as.numeric(metadata$Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Coprococcus_3)         
metadata$Tyzzerella = as.numeric(metadata$Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Tyzzerella)                
metadata$Bacteroides.pectinophilus_group = as.numeric(metadata$Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae..Bacteroides._pectinophilus_group)
metadata$Bifidobacterium = as.numeric(metadata$Bacteria.Actinobacteria.Actinobacteria.Bifidobacteriales.Bifidobacteriaceae.Bifidobacterium)
metadata$Lactobacillus = as.numeric(metadata$Bacteria.Firmicutes.Bacilli.Lactobacillales.Lactobacillaceae.Lactobacillus)



# Assuming  metadata dataset has bacterial taxa counts, lifestyle variables, and depression status.
# Bacterial taxa: Romboutsia, Bacteroides, etc.
# Lifestyle variables: qSmokingStatus, qEatFiber, qAlcoholUnitsPerWeek, etc.
# Depression outcome: qMedicalHistoryDepression (2 = Depression, 1 = No Depression)

# Define the SEM model
# Lifestyle variables impact bacterial taxa, and bacterial taxa impact depression
# Also I include direct effects estimation of lifestyle variables on depression


#list of trait associated with depression: 

Depression_metadata <- read.csv("/Users/allregression_multipletesting_Depression.csv") #%>% filter(BH_Depression <= 0.05) %>% select(DisplayName, estimate_Depression = estimate, BH_Depression)

#list of trait associated with bacteria associated with depression based on meta-analysis: 

Romboutsia_metadata <- read_csv("/Users/Bacteria/allregression2024_multipletesting_Romboutsia.csv") %>% dplyr::select (DisplayName, Romboutsia_estimate = estimate, BH_Romboutsia) %>% filter(BH_Romboutsia <= 0.05)
Clostridium.sensu_stricto_metadata <- read_csv("/Users/Bacteria/allregression2024_multipletesting_Clostridium.sensu_stricto.csv") %>% dplyr::select (DisplayName, Clostridium.sensu_stricto_estimate = estimate, BH_Clostridium.sensu_stricto) %>% filter(BH_Clostridium.sensu_stricto <= 0.05)
Agathobacter_metadata <- read_csv("/Users/Bacteria/allregression2024_multipletesting_Agathobacter.csv") %>% dplyr::select (DisplayName, Agathobacter_estimate = estimate, BH_Agathobacter) %>% filter(BH_Agathobacter <= 0.05)
Victivallis_metadata <- read_csv("/Users/Bacteria/allregression2024_multipletesting_Victivallis.csv") %>% dplyr::select (DisplayName, Victivallis_estimate = estimate, BH_Victivallis) %>% filter(BH_Victivallis <= 0.05)
Flavonifractor_metadata <- read_csv("/Users/Bacteria/allregression2024_multipletesting_Flavonifractor.csv") %>% dplyr::select (DisplayName, Flavonifractor_estimate = estimate, BH_Flavonifractor) %>% filter(BH_Flavonifractor <= 0.05)
Bacteroides_metadata <- read_csv("/Users/Bacteria/allregression2024_multipletesting_Bacteroides.csv") %>% dplyr::select (DisplayName, Bacteroides_estimate = estimate, BH_Bacteroides) %>% filter(BH_Bacteroides <= 0.05)
Coprococcus_metadata <- read_csv("/Users/Bacteria/allregression2024_multipletesting_Coprococcus.csv") %>% dplyr::select (DisplayName, Coprococcus_estimate = estimate, BH_Coprococcus) %>% filter(BH_Coprococcus <= 0.05)
Tyzzerella_metadata <- read_csv("/Users/Bacteria/allregression2024_multipletesting_Tyzzerella.csv") %>% dplyr::select (DisplayName, Tyzzerella_estimate = estimate, BH_Tyzzerella) %>% filter(BH_Tyzzerella <= 0.05)
Bacteroides.pectinophilus_group_metadata <- read_csv("/Users/Bacteria/allregression2024_multipletesting_Bacteroides.pectinophilus_group.csv") %>% dplyr::select (DisplayName, Bacteroides.pectinophilus_group_estimate = estimate, BH_Bacteroides.pectinophilus_group) %>% filter(BH_Bacteroides.pectinophilus_group <= 0.05)

# _______________----


#MA ----
view(metadata$qDrugHistoryProbioticsType_bifidobacterium) 
metadata$qMedicalHistoryDepression <- ifelse(metadata$qMedicalHistoryDepression == "2", 1, 0)
metadata$Healthyweight.status <- ifelse(metadata$Healthyweight.status == "2", 1, 0)
metadata$Obesity.status <- ifelse(metadata$Obesity.status == "2", 1, 0)

#food

metadata$qEatFiber <- ifelse(metadata$qEatFiber == "2", 1, 0)
metadata$qDrinkSweetBeverages <- ifelse(metadata$qDrinkSweetBeverages == "2", 1, 0)
metadata$qEatFoodRichInVitC <- ifelse(metadata$qEatFoodRichInVitC == "2", 1, 0)
metadata$qEatOliveOil <- ifelse(metadata$qEatOliveOil == "2", 1, 0)
metadata$FrVeg_less_5_a_day <- ifelse(metadata$FrVeg_less_5_a_day == "2", 1, 0)

#Health

metadata$qMedicalHistoryLactoseIntolerance <- ifelse(metadata$qMedicalHistoryLactoseIntolerance == "2", 1, 0)
metadata$qMedicalHistorySmallIntestinalBacterialOvergrowth <- ifelse(metadata$qMedicalHistorySmallIntestinalBacterialOvergrowth == "2", 1, 0)
metadata$qMedicalHistoryGastroesophagealRefluxDisease <- ifelse(metadata$qMedicalHistoryGastroesophagealRefluxDisease == "2", 1, 0)
metadata$qMedicalHistoryAtopicDermatitis <- ifelse(metadata$qMedicalHistoryAtopicDermatitis == "2", 1, 0)
metadata$qMedicalHistoryBloating <- ifelse(metadata$qMedicalHistoryBloating == "2", 1, 0)
metadata$qMedicalHistoryMigraine <- ifelse(metadata$qMedicalHistoryMigraine == "2", 1, 0)
metadata$qMedicalHistoryHypothyroidism <- ifelse(metadata$qMedicalHistoryHypothyroidism == "2", 1, 0)
metadata$qMedicalHistoryHighBloodPressure <- ifelse(metadata$qMedicalHistoryHighBloodPressure == "2", 1, 0)

#Lifestyle

metadata$qDrugHistoryProbioticsType_bifidobacterium <- ifelse(metadata$qDrugHistoryProbioticsType_bifidobacterium == "2", 1, 0)
metadata$qDrugHistoryProbioticsType_lactobacillus <- ifelse(metadata$qDrugHistoryProbioticsType_lactobacillus == "2", 1, 0)
metadata$qPArecreationalModeratePresent <- ifelse(metadata$qPArecreationalModeratePresent == "2", 1, 0)
metadata$qPArecreationalVigorousPresent <- ifelse(metadata$qPArecreationalVigorousPresent == "2", 1, 0)
metadata$qPAtravelBicycleUse <- ifelse(metadata$qPAtravelBicycleUse == "2", 1, 0)


#Romboutsia ----

#Microbal alpha-diversity "+"

mediator_model <- lm(alpha_diversity ~ Romboutsia + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Romboutsia + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "Romboutsia", mediator = "alpha_diversity", boot = TRUE)
summary(mediation_result)


#Agathobacter ----

#Microbal alpha-diversity "+"

mediator_model <- lm(alpha_diversity ~ Agathobacter + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Agathobacter + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "Agathobacter", mediator = "alpha_diversity", boot = TRUE)
summary(mediation_result)

mediator_model <- lm(Agathobacter ~ alpha_diversity + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Agathobacter + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "alpha_diversity", mediator = "Agathobacter", boot = TRUE)
summary(mediation_result)

#Victivallis ----

#Microbal alpha-diversity "+"

mediator_model <- lm(alpha_diversity ~ Victivallis + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Victivallis + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "Victivallis", mediator = "alpha_diversity", boot = TRUE)
summary(mediation_result)

mediator_model <- lm(Victivallis ~ alpha_diversity + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Victivallis + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "alpha_diversity", mediator = "Victivallis", boot = TRUE)
summary(mediation_result)

#BMI "-"

mediator_model <- lm(BMI ~ Victivallis + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ BMI + Victivallis + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "Victivallis", mediator = "BMI", boot = TRUE)
summary(mediation_result)


mediator_model <- lm(BMI ~ Victivallis + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ BMI + Victivallis + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "BMI", mediator = "Victivallis", boot = TRUE)
summary(mediation_result)

#Obesity "-"  - EXPOSURE

mediator_model <- lm(Victivallis ~ Obesity.status + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ Obesity.status + Victivallis + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "Obesity.status", mediator = "Victivallis", boot = TRUE)
summary(mediation_result)

# Obesity - outcome, Depression - EXPOSURE
mediator_model <- lm(Victivallis ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(Obesity.status ~ qMedicalHistoryDepression + Victivallis + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Victivallis", boot = TRUE)
summary(mediation_result)

#Healthy weight "+"  - EXPOSURE

mediator_model <- lm(Victivallis ~ Healthyweight.status + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ Healthyweight.status + Victivallis + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "Healthyweight.status", mediator = "Victivallis", boot = TRUE)
summary(mediation_result)

# Healthy weight, Depression - EXPOSURE

mediator_model <- lm(Victivallis ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(Healthyweight.status ~ qMedicalHistoryDepression + Victivallis + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Victivallis", boot = TRUE)
summary(mediation_result)

#Flavonifractor ----

#Microbal alpha-diversity "+"

mediator_model <- lm(alpha_diversity ~ Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "Flavonifractor", mediator = "alpha_diversity", boot = TRUE)
summary(mediation_result)

mediator_model <- lm(Flavonifractor ~ alpha_diversity + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "alpha_diversity", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

#Eat fibre rich foods "-"  - EXPOSURE

mediator_model <- lm(Flavonifractor ~ qEatFiber + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ qEatFiber + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qEatFiber", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

# Eat Fiber - OUTCOME, Depression - EXPOSURE

metadata_Fibre <- metadata %>% filter(!is.na(qEatFiber))

mediator_model <- lm(Flavonifractor ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre)
outcome_model <- glm(qEatFiber ~ qMedicalHistoryDepression + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

#Drinks sweet beverages regularly "+" - EXPOSURE

mediator_model <- lm(Flavonifractor ~ qDrinkSweetBeverages + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ qDrinkSweetBeverages + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qDrinkSweetBeverages", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

# Drink Sweet Beverages - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Flavonifractor ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qDrinkSweetBeverages ~ qMedicalHistoryDepression + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

#Physical activity - traveling by bicycle "-" - EXPOSURE

metadata_Bicycle <- metadata %>% filter(!is.na(qPAtravelBicycleUse))

mediator_model <- lm(Flavonifractor ~ qPAtravelBicycleUse + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Bicycle)
outcome_model <- glm(qMedicalHistoryDepression ~ qPAtravelBicycleUse + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Bicycle, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qPAtravelBicycleUse", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

# Physical activity - traveling by bicycle - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Flavonifractor ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Bicycle)
outcome_model <- glm(qPAtravelBicycleUse ~ qMedicalHistoryDepression + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Bicycle, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

#Had lactobacillus probiotics taken in the past "+" EXPOSURE

metadata_lactobacillus <- metadata %>% filter(!is.na(qDrugHistoryProbioticsType_lactobacillus))

mediator_model <- lm(Flavonifractor ~ qDrugHistoryProbioticsType_lactobacillus + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_lactobacillus)
outcome_model <- glm(qMedicalHistoryDepression ~ qDrugHistoryProbioticsType_lactobacillus + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_lactobacillus, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qDrugHistoryProbioticsType_lactobacillus", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

# lactobacillus probiotics - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Flavonifractor ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_lactobacillus)
outcome_model <- glm(qDrugHistoryProbioticsType_lactobacillus ~ qMedicalHistoryDepression + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_lactobacillus, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

#Had bifidobacterium probiotics taken in the past "+" EXPOSURE

metadata_bifido <- metadata %>% filter(!is.na(qDrugHistoryProbioticsType_bifidobacterium))
view(metadata$qDrugHistoryProbioticsType_bifidobacterium)

mediator_model <- lm(Flavonifractor ~ qDrugHistoryProbioticsType_bifidobacterium + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_bifido)
outcome_model <- glm(qMedicalHistoryDepression ~ qDrugHistoryProbioticsType_bifidobacterium + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_bifido, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qDrugHistoryProbioticsType_bifidobacterium", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

# bifidobacterium probiotics - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Flavonifractor ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_bifido)
outcome_model <- glm(qDrugHistoryProbioticsType_bifidobacterium ~ qMedicalHistoryDepression + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_bifido, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

#Physical activity - moderate, recreational "-" EXPOSURE

metadata_physio <- metadata %>% filter(!is.na(qPArecreationalModeratePresent))

mediator_model <- lm(Flavonifractor ~ qPArecreationalModeratePresent + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_physio)
outcome_model <- glm(qMedicalHistoryDepression ~ qPArecreationalModeratePresent + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_physio, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qPArecreationalModeratePresent", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

# Physical activity - moderate, recreational - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Flavonifractor ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_physio)
outcome_model <- glm(qPArecreationalModeratePresent ~ qMedicalHistoryDepression + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_physio, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

#Use olive oil regularly “-" EXPOSURE

metadata_oil <- metadata %>% filter(!is.na(qEatOliveOil))

mediator_model <- lm(Flavonifractor ~ qEatOliveOil + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_oil)
outcome_model <- glm(qMedicalHistoryDepression ~ qEatOliveOil + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_oil, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qEatOliveOil", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

# Use olive oil regularly - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Flavonifractor ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_oil)
outcome_model <- glm(qEatOliveOil ~ qMedicalHistoryDepression + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_oil, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

#Smoking “+” EXPOSURE

metadata_smoking <- metadata %>% filter(!is.na(qSmokingStatus))

mediator_model <- lm(Flavonifractor ~ qSmokingStatus + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_smoking)
outcome_model <- glm(qMedicalHistoryDepression ~ qSmokingStatus + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_smoking, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qSmokingStatus", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

# Smoking - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Flavonifractor ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_smoking)
outcome_model <- lm(qSmokingStatus ~ qMedicalHistoryDepression + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_smoking)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

#Adding salt “+” EXPOSURE

metadata_salt <- metadata %>% filter(!is.na(qMoreSalt))

mediator_model <- lm(Flavonifractor ~ qMoreSalt + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_salt)
outcome_model <- glm(qMedicalHistoryDepression ~ qMoreSalt + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_salt, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMoreSalt", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

# Adding salt - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Flavonifractor ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_salt)
outcome_model <- lm(qMoreSalt ~ qMedicalHistoryDepression + Flavonifractor + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_salt)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Flavonifractor", boot = TRUE)
summary(mediation_result)

#Bacteroides ----

#Microbal alpha-diversity "-" EXPOSURE

mediator_model <- lm(alpha_diversity ~ Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "Bacteroides", mediator = "alpha_diversity", boot = TRUE)
summary(mediation_result)

mediator_model <- lm(Bacteroides ~ alpha_diversity + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "alpha_diversity", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

#Physical activity - vigorous, recreational "-" is EXPOSURE

metadata_PAvigorous <- metadata %>% filter(!is.na(qPArecreationalVigorousPresent))

mediator_model <- lm(Bacteroides ~ qPArecreationalVigorousPresent + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_PAvigorous)
outcome_model <- glm(qMedicalHistoryDepression ~ qPArecreationalVigorousPresent + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_PAvigorous, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qPArecreationalVigorousPresent", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

# Physical activity - Vigorous, recreational - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_PAvigorous)
outcome_model <- glm(qPArecreationalVigorousPresent ~ qMedicalHistoryDepression + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_PAvigorous, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)


#Sleep from 7 to 9 hours "-" EXPOSURE

mediator_model <- lm(Bacteroides ~ normal_sleep + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ normal_sleep + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "normal_sleep", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

# Sleep from 7 to 9 hours - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(normal_sleep ~ qMedicalHistoryDepression + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)


#Obesity  "+"  EXPOSURE

mediator_model <- lm(Bacteroides ~ Obesity.status + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ Obesity.status + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "Obesity.status", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

# Obesity - outcome, Depression - EXPOSURE
mediator_model <- lm(Bacteroides ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(Obesity.status ~ qMedicalHistoryDepression + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)


#Sleep less than 7 hours "+"  EXPOSURE

mediator_model <- lm(Bacteroides ~ short_sleep + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ short_sleep + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "short_sleep", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

# Sleep less than  7 hours - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(normal_sleep ~ qMedicalHistoryDepression + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)


#Eat fibre rich foods "-" EXPOSURE

metadata_Fibre <- metadata %>% filter(!is.na(qEatFiber))

mediator_model <- lm(Bacteroides ~ qEatFiber + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre)
outcome_model <- glm(qMedicalHistoryDepression ~ qEatFiber + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qEatFiber", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

# Eat fibre rich foods - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre)
outcome_model <- glm(qEatFiber ~ qMedicalHistoryDepression + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)


#Physical activity - moderate, recreational "-" EXPOSURE

metadata_PAmoderate <- metadata %>% filter(!is.na(qPArecreationalModeratePresent))

mediator_model <- lm(Bacteroides ~ qPArecreationalModeratePresent + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_PAmoderate)
outcome_model <- glm(qMedicalHistoryDepression ~ qPArecreationalModeratePresent + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_PAmoderate, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qPArecreationalModeratePresent", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

#  Physical activity - moderate, recreational - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_PAmoderate)
outcome_model <- glm(qPArecreationalModeratePresent ~ qMedicalHistoryDepression + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_PAmoderate, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

#Healthy weight "-" EXPOSURE

mediator_model <- lm(Bacteroides ~ Healthyweight.status + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ Healthyweight.status + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "Healthyweight.status", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

#  Healthy weight - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(Healthyweight.status ~ qMedicalHistoryDepression + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

#Physical activity - traveling by bicycle "-"EXPOSURE

metadata_PAbicycle <- metadata %>% filter(!is.na(qPAtravelBicycleUse))

mediator_model <- lm(Bacteroides ~ qPAtravelBicycleUse + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_PAbicycle)
outcome_model <- glm(qMedicalHistoryDepression ~ qPAtravelBicycleUse + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_PAbicycle, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qPAtravelBicycleUse", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

#  travel Bicycle Use - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_PAbicycle)
outcome_model <- glm(qPAtravelBicycleUse ~ qMedicalHistoryDepression + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_PAbicycle, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

#Eat Vitamin C rich food "-" EXPOSURE

metadata_VitC <- metadata %>% filter(!is.na(qEatFoodRichInVitC))

mediator_model <- lm(Bacteroides ~ qEatFoodRichInVitC + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_VitC)
outcome_model <- glm(qMedicalHistoryDepression ~ qEatFoodRichInVitC + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_VitC, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qEatFoodRichInVitC", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

# Eat Food Rich In VitC - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_VitC)
outcome_model <- glm(qEatFoodRichInVitC ~ qMedicalHistoryDepression + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_VitC, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

#Drinks sweet beverages regularly "+"  EXPOSURE

metadata_SweetBev <- metadata %>% filter(!is.na(qDrinkSweetBeverages))

mediator_model <- lm(Bacteroides ~ qDrinkSweetBeverages + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_SweetBev)
outcome_model <- glm(qMedicalHistoryDepression ~ qDrinkSweetBeverages + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_SweetBev, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qDrinkSweetBeverages", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

# qDrinkSweetBeverages - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_SweetBev)
outcome_model <- glm(qDrinkSweetBeverages ~ qMedicalHistoryDepression + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_SweetBev, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

#Use olive oil regularly “-" EXPOSURE

metadata_oil <- metadata %>% filter(!is.na(qEatOliveOil))

mediator_model <- lm(Bacteroides ~ qEatOliveOil + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_oil)
outcome_model <- glm(qMedicalHistoryDepression ~ qEatOliveOil + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_oil, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qEatOliveOil", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

# Eat Olive Oil - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_oil)
outcome_model <- glm(qEatOliveOil ~ qMedicalHistoryDepression + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_oil, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

# Smoking "+" EXPOSURE

metadata_smoking <- metadata %>% filter(!is.na(qSmokingStatus))

mediator_model <- lm(Bacteroides ~ qSmokingStatus + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_smoking)
outcome_model <- glm(qMedicalHistoryDepression ~ qSmokingStatus + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_smoking, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qSmokingStatus", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)

# Smoking Status - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_smoking)
outcome_model <- lm(qSmokingStatus ~ qMedicalHistoryDepression + Bacteroides + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_smoking)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides", boot = TRUE)
summary(mediation_result)


#Bacteroides.pectinophilus_group ----

#Microbal alpha-diversity "+"  EXPOSURE

mediator_model <- lm(alpha_diversity ~ Bacteroides.pectinophilus_group + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Bacteroides.pectinophilus_group + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "Bacteroides.pectinophilus_group", mediator = "alpha_diversity", boot = TRUE)
summary(mediation_result)

mediator_model <- lm(Bacteroides.pectinophilus_group ~ alpha_diversity + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Bacteroides.pectinophilus_group + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "alpha_diversity", mediator = "Bacteroides.pectinophilus_group", boot = TRUE)
summary(mediation_result)

#Eat fibre rich foods "+"   EXPOSURE

mediator_model <- lm(Bacteroides.pectinophilus_group ~ qEatFiber + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre)
outcome_model <- glm(qMedicalHistoryDepression ~ qEatFiber + Bacteroides.pectinophilus_group + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qEatFiber", mediator = "Bacteroides.pectinophilus_group", boot = TRUE)
summary(mediation_result)

# qEatFiber - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides.pectinophilus_group ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre)
outcome_model <- glm(qEatFiber ~ qMedicalHistoryDepression + Bacteroides.pectinophilus_group + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides.pectinophilus_group", boot = TRUE)
summary(mediation_result)

#Smoking "-" EXPOSURE

mediator_model <- lm(Bacteroides.pectinophilus_group ~ qSmokingStatus + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_smoking)
outcome_model <- glm(qMedicalHistoryDepression ~ qSmokingStatus + Bacteroides.pectinophilus_group + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_smoking, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qSmokingStatus", mediator = "Bacteroides.pectinophilus_group", boot = TRUE)
summary(mediation_result)

# Smoking - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides.pectinophilus_group ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_smoking)
outcome_model <- lm(qSmokingStatus ~ qMedicalHistoryDepression + Bacteroides.pectinophilus_group + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_smoking)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides.pectinophilus_group", boot = TRUE)
summary(mediation_result)

#Eat less than five portions of fruits and vegetables "-" EXPOSURE

metadata_less_5 <- metadata %>% filter(!is.na(FrVeg_less_5_a_day))

mediator_model <- lm(Bacteroides.pectinophilus_group ~ FrVeg_less_5_a_day + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_less_5)
outcome_model <- glm(qMedicalHistoryDepression ~ FrVeg_less_5_a_day + Bacteroides.pectinophilus_group + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_less_5, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "FrVeg_less_5_a_day", mediator = "Bacteroides.pectinophilus_group", boot = TRUE)
summary(mediation_result)

# qEatFiber - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Bacteroides.pectinophilus_group ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_less_5)
outcome_model <- glm(FrVeg_less_5_a_day ~ qMedicalHistoryDepression + Bacteroides.pectinophilus_group + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_less_5, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Bacteroides.pectinophilus_group", boot = TRUE)
summary(mediation_result)

#Coprococcus ----

#Microbal alpha-diversity "+" 

mediator_model <- lm(alpha_diversity ~ Coprococcus + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Coprococcus + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "Coprococcus", mediator = "alpha_diversity", boot = TRUE)
summary(mediation_result)

mediator_model <- lm(Coprococcus ~ alpha_diversity + age + qGender + reg + weeks.since.testing.microbiome, data = metadata)
outcome_model <- glm(qMedicalHistoryDepression ~ alpha_diversity + Coprococcus + age + qGender + reg + weeks.since.testing.microbiome, data = metadata, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "alpha_diversity", mediator = "Coprococcus", boot = TRUE)
summary(mediation_result)

#Drink coffee “+” EXPOSURE

metadata_coffee <- metadata %>% filter(!is.na(qCoffeeFrequency))

mediator_model <- lm(Coprococcus ~ qCoffeeFrequency + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_coffee)
outcome_model <- glm(qMedicalHistoryDepression ~ qCoffeeFrequency + Coprococcus + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_coffee, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qCoffeeFrequency", mediator = "Coprococcus", boot = TRUE)
summary(mediation_result)

# Coffee Frequency - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Coprococcus ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_coffee)
outcome_model <- lm(qCoffeeFrequency ~ qMedicalHistoryDepression + Coprococcus + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_coffee)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Coprococcus", boot = TRUE)
summary(mediation_result)

#Eat meat “+” EXPOSURE

metadata_meat <- metadata %>% filter(!is.na(qEatMeatRegularly))

mediator_model <- lm(Coprococcus ~ qEatMeatRegularly + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_meat)
outcome_model <- glm(qMedicalHistoryDepression ~ qEatMeatRegularly + Coprococcus + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_meat, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qEatMeatRegularly", mediator = "Coprococcus", boot = TRUE)
summary(mediation_result)

# Eat Meat Regularly - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Coprococcus ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_meat)
outcome_model <- lm(qEatMeatRegularly ~ qMedicalHistoryDepression + Coprococcus + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_meat)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Coprococcus", boot = TRUE)
summary(mediation_result)

#Tyzzerella ----

#Eat fibre rich foods "-"    EXPOSURE

mediator_model <- lm(Tyzzerella ~ qEatFiber + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre)
outcome_model <- glm(qMedicalHistoryDepression ~ qEatFiber + Tyzzerella + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qEatFiber", mediator = "Tyzzerella", boot = TRUE)
summary(mediation_result)

# Eat fibre rich foods - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Tyzzerella ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre)
outcome_model <- glm(qEatFiber ~ qMedicalHistoryDepression + Tyzzerella + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_Fibre, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Tyzzerella", boot = TRUE)
summary(mediation_result)

# Eat less than five portions of fruits and vegetables "+" EXPOSURE

mediator_model <- lm(Tyzzerella ~ FrVeg_less_5_a_day + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_less_5)
outcome_model <- glm(qMedicalHistoryDepression ~ FrVeg_less_5_a_day + Tyzzerella + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_less_5, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "FrVeg_less_5_a_day", mediator = "Tyzzerella", boot = TRUE)
summary(mediation_result)

# Eat less than five portions of fruits and vegetables - OUTCOME, Depression - EXPOSURE

mediator_model <- lm(Tyzzerella ~ qMedicalHistoryDepression + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_less_5)
outcome_model <- glm(FrVeg_less_5_a_day ~ qMedicalHistoryDepression + Tyzzerella + age + qGender + reg + weeks.since.testing.microbiome, data = metadata_less_5, family = binomial)

mediation_result <- mediate(mediator_model, outcome_model, treat = "qMedicalHistoryDepression", mediator = "Tyzzerella", boot = TRUE)
summary(mediation_result)


# _______________----

#MR ----



## Tyzzerella outcome ----

#Prepare exposure data
# Rename the columns to match the expected format for MR analysis
depression_tyzzerella_snp <- tyzzerella_from_depression_SNPs %>%
  rename(
    SNP = MarkerName,
    beta.exposure = LogOR,
    se.exposure = StdErrLogOR,
    effect_allele.exposure = A1,
    other_allele.exposure = A2,
    pval.exposure = P
  )

# Add the required columns for identification
depression_tyzzerella_snp$id.exposure <- "Depression"
depression_tyzzerella_snp$exposure <- "Depression"

# #Prepare outcome data
# Rename the columns to match the expected format

# Find the common SNPs between the two datasets
common_snps <- intersect(depression_tyzzerella_snp$SNP, tyzzerella_depression_snp$SNP)
# Filter both datasets to retain only common SNPs
depression_tyzzerella_snp_filtered <- depression_tyzzerella_snp[depression_tyzzerella_snp$SNP %in% common_snps, ]
tyzzerella_depression_snp_filtered <- tyzzerella_depression_snp[tyzzerella_depression_snp$SNP %in% common_snps, ]



harmonized_data <- harmonise_data(
  exposure_dat = depression_tyzzerella_snp_filtered,  # Depression as exposure
  outcome_dat =  tyzzerella_depression_snp_filtered   # Tyzzerella as outcome
)

# Ensure that se.exposure and se.outcome are numeric
harmonized_data$beta.exposure <- as.numeric(harmonized_data$beta.exposure)
harmonized_data$se.exposure <- as.numeric(harmonized_data$se.exposure)
harmonized_data$beta.outcome <- as.numeric(harmonized_data$beta.outcome)
harmonized_data$se.outcome <- as.numeric(harmonized_data$se.outcome)

harmonized_data <- harmonized_data[!is.na(harmonized_data$se.exposure) & !is.na(harmonized_data$se.outcome), ]
str(harmonized_data)

# Perform MR analysis with the harmonized data
mr_results <- mr(harmonized_data, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
# Record MR results
Tyzzerella_outcome <- mr_results

## Victivallis outcome ----

#Prepare exposure data
# Rename the columns to match the expected format for MR analysis
depression_victivallis_snp <- victivallis_from_depression_SNPs %>%
  rename(
    SNP = MarkerName,
    beta.exposure = LogOR,
    se.exposure = StdErrLogOR,
    effect_allele.exposure = A1,
    other_allele.exposure = A2,
    pval.exposure = P
  )

# Add the required columns for identification
depression_victivallis_snp$id.exposure <- "Depression"
depression_victivallis_snp$exposure <- "Depression"

# #Prepare outcome data
# Rename the columns to match the expected format

# Find the common SNPs between the two datasets
common_snps <- intersect(depression_victivallis_snp$SNP, victivallis_depression_snp$SNP)
# Filter both datasets to retain only common SNPs
depression_victivallis_snp_filtered <- depression_victivallis_snp[depression_victivallis_snp$SNP %in% common_snps, ]
victivallis_depression_snp_filtered <- victivallis_depression_snp[victivallis_depression_snp$SNP %in% common_snps, ]

if (!"eaf.exposure" %in% colnames(depression_victivallis_snp_filtered)) {
  depression_victivallis_snp_filtered$eaf.exposure <- NA
}
if (!"eaf.outcome" %in% colnames(victivallis_depression_snp_filtered)) {
  victivallis_depression_snp_filtered$eaf.outcome <- NA
}

harmonized_data <- harmonise_data(
  exposure_dat = depression_victivallis_snp_filtered,  # Depression as exposure
  outcome_dat =  victivallis_depression_snp_filtered   # Victivallis as outcome
)

# Ensure that se.exposure and se.outcome are numeric
harmonized_data$beta.exposure <- as.numeric(harmonized_data$beta.exposure)
harmonized_data$se.exposure <- as.numeric(harmonized_data$se.exposure)
harmonized_data$beta.outcome <- as.numeric(harmonized_data$beta.outcome)
harmonized_data$se.outcome <- as.numeric(harmonized_data$se.outcome)

harmonized_data <- harmonized_data[!is.na(harmonized_data$se.exposure) & !is.na(harmonized_data$se.outcome), ]
str(harmonized_data)

# Perform MR analysis with the harmonized data
mr_results <- mr(harmonized_data, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
# Record MR results
Victivallis_outcome <- mr_results


## Bacteroides outcome ----

# #Prepare outcome data
# Rename the columns to match the expected format

# Find the common SNPs between the two datasets
common_snps <- intersect(depression_bacteroides_snp$SNP, bacteroides_depression_snp$SNP)
# Filter both datasets to retain only common SNPs
depression_bacteroides_snp_filtered <- depression_bacteroides_snp[depression_bacteroides_snp$SNP %in% common_snps, ]
bacteroides_depression_snp_filtered <- bacteroides_depression_snp[bacteroides_depression_snp$SNP %in% common_snps, ]

if (!"eaf.exposure" %in% colnames(depression_bacteroides_snp_filtered)) {
  depression_bacteroides_snp_filtered$eaf.exposure <- NA
}
if (!"eaf.outcome" %in% colnames(bacteroides_depression_snp_filtered)) {
  bacteroides_depression_snp_filtered$eaf.outcome <- NA
}

harmonized_data <- harmonise_data(
  exposure_dat = depression_bacteroides_snp_filtered,  # Depression as exposure
  outcome_dat =  bacteroides_depression_snp_filtered   # Bacteroides as outcome
)

# Ensure that se.exposure and se.outcome are numeric
harmonized_data$beta.exposure <- as.numeric(harmonized_data$beta.exposure)
harmonized_data$se.exposure <- as.numeric(harmonized_data$se.exposure)
harmonized_data$beta.outcome <- as.numeric(harmonized_data$beta.outcome)
harmonized_data$se.outcome <- as.numeric(harmonized_data$se.outcome)

harmonized_data <- harmonized_data[!is.na(harmonized_data$se.exposure) & !is.na(harmonized_data$se.outcome), ]
str(harmonized_data)

# Perform MR analysis with the harmonized data
mr_results <- mr(harmonized_data, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
# Record MR results
bacteroides_outcome <- mr_results


## Romboutsia outcome ----

#Prepare exposure data
# Rename the columns to match the expected format for MR analysis
depression_romboutsia_snp <- romboutsia_from_depression_SNPs %>%
  rename(
    SNP = MarkerName,
    beta.exposure = LogOR,
    se.exposure = StdErrLogOR,
    effect_allele.exposure = A1,
    other_allele.exposure = A2,
    pval.exposure = P
  )

# Add the required columns for identification
depression_romboutsia_snp$id.exposure <- "Depression"
depression_romboutsia_snp$exposure <- "Depression"

# #Prepare outcome data
# Rename the columns to match the expected format

# Find the common SNPs between the two datasets
common_snps <- intersect(depression_romboutsia_snp$SNP, romboutsia_depression_snp$SNP)
# Filter both datasets to retain only common SNPs
depression_romboutsia_snp_filtered <- depression_romboutsia_snp[depression_romboutsia_snp$SNP %in% common_snps, ]
romboutsia_depression_snp_filtered <- romboutsia_depression_snp[romboutsia_depression_snp$SNP %in% common_snps, ]

if (!"eaf.exposure" %in% colnames(depression_romboutsia_snp_filtered)) {
  depression_romboutsia_snp_filtered$eaf.exposure <- NA
}
if (!"eaf.outcome" %in% colnames(romboutsia_depression_snp_filtered)) {
  romboutsia_depression_snp_filtered$eaf.outcome <- NA
}

harmonized_data <- harmonise_data(
  exposure_dat = depression_romboutsia_snp_filtered,  # Depression as exposure
  outcome_dat =  romboutsia_depression_snp_filtered   # Romboutsia as outcome
)

# Ensure that se.exposure and se.outcome are numeric
harmonized_data$beta.exposure <- as.numeric(harmonized_data$beta.exposure)
harmonized_data$se.exposure <- as.numeric(harmonized_data$se.exposure)
harmonized_data$beta.outcome <- as.numeric(harmonized_data$beta.outcome)
harmonized_data$se.outcome <- as.numeric(harmonized_data$se.outcome)

harmonized_data <- harmonized_data[!is.na(harmonized_data$se.exposure) & !is.na(harmonized_data$se.outcome), ]
str(harmonized_data)

# Perform MR analysis with the harmonized data
mr_results <- mr(harmonized_data, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
# Record MR results
romboutsia_outcome <- mr_results


## Flavonifractor outcome ----

#Prepare exposure data
# Rename the columns to match the expected format for MR analysis
depression_flavonifractor_snp <- flavonifractor_from_depression_SNPs %>%
  rename(
    SNP = MarkerName,
    beta.exposure = LogOR,
    se.exposure = StdErrLogOR,
    effect_allele.exposure = A1,
    other_allele.exposure = A2,
    pval.exposure = P
  )

# Add the required columns for identification
depression_flavonifractor_snp$id.exposure <- "Depression"
depression_flavonifractor_snp$exposure <- "Depression"

# #Prepare outcome data
# Rename the columns to match the expected format

# Find the common SNPs between the two datasets
common_snps <- intersect(depression_flavonifractor_snp$SNP, flavonifractor_depression_snp$SNP)
# Filter both datasets to retain only common SNPs
depression_flavonifractor_snp_filtered <- depression_flavonifractor_snp[depression_flavonifractor_snp$SNP %in% common_snps, ]
flavonifractor_depression_snp_filtered <- flavonifractor_depression_snp[flavonifractor_depression_snp$SNP %in% common_snps, ]

if (!"eaf.exposure" %in% colnames(depression_flavonifractor_snp_filtered)) {
  depression_flavonifractor_snp_filtered$eaf.exposure <- NA
}
if (!"eaf.outcome" %in% colnames(flavonifractor_depression_snp_filtered)) {
  flavonifractor_depression_snp_filtered$eaf.outcome <- NA
}

harmonized_data <- harmonise_data(
  exposure_dat = depression_flavonifractor_snp_filtered,  # Depression as exposure
  outcome_dat =  flavonifractor_depression_snp_filtered   # Flavonifractor as outcome
)

# Ensure that se.exposure and se.outcome are numeric
harmonized_data$beta.exposure <- as.numeric(harmonized_data$beta.exposure)
harmonized_data$se.exposure <- as.numeric(harmonized_data$se.exposure)
harmonized_data$beta.outcome <- as.numeric(harmonized_data$beta.outcome)
harmonized_data$se.outcome <- as.numeric(harmonized_data$se.outcome)

harmonized_data <- harmonized_data[!is.na(harmonized_data$se.exposure) & !is.na(harmonized_data$se.outcome), ]
str(harmonized_data)

# Perform MR analysis with the harmonized data
mr_results <- mr(harmonized_data, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
# Record MR results
flavonifractor_outcome <- mr_results


## Coprococcus outcome ----

#Prepare exposure data
# Rename the columns to match the expected format for MR analysis
depression_coprococcus_snp <- coprococcus_from_depression_SNPs %>%
  rename(
    SNP = MarkerName,
    beta.exposure = LogOR,
    se.exposure = StdErrLogOR,
    effect_allele.exposure = A1,
    other_allele.exposure = A2,
    pval.exposure = P
  )

# Add the required columns for identification
depression_coprococcus_snp$id.exposure <- "Depression"
depression_coprococcus_snp$exposure <- "Depression"

# #Prepare outcome data
# Rename the columns to match the expected format


# Find the common SNPs between the two datasets
common_snps <- intersect(depression_coprococcus_snp$SNP, coprococcus_depression_snp$SNP)
# Filter both datasets to retain only common SNPs
depression_coprococcus_snp_filtered <- depression_coprococcus_snp[depression_coprococcus_snp$SNP %in% common_snps, ]
coprococcus_depression_snp_filtered <- coprococcus_depression_snp[coprococcus_depression_snp$SNP %in% common_snps, ]

if (!"eaf.exposure" %in% colnames(depression_coprococcus_snp_filtered)) {
  depression_coprococcus_snp_filtered$eaf.exposure <- NA
}
if (!"eaf.outcome" %in% colnames(coprococcus_depression_snp_filtered)) {
  coprococcus_depression_snp_filtered$eaf.outcome <- NA
}

harmonized_data <- harmonise_data(
  exposure_dat = depression_coprococcus_snp_filtered,  # Depression as exposure
  outcome_dat =  coprococcus_depression_snp_filtered   # Coprococcus as outcome
)

# Ensure that se.exposure and se.outcome are numeric
harmonized_data$beta.exposure <- as.numeric(harmonized_data$beta.exposure)
harmonized_data$se.exposure <- as.numeric(harmonized_data$se.exposure)
harmonized_data$beta.outcome <- as.numeric(harmonized_data$beta.outcome)
harmonized_data$se.outcome <- as.numeric(harmonized_data$se.outcome)

harmonized_data <- harmonized_data[!is.na(harmonized_data$se.exposure) & !is.na(harmonized_data$se.outcome), ]
str(harmonized_data)

# Perform MR analysis with the harmonized data
mr_results <- mr(harmonized_data, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
# Record MR results
coprococcus_outcome <- mr_results


## Clostridium outcome ----

#Prepare exposure data
# Rename the columns to match the expected format for MR analysis
depression_clostridium_snp <- clostridium_from_depression_SNPs %>%
  rename(
    SNP = MarkerName,
    beta.exposure = LogOR,
    se.exposure = StdErrLogOR,
    effect_allele.exposure = A1,
    other_allele.exposure = A2,
    pval.exposure = P
  )

# Add the required columns for identification
depression_clostridium_snp$id.exposure <- "Depression"
depression_clostridium_snp$exposure <- "Depression"

# #Prepare outcome data
# Rename the columns to match the expected format

# Find the common SNPs between the two datasets
common_snps <- intersect(depression_clostridium_snp$SNP, clostridium_depression_snp$SNP)
# Filter both datasets to retain only common SNPs
depression_clostridium_snp_filtered <- depression_clostridium_snp[depression_clostridium_snp$SNP %in% common_snps, ]
clostridium_depression_snp_filtered <- clostridium_depression_snp[clostridium_depression_snp$SNP %in% common_snps, ]

if (!"eaf.exposure" %in% colnames(depression_clostridium_snp_filtered)) {
  depression_clostridium_snp_filtered$eaf.exposure <- NA
}
if (!"eaf.outcome" %in% colnames(clostridium_depression_snp_filtered)) {
  clostridium_depression_snp_filtered$eaf.outcome <- NA
}

harmonized_data <- harmonise_data(
  exposure_dat = depression_clostridium_snp_filtered,  # Depression as exposure
  outcome_dat =  clostridium_depression_snp_filtered   # Clostridium as outcome
)

# Ensure that se.exposure and se.outcome are numeric
harmonized_data$beta.exposure <- as.numeric(harmonized_data$beta.exposure)
harmonized_data$se.exposure <- as.numeric(harmonized_data$se.exposure)
harmonized_data$beta.outcome <- as.numeric(harmonized_data$beta.outcome)
harmonized_data$se.outcome <- as.numeric(harmonized_data$se.outcome)

harmonized_data <- harmonized_data[!is.na(harmonized_data$se.exposure) & !is.na(harmonized_data$se.outcome), ]
str(harmonized_data)

# Perform MR analysis with the harmonized data
mr_results <- mr(harmonized_data, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
# Record MR results
clostridium_outcome <- mr_results



# _______________----

#MRS ----

effects_for_MRS <- read_csv("/Users/MRS_estimate_depression.csv") %>% t() %>% as.data.frame() %>% rownames_to_column("taxa") %>% rename(estimate = V1)
#data from ANCOM-BC2 and MaAsLin2 meta-analysis
#bacteriadata <- read.csv("/Users/multiple_tests_bacteriadata.csv") #2,025 
bacteriadata <- read.csv("/Users/differential_analysis_counts.csv") #7,506, 

bacteriadata[bacteriadata == 0] <- 1 #This is a step recommended by Author of ANCOM-BC2 method from the email 976 taxa
DS_barcodes <- read_csv("/Users/DS_barcodes.csv") #4,456
RS_barcodes <- read_csv("/Users/RS_barcodes.csv") #2,452
metadata_full <- read_csv("/Users/metadata.csv" ) %>% dplyr::select(-1) %>% select(barcode, age, qGender, reg, weight.status, weeks.since.testing.microbiome, BMI , qEatFiber , alpha_diversity ) %>% mutate(barcode = paste0("X",barcode)) #7,506 individuals,  238 questions

bacterianames <- colnames(effects_for_MRS)
bacterianames <- gsub("...",".",bacterianames,fixed = TRUE)
bacterianames <- gsub("..",".",bacterianames, fixed = TRUE) #swap .. to .
colnames(effects_for_MRS) <- bacterianames

bacterianames <- colnames(bacteriadata)
bacterianames <- gsub("...",".",bacterianames,fixed = TRUE)
bacterianames <- gsub("..",".",bacterianames, fixed = TRUE) #swap .. to .
colnames(bacteriadata) <- bacterianames


relevant_taxa <- c("Bacteria.Firmicutes.Clostridia.Clostridiales.Peptostreptococcaceae.Romboutsia",
                   "Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae_1.Clostridium_sensu_stricto_1",
                   "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Agathobacter",
                   "Bacteria.Lentisphaerae.Lentisphaeria.Victivallales.Victivallaceae.Victivallis",
                   "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Coprococcus_3",
                   "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Bacteroides._pectinophilus_group",
                   "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Flavonifractor",
                   "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Tyzzerella",
                   
                   #  "Bacteria.Actinobacteria.Actinobacteria.Bifidobacteriales.Bifidobacteriaceae.Bifidobacterium",  #Genus MaAsLin2 meta-aanalysis
                   #  "Bacteria.Actinobacteria.Coriobacteriia.Coriobacteriales.Eggerthellaceae.Eggerthella",
                   #  "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Clostridium._innocuum_group", #Genus MaAsLin2 meta-aanalysis
                   #  "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_UCG.001", #Genus MaAsLin2 meta-aanalysis
                   #  "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Ruminococcus_1", #Genus MaAsLin2 meta-aanalysis
                   #  "Bacteria.Actinobacteria.Coriobacteriia.Coriobacteriales.Atopobiaceae.Olsenella", #Genus MaAsLin2 meta-aanalysis
                   #  "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Eubacterium._ruminantium_group", #Genus MaAsLin2 meta-aanalysis
                   #   "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Other",
                   #   "Bacteria.Firmicutes.Other.Other.Other.Other",
                   #   "Bacteria.Actinobacteria.Coriobacteriia.Coriobacteriales.Other.Other",
                   #   "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Other.Other",
                   #   "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Other",
                   #   "Bacteria.Firmicutes.Clostridia.Clostridiales.Other.Other",
                   
                   "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Bacteroidaceae.Bacteroides")

# Filter bacteriadata to retain only the relevant taxa
filtered_bacteriadata <- bacteriadata %>% 
  select(c("barcode", relevant_taxa)) 

# Make sure the taxa names in effects_for_MRS match with bacteriadata column names
# Merge bacteriadata with effects_for_MRS on relevant taxa
bacteriadata_long <- gather(filtered_bacteriadata, key = "taxa", value = "abundance", -barcode)

# Join with the effect estimates
bacteriadata_long <- bacteriadata_long %>%
  left_join(effects_for_MRS, by = c("taxa")) 
bacteriadata_long$estimate <- as.numeric(bacteriadata_long$estimate)
bacteriadata_long$abundance <- as.numeric(bacteriadata_long$abundance)

# Calculate MRS for each individual by multiplying abundance by the effect size
bacteriadata_long <- bacteriadata_long %>%
  mutate(weighted_abundance = abundance * estimate)

# Summarize MRS for each individual by summing weighted abundances across all taxa
MRS <- bacteriadata_long %>%
  group_by(barcode) %>%
  summarise(Microbiome_Risk_Score_MDD = sum(weighted_abundance, na.rm = TRUE)) #7,506 MRS data 

write_csv(MRS, "/Users/MRS_depression.csv")#7,506 MRS data by barcode
#write_csv(MRS, "/Users/MRS_depression_multipletests.csv")#2,025 MRS data by barcode

## MRS distribution ----

metadata_MRS <- metadata %>% select(barcode, qMedicalHistoryDepression) #6,942
MRS <- MRS %>% left_join(metadata_MRS, by = c("barcode"))

NoDepressionData <- MRS %>% filter(is.na(MRS$qMedicalHistoryDepression)) %>% mutate(barcode = paste0("X",barcode)) %>% left_join(metadata_full, by = c("barcode")) %>% distinct() #564 individuals with no Depression data
MRS_Depression_Data <- MRS %>% filter(!is.na(MRS$qMedicalHistoryDepression)) %>% distinct() #6,942 individuals with Depression data
MRS_Depression_Data <- MRS_Depression_Data %>% mutate(barcode = paste0("X",barcode)) %>% 
  filter(!barcode == "650.243.191") %>%  #data error
  filter(!barcode == "278.739.698") #data error

ggplot(MRS_Depression_Data, aes(x = Microbiome_Risk_Score)) +
  geom_histogram(binwidth = 50, fill = "blue", color = "white") +
  labs(title = "Distribution of Microbiome Risk Score (MRS)", x = "MRS", y = "Frequency") +
  theme_minimal()

#START ANALYSIS  
DS <- MRS_Depression_Data %>% filter(barcode%in%DS_barcodes$barcode)%>% mutate(depression = ifelse(qMedicalHistoryDepression == 2, 1, 0)) %>% left_join(DS_barcodes, by = c("barcode")) #4,455
RS <- MRS_Depression_Data %>% filter(barcode%in%RS_barcodes$barcode)%>% mutate(depression = ifelse(qMedicalHistoryDepression == 2, 1, 0)) %>% left_join(RS_barcodes, by = c("barcode")) #2,451
MRS_Depression_Data <- MRS_Depression_Data %>% left_join(metadata_full, by = c("barcode")) %>% mutate(depression = ifelse(qMedicalHistoryDepression == 2, 1, 0)) #6,940

DS$Set <- "Discovery Set"
RS$Set <- "Replication Set"
combined_data <- bind_rows(DS, RS)

# Combine the datasets

ggplot(combined_data, aes(x = Microbiome_Risk_Score, fill = factor(depression))) +
  geom_histogram(binwidth = 0.5, position = "dodge", alpha = 0.7) +
  facet_wrap(~ Set, scales = "free_y") +
  scale_fill_manual(values = c("#1f78b4", "red"), 
                    name = "Depression Status", 
                    labels = c("No", "Yes")) +
  labs(title = "Distribution of Microbiome Risk Score (MRS) in Discovery and Replication Sets",
       x = "MRS",
       y = "Frequency") +
  theme_minimal()

# Fit logistic regression model using the discovery set
logit_model <- glm(depression ~ Microbiome_Risk_Score, data = DS, family = "binomial")
# Predict depression probabilities in the discovery set
DS$predicted_prob <- predict(logit_model, type = "response")
# Predict depression probabilities in the replication set
RS$predicted_prob <- predict(logit_model, newdata = RS, type = "response")

## ROC for Discovery and Replication Set ----
roc_discovery <- roc(DS$depression, DS$predicted_prob)
auc_discovery <- auc(roc_discovery)

# ROC for Replication Set
roc_replication <- roc(RS$depression, RS$predicted_prob)
auc_replication <- auc(roc_replication)

# Plot ROC curves
plot(roc_discovery, col = "blue", main = "ROC Curves for Discovery and Replication Sets")
lines(roc_replication, col = "red")
legend("bottomright", legend = c(paste("Discovery AUC =", round(auc_discovery, 2)),
                                 paste("Replication AUC =", round(auc_replication, 2))),
       col = c("blue", "red"), lwd = 2)

##Calculate AUC ----

# Calculate AUC for Discovery Set
roc_discovery <- roc(DS$depression, DS$predicted_prob)
auc_discovery <- auc(roc_discovery)
print(paste("AUC for Discovery Set:", round(auc_discovery, 2))) #[1] "AUC for Replication Set: 0.55"

# Calculate AUC for Replication Set
roc_replication <- roc(RS$depression, RS$predicted_prob)
auc_replication <- auc(roc_replication)
print(paste("AUC for Replication Set:", round(auc_replication, 2))) #[1] "AUC for Replication Set: 0.57"

# Calculate AUC for No Depression Data

NoDepressionData$predicted_prob <- predict(logit_model, newdata = NoDepressionData, type = "response")
NoDepressionData$dummy_depression <- sample(c(0, 1), nrow(NoDepressionData), replace = TRUE)
roc_no_records <- roc(NoDepressionData$dummy_depression, NoDepressionData$predicted_prob)
auc_no_records <- auc(roc_no_records)
print(paste("AUC for No Medical Records:", round(auc_no_records, 2))) #[1] "AUC for No Medical Records: 0.5"

#write_csv(DS, "/Users/Bacteria/DS_ROC_ALL_TAXA.csv")
#write_csv(RS, "/Users/Bacteria/RS_ROC_ALL_TAXA.csv")
#write_csv(NoDepressionData, "/Users/Bacteria/NoDepressionData_ROC_ALL_TAXA.csv")

DS <- read_csv("/Users/Bacteria/DS_ROC_ALL_TAXA.csv") %>% distinct() #4,456, 544 depressed
RS <- read_csv("/Users/Bacteria/RS_ROC_ALL_TAXA.csv") %>% distinct() #2,452, 285 depressed
NoDepressionData <- read_csv("/Users/Bacteria/NoDepressionData_ROC_ALL_TAXA.csv") %>% distinct() #564

DS %>% count(qMedicalHistoryDepression)

# Calculate ROC and AUC for Discovery Set
roc_discovery <- roc(DS$depression, DS$predicted_prob)
auc_discovery <- auc(roc_discovery)

# Calculate ROC and AUC for Replication Set
roc_replication <- roc(RS$depression, RS$predicted_prob)
auc_replication <- auc(roc_replication)

# Calculate ROC and AUC for No Depression Data
roc_no_records <- roc(NoDepressionData$dummy_depression, NoDepressionData$predicted_prob)
auc_no_records <- auc(roc_no_records)

# Plot ROC curves with pastel colors and add AUC values
plot(roc_discovery, col = "#66c2a5", lwd = 2, main = "ROC Curves for Discovery, Replication Sets, and No Records")
lines(roc_replication, col = "#e78ac3", lwd = 2)
lines(roc_no_records, col = "#1f78b4", lwd = 2)

# Add AUC values as text on the plot
legend(x = 0.6, y = 0.15,
       legend = c(paste("Discovery Set AUC =", round(auc_discovery, 2)),
                  paste("Replication Set AUC =", round(auc_replication, 2)),
                  paste("No Medical Records AUC =", round(auc_no_records, 2))),
       col = c("#66c2a5", "#e78ac3", "#1f78b4"),
       lwd = 2)

# Customize the axis labels and title
title(xlab = "Specificity", ylab = "Sensitivity")

## Adjustment for Gender + age ----

# Fit logistic regression model adjusted for age and gender using the discovery set
logit_model_DS_epid <- glm(depression ~ age + BMI , data = DS, family = "binomial")
logit_model_adjusted_DS <- glm(depression ~ Microbiome_Risk_Score + age + BMI, data = DS, family = "binomial")
logit_model_adjusted_NO <- glm(dummy_depression ~ Microbiome_Risk_Score, data = NoDepressionData, family = "binomial")


# Predict depression probabilities for the discovery set with epid factors
DS$predicted_prob_epid <- predict(logit_model_DS_epid, type = "response")

# Predict depression probabilities for the discovery set with MRS + epid factors
DS$predicted_prob_adjusted <- predict(logit_model_adjusted_DS, type = "response")

# Predict depression probabilities for the replication set
RS$predicted_prob_adjusted <- predict(logit_model_adjusted_DS, newdata = RS, type = "response")

# Predict depression probabilities for the replication set with epid factors
RS$predicted_prob_epid <- predict(logit_model_DS_epid, newdata = RS, type = "response")

# Predict probabilities for the individuals without depression data
NoDepressionData$predicted_prob_adjusted <- predict(logit_model_adjusted_NO, newdata = NoDepressionData, type = "response")



# Calculate ROC curves for the adjusted model
roc_discovery_epid <- roc(DS$depression, DS$predicted_prob_epid)
auc_discovery_epid <- auc(roc_discovery_epid)

roc_discovery_adjusted <- roc(DS$depression, DS$predicted_prob_adjusted)
auc_discovery_adjusted <- auc(roc_discovery_adjusted)

roc_replication_adjusted <- roc(RS$depression, RS$predicted_prob_adjusted)
auc_replication_adjusted <- auc(roc_replication_adjusted)

roc_replication_epid <- roc(RS$depression, RS$predicted_prob_epid)
auc_replication_epid <- auc(roc_replication_epid)

roc_no_records_adjusted <- roc(NoDepressionData$dummy_depression, NoDepressionData$predicted_prob_adjusted)
auc_no_records_adjusted <- auc(roc_no_records_adjusted)


# Plot ROC curves with pastel colors and add AUC values
plot(
  #  roc_discovery, col = "#66c2a5", lwd = 2, main = "ROC curves for discovery set")
  #lines(roc_discovery_adjusted, col = "#1f78b4", lwd = 2)
  #lines(roc_discovery_epid, col = "lightblue", lwd = 2)
  roc_replication, col = "#e78ac3", lwd = 2, main = "ROC curves for trained model performance on replication set")
lines(roc_replication_adjusted, col = "orange", lwd = 2)
lines(roc_replication_epid, col = "lightpink", lwd = 2)
lines(roc_no_records_adjusted, col = "grey", lwd = 2)

# Add AUC values as text on the plot
legend(x = 0.7, y = 0.2,
       legend = c(
         #         paste("Prediction with MRS + BMI+age, discovery AUC =", round(auc_discovery_adjusted, 2)),
         #                  paste("Prediction with MRS, discovery AUC =", round(auc_discovery, 2)),
         #                  paste("Prediction with BMI+age, discovery AUC =", round(auc_discovery_epid, 2)),
         paste("Prediction with MRS + BMI+age, replication AUC =", round(auc_replication_adjusted, 2)),
         paste("Prediction with MRS, replication AUC =", round(auc_replication, 2)),
         paste("Prediction with BMI+age, replication AUC =", round(auc_replication_epid, 2)),
         paste("No Medical Records AUC =", round(auc_no_records_adjusted, 2))),
       #       col = c("#66c2a5","#1f78b4","lightblue", "#e78ac3", "orange", "lightpink","grey" ),
       #       col = c("#1f78b4","#66c2a5","lightblue","grey"),
       col = c( "orange", "#e78ac3","lightpink","grey" ),
       
       lwd = 2)

# Customize the axis labels and title
title(xlab = "Specificity", ylab = "Sensitivity")


## Split by gender ----

DS_males <- DS %>% filter(qGender == 1)
DS_females <- DS %>% filter(qGender == 0)

RS_males <- RS %>% filter(qGender == 1)
RS_females <- RS %>% filter(qGender == 0)

NoDepressionData_males <- NoDepressionData %>% filter(qGender == 1)
NoDepressionData_females <- NoDepressionData %>% filter(qGender == 0)

# Logistic regression for males
logit_model_males <- glm(depression ~ Microbiome_Risk_Score + age, data = DS_males, family = "binomial")
DS_males$predicted_prob_males <- predict(logit_model_males, type = "response")
RS_males$predicted_prob_males <- predict(logit_model_males, newdata = RS_males, type = "response")
NoDepressionData_males$predicted_prob_males <- predict(logit_model_males, newdata = NoDepressionData_males, type = "response")

# Logistic regression for females
logit_model_females <- glm(depression ~ Microbiome_Risk_Score + age, data = DS_females, family = "binomial")
DS_females$predicted_prob_females <- predict(logit_model_females, type = "response")
RS_females$predicted_prob_females <- predict(logit_model_females, newdata = RS_females, type = "response")
NoDepressionData_females$predicted_prob_females <- predict(logit_model_females, newdata = NoDepressionData_females, type = "response")

# ROC curves for males
roc_discovery_males <- roc(DS_males$depression, DS_males$predicted_prob_males)
auc_discovery_males <- auc(roc_discovery_males)

roc_replication_males <- roc(RS_males$depression, RS_males$predicted_prob_males)
auc_replication_males <- auc(roc_replication_males)

roc_no_records_males <- roc(NoDepressionData_males$dummy_depression, NoDepressionData_males$predicted_prob_males)
auc_no_records_males <- auc(roc_no_records_males)

# ROC curves for females
roc_discovery_females <- roc(DS_females$depression, DS_females$predicted_prob_females)
auc_discovery_females <- auc(roc_discovery_females)

roc_replication_females <- roc(RS_females$depression, RS_females$predicted_prob_females)
auc_replication_females <- auc(roc_replication_females)

roc_no_records_females <- roc(NoDepressionData_females$dummy_depression, NoDepressionData_females$predicted_prob_females)
auc_no_records_females <- auc(roc_no_records_females)

# Plot ROC curves for males
par(mfrow = c(1, 2)) # Set up side-by-side plots
plot(roc_discovery_males, col = "#66c2a5", lwd = 2, main = "ROC Curves for Males")
lines(roc_replication_males, col = "#e78ac3", lwd = 2)
lines(roc_no_records_males, col = "#1f78b4", lwd = 2)
legend(x = 0.6, y = 0.15,
       legend = c(paste("Discovery Set AUC =", round(auc_discovery_males, 2)),
                  paste("Replication Set AUC =", round(auc_replication_males, 2)),
                  paste("No Medical Records AUC =", round(auc_no_records_males, 2))),
       col = c("#66c2a5", "#e78ac3", "#1f78b4"),
       lwd = 2)

# Plot ROC curves for females
plot(roc_discovery_females, col = "#66c2a5", lwd = 2, main = "ROC Curves for Females")
lines(roc_replication_females, col = "#e78ac3", lwd = 2)
lines(roc_no_records_females, col = "#1f78b4", lwd = 2)
legend(x = 0.6, y = 0.15,
       legend = c(paste("Discovery Set AUC =", round(auc_discovery_females, 2)),
                  paste("Replication Set AUC =", round(auc_replication_females, 2)),
                  paste("No Medical Records AUC =", round(auc_no_records_females, 2))),
       col = c("#66c2a5", "#e78ac3", "#1f78b4"),
       lwd = 2)

## Boxplot by gender ----

MRS_males <- MRS_Depression_Data %>% filter(qGender == "1") %>% distinct() #2,936
MRS_females <- MRS_Depression_Data %>% filter(qGender == "0")%>% distinct() #4,004
# Calculate the number of cases for each group
counts_males <- MRS_males %>%
  group_by(qMedicalHistoryDepression) %>%
  summarise(count = n())

counts_females <- MRS_females %>%
  group_by(qMedicalHistoryDepression) %>%
  summarise(count = n())

# Create boxplots for MRS distribution for males with p-value and number of cases
boxplot_males <- ggplot(MRS_males, aes(x = as.factor(qMedicalHistoryDepression), y = Microbiome_Risk_Score, fill = as.factor(qMedicalHistoryDepression))) +
  geom_boxplot() +
  scale_fill_manual(values = c("lightblue", "lightpink"), name = "Depression Status", labels = c("No", "Yes")) +
  labs(title = "Microbiome Risk Score for Males", x = "Depression Status", y = "Microbiome Risk Score") +
  theme_minimal() +
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     label.y = max(MRS_males$Microbiome_Risk_Score, na.rm = TRUE) - 1,
                     hjust = -0.55) + # Adjust p-value position with hjust
  geom_text(data = counts_males, aes(x = as.factor(qMedicalHistoryDepression), y = max(MRS_males$Microbiome_Risk_Score, na.rm = TRUE) + 0.1, 
                                     label = paste("N =", count)), 
            color = "black", vjust = -0.5)

# Create boxplots for MRS distribution for females with p-value and number of cases
boxplot_females <- ggplot(MRS_females, aes(x = as.factor(qMedicalHistoryDepression), y = Microbiome_Risk_Score, fill = as.factor(qMedicalHistoryDepression))) +
  geom_boxplot() +
  scale_fill_manual(values = c("lightblue", "lightpink"), name = "Depression Status", labels = c("No", "Yes")) +
  labs(title = "Microbiome Risk Score for Females", x = "Depression Status", y = "Microbiome Risk Score") +
  theme_minimal() +
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     label.y = max(MRS_females$Microbiome_Risk_Score, na.rm = TRUE) - 1,
                     hjust = -0.55) + # Adjust p-value position with hjust
  geom_text(data = counts_females, aes(x = as.factor(qMedicalHistoryDepression), y = max(MRS_females$Microbiome_Risk_Score, na.rm = TRUE) + 0.1, 
                                       label = paste("N =", count)), 
            color = "black", vjust = -0.5)


# Arrange the plots side by side

grid.arrange(boxplot_males, boxplot_females, ncol = 2)

# Density plot

ggplot(MRS_Depression_Data, aes(x = Microbiome_Risk_Score, fill = as.factor(qMedicalHistoryDepression))) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("lightblue", "lightpink"), name = "Depression Status", labels = c("No", "Yes")) +
  labs(title = "Density Plot of Microbiome Risk Score by Depression Status", x = "Microbiome Risk Score", y = "Density") +
  theme_minimal()

## Risk score nomination ----

# Normalize MRS values to a range of 0 to 1
MRS_normalisation <- MRS_Depression_Data %>%
  mutate(normalized_MRS = (Microbiome_Risk_Score - min(Microbiome_Risk_Score, na.rm = TRUE)) / 
           (max(Microbiome_Risk_Score, na.rm = TRUE) - min(Microbiome_Risk_Score, na.rm = TRUE)))

# Scale the normalized MRS values to a range of 1 to 10
MRS_normalisation <- MRS_normalisation %>%
  mutate(scaled_MRS = 1 + (normalized_MRS * 9))

# Round the scaled MRS values to one decimal place for better readability
MRS_normalisation <- MRS_normalisation %>%
  mutate(scaled_MRS = round(scaled_MRS, 1))

# Round the MRS values to the nearest integer
MRS_normalisation <- MRS_normalisation %>%
  mutate(rounded_MRS = round(scaled_MRS))

# Ensure the rounded scores are within the 0 to 10 range
MRS_normalisation <- MRS_normalisation %>%
  mutate(rounded_MRS = ifelse(rounded_MRS <= 0, 0, ifelse(rounded_MRS > 10, 10, rounded_MRS)))

# Count the number of individuals for each rounded MRS score
MRS_counts <- MRS_normalisation %>%
  group_by(rounded_MRS) %>%
  summarise(count = n())

# Plotting the distribution of the scaled MRS

ggplot(MRS_normalisation, aes(x = rounded_MRS, fill = as.factor(qMedicalHistoryDepression))) +
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
  scale_fill_manual(values = c("blue", "red"), name = "Depression Status", labels = c("No", "Yes")) +
  labs(title = "Distribution of Scaled Microbiome Risk Score (1-10)", x = "Scaled Microbiome Risk Score", y = "Count") +
  theme_minimal()


# Group the data by Rounded_MRS and calculate the range
MRS_ranges <- MRS_normalisation %>%
  group_by(rounded_MRS) %>%
  summarise(
    Min_Microbiome_Risk_Score = min(Microbiome_Risk_Score),
    Max_Microbiome_Risk_Score = max(Microbiome_Risk_Score)
  )

print(MRS_ranges)

## Depression probability ----

MRS_Depression_Data_prob <- MRS_normalisation

# Plotting the distribution of depression probability by scaled MRS

logit_model <- glm(depression ~ scaled_MRS, data = MRS_Depression_Data_prob, family = "binomial")

# Predict depression probability using the logistic model
MRS_Depression_Data_prob$predicted_prob <- predict(logit_model, type = "response")

# Calculate the percentage of individuals with depression for each rounded_MRS

depression_prob_data <- MRS_Depression_Data_prob %>%
  group_by(rounded_MRS) %>%
  summarise(depression_percentage = mean(depression) * 100) %>% filter(rounded_MRS < 10)

# Calculate the maximum count in the histogram
max_count <- max(table(MRS_Depression_Data_prob$rounded_MRS))

# Create the plot
ggplot() +
  # Histogram layer using MRS_Depression_Data_prob
  geom_histogram(data = MRS_Depression_Data_prob, aes(x = rounded_MRS, fill = as.factor(qMedicalHistoryDepression)),
                 binwidth = 1, position = "dodge", alpha = 0.7) +
  scale_fill_manual(values = c("#1f78b4", "#e78ac3"), name = "Depression Status", labels = c("No", "Yes")) +
  labs(title = "Distribution of Scaled Microbiome Risk Score (1-10) with Depression Probability",
       x = "Scaled Microbiome Risk Score",
       y = "Count") +
  theme_minimal() +
  # Smooth line using depression_prob_data and scaling to the max count
  geom_smooth(data = depression_prob_data, aes(x = rounded_MRS, y = depression_percentage * max_count / 100),
              method = "loess", color = "red", size = 1.5, se = FALSE) +
  # Secondary axis on the right with appropriate scaling
  scale_y_continuous(sec.axis = sec_axis(~ . / max_count * 100, name = "Depression Probability (%)"))

# _______________----


# REGRESSION Microbiome-METADATA ----

variables <- c( "qSmokingStatus",
                "qPAworkVigorousPresent", #PHYSICAL ACTIVITY
                "qPArecreationalVigorousPresent",
                "qPAtravelBicycleUse",
                "qMedicalHistoryDepression",
                "qPAworkModeratePresent",
                "qPAyoga",
                "qPArecreationalModeratePresent", #MEDICAL
                "qDrugHistoryProbiotics",
                "qDrugHistoryPrebiotics",
                "qDrugHistoryProbioticsType_other",
                "qDrugHistoryProbioticsType_bifidobacterium",
                "qDrugHistoryProbioticsType_lactobacillus",
                "qMedicalHistoryPolycysticOvarySyndrome",
                "qDrugHistoryAntibioticsLast3M",
                "qMedicalHistoryCoeliacDisease",
                "qMedicalHistoryMigraine",
                "qMedicalHistoryDiabetesType2",
                "qMedicalHistoryCrohnsDisease",
                "qMedicalHistoryBloating",
                "qMedicalHistorySmallIntestinalBacterialOvergrowth",
                "qMedicalHistoryAsthma",
                "qMedicalHistoryPsoriasis",
                "qMedicalHistoryGastroesophagealRefluxDisease",
                "qMedicalHistoryAtopicDermatitis",
                "qMedicalHistoryHiatusHernia",
                "qMedicalHistoryAppendectomy",
                "qMedicalHistoryHighBloodPressure",
                "qMedicalHistoryLactoseIntolerance",
                "qMedicalHistoryVaricoseVeins",
                "qMedicalHistoryHypothyroidism",
                "qMedicalHistoryDepression",
                "qMedicalHistoryUrolithiasis",
                "qMedicalHistoryGout",
                "qMedicalHistoryOtherAllergicDiseases",
                "qMedicalHistoryGallstones",
                "qEatMeatRegularly" ,#FOOD
                "qEatFoodRichInIron",
                "qEatNonOilyFishRegularly",
                "qAlcoholStatus",
                "qMoreSalt",
                "qEatFoodRichInVitD",
                "qEatNuts", 
                "qEatOliveOil",
                "qEatFoodRichInOmega3",
                "qEatWholeGrains",
                "qCoffeeFrequency",
                "qEatFoodRichInVitC",
                "qEatFoodRichInTransFat",
                "qEatFiber",
                "qTotalFluidIntake",
                "IBD",
                "qAlcoholUnitsPerWeek",
                "qEatBackedSweets",
                "qDrinkSweetBeverages",
                #"qNightSleepDuration",
                "qMedicalHistoryUlcerativeColitis",
                "BMI",
                "alpha_diversity",
                "cohort",
                "short_sleep",
                "long_sleep",
                "normal_sleep",
                "Obesity.status",
                "Underweight.status",
                "Healthyweight.status",
                "Overweight.status",
                "SCFA.percentage",
                "NoFrVeg",
                "FrVeg_less_5_a_day",
                "FrVeg_5_a_day_and_more",
                "FrVeg_10_a_day_and_more",
                "qCalciumVitDIntake",
                "daily_fruits",
                "daily_veg",
                "daily_veg_cooked",
                "qPregnancyStatus",
                "qMenopauseStatus",
                "qOralContraceptivesTaking",
                "Diarrhea",
                "Constipation",
                "Normal_stool"
)

# Romboutsia ----

regress_output <- list()
for ( i in 1:length(variables)){
  print(paste("Fitting ", i, variables[i])) 
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  mod <- as.formula(sprintf("Romboutsia ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = metadata ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  #let's add Odds Ratio and confidence interval to our loop and visualize it in a summary spreadsheet fit.summary
  #OR.CI <- exp(cbind(OR = coef(fit), confint(fit))) # calculate OR and confidence interval
  #we57rix bind two spreadsheets together into nice summary (cbind)
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL
write.csv(allregression,"/Users/Bacteria/allregression2024_Romboutsia.csv")

#### MULTIPLE TESTING -----
#Let's prepare spreadsheet for multiple testing Shannon diversity index
#remove Intersepr rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_Romboutsia <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)
#p-value correction 
allregression_Romboutsia$Bonferroni =
  p.adjust(allregression_Romboutsia$p.value,
           method = "bonferroni")
allregression_Romboutsia$BH_Romboutsia =
  p.adjust(allregression_Romboutsia$p.value,
           method = "BH")
#Add description
Description <- read.csv("/Users/Metadata_names.csv") %>% select(term = Name, DisplayName)
allregression_Romboutsia <- left_join(allregression_Romboutsia, Description, by = "term") %>% select(DisplayName, everything())
write.csv(allregression_Romboutsia,"/Users/Bacteria/allregression2024_multipletesting_Romboutsia.csv")

# Clostridium.sensu_stricto ----

regress_output <- list()
for ( i in 1:length(variables)){
  print(paste("Fitting ", i, variables[i])) 
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  mod <- as.formula(sprintf("Clostridium.sensu_stricto ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = metadata ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  #let's add Odds Ratio and confidence interval to our loop and visualize it in a summary spreadsheet fit.summary
  #OR.CI <- exp(cbind(OR = coef(fit), confint(fit))) # calculate OR and confidence interval
  #we57rix bind two spreadsheets together into nice summary (cbind)
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL
write.csv(allregression,"/Users/Bacteria/allregression2024_Clostridium.sensu_stricto.csv")

#### MULTIPLE TESTING -----
#Let's prepare spreadsheet for multiple testing Shannon diversity index
#remove Intersepr rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_Clostridium.sensu_stricto <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)
#p-value correction 
allregression_Clostridium.sensu_stricto$Bonferroni =
  p.adjust(allregression_Clostridium.sensu_stricto$p.value,
           method = "bonferroni")
allregression_Clostridium.sensu_stricto$BH_Clostridium.sensu_stricto =
  p.adjust(allregression_Clostridium.sensu_stricto$p.value,
           method = "BH")
#Add description
Description <- read.csv("/Users/Metadata_names.csv") %>% select(term = Name, DisplayName)
allregression_Clostridium.sensu_stricto <- left_join(allregression_Clostridium.sensu_stricto, Description, by = "term") %>% select(DisplayName, everything())
write.csv(allregression_Clostridium.sensu_stricto,"/Users/Bacteria/allregression2024_multipletesting_Clostridium.sensu_stricto.csv")
allregression_Clostridium.sensu_stricto_list <- allregression_Clostridium.sensu_stricto %>% filter(BH_Clostridium.sensu_stricto<= 0.051) %>% select(term, estimate)


# Agathobacter ----

regress_output <- list()
for ( i in 1:length(variables)){
  print(paste("Fitting ", i, variables[i])) 
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  mod <- as.formula(sprintf("Agathobacter ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = metadata ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  #let's add Odds Ratio and confidence interval to our loop and visualize it in a summary spreadsheet fit.summary
  #OR.CI <- exp(cbind(OR = coef(fit), confint(fit))) # calculate OR and confidence interval
  #we57rix bind two spreadsheets together into nice summary (cbind)
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL
write.csv(allregression,"/Users/Bacteria/allregression2024_Agathobacter.csv")

#### MULTIPLE TESTING -----
#Let's prepare spreadsheet for multiple testing Shannon diversity index
#remove Intersepr rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_Agathobacter <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)
#p-value correction 
allregression_Agathobacter$Bonferroni =
  p.adjust(allregression_Agathobacter$p.value,
           method = "bonferroni")
allregression_Agathobacter$BH_Agathobacter =
  p.adjust(allregression_Agathobacter$p.value,
           method = "BH")
#Add description
Description <- read.csv("/Users/Metadata_names.csv") %>% select(term = Name, DisplayName)
allregression_Agathobacter <- left_join(allregression_Agathobacter, Description, by = "term") %>% select(DisplayName, everything())
write.csv(allregression_Agathobacter,"/Users/Bacteria/allregression2024_multipletesting_Agathobacter.csv")
allregression_Agathobacter_list <- allregression_Agathobacter %>% filter(BH_Agathobacter<= 0.051) %>% select(term, estimate)


# Victivallis ----

regress_output <- list()
for ( i in 1:length(variables)){
  print(paste("Fitting ", i, variables[i])) 
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  mod <- as.formula(sprintf("Victivallis ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = metadata ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  #let's add Odds Ratio and confidence interval to our loop and visualize it in a summary spreadsheet fit.summary
  #OR.CI <- exp(cbind(OR = coef(fit), confint(fit))) # calculate OR and confidence interval
  #we57rix bind two spreadsheets together into nice summary (cbind)
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL
write.csv(allregression,"/Users/Bacteria/allregression2024_Victivallis.csv")

#### MULTIPLE TESTING -----
#Let's prepare spreadsheet for multiple testing Shannon diversity index
#remove Intersepr rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_Victivallis <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)
#p-value correction 
allregression_Victivallis$Bonferroni =
  p.adjust(allregression_Victivallis$p.value,
           method = "bonferroni")
allregression_Victivallis$BH_Victivallis =
  p.adjust(allregression_Victivallis$p.value,
           method = "BH")
#Add description
Description <- read.csv("/Users/Metadata_names.csv") %>% select(term = Name, DisplayName)
allregression_Victivallis <- left_join(allregression_Victivallis, Description, by = "term") %>% select(DisplayName, everything())
write.csv(allregression_Victivallis,"/Users/Bacteria/allregression2024_multipletesting_Victivallis.csv")
allregression_Victivallis_list <- allregression_Victivallis %>% filter(BH_Victivallis<= 0.051) %>% select(term, estimate)


# Flavonifractor ----

regress_output <- list()
for ( i in 1:length(variables)){
  print(paste("Fitting ", i, variables[i])) 
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  mod <- as.formula(sprintf("Flavonifractor ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = metadata ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  #let's add Odds Ratio and confidence interval to our loop and visualize it in a summary spreadsheet fit.summary
  #OR.CI <- exp(cbind(OR = coef(fit), confint(fit))) # calculate OR and confidence interval
  #we57rix bind two spreadsheets together into nice summary (cbind)
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL
write.csv(allregression,"/Users/Bacteria/allregression2024_Flavonifractor.csv")

#### MULTIPLE TESTING -----
#Let's prepare spreadsheet for multiple testing Shannon diversity index
#remove Intersepr rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_Flavonifractor <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)
#p-value correction 
allregression_Flavonifractor$Bonferroni =
  p.adjust(allregression_Flavonifractor$p.value,
           method = "bonferroni")
allregression_Flavonifractor$BH_Flavonifractor =
  p.adjust(allregression_Flavonifractor$p.value,
           method = "BH")
#Add description
Description <- read.csv("/Users/Metadata_names.csv") %>% select(term = Name, DisplayName)
allregression_Flavonifractor <- left_join(allregression_Flavonifractor, Description, by = "term") %>% select(DisplayName, everything())
write.csv(allregression_Flavonifractor,"/Users/Bacteria/allregression2024_multipletesting_Flavonifractor.csv")

allregression_Flavonifractor_list <- allregression_Flavonifractor %>% filter(BH_Flavonifractor<= 0.051) %>% select(term, estimate)

# Bacteroides ----

regress_output <- list()
for ( i in 1:length(variables)){
  print(paste("Fitting ", i, variables[i])) 
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  mod <- as.formula(sprintf("Bacteroides ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = metadata ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  #let's add Odds Ratio and confidence interval to our loop and visualize it in a summary spreadsheet fit.summary
  #OR.CI <- exp(cbind(OR = coef(fit), confint(fit))) # calculate OR and confidence interval
  #we57rix bind two spreadsheets together into nice summary (cbind)
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL
write.csv(allregression,"/Users/Bacteria/allregression2024_Bacteroides.csv")

#### MULTIPLE TESTING -----
#Let's prepare spreadsheet for multiple testing Shannon diversity index
#remove Intersepr rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_Bacteroides <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)
#p-value correction 
allregression_Bacteroides$Bonferroni =
  p.adjust(allregression_Bacteroides$p.value,
           method = "bonferroni")
allregression_Bacteroides$BH_Bacteroides =
  p.adjust(allregression_Bacteroides$p.value,
           method = "BH")
#Add description
Description <- read.csv("/Users/Metadata_names.csv") %>% select(term = Name, DisplayName)
allregression_Bacteroides <- left_join(allregression_Bacteroides, Description, by = "term") %>% select(DisplayName, everything())
write.csv(allregression_Bacteroides,"/Users/Bacteria/allregression2024_multipletesting_Bacteroides.csv")
allregression_Bacteroides_list <- allregression_Bacteroides %>% filter(BH_Bacteroides<= 0.051) %>% select(term, estimate)

# Coprococcus ----

regress_output <- list()
for ( i in 1:length(variables)){
  print(paste("Fitting ", i, variables[i])) 
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  mod <- as.formula(sprintf("Coprococcus ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = metadata ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  #let's add Odds Ratio and confidence interval to our loop and visualize it in a summary spreadsheet fit.summary
  #OR.CI <- exp(cbind(OR = coef(fit), confint(fit))) # calculate OR and confidence interval
  #we57rix bind two spreadsheets together into nice summary (cbind)
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL
write.csv(allregression,"/Users/Bacteria/allregression2024_Coprococcus.csv")

#### MULTIPLE TESTING -----
#Let's prepare spreadsheet for multiple testing Shannon diversity index
#remove Intersepr rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_Coprococcus <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)
#p-value correction 
allregression_Coprococcus$Bonferroni =
  p.adjust(allregression_Coprococcus$p.value,
           method = "bonferroni")
allregression_Coprococcus$BH_Coprococcus =
  p.adjust(allregression_Coprococcus$p.value,
           method = "BH")
#Add description
Description <- read.csv("/Users/Metadata_names.csv") %>% select(term = Name, DisplayName)
allregression_Coprococcus <- left_join(allregression_Coprococcus, Description, by = "term") %>% select(DisplayName, everything())
write.csv(allregression_Coprococcus,"/Users/Bacteria/allregression2024_multipletesting_Coprococcus.csv")

# Tyzzerella ----

regress_output <- list()
for ( i in 1:length(variables)){
  print(paste("Fitting ", i, variables[i])) 
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  mod <- as.formula(sprintf("Tyzzerella ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = metadata ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  #let's add Odds Ratio and confidence interval to our loop and visualize it in a summary spreadsheet fit.summary
  #OR.CI <- exp(cbind(OR = coef(fit), confint(fit))) # calculate OR and confidence interval
  #we57rix bind two spreadsheets together into nice summary (cbind)
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL
write.csv(allregression,"/Users/Bacteria/allregression2024_Tyzzerella.csv")

#### MULTIPLE TESTING -----
#Let's prepare spreadsheet for multiple testing Shannon diversity index
#remove Intersepr rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_Tyzzerella <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)
#p-value correction 
allregression_Tyzzerella$Bonferroni =
  p.adjust(allregression_Tyzzerella$p.value,
           method = "bonferroni")
allregression_Tyzzerella$BH_Tyzzerella =
  p.adjust(allregression_Tyzzerella$p.value,
           method = "BH")
#Add description
Description <- read.csv("/Users/Metadata_names.csv") %>% select(term = Name, DisplayName)
allregression_Tyzzerella <- left_join(allregression_Tyzzerella, Description, by = "term") %>% select(DisplayName, everything())
write.csv(allregression_Tyzzerella,"/Users/Bacteria/allregression2024_multipletesting_Tyzzerella.csv")

# Bacteroides.pectinophilus_group ----

regress_output <- list()
for ( i in 1:length(variables)){
  print(paste("Fitting ", i, variables[i])) 
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  mod <- as.formula(sprintf("Bacteroides.pectinophilus_group ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = metadata ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  #let's add Odds Ratio and confidence interval to our loop and visualize it in a summary spreadsheet fit.summary
  #OR.CI <- exp(cbind(OR = coef(fit), confint(fit))) # calculate OR and confidence interval
  #we57rix bind two spreadsheets together into nice summary (cbind)
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL
write.csv(allregression,"/Users/Bacteria/allregression2024_Bacteroides.pectinophilus_group.csv")

#### MULTIPLE TESTING -----
#Let's prepare spreadsheet for multiple testing Shannon diversity index
#remove Intersepr rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_Bacteroides.pectinophilus_group <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)
#p-value correction 
allregression_Bacteroides.pectinophilus_group$Bonferroni =
  p.adjust(allregression_Bacteroides.pectinophilus_group$p.value,
           method = "bonferroni")
allregression_Bacteroides.pectinophilus_group$BH_Bacteroides.pectinophilus_group =
  p.adjust(allregression_Bacteroides.pectinophilus_group$p.value,
           method = "BH")
#Add description
Description <- read.csv("/Users/Metadata_names.csv") %>% select(term = Name, DisplayName)
allregression_Bacteroides.pectinophilus_group <- left_join(allregression_Bacteroides.pectinophilus_group, Description, by = "term") %>% select(DisplayName, everything())
write.csv(allregression_Bacteroides.pectinophilus_group,"/Users/Bacteria/allregression2024_multipletesting_Bacteroides.pectinophilus_group.csv")

#Bifidobacterium ----

regress_output <- list()
for ( i in 1:length(variables)){
  print(paste("Fitting ", i, variables[i])) 
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  mod <- as.formula(sprintf("Bifidobacterium ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = metadata ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  #let's add Odds Ratio and confidence interval to our loop and visualize it in a summary spreadsheet fit.summary
  #OR.CI <- exp(cbind(OR = coef(fit), confint(fit))) # calculate OR and confidence interval
  #we57rix bind two spreadsheets together into nice summary (cbind)
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL
write.csv(allregression,"/Users/Bacteria/allregression2024_Bifidobacterium.csv")

#### MULTIPLE TESTING -----
#Let's prepare spreadsheet for multiple testing Shannon diversity index
#remove Intersepr rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_Bifidobacterium <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)
#p-value correction 
allregression_Bifidobacterium$Bonferroni =
  p.adjust(allregression_Bifidobacterium$p.value,
           method = "bonferroni")
allregression_Bifidobacterium$BH_Bifidobacterium =
  p.adjust(allregression_Bifidobacterium$p.value,
           method = "BH")
#Add description
Description <- read.csv("/Users/Metadata_names.csv") %>% dplyr::select(term = Name, DisplayName)
allregression_Bifidobacterium <- left_join(allregression_Bifidobacterium, Description, by = "term") %>% dplyr::select(DisplayName, everything())
write.csv(allregression_Bifidobacterium,"/Users/Bacteria/allregression2024_multipletesting_Bifidobacterium.csv")


#Lactobacillus ----


regress_output <- list()
for ( i in 1:length(variables)){
  print(paste("Fitting ", i, variables[i])) 
  #creating formula for each variable in a loop so we could get a list of associations and don't check association with every factor
  mod <- as.formula(sprintf("Lactobacillus ~ %s + qGender + age + reg + weeks.since.testing.microbiome + cohort", variables[i])) 
  #fitting regression model for each variable. Using formula created each variable
  fit <- glm(formula = mod, family=gaussian(link="identity"), data = metadata ) #regression
  #let's add previously created regression  summary using broom tidy to save data as table into the spreadsheet
  fit.summary.tidy <- broom::tidy(fit) #spreadsheet view
  #let's add Odds Ratio and confidence interval to our loop and visualize it in a summary spreadsheet fit.summary
  #OR.CI <- exp(cbind(OR = coef(fit), confint(fit))) # calculate OR and confidence interval
  #we57rix bind two spreadsheets together into nice summary (cbind)
  regress_output[[variables[i]]] <- cbind(fit.summary.tidy) # add it all together in one spreadsheet 
}
allregression <- bind_rows(regress_output)
rownames(allregression) <- NULL
write.csv(allregression,"/Users/Bacteria/allregression2024_Lactobacillus.csv")

#### MULTIPLE TESTING -----
#Let's prepare spreadsheet for multiple testing Shannon diversity index
#remove Intersepr rows and sorting p-value from small to large(as recommended for multiple testing https://rcompanion.org/rcompanion/f_01.html)
allregression_Lactobacillus <- allregression %>% filter(!term %in% c("(Intercept)", "regru", "age", "qGender", "weeks.since.testing.microbiome", "cohort")) %>% arrange(p.value)
#p-value correction 
allregression_Lactobacillus$Bonferroni =
  p.adjust(allregression_Lactobacillus$p.value,
           method = "bonferroni")
allregression_Lactobacillus$BH_Lactobacillus =
  p.adjust(allregression_Lactobacillus$p.value,
           method = "BH")
#Add description
Description <- read.csv("/Users/Metadata_names.csv") %>% dplyr::select(term = Name, DisplayName)
allregression_Lactobacillus <- left_join(allregression_Lactobacillus, Description, by = "term") %>% dplyr::select(DisplayName, everything())
allregression_Lactobacillus_bifidobacteriun <- left_join(allregression_Lactobacillus, allregression_Bifidobacterium, by = "term")
write.csv(allregression_Lactobacillus,"/Users/Bacteria/allregression2024_multipletesting_Lactobacillus.csv")
write.csv(allregression_Lactobacillus_bifidobacteriun,"/Users/Bacteria/Lactobacillus_bifidobacteriun.csv")

