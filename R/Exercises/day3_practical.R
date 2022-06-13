#### Course 22123: Computational precision medicine
#### Day 3 practical: clustering and feature selection for molecular subtyping
#### By Lars Ronn Olsen
#### 07/06/2022
#### Technical University of Denmark

## Load necessary packages
library("tidyverse")
library("e1071")
library("caTools")
library("class")


## Load expression data (NOTE: this is microarray data, not RNA-seq!)
# CIT_full contains the expression of all probes for all samples
# CIT_annot contains the annotation for the 375 probes used for classification
# CIT_class contains the subtypes
load("data/CIT_data.Rdata")


## Check out the distribution of the expression of a probe or two (or all)
df <- data.frame(t(CIT_full))

df %>% ggplot(aes(x = X1007_s_at)) +
  geom_density()


# QUESTION: Does it look different from RNA-seq? How? Why? Do you anticipate any problems comparing these data with RNA-seq data?
# IDK lav dag 2 og læs paperet fuuuuuck


## Subset to core probes (in CIT_annot) and do a leave-one-out classification of CIT training samples using distance to centroid (DTC)
CIT_sub <- CIT_full[rownames(CIT_full) %in% CIT_annot$Probe.Set.ID,]

predicted <- c()

for (i in 1:ncol(CIT_sub)) {
  training <- CIT_sub[,-i]
  training_class <- CIT_classes[-i]
  
  test <- CIT_sub[,i]
  test_class <- CIT_classes[i]
  
  distances <- c()
  
  for (class in unique(CIT_classes)) {
    curr_centroid <- rowMeans(training[, training_class == class])
    
    distances[[class]] <- dist(rbind(curr_centroid, test))
    
  }
  predicted[i] <- names(distances)[which.min(distances)]
  
}

hit_rate <- mean(predicted == CIT_classes)
hit_rate


# QUESTION: Does DTC classification work?
# Yes it works quite well. The hitrate I get is ~0.997, which is quite good. 


## Do a leave-one-out classification of CIT training samples using kNN
knn.cross <- tune.knn(x = t(CIT_sub), 
                      y = as.factor(CIT_classes), 
                      k = 1:20,
                      tunecontrol=tune.control(sampling = "cross"), 
                      cross=1)
summary(knn.cross)

predictions <- c()


for (i in 1:ncol(CIT_sub)) {
  training <- t(CIT_sub[,-i])
  training_class <- CIT_classes[-i]
  
  test <- t(CIT_sub[,i])
  test_class <- CIT_classes[i]
  
  prediction <- knn(train = training,
                    test = test,
                    cl = training_class,
                    k = 5)
  predictions[i] <- (prediction == test_class)
  
}
knn_hitrate <- mean(predictions)
knn_hitrate


# QUESTION: Does kNN classification work?
# Yes but not as good as DTC as it has a hitrate of ~0.972.


## Load test data set and try to classify these samples. Calculate some metrics that will help you get an idea of how well you did
load("data/test_data.Rdata")

## Load RNA-seq data for breast cancer from Xena browser ("TCGA Breast Cancer (BRCA)" -> "IlluminaHiSeq (n=1,218) TCGA Hub") + phenotype data ("Phenotypes")


## See if any of the following features are confounded with PAM50 subtype ("PAM50Call_RNAseq" in the phenotype matrix)

# Age_at_Initial_Pathologic_Diagnosis_nature2012
# AJCC_Stage_nature2012
# ER_Status_nature2012
# PR_Status_nature2012
# HER2_Final_Status_nature2012
# radiation_therapy



