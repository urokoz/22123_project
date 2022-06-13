# subset to CIT core probes (first column in CIT_annot)
CIT_sub <- CIT_full[rownames(CIT_full) %in% CIT_annot$Probe.Set.ID,]
# make empty prediction vector
pred_vector <- c()
# loop over samples (columns)
for(i in 1:ncol(CIT_sub)) {
  # make a training matrix consisting of all samples but the one you leave out
  training <- CIT_sub[,-i]
  # remove that samples class from the class vector
  training_classes <- CIT_classes[-i]
  # make a test vector consisting of only that sample
  test <- CIT_sub[,i]
  # get the class of that sample
  test_class <- CIT_classes[i]
  
  # make an empty centroid matrix
  centroids <- NULL
  # loop over each of the six classes
  for (class in unique(CIT_classes)) {
    # for each of these classes, subset the training matrix to samples belonging to that class, and calculate the mean expression of each gene in the class
    class_centroid <- rowMeans(training[,training_classes==class])
    # add the mean vector to the centroids matrix
    centroids <- cbind(centroids, class_centroid)
  }
  # add colnames to the centroid matrix
  colnames(centroids) <- unique(CIT_classes)
  # calculate the distance of the test sample to the centroids
  d <- as.matrix(dist(t(cbind(centroids, test))))
  # assign the class of the closest centroid
  class_pred <- names(which.min(d[1:6,7]))
  # check if you got it right and make a logical vector
  pred_vector <- c(pred_vector, test_class==class_pred)
}
# check how many of the cases you got it right
table(pred_vector)