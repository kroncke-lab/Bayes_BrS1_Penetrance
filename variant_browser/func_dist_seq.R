# B Kroncke
# 2019.11.15
# functions to calculate feature densities from 
# structures.

# weight_function == "sigmoid" or "sin"
funcdist <- function(resnum, var, deval, distances, function_name, weight_function, a){
  #resnum<-toString(resnum)
  sum_weights <- 0
  sum_function <- 0
  if (is.na(match(resnum, distances$V1))) {
    seq_dists<-seq(as.integer(resnum)-30, as.integer(resnum)+30,1)
    seq_dists<-seq_dists[seq_dists<1160 & seq_dists>0]
    #3.8*(NumberOfResiduesAway)^0.5 (30 residues before and 30 after) (ref. Folding & Design Vol 1 No 5 Dobson, 1996)
    for(dists in seq_dists){
      est_dist<-3.8*(abs(as.integer(resnum)-dists))^0.5
      for (flist in deval[deval$resnum==dists & deval$var != var, function_name]){
        if (!is.na(flist) & weight_function == "sigmoid"){
          sum_weights <- sum_weights + sigmoid(as.numeric(est_dist))
          sum_function <- sum_function + flist*sigmoid(as.numeric(est_dist))
        } else if (!is.na(flist) & weight_function == "sin") {
          sum_weights <- sum_weights + sin_dist(as.numeric(est_dist),a)
          sum_function <- sum_function + flist*sin_dist(as.numeric(est_dist),a)
        }
      }
    }
    return(c(sum_function/sum_weights, sum_weights))
  }
  neighbors<-distances[distances$V1==resnum, "V2"]
  dists<-distances[distances$V1==resnum, "V3"]
  i = 0
  for (neig in neighbors) {
    i = i + 1
    for (flist in deval[deval$resnum==neig & deval$var != var, function_name]){
      if (!is.na(flist) & weight_function == "sigmoid"){
        sum_weights <- sum_weights + sigmoid(as.numeric(dists[i]))
        sum_function <- sum_function + flist*sigmoid(as.numeric(dists[i]))
      } else if (!is.na(flist) & weight_function == "sin") {
        sum_weights <- sum_weights + sin_dist(as.numeric(dists[i]),a)
        sum_function <- sum_function + flist*sin_dist(as.numeric(dists[i]),a)
      }
    }
  }
  return(c(sum_function/sum_weights, sum_weights))
}

sigmoid <- function(x){
  return(1/(1+exp(x/2)))
}

sin_dist <- function(x,a){
  if (x<=a-3.14){
    return(1)
  } else if (a-3.14<x & x<=3.14*2+a-3.14) {
    return(-sin((x-a)/2)/2+0.5)
  } else {
    return(0)
  }
}

LASSO_CV <- function(d, k, features, EP_function){
  feat_dist <- paste(EP_function, "_dist", sep = "")
  feat_dist_weight <- paste(EP_function,"_dist_weight", sep = "")
  features <- c(features, EP_function, 
                feat_dist, feat_dist_weight)
  d2 <- d[, which(names(d) %in% c("resnum", "mutAA", features))]
  d2 <- d2[complete.cases(d2),]
  kfold=k
  d3<-createFolds(d2$pamscore, k=kfold)
  r2<-matrix()
  j=0
  
  for(i in d3){
    j = j+1
    d2[, feat_dist]<-NA
    d2[, feat_dist_weight]<-NA
    d4<-d2[-i,]
    for (rec in 1:nrow(d2)){
      d2[rec, c(feat_dist, feat_dist_weight)] <- funcdist(d2[rec, "resnum"], d2[rec, "mutAA"], d4, distances, EP_function, "sigmoid", 7)
    }
    d4[,feat_dist] <- d2[-i,feat_dist]
    d4[,feat_dist_weight] <- d2[-i,feat_dist_weight]
    d2 <- d2[complete.cases(d2),]
    d4 <- d4[complete.cases(d4),]
    d4 <- d4[, which(names(d4) %in% c(features))]
    d5 <- d2[, which(names(d2) %in% c(features))]
    cv.lass <- cv.glmnet( as.matrix(d4[,-which(names(d4)==EP_function)]), d4[, EP_function]) 
    fin.lass <- glmnet(as.matrix(d4[,-which(names(d4)==EP_function)]), d4[, EP_function], lambda = cv.lass$lambda.min)  
    print(names(d4[, c(fin.lass$beta@Dimnames[[1]][1+fin.lass$beta@i])]))
    tmpfit <- lm(d4[,EP_function] ~ ., data = d4[, c(fin.lass$beta@Dimnames[[1]][1+fin.lass$beta@i])])
    r2[j] = cor(x = predict(tmpfit, rows.in.a1.not.in.a2(d5,d4)), rows.in.a1.not.in.a2(d5,d4)[EP_function], method = "pearson")
    print(r2[j])
  }
  print(r2)
  return(print(paste("mean adj R2: ", 1-(1-mean(r2)^2)*(nrow(d5)-1)/(nrow(d5)-length(fin.lass$beta@Dimnames[[1]][1+fin.lass$beta@i])-1), " 95% CI: ", (mean(r2)-2*std.error(r2))^2, ", ", (mean(r2) + 2*std.error(r2))^2)))
}


a1 <- data.frame(a = 1:5, b=letters[1:5])
a2 <- data.frame(a = 1:3, b=letters[1:3])

rows.in.a1.not.in.a2  <- function(a1,a2){
  a1.vec <- apply(a1, 1, paste, collapse = "")
  a2.vec <- apply(a2, 1, paste, collapse = "")
  a1.without.a2.rows <- a1[!a1.vec %in% a2.vec,]
  return(a1.without.a2.rows)
}


