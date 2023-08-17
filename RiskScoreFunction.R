###########################
##  Risk Score Function  ##
###########################

library(caret)
library(glmnet)
library(RcppEigen)


RiskScoreFun <- function(DF,Ncov=5,pval.thres.vec,FDR=0.2,miss=FALSE,iseed,level=0.95){
  
  Nsam <- nrow(DF)
  
  set.seed(iseed) #edit for real data (can use anything)
  train <- sample(1:Nsam,floor(Nsam/2),replace=FALSE)
  test <- setdiff(1:Nsam,train)
 
  trainDF <- DF[train,]
  testDF <- DF[test,]
  
trap <- names(trainDF)[1]

zval <- abs(qnorm((1-level)/2))

#####################################
# elastic net (for imputed data only)
#####################################

if(miss==FALSE){
  print("Beginning elastic net")
  f <- as.formula(paste(trap,"~."))
  model <- train(
    f, data = trainDF, method = "glmnet",
    trControl = trainControl("cv", number = 10),
    tuneLength = 50
  )
  # get coefficients
  coefs <- coef(model$finalModel,model$bestTune$lambda)
  coefs.mets <- coefs[(Ncov+2):length(coefs)] #Ncov+2 accounts for intercept to capture only metabolite coefficients

  # how many total metabolites were captured
  coefs.nonzero <- which(coefs.mets!=0)
  captured.e <- length(coefs.nonzero)

  metnames.e <- colnames(testDF)[coefs.nonzero+(Ncov+1)]
  
  if(captured.e==0){print("No nonzero coefficients found with elastic net")}

  if(captured.e>0){
    # calculate risk scores
    test.mets <- testDF[,(Ncov+2):ncol(testDF)] #remove covariates and exposure
    risk.scores.enet <- as.matrix(test.mets) %*% coefs.mets
    
    # predictive model R square
    f2 <- as.formula(paste(trap,"~ risk.scores.enet +",paste(colnames(testDF)[2:(Ncov+1)],
                                                             collapse = "+")))
    enetmod <- summary(fastLm(f2,data=testDF))
    coefmat.e <- enetmod$coefficients[,c(1,2,4)]
    lowerCI <- coefmat.e[,1] - zval*(coefmat.e[,2])
    upperCI <- coefmat.e[,1] + zval*(coefmat.e[,2])
    Coefs.e <- cbind(coefmat.e,lowerCI,upperCI)

    rsquare.e <- enetmod$adj.r.squared

    Enetlist <- list(captured=captured.e, RiskScores=risk.scores.enet, Coefficients=Coefs.e,
                     Rsquare=rsquare.e, MetNames=metnames.e, RSCoefs=coefs.mets)
  }
  print("Elastic net complete")
} #end elastic net


  ##############
  # thresholding
  ##############
  
  print("Beginning thresholding")
  
  pvals <- apply(trainDF[,(Ncov+2):ncol(trainDF)],2,
                 function(x) summary(fastLm(x ~ as.matrix(trainDF[,1:(Ncov+1)])))$coefficients[,4])
  
  qvals <- p.adjust(pvals[2,],method="BH") #Benjamini-Hochberg correction
  
  coefs <- apply(trainDF[,(Ncov+2):ncol(trainDF)],2,
                 function(x) summary(fastLm(x ~ as.matrix(trainDF[,1:(Ncov+1)])))$coefficients[,1])
  
  
  
  p.list <- list() #set up output
  f4 <- as.formula(paste(trap," ~ risk.scores.p +",paste(colnames(testDF)[2:(Ncov+1)],collapse = "+")))
  
  # replace NA with 0 for matrix multiplication in loop below (equivalent of na.rm=TRUE)
  testDF0 <- testDF
  testDF0[is.na(testDF0)] <- 0
  
  
  for(i in 1:(length(pval.thres.vec)+1)){ #loop through p-value thresholds
    tryCatch(
      {
    if(i<(length(pval.thres.vec)+1)){
      p1 <- which(pvals[2,]<pval.thres.vec[i])+(Ncov+1)
    }
    if(i==(length(pval.thres.vec)+1)){
      p1 <- which(qvals<FDR)+(Ncov+1)
    }
        
    captured.p <- length(p1)
    metnames <- colnames(testDF0)[p1]
    if(captured.p == 0){
      print('No significant metabolites')
      break
    }
    
    RSCoefs <- coefs[2,p1-(Ncov+1)]
    risk.scores.p <- as.matrix(testDF0[,p1]) %*% coefs[2,p1-(Ncov+1)]

    pmod <- summary(fastLm(f4,data=testDF0))
    coefmat.p <- pmod$coefficients[,c(1,2,4)]
    LowerCI <- coefmat.p[,1] - zval*coefmat.p[,2]
    UpperCI <- coefmat.p[,1] + zval*coefmat.p[,2]
    Coefs.p <- cbind(coefmat.p,LowerCI,UpperCI)
    rsquare.p <- pmod$adj.r.squared

    out.p <- list(captured=captured.p,
                  RiskScores=risk.scores.p,
                  Coefficients=Coefs.p,
                  Rsquare=rsquare.p, 
                  MetNames=metnames,
                  RSCoefs=RSCoefs)
    
    if(i<(length(pval.thres.vec)+1)){
      p.list[[as.character(pval.thres.vec[i])]] <- out.p
    }
    if(i==(length(pval.thres.vec)+1)){
      p.list[["BH"]] <- out.p
    }
    
        }, error = function(e){})
    
  } # end loop through p-value thresholds
  

  if(miss==TRUE){
    p.list[["test.DF"]] <- testDF0
    return(p.list)
  }
  
  if(miss==FALSE){
    p.list[["test.DF"]] <- testDF0
    p.list[["Enet"]] <- Enetlist
    return(p.list)
  }
  
}
  
  
    
