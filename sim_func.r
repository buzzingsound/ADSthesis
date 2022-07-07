### Data generation, Amputation and Imputation ###

simData <- function(parms) {
  ## Generate predictor covariance matrix:
  p        <- length(parms$varNames) - 1
  sX       <- matrix(parms$covX, p, p)
  diag(sX) <- 1.0
  
  ## Simulate predictor data:
  X <- rmvnorm(n = parms$nObs, mean = rep(0, p), sigma = sX)
  
  ## Extract coefficients:
  alpha  <- parms$coefs[1, ]
  beta   <- matrix(parms$coefs[-1, ])
  
  ## Define residual variance:
  signal <- t(beta) %*% cov(X) %*% beta
  sY     <- (signal / parms$r2) - signal
  
  ## Simulate outcome data:
  y <- alpha + X %*% beta + rnorm(parms$nObs, 0, sqrt(sY))
  
  outData           <- data.frame(y, X)
  colnames(outData) <- parms$varNames
  
  outData
}


imposeMissing <- function(data, parms) {
  data_comp <- data
  for(v in 1 : length(parms$incompVars)){
    if(parms$missingtype == "pMNAR"){
      ## Generate nonresponse vector:
      rVec <- with(parms,
                   simMissingness(pm    = pm,
                                  data  = data,
                                  preds = incompVars[v],
                                  type  = missType[v])
      )
    }
    else{
      ## Generate nonresponse vector:
      rVec <- with(parms,
                   simMissingness(pm    = pm,
                                  data  = data_comp,
                                  preds = pred,
                                  type  = missType[v])
      )
    }
    data[rVec, parms$incompVars[v]] <- NA
  }
  data
}

#data <- derivation

fitModels <- function(data, parms) {
  if(is.data.frame(data)) { # We're analyzing a single dataset
    lmFit <- try(
      lm(parms$model, data = data, na.action = "na.omit"),
      silent = TRUE
    )
  }
  else { # We're analyzing multiply imputed data
    fitList <- try(
      lapply(data,
             FUN = function(miData, parms) lm(parms$model, data = miData),
             parms = parms),
      silent = TRUE
    )
    
    lmFit <- fitList
  }
  lmFit
}


runImp <- function(parms, missData) {
  
  ## Impute the missing data:
  miceOut <- with(parms,
                  mice(data      = missData,
                       m         = nImps,
                       maxit     = miceIters,
                       method    = impmethod,
                       printFlag = verbose)
  )
  #saveRDS(miceOut, "miceOut.rds") 
  
  complete(miceOut, "all")
  #impLists <- complete(miceOut, "all") 
  
}



.dvName <- function(x) all.vars(x$terms)[1]

.lambda <- function(x) (x + 1) / (x + 3)

.miDf <- function(m, b, t, dfCom) {
  fmi   <- .fmi(m, b, t)
  df0   <- (m - 1) * (1 / fmi^2)
  dfObs <- .lambda(dfCom) * dfCom * (1 - fmi)
  
  df0 / (1 + (df0 / dfObs))
}

.fmi <- function(m, b, t) (1 + (1 / m)) * (b / t)

## Generate predictions from MI data:
predictMi <- function(fits, newData, interval = NA, pooled = TRUE) {
  ## Generate imputation-specific predictions:
  predStuff <- list()
  for(m in 1 : length(fits))
    predStuff[[m]] <-
      predict(fits[[m]], newdata = newData, se.fit = TRUE)
  
  ## Extract components from prediction objects:
  preds  <- lapply(predStuff, "[[", x = "fit")
  seMat  <- sapply(predStuff, "[[", x = "se.fit")
  scales <- sapply(predStuff, "[[", x = "residual.scale")
  
  ## Return early when imputation-specific predictions are requested:
  if(!pooled)
    return(preds)
  else
    predMat <- do.call(cbind, preds)
  
  ## Compute pooled predictions:
  preds <- rowMeans(predMat)
  
  ## Compute within-imputation variance:
  if(is.na(interval) | interval == "confidence")
    w <- rowMeans(seMat^2)
  else if(interval == "prediction")
    w <- colMeans(
      apply(seMat^2, 1, function(x, y) x + y, y = scales^2)
    )
  else
    stop(paste0("'",
                interval,
                "' is not a valid argument for 'interval'. Should be either 'confidence' or 'prediction'."
    )
    )
  
  ## Compute between-imputation variance:
  b <- apply(predMat, 1, var)
  
  ## Compute total variance:
  t <- w + b + (b / length(fits))
  
  ## Compute df:
  df <- .miDf(length(fits), b, t, predStuff[[1]]$df)
  
  ## Compute the FMIs:
  out <- cbind(est = preds, fmi = .fmi(length(fits), b, t))
  
  ## Compute CIs and aggregate output:
  if(!is.na(interval)) {
    moe <- qt(0.975, df) * sqrt(t)
    out <- cbind(out, lwr = preds - moe, upr = preds + moe)
  }
  
  out
}


## Compute Mean Squared Errors from MI data:
mseMi <- function(fits, newData) {
  ## Generate imputation-specific predictions:
  preds <- predictMi(fits = fits, newData = newData, pooled = FALSE)
  
  ## Compute imputation-specific MSEs:
  mse <- c()
  for(m in 1 : length(fits))
    mse[m] <- MSE(y_pred = preds[[m]],
                  y_true = newData[[m]][[.dvName(fits[[1]])]]
    )
  
  ## Return the aggregate MSE:
  mean(mse)
}


getStats <- function(lmOut, valData, parms, mi = FALSE) {
  if(class(lmOut) != "try-error") {
    if(mi){
      PIpred <- predictMi(lmOut, valData, interval = "prediction")[, c('est', 'lwr','upr')]
      CIpred <- predictMi(lmOut, valData, interval = "confidence")[, c('lwr','upr')]
      
      outvec <- list()
      outvec$y <- valData$y ### KML: This is the wrong Y value. Don't average the imputed data. 
      ### BVV: changed to separate val set for everything to evaluate onc
      outvec$CI <- CIpred
      outvec$PI <- PIpred
      outvec$pm <- parms$pm
      
    }
    else{
      PIpred <- predict(lmOut, newdata = valData, interval = "prediction")
      CIpred <- predict(lmOut, newdata = valData, interval = "confidence")[, c('lwr','upr')]
      
      outvec <- list()
      outvec$y <-  valData$y
      outvec$CI <- CIpred
      outvec$PI <- PIpred
      outvec$pm <- parms$pm
    }
  }
  else
    outvec <- list("LM_CONVERGENCE_FAILURE", lmOut)
  
  outvec
}# END getStats()


runCell <- function(rp, compData, missData, valData, impLists, parms) {
  
  ## Create a condition tag to label output objects:
  tag1 <- with(parms,
               paste0( "_rs", 100 * r2)
  )
  
  tag2 <- with(parms,
               paste0(tag1, "_missingtype_", missingtype)
  )
  
  tag3 <- with(parms,
               paste0(tag2, "_impmethod_", impmethod)
  )
  
  ## Save the current parameter set:
  if(rp == 1)
    saveRDS(parms,
            file = paste0(parms$outDir,
                          "parms",
                          tag3,
                          ".rds")
    )
  
  ### Save the results ###
  
  ## Fit complete data model:
  compFit <- fitModels(compData, parms) 
  compOut <- getStats(compFit, valData, parms) 
  
  saveRDS(compOut,
          file = paste0(parms$outDir,
                        "compOut",
                        tag1,
                        "_rep",
                        rp,
                        ".rds")
  )
  
  ## Fit listwise deleted models:
  ldFit <- fitModels(missData, parms) ### KML
  ldOut <- getStats(ldFit, valData, parms) ### KML
  
  saveRDS(ldOut,
          file = paste0(parms$outDir,
                        "ldOut",
                        tag2,
                        "_rep",
                        rp,
                        ".rds")
  )
  
  ## Fit SI models
  siFit <- fitModels(impLists[[1]], parms)
  siOut <- getStats(siFit, valData, parms)
  
  saveRDS(siOut,
          file = paste0(parms$outDir,
                        "siOut",
                        tag3,
                        "_rep",
                        rp,
                        ".rds")
  )
  
  ## Fit MI models
  miFit <- fitModels(impLists, parms)
  miOut <- getStats(miFit, valData, parms, mi = TRUE)
  
  saveRDS(miOut,
          file = paste0(parms$outDir,
                        "miOut",
                        tag3,
                        "_rep",
                        rp,
                        ".rds")
  )
  
}# END runCell()

doRep <- function(rp, conds, parms) {
  ## Setup the PRNG:
  .lec.SetPackageSeed(rep(parms$seed, 6))
  if(!rp %in% .lec.GetStreams())
    .lec.CreateStream(c(1 : parms$nStreams))
  .lec.CurrentStream(rp)
  
  ## Sample a value for proportion of missing data:
  parms$pm <- runif(1, min = 0.1, max = 0.3)
  
  ## Loop over conditions:
  for(i in 1 : nrow(conds)) {
    ## Save the current r2
    r2 <- parms$r2
    
    ## Update the values of imputation method and r2:
    parms$r2   <- conds[i, "r2"]
    parms$impmethod <- conds[i, "impmethod"]
    
    ## Simulate new complete and validation data, if r2 has changed:
    ## Add separate testing set
    check <-
      (is.null(r2)) || (r2 != parms$r2)
    if(check){
      compData <- simData(parms)
      
      valparms <- parms
      valparms$nObs <- parms$nVal
      
      val <- simData(valparms)
    }
    
    ## Define the missing data-related design parameters:
    parms$missingtype <- conds[i, "missingtype"]
    
    ## Generate missing data:
    if(parms$missingtype == "MAR"){
      parms$pred <- parms$MARpred
      missData <- imposeMissing(compData, parms)
    }
    else if(parms$missingtype == "MNAR" || parms$missingtype == "pMNAR"){
      parms$pred <- parms$MNARpred
      missData <- imposeMissing(compData, parms)
    }
    else{
      missData <- compData
      
      ## Response vector for MCAR missingness
      tmp1 <- with(parms,
                   rbinom(n = nrow(missData) * length(incompVars),
                          size = 1,
                          prob = pm)
      )
      
      ## Response vector for observed variables
      tmp2 <- rep(FALSE, prod(dim(missData)) - length(tmp1))
      
      ## Response matrix (assumes column-major wrapping)
      rMat <- matrix(c(as.logical(tmp1), tmp2), nrow = nrow(missData))
      
      ## Punch holes in the data:
      missData[rMat] <- NA
    }
    
    ## Impute missdatasets:
    impLists <-
      runImp(parms = parms, missData = missData)
  
    
    ## Run the computations for the current condition:
    
    runCell(rp       = rp,
            compData = compData,
            missData = missData,
            valData  = val,
            impList  = impLists,
            parms    = parms)
  }
  
  rp # return rep index
}


applyLib <- function(pkgList)
  lapply(pkgList, library, character.only = TRUE, logical = TRUE)