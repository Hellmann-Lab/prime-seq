
# Slide Theme -------------------------------------------------------------

theme_slide <- function (base_size = 11, base_family = "") {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(plot.background = element_rect(colour = "black", fill=NA, size=1.25),
          panel.grid.major = element_line(colour = "grey75"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey25", size = 1),
          strip.background = element_rect(fill = NA, colour = NA),
          strip.text = element_text(colour = "black", size = rel(1.25), face = 'bold'),
          plot.title = element_text(size = rel(1.25), face = 'bold', colour = 'black',
                                    margin = margin(t = 5, r = 0,  b = 15, l = 0)),
          axis.text = element_text(colour = "black", size = rel(1.25)),
          axis.title = element_text(colour = "black", face="bold", size = rel(1.25)),
          axis.title.y = element_text(margin = margin(t = 0, r = 15,  b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15,  r = 0, b = 0, l = 0)),
          legend.title = element_text(colour = "black", size = rel(1), face = "bold"),
          legend.key.size = unit(0.75, "lines"),
          legend.text = element_text(size = rel(0.9), colour = "black"),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_rect(colour = NA, fill = NA))
}


# Ambient Profiles --------------------------------------------------------

# estimate ambient profile
get_ambient_mle<-function(countData, AmbProp){
  
  sample_number  <- ncol(countData)
  feature_number <- nrow(countData)
  
  if( length(AmbProp) == 1){
    message("assume equal absolute amounts of contaminating molecules.")
    ssize <- rep(sum(countData)*AmbProp/sample_number, sample_number)
  }else if( length(AmbProp) == ncol(countData)){
    message("use provided different ambient proportions.")
    ssize <- sapply(1:sample_number, function(i){ sum(countData[,i])*AmbProp[i] })
  }else{ 
    stop(message("AmbProp not equal to number of samples."))
  }
  est <- sapply(1:sample_number, 
                function(i){ ssize[i]*rowSums(countData)/sum(countData)} )
  rownames(est) <- rownames(countData)
  colnames(est) <- colnames(countData)
  return(est)
}

# get_ambient_mle_IH <- function(y, amb=0.05){
#   
#   sample_number  <- ncol(y)
#   feature_number <- nrow(y)
#   
#   if( length(amb)==1){
#     print("Assume equal absolute amounts of contaminating molecules.")
#     ssize <- rep(sum(y)*amb/sample_number, sample_number )
#   }else if( length(amb)== sample_number ){
#     print("use provided different ambient proportions.")
#     ssize<-sapply(1:sample_number,function(i){ sum(y[,i])*amb[i] })
#   }else{ 
#     print("ERROR. amb not equal number of samples.")
#   }
#   sapply(1:sample_number , function(i){ ssize*rowSums(y)/sum(y)} )
# }

# sample ambient profile
get_ambient_sample<-function(countData, AmbProp){
  
  sample_number  <- ncol(countData)
  feature_number <- nrow(countData)
  mtx<-matrix(0, nrow = feature_number, ncol = sample_number)
  if( length(AmbProp) == 1){
    message("assume equal absolute amounts of contaminating molecules.")
    ssize <- rep(sum(countData) * AmbProp /sample_number, sample_number)
  }else if( length(AmbProp) == ncol(countData)){
    message("use provided different ambient proportions.")
    ssize <- sapply(1:sample_number, 
                    function(i){ sum(countData[,i]) * AmbProp[i] })
  }else{ 
    stop(message("AmbProp not equal to number of samples."))
  }

  for(i in 1:sample_number ) {
    xx<-table(sample(1:feature_number,
                     replace = T,
                     size= ssize[i], prob=rowSums(countData)/sum(countData))) 
    mtx[as.numeric(names(xx)),i]<-xx
  }
  rownames(mtx) <- rownames(countData)
  colnames(mtx) <- colnames(countData)
  return(mtx)
}

# get_ambient_sample_IH<-function(y, amb=0.05){
#   
#   sample_number  <- ncol(y)
#   feature_number <- nrow(y)
#   mtx<-matrix(0, nrow = feature_number, ncol = sample_number)
#   
#   if( length(amb)==1){
#     print("Assume equal absolute amounts of contaminating molecules.")
#     ssize <- rep(sum(y)*amb/sample_number, sample_number )
#   }else if( length(amb)== sample_number){
#     print("use provided different ambient proportions.")
#     ssize<-sapply(1:sample_number,function(i){ sum(y[,i])*amb[i] })
#   }else{ 
#     print("ERROR. amb not equal number of samples.")
#   }
#   
#   for(i in 1:sample_number ) {
#     xx<-sample(1:feature_number,replace = T, size= ssize[i], prob=rowSums(y)/sum(y)) %>% table 
#     mtx[as.numeric(names(xx)),i]<-xx
#   }
#   
#   return(mtx)
# }


run_simulate <- function(nsamples, ngenes, mu = NULL, disp = NULL, sf = NULL, EstRes = NULL){
  # MEAN
  if(is.null(mu)){
    if(!is.null(EstRes)){
      est.mu <- EstRes$Estimates$MoM$mean
      est.mu <- na.omit(ifelse(is.infinite(est.mu), NA, est.mu))
      est.mu <- est.mu[est.mu>0]
      if(ngenes > length(est.mu)){
        muvec <- sample(x = est.mu, size = ngenes, replace = T)
      }
      if(ngenes <= length(est.mu)){
        muvec <- sample(x = est.mu, size = ngenes, replace = F)
      }
    } else {
      stop(message('You have not defined a mean expression vector nor provided estimates. Aborting.'))
    }
  } else {
    if(length(mu) == ngenes){
      muvec <- mu
    } else {
      stop(message('The vector containing mean expression values is not the same length as the number of genes requested for simulation. Aborting.'))
    } 
  }
  
  lmu <- log2(muvec)
  
  # Size (1/DISPERSION)
  if(is.null(disp)){
    if(!is.null(EstRes)){
      fit.obj <- EstRes$MeanvsDisp$MoM
      preddisp.mean = suppressWarnings(approx(fit.obj$x,
                                              fit.obj$y,
                                              xout = lmu, rule=2)$y)
      preddisp.sd = suppressWarnings(approx(fit.obj$x,
                                            fit.obj$sd,
                                            xout = lmu, rule=2)$y)
      ldisp = stats::rnorm(n = length(lmu),
                           mean = preddisp.mean,
                           sd = preddisp.sd)
      dispvec <- 2^ldisp
      sizevec <- 1/dispvec
    } else {
      stop(message('You have not defined a dispersion vector nor provided estimates. Aborting.'))
    }
  } else {
    if(length(disp) == ngenes){
      dispvec <- disp
      sizevec <- 1/disp
    } else {
      stop(message('The vector containing dispersion values is not the same length as the number of genes requested for simulation. Aborting.'))
    } 
  }
  
  # library size
  if(is.null(sf)){
    if(!is.null(EstRes)){
      est.sf <- EstRes$Estimates$SF
      if(nsamples > length(est.sf)){
        sfvec <- sample(x = est.sf, size = nsamples, replace = T)
      }
      if(nsamples <= length(est.sf)){
        sfvec <- sample(x = est.sf, size = nsamples, replace = F)
      }
    } else {
      stop(message('You have not defined a library size factor vector nor provided estimates. Aborting.'))
    }
  } else {
    if(length(sf) == nsamples){
      sfvec <- sf
    } else {
      stop(message('The vector containing library size factor values is not the same length as the number of samples requested for simulation. Aborting.'))
    } 
  }
  
  mumat <- outer(muvec , sfvec, "*")
  
  # simulate NB counts
  sim.cnts <- matrix(rnbinom(n = length(mumat), 
                             mu = mumat, 
                             size = sizevec), 
                     ngenes, nsamples)
  
  return(sim.cnts)
}


# EdgeR estimation --------------------------------------------------------

run_estimate <- function(countData, 
                      Ambient = FALSE, 
                      AmbProp = 0.025, 
                      normalisation = 'MR', 
                      downsample = T, 
                      nsigma = 1.96){
  
  # kick out empty genes and samples!
  keep.genes <- rowSums(countData) > 0
  keep.samples <- colSums(countData) > 0
  countData <- countData[keep.genes, keep.samples]
  
  if(isTRUE(Ambient)){
    ambData <- get_ambient_mle(countData = countData, AmbProp = AmbProp)
    countData <- countData - round(ambData)
    countData[countData<0] <- 0
  } else {
    ambData <- NULL
  }
  
  # normalise data
  if(normalisation == "MR"){
    normData <- .MR.calc(countData = countData)
  }
  
  if(normalisation == "TMM"){
    normData <- .TMM.calc(countData = countData)
  }
  if(normalisation == "depth"){
    normData <- .depth.calc(countData = countData)
  }
    
  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))
  
  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = NULL,
                        remove.zeros = FALSE)
  
  # apply edgeR estimation for mean and dispersion
  dge <- suppressMessages(edgeR::estimateDisp(dge))
  
  # method of moments estimates of mean, dispersion, dropout
  mu <- rowMeans(normData$NormCounts)
  lmean <- log2(mu)
  lmu <- rowMeans(log2(normData$NormCounts+1))
  vars <- matrixStats::rowVars(normData$NormCounts)
  size <- mu^2 / (vars -  mu + 1e-04)
  size[size<0] <-  NA
  lsize <- log2(size)
  disp <- 1/size
  ldisp <- log2(disp)
  countData0 = countData == 0
  nsamples = ncol(countData)
  ng0 = rowSums(!countData0)
  g0 = (nsamples - ng0)/nsamples
  
  # estimates from edgeR
  fitestimates <- data.frame("mean" = mu, 
                             'lmean' = lmean,
                             'lmu' = dge$AveLogCPM,
                             'dispersion' = dge$tagwise.dispersion, 
                             'ldispersion' = log2(dge$tagwise.dispersion),
                             "size" = 1/dge$tagwise.dispersion,
                             'lsize' = log2(1/dge$tagwise.dispersion),
                             'dropout' = g0)
  # method of moments estimates
  momestimates <- data.frame("mean" = mu, 
                             'lmean' = lmean,
                             'lmu' = lmu,
                             'dispersion' = disp, 
                             'ldispersion' = ldisp,
                             "size" = size,
                             'lsize' = lsize,
                             'dropout' = g0)
  
  # reduce to available estimates
  fitestimates <- do.call(data.frame,lapply(data.table::data.table(fitestimates), 
                                            function(x) replace(x, is.infinite(x),NA)))
  fitestimates$mean <- ifelse(fitestimates$mean>0, fitestimates$mean, NA)
  # fitestimates <- na.omit(fitestimates)
  momestimates <- do.call(data.frame,lapply(data.table::data.table(momestimates), 
                                            function(x) replace(x, is.infinite(x),NA)))
  momestimates$mean <- ifelse(momestimates$mean>0, momestimates$mean, NA)
  momestimates <- na.omit(momestimates)
  
  # downsample
  if(isTRUE(downsample)){
    fgid <- sample(1:nrow(fitestimates), size = 10000, replace = F)
    fitestimates.red <- fitestimates[fgid,]
    mgid <- sample(1:nrow(momestimates), size = 10000, replace = F)
    momestimates.red <- momestimates[mgid,]
  } else {
    momestimates.red <- momestimates
    fitestimates.red <- fitestimates
  }
  
  # fit mean-size relation
  fit.mvs.fit <- msir::loess.sd(x = fitestimates.red$lmu, y = fitestimates.red$size,
                                nsigma = nsigma)
  mom.mvs.fit <- msir::loess.sd(x = momestimates.red$lmean, y = momestimates.red$lsize,
                                nsigma = nsigma)
  # fit mean-disp relation
  fit.mvd.fit <- msir::loess.sd(x = fitestimates.red$lmu, y = fitestimates.red$dispersion,
                                nsigma = nsigma)
  mom.mvd.fit <- msir::loess.sd(x = momestimates.red$lmean, y = momestimates.red$ldispersion,
                                nsigma = nsigma)
  
  # return object
  res <- list("Estimates" = list(edgeR = fitestimates,
                                 MoM = momestimates,
                                 SF = sf),
              "MeanvsDisp" = list(edgeR = fit.mvd.fit,
                            MoM = mom.mvd.fit),
              "MeanvsSize" = list(edgeR = fit.mvs.fit,
                                  MoM = mom.mvs.fit),
              "AmbientProfile" = ambData,
              "edgeRObj" = dge)
  
  invisible(gc())
  
  return(res)
  
}


# Normalisation -----------------------------------------------------------

.MR.calc <- function(countData) {
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                                         colData=data.frame(group=rep('A', ncol(countData))),
                                                         design=~1))
  dds <- DESeq2::estimateSizeFactors(dds, type='ratio')
  
  sf <- DESeq2::sizeFactors(dds)
  names(sf) <- colnames(countData)
  norm.counts <- DESeq2::counts(dds, normalized=TRUE)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "MR"
  return(res)
}

.TMM.calc <- function(countData) {
  norm.factors <- edgeR::calcNormFactors(object=countData, method='TMM')
  sf <- norm.factors * colSums(countData)
  sf <- sf/exp(mean(log(sf)))
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "TMM"
  return(res)
}

.depth.calc <- function(countData) {
  sf <- colSums(countData) / mean(colSums(countData))
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "depth"
  return(res)
}