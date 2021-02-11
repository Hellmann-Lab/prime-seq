### calculate number of spike in molecules----

Spike_calc<- function(Mix, Dilution, Volume){
  ERCC<- read.table(paste0("/data/share/htp/prime-seq_Paper/genomes/ERCC/ERCC_Mix",Mix,".txt"), header=T, sep="\t")
  colnames(ERCC)<- c("SpikeID","molecules" )
  ERCC$SpikeInput<- Volume*(ERCC$molecules/Dilution)
  SpikeInfo<- ERCC[, c(1,3)]
  rownames(SpikeInfo)<-ERCC[,"SpikeID"]
  return(SpikeInfo)
}

## extract dataframes for plotting

get_Spike_plot_dfs<-function(estSpike,cond,spike_info){
  
  df_list<-list(cal.info.dat=NULL,capture.dat=NULL)
  # calibration curve data
  cal.dat <- reshape2::melt(estSpike$normCounts)
  names(cal.dat) <- c("SpikeID", "SampleID", "normCounts")
  df_list$cal.info.dat <- cal.dat %>%
    dplyr::left_join(spike_info, by="SpikeID") %>%
    dplyr::mutate(FSpike = factor(.data$SpikeInput)) %>%
    dplyr::group_by(.data$SpikeID) %>% ## why not by SpikeID here?
    dplyr::mutate(Expectation=mean(.data$normCounts),
                  LExpectation=mean(log10(.data$normCounts)),
                  Deviation=stats::sd(.data$normCounts),
                  Error=stats::sd(.data$normCounts)/sqrt(dplyr::n()),
                  LError=stats::sd(log10(.data$normCounts))/sqrt(dplyr::n()),
                  CV = Deviation/Expectation ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(LSpikeInput = log10(.data$SpikeInput),
                  LExpectation = log10(.data$Expectation),## why log 10 a second time?
                  cond=as.character(cond))
  # capture efficiency data
  df_list$capture.dat <- estSpike$CaptureEfficiency$`Spike-In` %>%
    tibble::rownames_to_column(var = "SpikeID") %>%
    dplyr::select(.data$SpikeID, .data$p_success,
                  .data$hat_p_success_cilower, .data$hat_p_success_ciupper) %>%
    dplyr::left_join(spike_info, by="SpikeID") %>%
    dplyr::mutate(LSpikeInput = log10(.data$SpikeInput),
                  cond=as.character(cond))
  return(df_list)
}

### Main function

estimateSpike2 <- function(spikeData,
                          spikeInfo,
                          MeanFragLengths = NULL,
                          batchData = NULL,
                          Normalisation=c('depth','none'),
                          SampleFilter = 3,
                          RNAseq = c('bulk', 'singlecell'),
                          Protocol = c('UMI', 'Read'),
                          verbose = TRUE,
                          Theta_start=0.01) {
  # check provided input
  if(!is.null(batchData)) {
    if(!nrow(batchData) == ncol(spikeData)) {
      stop(message(paste0("batchData and spikeData is not the same length!")))
    }
    if(is.null(rownames(batchData))) {
      stop(message(paste0("batchData has no sample annotation as row names!")))
    }
    if(is.null(colnames(spikeData))) {
      stop(message(paste0("spikeData has no sample annotation as column names!")))
    }
  }
  if(is.null(rownames(spikeData)) || is.null(rownames(spikeInfo))) {
    stop(message(paste0("The input data frames of spikeInfo and spikeData need row names!")))
  }
  if(is.null(colnames(spikeInfo)) || !any(colnames(spikeInfo) %in% "SpikeInput")) {
    stop(message(paste0("The input data frame of spikeInfo has no column labelled SpikeInput.")))
  }
  if (!Normalisation %in% c('depth', 'none')) {
    stop(message(paste0("For spike-in count data normalisation, only depth normalization or providing prenormalized data is implemented.")))
  }
  if (Normalisation=='none' && all((spikeData - round(spikeData)) == 0)) {
    stop(message(paste0("Skipping the normalisation should only be done with pre-normalized data, eg. RSEM output, but the provided input matrix contains only integer values!")))
  }
  
  
  # clean out extreme samples and spike-ins prior to normalisation
  if(!is.null(batchData)){
    bData <- as.vector(batchData[, 1])
  }
  if(is.null(batchData)){
    bData <- NULL
  }
  
  # define outliers as determined by SampleFilter using sequencing depth, detected features
  totCounts <- colSums(spikeData)
  libsize.drop <- scater::isOutlier(totCounts, nmads=SampleFilter, type="both",
                                    log=TRUE, batch = bData)
  totFeatures <- colSums(spikeData>0)
  feature.drop <- scater::isOutlier(totFeatures, nmads=SampleFilter, type="both",
                                    log=TRUE, batch = bData)
  minexpr.drop <- colSums(spikeData)<100
  # kick out spike-ins with no expression values
  spike.keep <- rowSums(spikeData )>0
  
  # remove outliers from spikeData
  spikeData.red <- spikeData[spike.keep, c(!libsize.drop | !feature.drop | !minexpr.drop)]
  
  # match spike Info to reduced spikeData
  spikeInfo.red <- spikeInfo[rownames(spikeInfo) %in% rownames(spikeData.red),]
  spikeInfo.red <- spikeInfo.red[match(rownames(spikeData.red), rownames(spikeInfo.red)),]
  
  # Mean fragment lengths
  if(!is.null(MeanFragLengths)) {
    cell.id <- colnames(spikeData.red)
    MeanFragLengths.red <- MeanFragLengths[match(cell.id,names(MeanFragLengths))]
  }
  if(is.null(MeanFragLengths)) {
    MeanFragLengths.red <- NULL
  }
  
  # batches
  # define batches
  if(!is.null(batchData)){
    batchData.red <- batchData[rownames(batchData) %in% colnames(spikeData.red), , drop = FALSE]
  }
  if(is.null(batchData)){
    batchData.red <- NULL
  }
  
  if(verbose){
    message(paste0(nrow(spikeData.red), " spike-ins have been detected in ",
                   ncol(spikeData.red), " cells."))
  }
  
  if(nrow(spikeInfo.red)<10) {
    stop(message(paste0("Not enough spike-ins detected to allow reliable variance estimation.
                          Please proceed without spike-in estimation.")))
  }
  
  if(!is.null(batchData) && !is.null(rownames(batchData))) {
    batch <- stats::setNames(as.character(row.names(batchData.red)), batchData.red[,1])
    # calculate normalized spike-ins per batch
    normspikeData <- sapply(unique(names(batch)), function(b){
      tmpData <- spikeData.red[,grepl(pattern = b, names(batch))]
      normData <- .norm.calc(Normalisation = Normalisation,
                             sf = NULL,
                             countData = tmpData,
                             spikeData = NULL,
                             spikeInfo = NULL,
                             batchData = NULL,
                             Lengths = NULL,
                             MeanFragLengths = NULL,
                             PreclustNumber = NULL,
                             Label = 'none',
                             Step = 'Estimation',
                             Protocol = Protocol,
                             NCores = NULL,
                             verbose = verbose)
      normData
    }, simplify=FALSE)
    
    # estimate gamma theta per batch
    gammaThetaEstimate <- sapply(unique(names(batch)), function(b){
      tmpData <- normspikeData[[b]]$NormCounts
      tmpSF <- normspikeData[[b]]$size.factors
      if(!is.null(spikeInfo.red$Lengths) && !is.null(MeanFragLengths.red)) {
        tmpSFmat <- .repmat(t(as.matrix(tmpSF)), nrow(spikeInfo.red), 1) *
          .repmat(as.matrix( (spikeInfo.red$Lengths - MeanFragLengths.red + 1) / 10^3), 1, ncol(tmpData))
      }
      if(!is.null(spikeInfo.red$Lengths) && is.null(MeanFragLengths.red)) {
        tmpSFmat <- .repmat(t(as.matrix(tmpSF)), nrow(spikeInfo.red), 1) *
          .repmat(as.matrix( spikeInfo.red$Lengths / 10^3), 1, ncol(tmpData))
      }
      if(is.null(spikeInfo.red$Lengths) && is.null(MeanFragLengths.red)) {
        tmpSFmat <- .repmat(t(as.matrix(tmpSF)), nrow(spikeInfo.red), 1)
      }
      est <- .estimateGammaTheta2(nCountSpikes = tmpData,
                                 numberSpikes = spikeInfo.red[,"SpikeInput", drop = FALSE],
                                 sizeFactorMatrix = tmpSFmat,
                                 Theta_start=Theta_start)
      est
    }, simplify=FALSE)
    # correct for batch effect by dividing norm counts with gamma theta
    normCountsL <- sapply(unique(names(batch)), function(b){
      tmpData <- normspikeData[[b]]$NormCounts
      tmpgammaTheta <- gammaThetaEstimate[[b]]$gammaTheta[[1]]
      tmpData / tmpgammaTheta
    }, simplify=FALSE)
    
    # return normalized, batch-corrected counts and size factors for spike-ins
    normCounts <- do.call('cbind', normCountsL)
    sizeFactor <- as.vector(unlist(sapply(normspikeData, function(x) x$size.factors)))
    names(sizeFactor) <- colnames(spikeData.red)
    seqDepth <- colSums(spikeData.red)
    names(seqDepth) <- colnames(spikeData.red)
  }
  
  if(is.null(batchData)) {
    # calculate normalized spike-ins
    normspikeData <- .norm.calc(Normalisation = Normalisation,
                                sf = NULL,
                                countData = spikeData.red,
                                spikeData = NULL,
                                spikeInfo = NULL,
                                batchData = NULL,
                                Lengths = NULL,
                                MeanFragLengths = NULL,
                                PreclustNumber = NULL,
                                Label = 'none',
                                Step = 'Estimation',
                                Protocol = Protocol,
                                NCores = NULL,
                                verbose = verbose)
    # return normalized, batch-corrected counts and size factors for spike-ins
    normCounts <- normspikeData$NormCounts
    sizeFactor <- normspikeData$size.factors
    seqDepth <- colSums(spikeData.red)
    names(seqDepth) <- colnames(spikeData.red)
  }
  
  # technical noise fit
  if(!is.null(spikeInfo.red$Lengths) && !is.null(MeanFragLengths.red)) {
    sizeFactorMatrix <- .repmat(t(as.matrix(sizeFactor)), nrow(spikeInfo.red), 1) *
      .repmat(as.matrix( (spikeInfo.red$Lengths - MeanFragLengths.red + 1) / 10^3), 1, ncol(normCounts))
  }
  if(!is.null(spikeInfo.red$Lengths) && is.null(MeanFragLengths.red)) {
    sizeFactorMatrix <- .repmat(t(as.matrix(sizeFactor)), nrow(spikeInfo.red), 1) *
      .repmat(as.matrix( spikeInfo.red$Lengths / 10^3), 1, ncol(normCounts))
  }
  if(is.null(spikeInfo.red$Lengths) && is.null(MeanFragLengths.red)) {
    sizeFactorMatrix <- .repmat(t(as.matrix(sizeFactor)), nrow(spikeInfo.red), 1)
  }
  EVGammaThetaEstimate <- .estimateEVGammaTheta2(nCountSpikes = normCounts,
                                                numberSpikes = spikeInfo.red[,"SpikeInput", drop = FALSE],
                                                sizeFactorMatrix = sizeFactorMatrix,
                                                Theta_start= Theta_start)
  
  ## gene capture efficiency
  # sum of complete failures (Count=0) and successes= total- failures over all libraries
  spikecapture.dat <- data.frame(detect_fail=apply(spikeData.red==0,1,sum),
                                 detect_success=apply(spikeData.red>0,1,sum),
                                 row.names = rownames(spikeData.red))
  noobs <- ncol(spikeData.red)
  spikecapture.dat[,"p_fail"] <- spikecapture.dat["detect_fail"]/noobs
  spikecapture.dat[,"p_success"] <- spikecapture.dat["detect_success"]/noobs
  
  # confidence intervals of MLE p (applying Wilson interval which is a score based method as normal approximation can be biased
  # if n*p is not sufficiently large which is the case for some of the spike ins (very low p))
  spikecapture.dat [, "hat_p_success_cilower"] <- Hmisc::binconf(spikecapture.dat$detect_success, noobs, method="wilson")[,2]
  spikecapture.dat [, "hat_p_success_ciupper"] <- Hmisc::binconf(spikecapture.dat$detect_success, noobs, method="wilson")[,3]
  spikecapture.dat [, "hat_p_fail_cilower"] <- Hmisc::binconf(spikecapture.dat$detect_fail, noobs, method="wilson")[,2]
  spikecapture.dat [, "hat_p_fail_ciupper"] <- Hmisc::binconf(spikecapture.dat$detect_fail, noobs, method="wilson")[,3]
  
  ## cell capture efficiency
  total_ercc_molecules <- sum(spikeInfo.red$SpikeInput)
  cellcapture.dat <- apply(spikeData.red, 2, function(i) {
    sum(i) / sum(spikeInfo.red$SpikeInput)
  })
  
  ## estimated moments of normalized counts
  normParams <- data.frame(Log2Means=rowMeans(log2(as.matrix(normCounts)+1)),
                           Log2Sd=matrixStats::rowSds(log2(as.matrix(normCounts)+1)),
                           row.names = row.names(normCounts))
  
  # return object
  res <- list("normCounts" = normCounts,
              "normParams" = normParams,
              "CaptureEfficiency"=list("Spike-In"=spikecapture.dat, "Cell" =cellcapture.dat),
              "size.factors" = sizeFactor,
              "seqDepth" = seqDepth,
              "EVGammaThetaEstimates" = EVGammaThetaEstimate,
              "FilteredInput" = list("spikeData"=spikeData.red,
                                     "spikeInfo"=spikeInfo.red,
                                     "batchData"=batchData.red,
                                     "MeanFragLengths"=MeanFragLengths.red),
              "Settings"=list("normFramework"=Normalisation))
  attr(res, 'Protocol') <- Protocol
  attr(res, 'RNAseq') <- RNAseq
  
  return(res)
}


### Helper functions

.repmat <- function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}


.depth.calc <- function(countData, verbose) {
  sf <- colSums(countData) / mean(colSums(countData))
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "depth"
  return(res)
}



.norm.calc <- function(Normalisation,
                       sf,
                       countData,
                       spikeData,
                       spikeInfo,
                       batchData,
                       Lengths,
                       MeanFragLengths,
                       PreclustNumber,
                       Label,
                       Step,
                       Protocol,
                       NCores,
                       verbose) {
  if(Normalisation=='TMM') {
    NormData <- .TMM.calc(countData = countData,
                          verbose = verbose)
  }
  if(Normalisation=='UQ') {
    NormData <- .UQ.calc(countData = countData,
                         verbose = verbose)
  }
  if(Normalisation=='MR') {
    NormData <- .MR.calc(countData = countData,
                         spikeData = spikeData,
                         verbose = verbose)
  }
  if(Normalisation=='PosCounts') {
    NormData <- .PosCounts.calc(countData = countData,
                                spikeData = spikeData,
                                verbose = verbose)
  }
  if(Normalisation=='Linnorm') {
    NormData <- .linnormnorm.calc(countData = countData,
                                  spikeData = spikeData,
                                  verbose = verbose)
  }
  if(Normalisation=='scran' && Label == "none") {
    NormData <- .scran.calc(countData = countData,
                            spikeData = spikeData,
                            verbose = verbose)
  }
  if(Normalisation=='scran' && Label == "known") {
    if(Step=="Estimation"){
      NormData <- .scran.calc(countData = countData,
                              spikeData = spikeData,
                              verbose = verbose)}
    if(Step=="Simulation"){
      NormData <- .scrangroups.calc(countData = countData,
                                    batchData = batchData,
                                    verbose = verbose)}
  }
  if(Normalisation=='scran' && Label == "clustering") {
    NormData <- suppressWarnings(
      .scranclust.calc(countData = countData,
                       PreclustNumber = PreclustNumber,
                       verbose = verbose))
  }
  if(Normalisation=='SCnorm' && Step == "Estimation"){
    NormData <- .SCnorm.calc(countData = countData,
                             spikeData = spikeData,
                             batchData = NULL,
                             NCores = NCores,
                             verbose = verbose)
  }
  if(Normalisation=='SCnorm' && Label == "none" && Step == "Simulation"){
    NormData <- .SCnorm.calc(countData = countData,
                             spikeData = spikeData,
                             batchData = NULL,
                             NCores = NCores,
                             verbose = verbose)
  }
  if(Normalisation=='SCnorm' && Label == "known" && Step == "Simulation"){
    NormData <- .SCnorm.calc(countData = countData,
                             spikeData = spikeData,
                             batchData = batchData,
                             NCores = NCores,
                             verbose = verbose)
  }
  if(Normalisation=='SCnorm' && Label == "clustering" && Step == "Simulation"){
    NormData <- .SCnormclust.calc(countData = countData,
                                  spikeData = spikeData,
                                  PreclustNumber=PreclustNumber,
                                  NCores = NCores,
                                  verbose = verbose)
  }
  if(Normalisation=='sctransform') {
    NormData <- .sctransform.calc(countData = countData,
                                  batchData = batchData,
                                  Step = Step,
                                  NCores = NCores,
                                  verbose = verbose)
  }
  if(Normalisation=="bayNorm") {
    NormData <- .baynorm.calc(countData = countData,
                              batchData = batchData,
                              spikeData = spikeData,
                              spikeInfo = spikeInfo,
                              NCores  = NCores,
                              Step = Step,
                              Protocol = Protocol,
                              verbose = verbose)
  }
  if(Normalisation=='Census') {
    NormData <- .Census.calc(countData = countData,
                             Lengths = Lengths,
                             MeanFragLengths = MeanFragLengths,
                             spikeData = spikeData,
                             spikeInfo = spikeInfo,
                             Protocol = Protocol,
                             NCores = NCores,
                             verbose = verbose)
  }
  if(Normalisation=='depth') {
    NormData <- .depth.calc(countData = countData,
                            verbose = verbose)
  }
  if(Normalisation=='SF') {
    NormData <- .sf.calc(countData = countData,
                         sf = sf,
                         verbose = verbose)
  }
  if(Normalisation=='none') {
    NormData <- .none.calc(countData = countData,
                           verbose = verbose)
  }
  
  
  # inform the user of unusual / problematic normalisation results
  
  if(attr(NormData, "normFramework")=="scran" && any(NormData$size.factors<0)){
    if (verbose) { message(paste0("Negative size factors estimated.
                                  Apply stronger gene and sample filtering for parameter estimation.")) }
  }
  
  return(NormData)
}

.estimateEVGammaTheta2 <- function(nCountSpikes, numberSpikes, sizeFactorMatrix,Theta_start) {
  gammaThetaEstimate = .estimateGammaTheta2(nCountSpikes, numberSpikes, sizeFactorMatrix,Theta_start = Theta_start)
  EGamma = gammaThetaEstimate$gammaTheta[[1]] / gammaThetaEstimate$Theta[[1]]
  ETheta = gammaThetaEstimate$Theta[[1]]
  E2Gamma = EGamma^2
  E2Theta = ETheta^2
  
  varianceGammaThetaEstimate = powsimR:::.estimateVGammaThetaInitial(nCountSpikes, numberSpikes, sizeFactorMatrix)
  VGamma = varianceGammaThetaEstimate[1]
  VTheta = varianceGammaThetaEstimate[2]
  
  fitData = data.frame(Xi=numberSpikes[,1],
                       Y=apply(nCountSpikes,1,stats::var)/rowMeans(nCountSpikes),
                       Bi=rowMeans(1/sizeFactorMatrix))
  invisible(capture.output(
    fit <- suppressMessages(
      minpack.lm::nlsLM(Y ~ Bi + (VGamma+E2Gamma)/EGamma*(1-(E2Theta+VTheta)/ETheta) +
                          (VGamma+E2Gamma)*VTheta*Xi/(EGamma*ETheta) + (ETheta*VGamma)*Xi/EGamma, data=fitData,
                        start=list(VGamma=VGamma, VTheta=VTheta), lower=c(0,0))
    )
  ))
  VGamma=coefficients(fit)[[1]]
  VTheta=coefficients(fit)[[2]]
  # second optimization for robust estimates
  invisible(capture.output(
    fit <- suppressMessages(
      minpack.lm::nlsLM(Y ~ Bi + (VGamma+E2Gamma)/EGamma*(1-(E2Theta+VTheta)/ETheta) +
                          (VGamma+E2Gamma)*VTheta*Xi/(EGamma*ETheta) + (ETheta*VGamma)*Xi/EGamma, data=fitData,
                        start=list(VGamma=VGamma, VTheta=VTheta), lower=c(0,0))
    )
  ))
  VGamma=coefficients(fit)[[1]]
  VTheta=coefficients(fit)[[2]]
  list(EGamma=EGamma, ETheta=ETheta, E2Gamma=E2Gamma, E2Theta=E2Theta, VGamma=VGamma, VTheta=VTheta)
}

.estimateGammaTheta2 <- function(nCountSpikes, numberSpikes, sizeFactorMatrix,Theta_start) {
  fitData = data.frame(Xi=numberSpikes[,1], Ki=rowMeans(nCountSpikes))
  fit1 <- stats::lm(Ki ~ 0 + Xi, data=fitData)
  gammaTheta = stats::coefficients(fit1)
  PropCell = sapply(1:nrow(nCountSpikes), function(x) {
    sum(nCountSpikes[x,]>0)/ncol(nCountSpikes)
  })
  fitData = data.frame(Xi=numberSpikes[,1], Y=PropCell, Ai=rowMeans(sizeFactorMatrix))
  initialTheta = mean(fitData$Y[fitData$Xi>0.5 & fitData$Xi<5])
  if (initialTheta==0 | is.na(initialTheta)) {
    initialTheta = Theta_start
  }
  invisible(capture.output(
    fit2 <- suppressMessages(
      minpack.lm::nlsLM(Y ~ 1 - (1 - Theta + Theta*exp(-(gammaTheta/Theta)*Ai))^Xi,
                        data=fitData, start=list(Theta=initialTheta),
                        lower=0, upper=1, control=stats::nls.control(warnOnly=TRUE))
    )
  ))
  invisible(capture.output(
    Theta <- suppressMessages(
      stats::coefficients(fit2)
    )
  ))
  
  list(gammaTheta=gammaTheta, Theta=Theta)
}

#### evaluate DE

evaluateDE2<-function (simRes, alpha.type = c("adjusted", "raw"), MTC = c("BY", 
                                                             "BH", "holm", "hochberg", "hommel", "bonferroni", "Storey", 
                                                             "IHW"), alpha.nominal = 0.1, stratify.by = c("mean", "dispersion", 
                                                                                                          "dropout", "lfc","lfc_abs"), strata.probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 
                                                                                                                                              0.6, 0.7, 0.8, 0.9), filter.by = c("none", "mean", "dispersion", 
                                                                                                                                                                                 "dropout"), strata.filtered = 1, target.by = c("lfc", "effectsize"), 
          delta = 0, Table = TRUE) 
{
  alpha.type = match.arg(alpha.type)
  MTC = match.arg(MTC)
  stratify.by = match.arg(stratify.by)
  filter.by = match.arg(filter.by)
  target.by = match.arg(target.by)
  Nreps1 = simRes$SimSetup$n1
  Nreps2 = simRes$SimSetup$n2
  ngenes = simRes$DESetup$ngenes
  DEids = simRes$DESetup$DEid
  lfcs = simRes$DESetup$pLFC
  tlfcs = lapply(1:length(lfcs), function(i) {
    lfcs[[i]]
  })
  nsims = simRes$DESetup$nsims
  estmeans = simRes$estParamRes$Parameters[[simRes$DESetup$Draw$MoM]]$means
  estdisps = simRes$estParamRes$Parameters[[simRes$DESetup$Draw$MoM]]$dispersion
  estdropout = simRes$estParamRes$Parameters[[simRes$DESetup$Draw$MoM]]$gene.dropout
  mu = simRes$SimulateRes$mu
  disp = simRes$SimulateRes$disp
  dropout = simRes$SimulateRes$dropout
  elfc = simRes$SimulateRes$elfc
  DEmethod = simRes$Pipeline$DEmethod
  pvalue = simRes$SimulateRes$pvalue
  fdr = simRes$SimulateRes$fdr
  tmp.ecdf.mean = stats::ecdf(log2(estmeans + 1))
  tmp.quantile.mean = stats::quantile(tmp.ecdf.mean, probs = strata.probs)
  strata.mean = unique(c(0, unname(tmp.quantile.mean), Inf))
  strata.mean = unique(round(strata.mean, digits = 2))
  tmp.ecdf.disps = stats::ecdf(log2(estdisps + 1))
  tmp.quantile.disps = stats::quantile(tmp.ecdf.disps, probs = strata.probs)
  strata.disps = unique(c(0, unname(tmp.quantile.disps), Inf))
  strata.disps = unique(round(strata.disps, digits = 2))
  tmp.ecdf.drop = stats::ecdf(estdropout)
  tmp.quantile.drop = stats::quantile(tmp.ecdf.drop, probs = strata.probs)
  strata.drop = unique(c(0, unname(tmp.quantile.drop), 1))
  strata.drop = unique(round(strata.drop, digits = 2))
  tmp.ecdf.lfc = stats::ecdf(unique(unlist(tlfcs)))
  tmp.quantile.lfc = stats::quantile(tmp.ecdf.lfc, probs = strata.probs)
  strata.lfc = unique(c(-Inf, unname(tmp.quantile.lfc), Inf))
  strata.lfc = unique(round(strata.lfc, digits = 2))
  tmp.ecdf.lfc.abs = stats::ecdf(unique(abs(unlist(tlfcs))))
  tmp.quantile.lfc.abs = stats::quantile(tmp.ecdf.lfc.abs, probs = strata.probs)
  strata.lfc.abs = unique(c(0, unname(tmp.quantile.lfc.abs), Inf))
  strata.lfc.abs = unique(round(strata.lfc.abs, digits = 2))
  
  if (stratify.by == "mean") {
    nr = length(strata.mean) - 1
  }
  if (stratify.by == "dispersion") {
    nr = length(strata.disps) - 1
  }
  if (stratify.by == "dropout") {
    nr = length(strata.drop) - 1
  }
  if (stratify.by == "lfc") {
    nr = length(strata.lfc) - 1
  }
  if (stratify.by == "lfc_abs") {
    nr = length(strata.lfc.abs) - 1
  }
  if (filter.by %in% c("mean", "dispersion", "dropout")) {
    if (filter.by == stratify.by) {
      nstrata = nr - strata.filtered
    }
    if (!filter.by == stratify.by) {
      nstrata = nr
    }
  }
  if (filter.by == "none") {
    nstrata = nr
  }
  TP = TN = FP = FN = TPR = TNR = FPR = FNR = FDR = xgrl = xgrld = array(NA, 
                                                                         dim = c(nstrata, length(Nreps1), nsims))
  TP.marginal = TN.marginal = FP.marginal = FN.marginal = TPR.marginal = TNR.marginal = FPR.marginal = FNR.marginal = FDR.marginal = matrix(NA, 
                                                                                                                                            length(Nreps1), nsims)
  for (i in 1:nsims) {
    for (j in seq(along = Nreps1)) {
      Nrep1 = Nreps1[j]
      Nrep2 = Nreps2[j]
      DEid = DEids[[i]]
      lfc = lfcs[[i]]
      lfc_abs = abs(lfcs[[i]])
      Zg = Zg2 = rep(0, ngenes)
      Zg[DEid] = 1
      if (delta == 0) {
        Zg2 = Zg
      }
      if (!delta == 0) {
        if (target.by == "lfc") {
          ix = abs(lfc) > delta
        }
        else if (target.by == "effectsize") {
          effectsize = lfc/sqrt(1/((log2(mu[, j, i] + 
                                           1) + log2(disp[, j, i] + 1))))
          ix = abs(effectsize) > delta
        }
        Zg2[ix] = 1
      }
      X.bar1 = mu[, j, i]
      ix.keep.mean = which(!is.na(X.bar1))
      xgr.mean = cut(log2(X.bar1[ix.keep.mean] + 1), strata.mean)
      xgrd.mean = cut(log2(X.bar1[DEid] + 1), strata.mean)
      X.disp1 = disp[, j, i]
      ix.keep.disps = which(!is.na(X.disp1))
      xgr.disps = cut(log2(X.disp1[ix.keep.disps] + 1), 
                      strata.disps)
      xgrd.disps = cut(log2(X.disp1[DEid] + 1), strata.disps)
      X.drop1 = dropout[, j, i]
      ix.keep.drop = which(!is.na(X.drop1))
      xgr.drop = cut(X.drop1[ix.keep.drop], strata.drop)
      xgrd.drop = cut(X.drop1[DEid], strata.drop)
      X.lfc1 = elfc[, j, i]
      ix.keep.lfc = which(!is.na(X.lfc1))
      xgr.lfc = cut(X.lfc1[ix.keep.lfc], strata.lfc)
      xgrd.lfc = cut(X.lfc1[DEid], strata.lfc)
      X.lfc_abs1 = abs(elfc[, j, i])
      ix.keep.lfc_abs = which(!is.na(X.lfc_abs1))
      xgr.lfc_abs = cut(X.lfc_abs1[ix.keep.lfc_abs], strata.lfc.abs)
      xgrd.lfc_abs = cut(X.lfc_abs1[DEid], strata.lfc.abs)
      if (stratify.by == "mean") {
        if (filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          fix.keep.mean = ix.keep.mean[!(xgr.mean %in% 
                                           lev.mean[strata.filt.mean])]
          fix.dekeep.mean = intersect(fix.keep.mean, 
                                      DEid)
          fxgr.mean = cut(log2(X.bar1[fix.keep.mean] + 
                                 1), strata.mean[-strata.filt.mean])
          fxgrd.mean = cut(log2(X.bar1[fix.dekeep.mean] + 
                                  1), strata.mean[-strata.filt.mean])
        }
        if (filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((length(lev.disps) - 
                                   strata.filtered + 1):length(lev.disps))
          fix.keep.mean = ix.keep.mean[!(xgr.disps %in% 
                                           lev.disps[strata.filt.disps])]
          fix.dekeep.mean = intersect(fix.keep.mean, 
                                      DEid)
          fxgr.mean = cut(log2(X.bar1[fix.keep.mean] + 
                                 1), strata.mean)
          fxgrd.mean = cut(log2(X.bar1[fix.dekeep.mean] + 
                                  1), strata.mean)
        }
        if (filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((length(lev.drop) - strata.filtered + 
                                  1):length(lev.drop))
          fix.keep.mean = ix.keep.mean[!(xgr.drop %in% 
                                           lev.drop[strata.filt.drop])]
          fix.dekeep.mean = intersect(fix.keep.mean, 
                                      DEid)
          fxgr.mean = cut(log2(X.bar1[fix.keep.mean] + 
                                 1), strata.mean)
          fxgrd.mean = cut(log2(X.bar1[fix.dekeep.mean] + 
                                  1), strata.mean)
        }
        if (filter.by == "none") {
          fix.keep.mean = ix.keep.mean
          fxgr.mean = xgr.mean
          fxgrd.mean = xgrd.mean
        }
      }
      if (stratify.by == "dispersion") {
        if (filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          fix.keep.disps = ix.keep.disps[!(xgr.mean %in% 
                                             lev.mean[strata.filt.mean])]
          fix.dekeep.disps = intersect(fix.keep.disps, 
                                       DEid)
          fxgr.disps = cut(log2(X.disp1[fix.keep.disps] + 
                                  1), strata.disps)
          fxgrd.disps = cut(log2(X.disp1[fix.dekeep.disps] + 
                                   1), strata.disps)
        }
        if (filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          lev.filt.disps = c((length(lev.disps) - strata.filtered + 
                                1):length(lev.disps))
          fix.keep.disps = ix.keep.disps[!(xgr.disps %in% 
                                             lev.disps[lev.filt.disps])]
          fix.dekeep.disps = intersect(fix.keep.disps, 
                                       DEid)
          strata.filt.disps = length(strata.disps) - 
            strata.filtered
          fxgr.disps = cut(log2(X.disp1[fix.keep.disps] + 
                                  1), strata.disps[1:strata.filt.disps])
          fxgrd.disps = cut(log2(X.disp1[fix.dekeep.disps] + 
                                   1), strata.disps[1:strata.filt.disps])
        }
        if (filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((length(lev.drop) - strata.filtered + 
                                  1):length(lev.drop))
          fix.keep.disps = ix.keep.disps[!(xgr.drop %in% 
                                             lev.drop[strata.filt.drop])]
          fix.dekeep.disps = intersect(fix.keep.disps, 
                                       DEid)
          fxgr.disps = cut(X.disp1[fix.keep.disps], strata.disps)
          fxgrd.disps = cut(X.disp1[fix.dekeep.disps], 
                            strata.disps)
        }
        if (filter.by == "none") {
          fix.keep.disps = ix.keep.disps
          fxgr.disps = xgr.disps
          fxgrd.disps = xgrd.disps
        }
      }
      if (stratify.by == "dropout") {
        if (filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          fix.keep.drop = ix.keep.drop[!(xgr.mean %in% 
                                           lev.mean[strata.filt.mean])]
          fix.dekeep.drop = intersect(fix.keep.drop, 
                                      DEid)
          fxgr.drop = cut(X.drop1[fix.keep.drop], strata.drop)
          fxgrd.drop = cut(X.drop1[fix.dekeep.drop], 
                           strata.drop)
        }
        if (filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((length(lev.disps) - 
                                   strata.filtered + 1):length(lev.disps))
          fix.keep.drop = ix.keep.drop[!(xgr.disps %in% 
                                           lev.disps[strata.filt.disps])]
          fix.dekeep.drop = intersect(fix.keep.drop, 
                                      DEid)
          fxgr.drop = cut(X.drop1[fix.keep.drop], strata.drop)
          fxgrd.drop = cut(X.drop1[fix.dekeep.drop], 
                           strata.drop)
        }
        if (filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((length(lev.drop) - strata.filtered + 
                                  1):length(lev.drop))
          fix.keep.drop = ix.keep.drop[!(xgr.drop %in% 
                                           lev.drop[strata.filt.drop])]
          fix.dekeep.drop = intersect(fix.keep.drop, 
                                      DEid)
          strata.filt.drop = length(strata.drop) - strata.filtered
          fxgr.drop = cut(X.drop1[fix.keep.drop], strata.drop[1:strata.filt.drop])
          fxgrd.drop = cut(X.drop1[fix.dekeep.drop], 
                           strata.drop[1:strata.filt.drop])
        }
        if (filter.by == "none") {
          fix.keep.drop = ix.keep.drop
          fxgr.drop = xgr.drop
          fxgrd.drop = xgrd.drop
        }
      }
      if (stratify.by == "lfc") {
        if (filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          fix.keep.lfc = ix.keep.lfc[!(xgr.mean %in% 
                                         lev.mean[strata.filt.mean])]
          fix.dekeep.lfc = intersect(fix.keep.lfc, DEid)
          fxgr.lfc = cut(X.lfc1[fix.keep.lfc], strata.lfc)
          fxgrd.lfc = cut(X.lfc1[fix.dekeep.lfc], strata.lfc)
        }
        if (filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((length(lev.disps) - 
                                   strata.filtered + 1):length(lev.disps))
          fix.keep.lfc = ix.keep.lfc[!(xgr.disps %in% 
                                         lev.mean[strata.filt.disps])]
          fix.dekeep.lfc = intersect(fix.keep.lfc, DEid)
          fxgr.lfc = cut(X.lfc1[fix.keep.lfc], strata.lfc)
          fxgrd.lfc = cut(X.lfc1[fix.dekeep.lfc], strata.lfc)
        }
        if (filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((length(lev.drop) - strata.filtered + 
                                  1):length(lev.drop))
          fix.keep.lfc = ix.keep.lfc[!(xgr.drop %in% 
                                         lev.drop[strata.filt.drop])]
          fix.dekeep.lfc = intersect(fix.keep.lfc, DEid)
          fxgr.lfc = cut(X.lfc1[fix.keep.lfc], strata.lfc)
          fxgrd.lfc = cut(X.lfc1[fix.dekeep.lfc], strata.lfc)
        }
        if (filter.by == "none") {
          fix.keep.lfc = ix.keep.lfc
          fxgr.lfc = xgr.lfc
          fxgr.lfc = xgr.lfc
          fxgrd.lfc = xgrd.lfc
        }
      }
      if (stratify.by == "lfc_abs") {
        if (filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          fix.keep.lfc_abs = ix.keep.lfc_abs[!(xgr.mean %in% 
                                         lev.mean[strata.filt.mean])]
          fix.dekeep.lfc_abs = intersect(fix.keep.lfc_abs, DEid)
          fxgr.lfc_abs = cut(X.lfc_abs1[fix.keep.lfc_abs], strata.lfc.abs)
          fxgrd.lfc_abs = cut(X.lfc_abs1[fix.dekeep.lfc_abs], strata.lfc.abs)
        }
        if (filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((length(lev.disps) - 
                                   strata.filtered + 1):length(lev.disps))
          fix.keep.lfc_abs = ix.keep.lfc_abs[!(xgr.disps %in% 
                                         lev.mean[strata.filt.disps])]
          fix.dekeep.lfc_abs = intersect(fix.keep.lfc_abs, DEid)
          fxgr.lfc_abs = cut(X.lfc_abs1[fix.keep.lfc_abs], strata.lfc.abs)
          fxgrd.lfc_abs = cut(X.lfc_abs1[fix.dekeep.lfc_abs], strata.lfc.abs)
        }
        if (filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((length(lev.drop) - strata.filtered + 
                                  1):length(lev.drop))
          fix.keep.lfc_abs = ix.keep.lfc_abs[!(xgr.drop %in% 
                                         lev.drop[strata.filt.drop])]
          fix.dekeep.lfc_abs = intersect(fix.keep.lfc_abs, DEid)
          fxgr.lfc_abs = cut(X.lfc_abs1[fix.keep.lfc_abs], strata.lfc.abs)
          fxgrd.lfc_abs = cut(X.lfc_abs1[fix.dekeep.lfc_abs], strata.lfc.abs)
        }
        if (filter.by == "none") {
          fix.keep.lfc_abs = ix.keep.lfc_abs
          fxgr.lfc_abs = xgr.lfc_abs
          fxgr.lfc_abs = xgr.lfc_abs
          fxgrd.lfc_abs = xgrd.lfc_abs
        }
      }
      if (stratify.by == "mean") {
        strata = strata.mean
        xgr = fxgr.mean
        xgrd = fxgrd.mean
        ix.keep = fix.keep.mean
      }
      if (stratify.by == "dispersion") {
        strata = strata.disps
        xgr = fxgr.disps
        xgrd = fxgrd.disps
        ix.keep = fix.keep.disps
      }
      if (stratify.by == "dropout") {
        strata = strata.drop
        xgr = fxgr.drop
        xgrd = fxgrd.drop
        ix.keep = fix.keep.drop
      }
      if (stratify.by == "lfc") {
        strata = strata.lfc
        xgr = fxgr.lfc
        xgrd = fxgrd.lfc
        ix.keep = fix.keep.lfc
      }
      if (stratify.by == "lfc_abs") {
        strata = strata.lfc.abs
        xgr = fxgr.lfc_abs
        xgrd = fxgrd.lfc_abs
        ix.keep = fix.keep.lfc_abs
      }
      if (alpha.type == "raw") {
        if (DEmethod %in% c("edgeR-QL", "edgeR-LRT", 
                            "limma-voom", "limma-trend", "NBPSeq", "T-Test", 
                            "DESeq2", "ROTS", "MAST", "scde", "BPSC", "scDD", 
                            "monocle", "DECENT", "edgeR-zingeR", "edgeR-ZINB-WaVE", 
                            "DESeq2-zingeR", "DESeq2-ZINB-WaVE")) {
          x = pvalue[ix.keep, j, i]
          x[is.na(x)] = 1
        }
        if (DEmethod %in% c("baySeq", "NOISeq", "EBSeq")) {
          message(paste0("The DE method ", DEmethod, 
                         " only provides adjusted p-values."))
          x = fdr[ix.keep, j, i]
          x[is.na(x)] = 1
        }
      }
      if (alpha.type == "adjusted") {
        if (DEmethod %in% c("edgeR-QL", "edgeR-LRT", 
                            "limma-voom", "limma-trend", "NBPSeq", "T-Test", 
                            "DESeq2", "ROTS", "MAST", "scde", "BPSC", "scDD", 
                            "monocle", "DECENT", "edgeR-zingeR", "edgeR-ZINB-WaVE", 
                            "DESeq2-zingeR", "DESeq2-ZINB-WaVE")) {
          pval = pvalue[ix.keep, j, i]
          meanexpr = log2(mu[ix.keep, j, i] + 1)
          if (MTC %in% stats::p.adjust.methods) {
            x = stats::p.adjust(pval, method = MTC)
            x[is.na(x)] = 1
          }
          if (MTC %in% "Storey") {
            tmp.p = pval[!is.na(pval)]
            tmp.q = qvalue::qvalue(p = tmp.p)$qvalues
            x = rep(NA, length(pval))
            x[!is.na(pval)] = tmp.q
            x[is.na(x)] = 1
          }
          if (MTC %in% "IHW") {
            in.dat = data.frame(pvalue = pval, meanexpr = meanexpr)
            tmp = IHW::ihw(pvalue ~ meanexpr, data = in.dat, 
                           alpha = alpha.nominal)
            x = IHW::adj_pvalues(tmp)
            x[is.na(x)] = 1
          }
        }
        if (DEmethod %in% c("baySeq", "NOISeq", "EBSeq")) {
          message(paste0("The DE method ", DEmethod, 
                         " only provides adjusted p-values."))
          x = fdr[ix.keep, j, i]
          x[is.na(x)] = 1
        }
      }
      Zg = Zg[ix.keep]
      Zg2 = Zg2[ix.keep]
      xgrl[, j, i] = table(xgr)
      xgrld[, j, i] = table(xgrd)
      error.mat = powsimR:::.error.matrix(p = x, p.crit = alpha.nominal, 
                                Zg = Zg, Zg2 = Zg2, xgr = xgr)
      TP[, j, i] = error.mat$TP
      TN[, j, i] = error.mat$TN
      FP[, j, i] = error.mat$FP
      FN[, j, i] = error.mat$FN
      TP.marginal[j, i] = error.mat$TP.marginal
      TN.marginal[j, i] = error.mat$TN.marginal
      FP.marginal[j, i] = error.mat$FP.marginal
      FN.marginal[j, i] = error.mat$FN.marginal
      TPR[, j, i] = error.mat$TPR
      TNR[, j, i] = error.mat$TNR
      FPR[, j, i] = error.mat$FPR
      FNR[, j, i] = error.mat$FNR
      FDR[, j, i] = error.mat$FDR
      TPR.marginal[j, i] = error.mat$TPR.marginal
      TNR.marginal[j, i] = error.mat$TNR.marginal
      FPR.marginal[j, i] = error.mat$FPR.marginal
      FNR.marginal[j, i] = error.mat$FNR.marginal
      FDR.marginal[j, i] = error.mat$FDR.marginal
    }
  }
  output <- list(stratagenes = xgrl, stratadiffgenes = xgrld, 
                 TN = TN, TP = TP, FP = FP, FN = FN, TN.marginal = TN.marginal, 
                 TP.marginal = TP.marginal, FP.marginal = FP.marginal, 
                 FN.marginal = FN.marginal, TNR = TNR, TPR = TPR, FPR = FPR, 
                 FNR = FNR, FDR = FDR, TNR.marginal = TNR.marginal, TPR.marginal = TPR.marginal, 
                 FPR.marginal = FPR.marginal, FNR.marginal = FNR.marginal, 
                 FDR.marginal = FDR.marginal, alpha.type = alpha.type, 
                 MTC = ifelse(alpha.type == "adjusted", MTC, "not applicable"), 
                 alpha.nominal = alpha.nominal, stratify.by = stratify.by, 
                 strata = strata, strata.levels = levels(xgr), target.by = target.by, 
                 n1 = Nreps1, n2 = Nreps2, delta = delta)
  if (Table) {
    printEvalDE(evalRes = output)
  }
  return(output)
}

