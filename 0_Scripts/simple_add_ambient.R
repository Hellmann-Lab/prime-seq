###
# real data
tmp<- readRDS("/data/share/htp/prime-seq_Paper/Fig_CrossCont/zUMIs/count_multi/zUMIs_output/expression/CrossCont_50bp_nomult.dgecounts.rds")
mtx <- as.matrix( tmp$umicount$exon$all)

mtxh <- mtx[ grepl("ENSG",rownames(mtx)), ]
mtxh <- mtx[rowSums(mtx)>10,]
mtxh <- mtx[,colSums(mtx)>1e5]

tmp<-mtxh[]

######### Functions -------------------
get_ambient_mle<-function(y, amb=0.05){
  
  if( length(amb)==1){
    amb<-rep(amb,ncol(y))
  }else if( length(amb)== ncol(y)){
    print("use provided different ambient proportions.")
  }else{ 
    print("ERROR. amb not equal number of samples.")
  }
  sapply(1:ncol(y), function(i){ sum(y[,i])*amb[i]*rowSums(y)/sum(y)} )
}

#####
get_ambient_sample<-function(y, amb=0.05){
  
  if( length(amb)==1){
     amb<-rep(amb,ncol(y))
  }else if( length(amb)== ncol(y)){
    print("use provided different ambient proportions.")
  }else{ 
    print("ERROR. amb not equal number of samples.")
  }
  
  mtx<-matrix(0, nrow = nrow(y), ncol = ncol(y))
  for(i in 1:ncol(y) ) {
    xx<-sample(1:nrow(y),replace = T, size= sum(y[,i])*amb[i], prob=rowSums(y)/sum(y)) %>% table 
    mtx[as.numeric(names(xx)),i]<-xx
  }
  
  return(mtx)
}
#####
size_est<-function(mat_norm_j,lowNA=T){
  mean_j <- Matrix::rowMeans(mat_norm_j)
  var_j <- sparseMatrixStats::rowVars(mat_norm_j)
  size_j <- mean_j^2 / (var_j -  mean_j + 1e-04)
  if(lowNA){
    size_j <- ifelse(size_j > 0, size_j, NA)
  }else{
    msj<-min(size_j[size_j>0])
    size_j <- ifelse(size_j > 0, size_j, msj)
   }
  size_j
  }

depth.calc <- function(countData) {
  sf <- colSums(countData) / mean(colSums(countData))
  names(sf) <- colnames(countData)
  norm.counts <- t(t(countData)/sf)
  res <- list(NormCounts=norm.counts,
              RoundNormCounts=round(norm.counts),
              size.factors=sf)
  attr(res, 'normFramework') <- "depth"
  return(res)
}

######## Test Simulation

# initial matrix of 100x10 
library(gamlss.dist)
sf <- rep(1, 10)
muvec <- rIGAMMA(n = 100, mu = 5, sigma = 0.5)
mumat <- outer(muvec , sf, "*")
bcv = 0.2
sizevec <- 1/bcv^2

y <- matrix(rnbinom(n = length(mumat), 
                    mu = mumat, size = sizevec), 
            100, 10)
# assume 5% ambient add to y
# so this would be the observed count matrix where we have endogenous plus ambient counts
y_amb<- y + get_ambient_sample(y)

# with real data we would start here:
normalisation <- TRUE
if(normalisation){
  norm_data <- depth.calc(countData = y_amb)
  y_amb <- norm_data$NormCounts
}

y_amb_corr<- y_amb - get_ambient_mle(y_amb)
y_amb_corr[y_amb_corr<0]=0

#simulate on the corrected matrix
mumat <- outer(rowMeans(y_amb_corr) , sf, "*")
y2 <- matrix(rnbinom(n = length(mumat), 
                     mu = mumat,
                     size = size_est(y_amb_corr,lowNA = F)), 
            100, 10)
#add ambient
y2_amb<-y2+get_ambient_sample(y2)

data.frame( `start clean`   = size_est(y),
            `start+ambient` = size_est(y_amb),
            `ambient corr`  = size_est(y_amb_corr),
            `simulated` = size_est(y2),
            `simulated+amb` = size_est(y2_amb)) %>% 
  pivot_longer(cols = 1:5,names_to = "type", values_to = "size") %>% 
  ggplot(aes(y=1/size,x =type, col=type,fill=type)) + geom_violin()



# DropletUtils ------------------------------------------------------------

# BiocManager::install('DropletUtils')
library(DropletUtils)

# function removeAmbience()
# Estimate and remove the ambient profile from a count matrix, given pre-existing groupings of similar cells.
# This function will aggregate counts from each group of related cells into an average profile. For each group, we estimate the contribution of the ambient profile and subtract it from the average. By default, this is done with ambientContribMaximum, but if enough is known about the biological system, users can specify feaures to use ambientContribNegative instead. We then perform quantile-quantile mapping of counts in y from the old to new averages. This approach preserves the mean-variance relationship and improves the precision of estimate of the ambient contribution, but relies on a sensible grouping of similar cells, e.g., unsupervised clusters or cell type annotations. As such, this function is best used at the end of the analysis to clean up expression matrices prior to visualization.


# Making up some data.
ngenes <- 1000
ambient <- runif(ngenes, 0, 0.1)
cells <- c(runif(100) * 10, integer(900)) # ie first 100 genes are expressed, remaining are empty
y <- matrix(rpois(n = ngenes * 100, lambda = cells + ambient), nrow=ngenes)
y[1:10,1:5]
# Pretending that all cells are in one group, in this example.
removed <- removeAmbience(y, ambient, groups=rep(1, ncol(y)))
summary(rowMeans(y[1:100,]))
summary(rowMeans(removed[1:100,]))
summary(rowMeans(removed[101:1000,]))
summary(rowMeans(y[101:1000,]))
