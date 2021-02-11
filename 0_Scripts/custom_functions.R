### This is a collection of custom functions used in the analysis of the prime-seq paper data
### general functions
### PCA from Klaus_pcafunction.R
pcaFunction12<-function(mat, inf, genes, ngenes, alpha, col, label="sid",size) { 

if(!missing(genes)){
  vsdMat <- mat[genes,rownames(inf)]
  pc<-prcomp(t(vsdMat),scale=F)
  pc.sum<-summary(pc)$importance
  varExp<-round(pc.sum[2,]*100,2)
  pcs<-data.frame(pc$x, inf)
  pcs$sid<-rownames(pcs)
}
  
rowVar<-apply(mat,1,var)

if(!missing(ngenes)){
  mvg<-order(rowVar,decreasing = T)[1:ngenes]
  vsdMat<-mat[mvg,rownames(inf)]
  pc<-prcomp(t(vsdMat),scale=F)
  pc.sum<-summary(pc)$importance
  varExp<-round(pc.sum[2,]*100,2)
  pcs<-data.frame(pc$x, inf)
  pcs$sid<-rownames(pcs)
}

if(!missing(genes)&!missing(ngenes)){
  vsdMat<- mat[,rownames(inf)]
  pc<-prcomp(t(vsdMat),scale=F)
  pc.sum<-summary(pc)$importance
  varExp<-round(pc.sum[2,]*100,2)
  pcs<-data.frame(pc$x, inf)
  pcs$sid<-rownames(pcs)
}

if(missing(size)){
  size<-3
}



pp<-ggplot(pcs,aes_q(x=quote(PC1),y=quote(PC2),col=as.name(col)))
if(!missing(label)){
  pp<-pp+geom_text(aes(label=as.name(label)),size=size)
}else{
    pp<-pp+geom_point(alpha=alpha,size=size)
  }

pp+xlab(paste("PC1 (",varExp[1],"%)"))+
  ylab(paste("PC2 (",varExp[2],"%)"))
}

pcaFunction23<-function(mat, inf, genes, ngenes, alpha, col="treatment", label="sid",shape,size ) { 
  
  if(!missing(genes)){
    vsdMat <- mat[genes,rownames(inf)]
  }else{
    vsdMat<- mat[,rownames(inf)]
  }
  rowVar<-apply(mat,1,var)
  
  if(!missing(ngenes)){
    mv500<-order(rowVar,decreasing = T)[1:ngenes]
  }else{
    mv500<-(rowVar>0)
  }
  if(missing(shape)){
    shape=col
  }
  
  pc<-prcomp(t(vsdMat[mv500,]),scale=T)
  pc.sum<-summary(pc)$importance
  varExp<-round(pc.sum[2,]*100,2)
  pcs<-data.frame(pc$x, inf)
  pcs$sid<-rownames(pcs)
  pp<-ggplot(pcs,aes_q(x=quote(PC2),y=quote(PC3),col=as.name(col)))
  if(!missing(label)){
    pp<-pp+geom_text(aes(label=as.name(label)),size=size)
  }else{
    if(missing(shape)){
      pp<-pp+geom_point(aes(shape=as.name(shape)),alpha=alpha,size=size)
    }else{
      
      pp<-pp+geom_point(alpha=alpha,size=size)
    }
  }
  pp+xlab(paste("PC2 (",varExp[2],"%)"))+
    ylab(paste("PC3 (",varExp[3],"%)"))
}


pcaFunction13<-function(mat, inf, genes, ngenes, alpha, col="treatment", label="sid",shape ) { 
  
  if(!missing(genes)){
    vsdMat <- mat[genes,rownames(inf)]
  }else{
    vsdMat<- mat[,rownames(inf)]
  }
  rowVar<-apply(mat,1,var)
  
  if(!missing(ngenes)){
    mv500<-order(rowVar,decreasing = T)[1:ngenes]
  }else{
    mv500<-(rowVar>0)
  }
  if(missing(shape)){
    shape=col
  }
  
  pc<-prcomp(t(vsdMat[mv500,]),scale=T)
  pc.sum<-summary(pc)$importance
  varExp<-round(pc.sum[2,]*100,2)
  pcs<-data.frame(pc$x, inf)
  pcs$sid<-rownames(pcs)
  pp<-ggplot(pcs,aes_q(x=quote(PC1),y=quote(PC3),col=as.name(col)))
  if(!missing(label)){
    pp<-pp+geom_text(aes(label=as.name(label)),size=3)
  }else{
    if(missing(shape)){
      pp<-pp+geom_point(aes(shape=as.name(shape)),alpha=alpha,size=3)
    }else{
      
      pp<-pp+geom_point(alpha=alpha,size=4)
    }
  }
  pp+xlab(paste("PC1 (",varExp[1],"%)"))+
    ylab(paste("PC3 (",varExp[3],"%)"))
}

## PCA function form Johannes

plotPCA2<-function(object, intgroup="condition", ntop=numgenes, returnData=FALSE)
{
  library(matrixStats)
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3],PC4=pca$x[,4],group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    d<-list(PCA_df=d, percentVar= percentVar[1:4])
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed()
}


## boxplot of normalized counts to evaluate success of normalisation

normBoxplot<-function(mm,title){
  df<- data.frame(mm) %>% gather(key=sample_name,value="raw",1:dim(mm)[2])
  box<-ggplot(df,aes(x=sample_name,y=raw))+
    geom_boxplot()+ylab(expression(log[2](counts)))+ 
    coord_flip()+
    theme_grey()+
    theme(axis.title.y=element_blank())
  if(!missing(title)){
    box+ggtitle(title)
  }else{
    box
  }
}

## Gene filtering
#Function for selecting genes expressed in at least x% of samples
whichgenes_reproducible <- function(exprtable,exprcutoff,reproducecutoff){
  expressedgenes <- row.names(exprtable[which(rowSums(exprtable)>=exprcutoff),])
  tmp <- exprtable[expressedgenes,]
  ndetections <- rowSums(exprtable[expressedgenes,]>=1)
  outgenes <- names(ndetections[which(ndetections/ncol(exprtable)>=reproducecutoff)])
  return(outgenes)
}

## Mean Gene Expression for filtering
mean_gene_exp <- function(data, genes, UMI)
{

if (UMI){
  data <- data$umicount$inex$all
  data <- data[genes, ]
  data <- as.matrix(data)
  data_10 <- data[whichgenes_reproducible(data,1,0.10), ]
  data_25 <- data[whichgenes_reproducible(data,1,0.25), ]

  plot_mean <- ggplot()+
  geom_density(aes(x= rowMeans(data), fill="Unfiltered"), alpha=0.7)+
  geom_density(aes(x= rowMeans(data_10), fill="10% Filtered"), alpha=0.7)+
  geom_density(aes(x= rowMeans(data_25), fill="25% Filtered"), alpha=0.7)+
  scale_fill_manual(values = c("#003F5A","#DE6600", "#696464"))+
  scale_x_log10()+
  xlab("Mean Gene Expression (log10)")+
  theme_bw() + 
  theme(
     plot.title = element_text(hjust = 0.5, size=18, face="bold"),
     axis.text = element_text(colour="black", size=14), 
     axis.title=element_text(size=16,face="bold"), 
     legend.text=element_text(size=14),
     legend.title = element_blank(),
     legend.position="right",
     axis.line.x = element_line(colour = "black"), 
     axis.line.y = element_line(colour = "black"),
     strip.background=element_blank(), 
     strip.text=element_text(size=16))
} else{
  data <- data$readcount$inex$all
  data <- data[genes, ]
  data <- as.matrix(data)
  data_10 <- data[whichgenes_reproducible(data,1,0.10), ]
  data_25 <- data[whichgenes_reproducible(data,1,0.25), ]

  plot_mean <- ggplot()+
  geom_density(aes(x= rowMeans(data), fill="Unfiltered"), alpha=0.7)+
  geom_density(aes(x= rowMeans(data_10), fill="10% Filtered"), alpha=0.7)+
  geom_density(aes(x= rowMeans(data_25), fill="25% Filtered"), alpha=0.7)+
  scale_fill_manual(values = c("#003F5A","#DE6600", "#696464"))+
  scale_x_log10()+
  xlab("Mean Gene Expression (log10)")+
  theme_bw() + 
  theme(
     plot.title = element_text(hjust = 0.5, size=18, face="bold"),
     axis.text = element_text(colour="black", size=14), 
     axis.title=element_text(size=16,face="bold"), 
     legend.text=element_text(size=14),
     legend.title = element_blank(),
     legend.position="right",
     axis.line.x = element_line(colour = "black"), 
     axis.line.y = element_line(colour = "black"),
     strip.background=element_blank(), 
     strip.text=element_text(size=16))
    }
}

## same for exonic reads
mean_gene_exp_ex <- function(data, genes, UMI)
{
  
  if (UMI){
    data <- data$umicount$exon$all
    data <- data[genes, ]
    data <- as.matrix(data)
    data_10 <- data[whichgenes_reproducible(data,1,0.10), ]
    data_25 <- data[whichgenes_reproducible(data,1,0.25), ]
    
    plot_mean <- ggplot()+
      geom_density(aes(x= rowMeans(data), fill="Unfiltered"), alpha=0.7)+
      geom_density(aes(x= rowMeans(data_10), fill="10% Filtered"), alpha=0.7)+
      geom_density(aes(x= rowMeans(data_25), fill="25% Filtered"), alpha=0.7)+
      scale_fill_manual(values = c("#003F5A","#DE6600", "#696464"))+
      scale_x_log10()+
      xlab("Mean Gene Expression (log10)")+
      theme_bw() + 
      theme(
        plot.title = element_text(hjust = 0.5, size=18, face="bold"),
        axis.text = element_text(colour="black", size=14), 
        axis.title=element_text(size=16,face="bold"), 
        legend.text=element_text(size=14),
        legend.title = element_blank(),
        legend.position="right",
        axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black"),
        strip.background=element_blank(), 
        strip.text=element_text(size=16))
  } else{
    data <- data$readcount$exon$all
    data <- data[genes, ]
    data <- as.matrix(data)
    data_10 <- data[whichgenes_reproducible(data,1,0.10), ]
    data_25 <- data[whichgenes_reproducible(data,1,0.25), ]
    
    plot_mean <- ggplot()+
      geom_density(aes(x= rowMeans(data), fill="Unfiltered"), alpha=0.7)+
      geom_density(aes(x= rowMeans(data_10), fill="10% Filtered"), alpha=0.7)+
      geom_density(aes(x= rowMeans(data_25), fill="25% Filtered"), alpha=0.7)+
      scale_fill_manual(values = c("#003F5A","#DE6600", "#696464"))+
      scale_x_log10()+
      xlab("Mean Gene Expression (log10)")+
      theme_bw() + 
      theme(
        plot.title = element_text(hjust = 0.5, size=18, face="bold"),
        axis.text = element_text(colour="black", size=14), 
        axis.title=element_text(size=16,face="bold"), 
        legend.text=element_text(size=14),
        legend.title = element_blank(),
        legend.position="right",
        axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black"),
        strip.background=element_blank(), 
        strip.text=element_text(size=16))
  }
}


## Get CV, p0, mean expression 
generate_stats <- function(exprdf,conds){
  require(matrixStats)
  if(length(conds)!=ncol(exprdf)){
    print("Erorr in Conditions given")
    stop()
  }
  conds_vec <- unique(conds)
  cvfrommatrix <- function(exprmatrix){
    exprtmp <- as.matrix(exprmatrix)
    cv <- (matrixStats::rowSds(exprtmp)/rowMeans(exprtmp))
    return(cv)
  }
  calcp0 <- function(exprtable){
    nsamples = dim(exprtable)[2]
    nn0 = rowSums(!exprtable)
    p0 = 1 - ((nsamples - nn0)/nsamples)
    return(p0)
  }
  output_df_collect <- data.frame(GeneID=NA,mu=NA,cv=NA,p0=NA,cond=NA,stringsAsFactors = F)
  for(i in conds_vec){
    tmp_expr <- exprdf[,grep(i,conds)]
    tmp_expr <- tmp_expr[which(rowSums(tmp_expr)>0),]
    output_df_tmp <- data.frame(GeneID=row.names(tmp_expr),
                                mu=rowMeans(tmp_expr),
                                cv=cvfrommatrix(tmp_expr),
                                p0=calcp0(tmp_expr),
                                cond=rep(i,nrow(tmp_expr)),
                                stringsAsFactors = F)
    output_df_collect <- rbind.data.frame(output_df_collect,output_df_tmp) 
  }
  output_df_collect <- output_df_collect[-1,]
  output_df_collect$EPV <- output_df_collect$cv - (sqrt(output_df_collect$mu)/output_df_collect$mu)
  output_df_collect$poisson_cv <- (sqrt(output_df_collect$mu)/output_df_collect$mu)
  return(output_df_collect)
}

## remove Gencode Gene name extension

remove_Geneversion<-function(countmatrix){
rownames(countmatrix) <-str_split_fixed(rownames(countmatrix), pattern="[.]",n=2)[,1]

if (length(unique((rownames(countmatrix))))!=length(rownames(countmatrix))){
  
  countmatrix<-countmatrix %>% 
    as.data.frame() %>% 
    rownames_to_column(var="Gene_ID") %>% 
    mutate(Gene_ID=str_split_fixed(Gene_ID, pattern="[.]",n=2)[,1]) %>% 
    group_by(Gene_ID) %>% 
    dplyr::summarize(across(where(is.double),sum)) %>% 
    ungroup() %>% 
    column_to_rownames(var="Gene_ID") %>% 
    as.matrix()
  return(countmatrix) 
}else {
 return(countmatrix) 
}

}

## Get Biotype

getbiotype <- function(dataset,species){
  require(biomaRt)
  
  if(species=="human"){
  ensembl <- useMart("ensembl", dataset =dataset, host="uswest.ensembl.org") #use west coast servers ince european server was down
  mito <- getBM(attributes=c("ensembl_gene_id","ensembl_gene_id_version"),
                filter = "chromosome_name",
                values = "MT",
                mart=ensembl )
  rrna <- getBM(attributes=c("ensembl_gene_id","ensembl_gene_id_version"),
                filter="biotype",
                values = "rRNA",
                mart=ensembl )
  lnc <- getBM(attributes=c("ensembl_gene_id","ensembl_gene_id_version"),
               filter = "biotype",
               values = "lncRNA",
               mart=ensembl )
  
  tmp <- bind_rows( data.frame(type = "rrna",
                               ENSEMBL=rrna$ensembl_gene_id,
                               Gencode=rrna$ensembl_gene_id_version),
                    data.frame(type = "mito",
                               ENSEMBL=mito$ensembl_gene_id,
                               Gencode=mito$ensembl_gene_id_version),
                    data.frame(type = "lnc",
                               ENSEMBL=lnc$ensembl_gene_id,
                               Gencode=lnc$ensembl_gene_id_version) )
  rm(ensembl)
  return(tmp)
}

if (species=="mouse"){
  
  ensembl <- useMart("ensembl", dataset =dataset, host="uswest.ensembl.org") #use west coast servers ince european server was down
  mito <- getBM(attributes=c("ensembl_gene_id","ensembl_gene_id_version"),
                filter = "chromosome_name",
                values = "mt",
                mart=ensembl )
  rrna <- getBM(attributes=c("ensembl_gene_id","ensembl_gene_id_version"),
                filter="biotype",
                values = "rRNA",
                mart=ensembl )
  lnc <- getBM(attributes=c("ensembl_gene_id","ensembl_gene_id_version"),
               filter = "biotype",
               values = "lncRNA",
               mart=ensembl )
  
  tmp <- bind_rows( data.frame(type = "rrna",
                               ENSEMBL=rrna$ensembl_gene_id,
                               Gencode=rrna$ensembl_gene_id_version),
                    data.frame(type = "mito",
                               ENSEMBL=mito$ensembl_gene_id,
                               Gencode=mito$ensembl_gene_id_version),
                    data.frame(type = "lnc",
                               ENSEMBL=lnc$ensembl_gene_id,
                               Gencode=lnc$ensembl_gene_id_version) )
  rm(ensembl)
  return(tmp)
}
}

## Collapse Downsampled Counts

collapse_downsampled_counts<-function(zumismat,type,frac.samples,genes,umi){
  if(missing(genes)){
    genes<-rownames(zumismat$umicount$inex$all)
  }
  
  if(umi){
  
  if (type=="exon")
  {
    dsDF <- data.frame() #initialise output
    for(i in 1:length(zumismat$umicount$exon$downsampling)){
      tmp <- as.matrix(zumismat$umicount$exon$downsampling[[i]])
      if (ncol(zumismat$umicount$exon$downsampling[[i]])>1)
        {
        tmp<-tmp[whichgenes_reproducible(tmp,1,reproducecutoff=frac.samples),]
        tmp<-tmp[rownames(tmp)%in%genes,]
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$umicount$exon$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=colSums(tmp),
                             Genes=colSums(tmp>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F)}
      else 
        {tmp2<-as.matrix(tmp[rownames(tmp)%in%genes,])
      tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$umicount$exon$downsampling)[i],
                                                                pattern="downsampled_"))),
                           UMIs=sum(tmp2),
                           Genes=sum(tmp2>0),
                           XC=colnames(tmp),
                           stringsAsFactors = F)} #calculate Genes/UMIs detected
      dsDF <- rbind.data.frame(dsDF,tmp_df) #collect output
    }
    dsDF$type <-"exon"
    return(dsDF)
  }
  else if (type=="intron"){
    dsDF <- data.frame() #initialise output
    for(i in 1:length(zumismat$umicount$intron$downsampling)){ #for each downsampling depth, get the UMI counts and Gene counts
      tmp <- as.matrix(zumismat$umicount$intron$downsampling[[i]])
      if (ncol(zumismat$umicount$exon$downsampling[[i]])>1){##Do filtering only if there's more than one sample/BC left
        tmp<-tmp[whichgenes_reproducible(tmp,1,reproducecutoff=frac.samples),]
        tmp<-tmp[rownames(tmp)%in%genes,]
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$umicount$intron$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=colSums(tmp),
                             Genes=colSums(tmp>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F)}
      else{tmp2<-as.matrix(tmp[rownames(tmp)%in%genes,])
      tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$umicount$intron$downsampling)[i],
                                                                pattern="downsampled_"))),
                           UMIs=sum(tmp2),
                           Genes=sum(tmp2>0),
                           XC=colnames(tmp),
                           stringsAsFactors = F)}
      
      dsDF <- rbind.data.frame(dsDF,tmp_df) #collect output
    }
    dsDF$type <-"intron"
    return(dsDF)
  }
  else if (type=="inex"){
    dsDF <- data.frame() #initialise output
    for(i in 1:length(zumismat$umicount$inex$downsampling)){ #for each downsampling depth, get the UMI counts and Gene counts
      tmp <- as.matrix(zumismat$umicount$inex$downsampling[[i]])
      if (ncol(zumismat$umicount$exon$downsampling[[i]])>1){##Do filtering only if there's more than one sample/BC left
        tmp<-tmp[whichgenes_reproducible(tmp,1,reproducecutoff=frac.samples),]
        tmp<-tmp[rownames(tmp)%in%genes,]
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$umicount$inex$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=colSums(tmp),
                             Genes=colSums(tmp>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F) }
        else{tmp2<-as.matrix(tmp[rownames(tmp)%in%genes,])
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$umicount$inex$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=sum(tmp2),
                             Genes=sum(tmp2>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F)} 
      dsDF <- rbind.data.frame(dsDF,tmp_df) #collect output
    }
    dsDF$type <-"inex"
    return(dsDF)
  }
  else {
    dsDF_ex <- data.frame() #initialise output
    for(i in 1:length(zumismat$umicount$exon$downsampling)){ #for each downsampling depth, get the UMI counts and Gene counts
      tmp <- as.matrix(zumismat$umicount$exon$downsampling[[i]])
      if (ncol(zumismat$umicount$exon$downsampling[[i]])>1){##Do filtering only if there's more than one sample/BC left
        tmp<-tmp[whichgenes_reproducible(tmp,1,reproducecutoff=frac.samples),]
        tmp<-tmp[rownames(tmp)%in%genes,]
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$umicount$exon$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=colSums(tmp),
                             Genes=colSums(tmp>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F)}
        else{tmp2<-as.matrix(tmp[rownames(tmp)%in%genes,])
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$umicount$exon$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=sum(tmp2),
                             Genes=sum(tmp2>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F)}
      dsDF_ex <- rbind.data.frame(dsDF_ex,tmp_df) #collect output
    }
    dsDF_ex$type <-"exon"
    dsDF_in <- data.frame() #initialise output
    for(i in 1:length(zumismat$umicount$intron$downsampling)){ #for each downsampling depth, get the UMI counts and Gene counts
      tmp <- as.matrix(zumismat$umicount$intron$downsampling[[i]])
      if (ncol(zumismat$umicount$intron$downsampling[[i]])>1){##Do filtering only if there's more than one sample/BC left
        tmp<-tmp[whichgenes_reproducible(tmp,1,reproducecutoff=frac.samples),]
        tmp<-tmp[rownames(tmp)%in%genes,]
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$umicount$intron$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=colSums(tmp),
                             Genes=colSums(tmp>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F) }
        else{tmp2<-as.matrix(tmp[rownames(tmp)%in%genes,])
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$umicount$intron$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=sum(tmp2),
                             Genes=sum(tmp2>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F)}
      dsDF_in <- rbind.data.frame(dsDF_in,tmp_df) #collect output
    } 
    dsDF_in$type <-"intron"
    dsDF_inex <- data.frame() #initialise output
    for(i in 1:length(zumismat$umicount$inex$downsampling)){ #for each downsampling depth, get the UMI counts and Gene counts
      tmp <- as.matrix(zumismat$umicount$inex$downsampling[[i]])
      if (ncol(zumismat$umicount$inex$downsampling[[i]])>1){##Do filtering only if there's more than one sample/BC left
        tmp<-tmp[whichgenes_reproducible(tmp,1,reproducecutoff=frac.samples),]
        tmp<-tmp[rownames(tmp)%in%genes,]
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$umicount$inex$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=colSums(tmp),
                             Genes=colSums(tmp>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F) }
        else{tmp2<-as.matrix(tmp[rownames(tmp)%in%genes,])
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$umicount$inex$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=sum(tmp2),
                             Genes=sum(tmp2>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F)}
      dsDF_inex <- rbind.data.frame(dsDF_inex,tmp_df) #collect output
    }
    dsDF_inex$type <-"inex"
    dsDF_all<-rbind.data.frame(dsDF_ex,dsDF_in,dsDF_inex)
    return(dsDF_all)
  }
  }
  else{
  if (type=="exon")
  {
    dsDF <- data.frame() #initialise output
    for(i in 1:length(zumismat$readcount$exon$downsampling)){
      tmp <- as.matrix(zumismat$readcount$exon$downsampling[[i]])
      if (ncol(zumismat$readcount$exon$downsampling[[i]])>1)
      {
        tmp<-tmp[whichgenes_reproducible(tmp,1,reproducecutoff=frac.samples),]
        tmp<-tmp[rownames(tmp)%in%genes,]
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$readcount$exon$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=colSums(tmp),
                             Genes=colSums(tmp>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F)}
      else 
      {tmp2<-as.matrix(tmp[rownames(tmp)%in%genes,])
      tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$readcount$exon$downsampling)[i],
                                                                pattern="downsampled_"))),
                           UMIs=sum(tmp2),
                           Genes=sum(tmp2>0),
                           XC=colnames(tmp),
                           stringsAsFactors = F)} #calculate Genes/UMIs detected
      dsDF <- rbind.data.frame(dsDF,tmp_df) #collect output
    }
    dsDF$type <-"exon"
    return(dsDF)
  }
  else if (type=="intron"){
    dsDF <- data.frame() #initialise output
    for(i in 1:length(zumismat$readcount$intron$downsampling)){ #for each downsampling depth, get the UMI counts and Gene counts
      tmp <- as.matrix(zumismat$readcount$intron$downsampling[[i]])
      if (ncol(zumismat$readcount$exon$downsampling[[i]])>1){##Do filtering only if there's more than one sample/BC left
        tmp<-tmp[whichgenes_reproducible(tmp,1,reproducecutoff=frac.samples),]
        tmp<-tmp[rownames(tmp)%in%genes,]
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$readcount$intron$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=colSums(tmp),
                             Genes=colSums(tmp>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F)}
      else{tmp2<-as.matrix(tmp[rownames(tmp)%in%genes,])
      tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$readcount$intron$downsampling)[i],
                                                                pattern="downsampled_"))),
                           UMIs=sum(tmp2),
                           Genes=sum(tmp2>0),
                           XC=colnames(tmp),
                           stringsAsFactors = F)}
      
      dsDF <- rbind.data.frame(dsDF,tmp_df) #collect output
    }
    dsDF$type <-"intron"
    return(dsDF)
  }
  else if (type=="inex"){
    dsDF <- data.frame() #initialise output
    for(i in 1:length(zumismat$readcount$inex$downsampling)){ #for each downsampling depth, get the UMI counts and Gene counts
      tmp <- as.matrix(zumismat$readcount$inex$downsampling[[i]])
      if (ncol(zumismat$readcount$exon$downsampling[[i]])>1){##Do filtering only if there's more than one sample/BC left
        tmp<-tmp[whichgenes_reproducible(tmp,1,reproducecutoff=frac.samples),]
        tmp<-tmp[rownames(tmp)%in%genes,]
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$readcount$inex$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=colSums(tmp),
                             Genes=colSums(tmp>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F) }
      else{tmp2<-as.matrix(tmp[rownames(tmp)%in%genes,])
      tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$readcount$inex$downsampling)[i],
                                                                pattern="downsampled_"))),
                           UMIs=sum(tmp2),
                           Genes=sum(tmp2>0),
                           XC=colnames(tmp),
                           stringsAsFactors = F)} 
      dsDF <- rbind.data.frame(dsDF,tmp_df) #collect output
    }
    dsDF$type <-"inex"
    return(dsDF)
  }
  else {
    dsDF_ex <- data.frame() #initialise output
    for(i in 1:length(zumismat$readcount$exon$downsampling)){ #for each downsampling depth, get the UMI counts and Gene counts
      tmp <- as.matrix(zumismat$readcount$exon$downsampling[[i]])
      if (ncol(zumismat$readcount$exon$downsampling[[i]])>1){##Do filtering only if there's more than one sample/BC left
        tmp<-tmp[whichgenes_reproducible(tmp,1,reproducecutoff=frac.samples),]
        tmp<-tmp[rownames(tmp)%in%genes,]
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$readcount$exon$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=colSums(tmp),
                             Genes=colSums(tmp>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F)}
      else{tmp2<-as.matrix(tmp[rownames(tmp)%in%genes,])
      tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$readcount$exon$downsampling)[i],
                                                                pattern="downsampled_"))),
                           UMIs=sum(tmp2),
                           Genes=sum(tmp2>0),
                           XC=colnames(tmp),
                           stringsAsFactors = F)}
      dsDF_ex <- rbind.data.frame(dsDF_ex,tmp_df) #collect output
    }
    dsDF_ex$type <-"exon"
    dsDF_in <- data.frame() #initialise output
    for(i in 1:length(zumismat$readcount$intron$downsampling)){ #for each downsampling depth, get the UMI counts and Gene counts
      tmp <- as.matrix(zumismat$readcount$intron$downsampling[[i]])
      if (ncol(zumismat$readcount$intron$downsampling[[i]])>1){##Do filtering only if there's more than one sample/BC left
        tmp<-tmp[whichgenes_reproducible(tmp,1,reproducecutoff=frac.samples),]
        tmp<-tmp[rownames(tmp)%in%genes,]
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$readcount$intron$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=colSums(tmp),
                             Genes=colSums(tmp>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F) }
      else{tmp2<-as.matrix(tmp[rownames(tmp)%in%genes,])
      tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$readcount$intron$downsampling)[i],
                                                                pattern="downsampled_"))),
                           UMIs=sum(tmp2),
                           Genes=sum(tmp2>0),
                           XC=colnames(tmp),
                           stringsAsFactors = F)}
      dsDF_in <- rbind.data.frame(dsDF_in,tmp_df) #collect output
    } 
    dsDF_in$type <-"intron"
    dsDF_inex <- data.frame() #initialise output
    for(i in 1:length(zumismat$readcount$inex$downsampling)){ #for each downsampling depth, get the UMI counts and Gene counts
      tmp <- as.matrix(zumismat$readcount$inex$downsampling[[i]])
      if (ncol(zumismat$readcount$inex$downsampling[[i]])>1){##Do filtering only if there's more than one sample/BC left
        tmp<-tmp[whichgenes_reproducible(tmp,1,reproducecutoff=frac.samples),]
        tmp<-tmp[rownames(tmp)%in%genes,]
        tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$readcount$inex$downsampling)[i],
                                                                  pattern="downsampled_"))),
                             UMIs=colSums(tmp),
                             Genes=colSums(tmp>0),
                             XC=colnames(tmp),
                             stringsAsFactors = F) }
      else{tmp2<-as.matrix(tmp[rownames(tmp)%in%genes,])
      tmp_df <- data.frame(depth=as.integer(str_trim(str_remove(string =names(zumismat$readcount$inex$downsampling)[i],
                                                                pattern="downsampled_"))),
                           UMIs=sum(tmp2),
                           Genes=sum(tmp2>0),
                           XC=colnames(tmp),
                           stringsAsFactors = F)}
      dsDF_inex <- rbind.data.frame(dsDF_inex,tmp_df) #collect output
    }
    dsDF_inex$type <-"inex"
    dsDF_all<-rbind.data.frame(dsDF_ex,dsDF_in,dsDF_inex)
    return(dsDF_all)
  }
}
}

#### PowsimR helper functions to extract results for plotting

getEvalDE_df_marginal<-function(evalres,method=c("prime-seq","tru-seq"),count_type=c("exon","inex")){
  dat.marginal <- evalres[grep('*R.marginal', names(evalres))]
  
  names(dat.marginal) <- substr(x = names(dat.marginal), start = 1, stop = 3)
  dat.marginal <- lapply(dat.marginal, "rownames<-", paste0(evalres[['n1']], " vs ", evalres[['n2']]))
  dat.marginal.long <- reshape2::melt(dat.marginal)
  dat.marginal.long$L1 <- factor(dat.marginal.long$L1,
                                 levels = c("TPR", "FNR", "FPR", "TNR", "FDR"))
  
  dat.marginal.long$method<-method
  
  dat.marginal.long$count_type<-count_type
  return(dat.marginal.long)
}


## conditional rates
getEvalDE_list_conditional<-function(evalRes,method=c("prime-seq","tru-seq"),count_type=c("exon","inex")){
  plot_list<-list(dat.genes.calc=NULL,dat.stratified.long=NULL,stratum.name=NULL)
  plot_list$stratum.name <- dplyr::case_when(evalRes$stratify.by == "mean" ~ "(Log2 Mean Expression)",
                                   evalRes$stratify.by == "dispersion" ~ "(Log2 Dispersion)",
                                   evalRes$stratify.by == "lfc" ~ "(Log2 Fold Change)",
                                   evalRes$stratify.by == "lfc_abs" ~ "absolute (Log2 Fold Change)",
                                   evalRes$stratify.by == "dropout" ~ "(Gene Dropout Rate)")
  # strata genes
  strata <- evalRes$strata.levels
  N <- length(evalRes$n1)
  dat.genes <- list("Ngenes"=evalRes$stratagenes[,N,],
                    'DEgenes'=evalRes$stratadiffgenes[,N,])
  dat.genes <- lapply(dat.genes, "rownames<-", strata)
  dat.genes.long <- reshape2::melt(dat.genes)
  dat.genes.calc <- dat.genes.long %>%
    dplyr::group_by(.data$Var1, .data$L1) %>%
    dplyr::summarise(Expectation=mean(.data$value),
                     Deviation=stats::sd(.data$value),
                     Error=stats::sd(.data$value)/sqrt(dplyr::n())) %>%
    dplyr::ungroup()
  dat.genes.calc$method<-method
  dat.genes.calc$count_type<-count_type
  plot_list$dat.genes.calc<-dat.genes.calc
  refval <- data.frame(L1 = c("FDR", "TPR"),
                       ref = c(evalRes$alpha.nominal, 0.8))
  
  
  dat.stratified <- evalRes[grep('*R$', names(evalRes))]
  
  dat.stratified <- lapply(dat.stratified, "dimnames<-",
                           list(strata, paste0(evalRes[['n1']], " vs ", evalRes[['n2']]),
                                NULL))
  dat.stratified.long <- reshape2::melt(dat.stratified)
  dat.stratified.long$L1 <- factor(dat.stratified.long$L1,
                                   levels = c("TPR", "FNR", "FPR", "TNR", "FDR"))
  
  dat.stratified.long$method<-method
  
  dat.stratified.long$count_type<-count_type
  
  plot_list$dat.stratified.long<-dat.stratified.long
  return(plot_list)
}
