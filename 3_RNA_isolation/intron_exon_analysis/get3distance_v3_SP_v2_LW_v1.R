#!/usr/bin/env Rscript
require(optparse)

########### READ options #######
option_list = list(
  make_option(c("--bamfile"), type="character", default=NULL, 
              help="bamfile file name, note you also need a bai", metavar="character"),
  make_option(c("--gtf"), type="character", default="nf", 
              help="gtf file containing the transcript info", metavar="character"),
  make_option(c("--txdb"), type="character", default="nf", 
              help="txdb created from gtf", metavar="character"),
  make_option(c("--genes"), type="character", default="all", 
              help="text file with one ENSG per line", metavar="character"),
  make_option(c("--maxTdist"), type="integer", default=10000, 
              help="max dist. to 3' end counted", metavar="integer"),
  make_option(c("--out"), type="character", default="out", 
              help="name base for plot and outtable", metavar="character"),
  make_option(c("--g2tmapping"), type="character", default="nf", 
              help="RData file with a list of gene to transcript mappings", metavar="character"),
  make_option(c("--ftype"), type="character", default="exon", 
              help="intron or exon", metavar="character"),
  make_option(c("--minTlength"), type="integer", default=1, 
              help="minimal length of transcripts considered", metavar="integer"),
  make_option(c("--canonical"), type="logical", default=T, 
              help="if T only canonical chromosomes will be used", metavar="logical"),
  make_option(c("--plot"), type="logical", default=T, 
              help="if T plots will be made", metavar="logical"),
  make_option(c("--max_exon_length"), type="integer", default=2000, 
              help="total length of coding sequence", metavar="integer"),
  make_option(c("--ngenes"), type="integer", default=2000, 
              help="maximum number of genes to be used", metavar="integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#opt<-list(bamfile= paste0(bam_path, i),gtf="nf",txdb="/data/share/htp/prime-seq_Paper/genomes/Hsap/gencode.v35.primary_assembly.sqlite",genes="all",maxTdist=maxTdist,out=paste0(out,"_", j),g2tmapping=g2t,ftype="exon",minTlength=minTlength,canonical=T,plot=F,max_exon_length=5000,ngenes=5000)

#################################################
suppressMessages(library(GenomicRanges, quietly=T, warn.conflicts = F))
suppressMessages(library(GenomicFeatures, quietly=T, warn.conflicts = F))
suppressMessages(library(GenomicAlignments, quietly=T, warn.conflicts = F))
suppressMessages(library(ggplot2, quietly=T, warn.conflicts = F))
suppressMessages(library(parallel, quietly=T, warn.conflicts = F))

# Function to get distances
get3distanceTx<-function(tx, bamfile,gid,mtl=1){
  coo<-reduce(unlist(tx))
  param<-ScanBamParam(which=coo)
  grl.ali<-readGAlignments(bamfile, param=param)
  if(length(grl.ali)>0){ 
    grl<- GRanges( seqnames= seqnames(grl.ali), ranges = ranges(grl.ali))
    seqlevels(grl,pruning.mode="coarse")<-seqlevels(tx)
    nr <- length(grl) 
    tpos <- mclapply(1:length(tx),function(i){
      tlength<-sum(width(tx[[i]]))
      if( tlength>mtl ){
        tmap   <- start(mapToTranscripts(grl, tx[i] ))
        tmap   <- tlength-tmap[tmap>0]
        tmap
      }else{
        0 
      }
    })
    ## correct for reads counted for multiple transcripts
    double<-nr/sum(unlist(tpos)>0)
    table(unlist(tpos))*double
  }
}
##########################################
#Load database

if(file.exists(opt$txdb)){
  print(paste( "loading", opt$txdb) )
  txdb<-loadDb(opt$txdb)
}else if( file.exists(opt$gtf)){
  print(paste("loading",opt$gtf))
  txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")
}else{
  return("NEED TRANSCRIPT ANNOTATION")
}

if(opt$canonical){
  canonical_chr<-seqlevels(txdb)[1:25]
  seqlevels(txdb)<-canonical_chr
}

seqlevels(txdb)

if(file.exists(opt$g2tmapping) ){
  load(opt$g2tmapping)
}else{ 
  tmp<-transcriptsBy(txdb, by="gene")
  gene2transcript<-lapply( tmp, function(x){  mcols(x)$tx_name })
}

geneExon<-exonsBy(txdb, by="gene")
ebglen <- sum(width(geneExon))

selectedgenes <- ebglen[ebglen<=opt$max_exon_length]
if(opt$ftype=="exon"){
  transcr<-exonsBy(txdb, by="tx", use.names=TRUE)
}else{
  transcr<-intronsByTranscript(txdb, use.names=TRUE)
}

print("Loaded all the data!")
######################################################
### Read Genes
if( file.exists(opt$genes)){
  genes<-scan(opt$genes,what="character")
if ( length(genes)> 1000 & length(genes)< opt$ngenes){
  opt$ngenes<-length(genes)
}
}else{
  genes<-names(gene2transcript)
}
print( paste("Found",length(selectedgenes),"Genes. Analyzing",opt$ngenes, "Genes"))
################################################

positionCount<-rep(0,opt$maxTdist)
names(positionCount)<-1:opt$maxTdist

print(opt$ftype)
if (length(selectedgenes)>1000) {
  for ( gid in names(selectedgenes)[sample(1:length(selectedgenes),size = opt$ngenes)] ){
    #print(gid="ENSG00000232482.2")
    ntx<-gene2transcript[[ gid ]]
    dd<- get3distanceTx(tx=transcr[ ntx ], bamfile=opt$bamfile,gid=gid,mtl=opt$minTlength)
    
    for( i in names(dd) ){
      positionCount[i]<-positionCount[i]+dd[i]
    }
  }
}else {
  for ( gid in names(selectedgenes)){
    print(gid)
    ntx<-gene2transcript[[ gid ]]
    dd<- get3distanceTx(tx=transcr[ ntx ], bamfile=opt$bamfile,gid=gid,mtl=opt$minTlength)
    
    for( i in names(dd) ){
      positionCount[i]<-positionCount[i]+dd[i]
    }
  }
}

df<-data.frame(  distance_to_3_end = as.numeric(names(positionCount)),
                 read_number       = positionCount )


write.table(df, file=paste(opt$out,"txt",sep="."), quote=F, sep="\t", col.names=T,row.names=F )

###### Plotting ########
### BARPLOT 
barplot=F
if(barplot){
  df.bin<-df %>% group_by( round(distance_to_3_end/10)  ) %>%
    summarise( distance_to_3_end= mean(distance_to_3_end),
               read_number=sum(read_number))
  
  ggplot( df.bin[df.bin$read_number>0  , ], aes(x=distance_to_3_end, y=read_number) ) +
    geom_bar(stat = "identity") + scale_y_log10() +
    scale_x_reverse() +
    theme_minimal()
}
#### SMOOTH PLOT


p<-ggplot(df[df$read_number>0,], aes(x=distance_to_3_end, y=read_number)) + 
  stat_smooth()+
  #scale_y_log10() +
  scale_x_reverse() +
  #ggtitle(opt$bamfile)+
  theme_minimal()
ggsave(p,filename = paste(opt$out,"pdf",sep="."),width=7,height=6)

q()
