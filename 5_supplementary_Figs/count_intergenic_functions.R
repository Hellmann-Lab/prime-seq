## intergenic counting custom functions modified from zUMIs version 2.4.5b
## Featurecounts
.runFeatureCount<-function(abamfile,RG,saf,strand,type,primaryOnly,cpu,mem,outdir){
  print(paste0("Assigning reads to features (",type,")"))
  fc.stat<-Rsubread::featureCounts(files=abamfile,
                                   annot.ext=saf,
                                   isGTFAnnotationFile=F,
                                   primaryOnly=primaryOnly,
                                   countMultiMappingReads=primaryOnly,
                                   nthreads=cpu,
                                   reportReads="BAM",
                                   reportReadsPath = outdir,
                                   strandSpecific=strand,
                                   isPairedEnd=F,
                                   countChimericFragments=F)$stat
  fn<-paste0(abamfile,".featureCounts.bam")
  nfn<-paste0(abamfile,".",type,".featureCounts.bam")
  
  system(paste0("mv ",fn," ",nfn,".tmp"))
  
  invisible(suppressWarnings(suppressMessages(gc(verbose=F))))
  return(nfn)
}

splitRG<-function(bccount,mem,read_layout){
  if(is.null(mem) || mem==0){
    maxR<- Inf
  }else{
    maxR<- floor( mem*1000 * 4500 )
  }
  if( (maxR > 2e+09 & read_layout == "SE") | (maxR > 1e+09 & read_layout == "PE") ){
    maxR <- ifelse(read_layout == "SE",2e+09,1e+09)
  } 
  print(paste(maxR,"Reads per chunk"))
  nc<-nrow(bccount)
  cs=0
  chunkID=1
  bccount[,chunkID:=0]
  for(i in 1:nc){
    cs=cs+bccount[i]$n
    if(bccount[i]$n>maxR){
      print(paste("Warning: Barcode",bccount[i]$XC,"has more reads than allowed for the memory limit!
                  Proceeding anyway..."))
    }
    if(cs>=maxR){
      chunkID=chunkID+1
      cs=bccount[i][,"n"]
    }
    bccount[i][,"chunkID"]=chunkID
  }
  return(bccount)
}

BCbin <- function(bccount_file, bc_detected,nReadsperCell,nthread,BarcodeBinning) {
  true_BCs <- bc_detected[,XC]
  nocell_bccount<-data.table::fread( bccount_file, col.names = c("XC","n"))[
    ,list(n=sum(n)),by=XC][
      n>=nReadsperCell][
        order(-n)][
          !( XC %in% true_BCs )   ]
  nocell_BCs <- nocell_bccount[,XC]
  
  #break up in pieces of 1000 real BCs in case the hamming distance calculation gets too large!
  true_chunks <- split(true_BCs, ceiling(seq_along(true_BCs)/1000))
  for(i in 1:length(true_chunks)){
    dists <- stringdist::stringdistmatrix(true_chunks[[i]],nocell_BCs,method="hamming", nthread = nthread)
    dists <- setDT(data.frame(dists))
    colnames(dists) <- nocell_BCs
    dists <- suppressWarnings(data.table::melt(dists,variable.factor = F,variable.name="falseBC", value.name="hamming"))
    dists <- dists[, trueBC := rep(true_chunks[[i]],length(nocell_BCs))][
      hamming<=BarcodeBinning,]
    if(i==1){
      binmap <- dists
    }else{
      binmap <- rbind(binmap,dists)
    }
  }
  #remove unused BCs that fit equally well to two true parent BCs
  binmap[    , min_ham :=  min(hamming), by = falseBC][
    , n_false :=  length(hamming), by = falseBC][
      , n_min := sum(hamming==min_ham), by =  falseBC]
  binmap <- binmap[n_min==1         ,][
    hamming==min_ham ,][
      , min_ham := NULL][
        , n_false := NULL][
          , n_min := NULL][
            , n := nocell_bccount[match(falseBC,nocell_bccount$XC),n]]
  
  print(paste("Found",nrow(binmap),"daughter barcodes that can be binned into",length(unique(binmap[,trueBC])),"parent barcodes."))
  print(paste("Binned barcodes correspond to",sum(binmap[,n]),"reads."))
  return(binmap)
}

prep_samtools <- function(featfiles,bccount,cores,project,BarcodeBinning,outdir){
  print("Extracting reads from bam file(s)...")
  nfiles=length(featfiles)
  nchunks <- length(unique(bccount$chunkID))
  all_rgfiles <- paste0(outdir,project,".RGgroup.",1:nchunks,".txt")
  
  
  for(i in unique(bccount$chunkID)){
    rgfile <- all_rgfiles[i]
    chunks <- bccount[chunkID==i]$XC
    if(BarcodeBinning > 0){
      write.table(file=rgfile,c(chunks,binmap[,falseBC]),col.names = F,quote = F,row.names = F)
    }else{
      write.table(file=rgfile,chunks,col.names = F,quote = F,row.names = F)
    }
  }
  
  headerXX <- paste( c(paste0("V",1:3)) ,collapse="\t")
  write(headerXX,paste0(outdir,"freadHeader"))
  
  headercommand <- paste0("cat ",outdir,"freadHeader > ")
  samcommand <- paste("samtools view -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x XS -@")
  grepcommand <- " | cut -f12,13,14 | sed 's/BC:Z://' | sed 's/UB:Z://' | sed 's/XT:Z://' | grep -F -f "
  
  outfiles_ex <- paste0(outdir,project,".ex.",1:nchunks,".txt")
  system(paste(headercommand,outfiles_ex,collapse = "; "))
  
  if(length(featfiles)==1){
    cpusperchunk <- round(cores/nchunks,0)
    ex_cmd <- paste(samcommand,cpusperchunk,featfiles[1],grepcommand,all_rgfiles,">>",outfiles_ex," & ",collapse = " ")
    
    system(paste(ex_cmd,"wait"))
  }else{
    cpusperchunk <- round(cores/(2*nchunks),0)
    ex_cmd <- paste(samcommand,cpusperchunk,featfiles[1],grepcommand,all_rgfiles,">>",outfiles_ex," & ",collapse = " ")
    
    outfiles_in <- paste0(outdir,project,".in.",1:nchunks,".txt")
    system(paste(headercommand,outfiles_in,collapse = "; "))
    
    in_cmd <- paste(samcommand,cpusperchunk,featfiles[2],grepcommand,all_rgfiles,">>",outfiles_in," & ",collapse = " ")
    
    outfiles_inter <- paste0(outdir,project,".inter.",1:nchunks,".txt")
    system(paste(headercommand,outfiles_inter,collapse = "; "))
    
    inter_cmd <- paste(samcommand,cpusperchunk,featfiles[3],grepcommand,all_rgfiles,">>",outfiles_inter," & ",collapse = " ")
    
    
    system(paste(ex_cmd,in_cmd,inter_cmd,"wait"))
  }
  system(paste0("rm ",outdir,"freadHeader"))
  system(paste("rm",all_rgfiles))
  
  return(outfiles_ex)
}



#### stats functions

sumstatBAM <-
  function(featfiles,
           cores,
           outdir,
           user_seq,
           spec_seq,
           bc,
           outfile,
           BCstats,
           nReadsperCell,
           mem_limit,
           bccount,
           project,
           BarcodeBinning,
           mem,
           read_layout,
           binmap) {
    require(data.table)
    # chunk barcodes
    bccount_file <- BCstats
    bccount <-
      data.table::fread(bccount_file , col.names = c("XC", "n"))[, list(n = sum(n)), by =
                                                                   XC][n >= nReadsperCell][order(-n)]
    if (sum(bccount$n) > 1e+09) {
      #for huge datasets, don't summarise "bad BCs"
      bccount <- bccount[XC %in% bc$XC]
    }
    bccount <- splitRG_stats(bccount = bccount, mem = mem_limit,read_layout = read_layout)
    
    samouts <- prep_samtools_stats(
      featfiles = featfiles,
      bccount   = bccount,
      cores     = cores,
      outdir = outdir,
      project = project,
      BarcodeBinning = BarcodeBinning,
      binmap=binmap
    )
    
    
                    for (i in unique(bccount$chunkID)) {
      print(paste("Working on chunk", i))
      
      samfile_ex <- paste0(outdir, project, ".ex.", i, ".txt")
      if (grepl(pattern = ".filtered.tagged.Aligned.out.bam.in.featureCounts.bam", featfiles[2])) {
        samfile_in <- paste0(outdir, project, ".in.", i, ".txt")
        samfile_inter <- paste0(outdir, project, ".inter.", i, ".txt")
      } else{
        samfile_in <- samfile_ex
      }
      
      
      tmp <- data.table::fread(
        samfile_ex,
        na.strings = c(""),
        select = c(1, 2, 3),
        header = T,
        fill = T,
        colClasses = "character" ,
        col.names = c("RG", "XS", "GE")
      )[, "GEin" := fread(
          samfile_in,
          select = 3,
          header = T,
          fill = T,
          na.strings = c(""),
          colClasses = "character"
        )][, "GEinter" := fread(
          samfile_inter,
          select = 3,
          header = T,
          fill = T,
          na.strings = c(""),
          colClasses = "character"
        )][, "ftype" := "NA"][
          is.na(GEinter) == F, ftype := "Intergenic"][
          is.na(GEin) == F, ftype := "Intron"][
            is.na(GE) == F  , ftype :="Exon"][
              is.na(GE)     , GE := GEin][
                is.na(GE)     , GE := GEinter][
                ftype != "NA",   XS := ftype][
                  GE %in% user_seq[, V1], XS :="User"][
                    GE %in% spec_seq[species=="mouse", V1], XS :="Intergenic_mouse"][
                      GE %in% spec_seq[species!="mouse", V1], XS :="Intergenic_human"][
                    , RG := as.character(RG)][
                      !(RG %in% bc[, XC]), RG := "bad"][
                        , c("GEin","GEinter", "ftype") :=NULL][
                          , list(.N), by = c("RG", "XS")][
                            , type := .rmUnassigned(XS)][
                              type=="NoFeatures",type:="Intergenic"][
                                , XS := NULL]
      
      if (i == 1) {
        mapCount <- tmp
      } else{
        mapCount <- rbind(mapCount, tmp)
      }
      system(paste("rm ", samfile_ex, samfile_in,samfile_inter))
    }
    saveRDS(mapCount, file = outfile)
    return(mapCount)
  }

getUserSeq <-function(gtf ) {
  if(!is.null(gtf)){
    user_seq<-fread(gtf,select=1:2,header = F)[V2=="User","V1"]
  }else{
    user_seq=NULL
  }
  return(user_seq)
}

getSpecies <-function(gtf ) {
  if(!is.null(gtf)){
    spec_seq<-fread(gtf,
                    select=1:2,header = F,
                    colClasses = "character")[V2!="User","V1"
                                              ][!duplicated(V1),"V1"
                                               ][,"species":="NA"
                                                 ][grep(V1,pattern="*_M"),species:="mouse"
                                                   ][grep(V1,pattern="*_M",invert = T),species:="human"]
  }else{
    spec_seq=NULL
  }
  return(spec_seq)
}

prep_samtools_stats <- function(featfiles,bccount,cores,outdir,project,BarcodeBinning,binmap){
  print("Extracting reads from bam file(s)...")
  nfiles=length(featfiles)
  nchunks <- length(unique(bccount$chunkID))
  all_rgfiles <- paste0(outdir,project,".RGgroup.",1:nchunks,".txt")
  
  
  for(i in unique(bccount$chunkID)){
    rgfile <- all_rgfiles[i]
    chunks <- bccount[chunkID==i]$XC
    if(BarcodeBinning > 0){
      write.table(file=rgfile,c(chunks,binmap[,falseBC]),col.names = F,quote = F,row.names = F)
    }else{
      write.table(file=rgfile,chunks,col.names = F,quote = F,row.names = F)
    }
  }
  
  headerXX <- paste( c(paste0("V",1:3)) ,collapse="\t")
  write(headerXX,paste0(outdir,"freadHeader"))
  
  headercommand <- paste0("cat ",outdir,"freadHeader > ")
  samcommand <- paste("samtools view -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x UB -@")
  grepcommand <- " | cut -f12,13,14 | sed 's/BC:Z://' | sed 's/XS:Z://' | sed 's/XT:Z://' | grep -F -f "
  
  outfiles_ex <- paste0(outdir,project,".ex.",1:nchunks,".txt")
  system(paste(headercommand,outfiles_ex,collapse = "; "))
  
  if(length(featfiles)==1){
    cpusperchunk <- round(cores/nchunks,0)
    ex_cmd <- paste(samcommand,cpusperchunk,featfiles[1],grepcommand,all_rgfiles,">>",outfiles_ex," & ",collapse = " ")
    
    system(paste(ex_cmd,"wait"))
  }else{
    cpusperchunk <- round(cores/(2*nchunks),0)
    ex_cmd <- paste(samcommand,cpusperchunk,featfiles[1],grepcommand,all_rgfiles,">>",outfiles_ex," & ",collapse = " ")
    
    outfiles_in <- paste0(outdir,project,".in.",1:nchunks,".txt")
    system(paste(headercommand,outfiles_in,collapse = "; "))
    
    in_cmd <- paste(samcommand,cpusperchunk,featfiles[2],grepcommand,all_rgfiles,">>",outfiles_in," & ",collapse = " ")
    
    outfiles_inter <- paste0(outdir,project,".inter.",1:nchunks,".txt")
    system(paste(headercommand,outfiles_inter,collapse = "; "))
    
    inter_cmd <- paste(samcommand,cpusperchunk,featfiles[3],grepcommand,all_rgfiles,">>",outfiles_inter," & ",collapse = " ")
    
    
    system(paste(ex_cmd,in_cmd,inter_cmd,"wait"))
  }
  system(paste0("rm ",outdir,"freadHeader"))
  system(paste("rm",all_rgfiles))
  
  return(outfiles_ex)
}


splitRG_stats<-function(bccount,mem,read_layout){
  if(is.null(mem) || mem==0){
    maxR<- Inf
  }else{
    maxR<- floor( mem*1000 * 4500 )
  }
  if( (maxR > 1e+09 & read_layout == "SE") | (maxR > 5e+08 & read_layout == "PE") ){
    maxR <- ifelse(read_layout == "SE",1e+09,5e+08)
  } 
  print(paste(maxR,"Reads per chunk"))
  nc<-nrow(bccount)
  cs=0
  chunkID=1
  bccount[,chunkID:=0]
  for(i in 1:nc){
    cs=cs+bccount[i]$n
    if(bccount[i]$n>maxR){
      print(paste("Warning: Barcode",bccount[i]$XC,"has more reads than allowed for the memory limit!
                  Proceeding anyway..."))
    }
    if(cs>=maxR){
      chunkID=chunkID+1
      cs=bccount[i][,"n"]
    }
    bccount[i][,"chunkID"]=chunkID
  }
  return(bccount)
}


.rmUnassigned<-function(b){ gsub("Unassigned_","",b)}