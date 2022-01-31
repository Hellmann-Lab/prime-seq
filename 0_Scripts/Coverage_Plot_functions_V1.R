makeGviz_geneModel_from_gencode<- function(gtf){
  gtf %>% as_tibble %>% 
    filter(type %in% c("exon","UTR","CDS")) %>% 
    transmute(chromosome =seqnames, 
              start, end, width, strand, type,
              gene = gene_id,
              exon = exon_id,
              transcript = transcript_id,
              symbol = gene_name) %>% 
    group_by(exon) %>% 
    dplyr::mutate( feature = case_when(type %in% "CDS" ~ "protein_coding",
                                       type %in% "UTR" ~ "utr",
                                       T ~ as.character(type))) %>% 
    dplyr::select(-type) %>% ungroup
  
}


### The plotting function

Plot_Coverage <- function(PRange,
                         annot.grs=NULL, 
                         annot.colors,
                         b.inf, 
                         gene.models,
                         figure.name,
                         gen="hg38",
                         type=c("ATAC","RNA"),
                         pdf=T,
                         title=NA,
                         extend.range=1e4,
                         ymax=1000){

  #annot.gr list of gr objects, the first one determines the plotted range and the highlighted sections
  #b.inf a data.frame with the columns cond1, cond2 and bamfile, the bamfile should give the full paths to the files for coverage plots
  #gtf path to a gtf to visualize the gene models
  # color.scheme list of colors to be used in tracks
  options(ucscChromosomeNames=FALSE)
  chr <- as.character(unique(seqnames(PRange)))
  axis   <- GenomeAxisTrack(genome = gen)
  
  # itrack <- IdeogramTrack(genome = gen, chromosome = chr)
  
  
  if(type=="ATAC"){
    # gene_models
    gencode.track<-GeneRegionTrack(gene.models, 
                                   chromosome =chr,
                                   start = ifelse(start(PRange)>1e6,start(PRange)-1e6,1),
                                   end = end(PRange)+1e6,
                                   genome = gen,
                                   transcriptAnnotation = "gene_name",
                                   showId=TRUE, 
                                   geneSymbol=TRUE, 
                                   fill= "darkblue", col = "darkblue", col.axis="black", 
                                   col.title="black", cex.title=0.8, name= "Gencode")
    
  # annotation tracks
  annot.tracks <- lapply( names(annot.grs), function(nn){
    annot.grs[[nn]]<-keepStandardChromosomes(annot.grs[[nn]])
    tr<-AnnotationTrack(annot.grs[[nn]],
                        name=nn , 
                        chromosome = chr,
                        fill = as.character(unname(annot.colors[nn])),
                        col = as.character(unname(annot.colors[nn])),
                        col.axis ="black", col.title="black",
                        col.title="black", cex.title=0.8) })
  
  
  # ATAC-seq coverage 
  atac.tracks<- lapply( 1:dim(b.inf)[1] , function(i){
    DataTrack( range =b.inf$bamfile[i], 
               type = 'h', 
               chromosome = chr,
               col = as.character(unname(annot.colors[  b.inf$cond1[i]])) ,
               name = paste(b.inf$cond2[i], b.inf$cond1[i]),
               window = -1, windowSize = 100, genome = gen,
               col.axis="black", col.title="black", cex.title=0.5) })
  
  ht <- HighlightTrack( trackList = c(annot.tracks, atac.tracks) ,
                        start = start(annot.grs[[1]]), 
                        end   = end(annot.grs[[1]]),
                        chromosome = chr, 
                        fill = "grey",
                        col  = "grey" )
  
  track.list<- c(axis, gencode.track, ht)
  
  }else{
    # gene_models
    gencode.track<-GeneRegionTrack(gene.models, 
                                   chromosome =chr,
                                   start = ifelse(start(PRange)>1e6,start(PRange)-1e6,1),
                                   end = end(PRange)+1e6,
                                   genome = gen,
                                   transcriptAnnotation = "gene_name",
                                   showId=TRUE, 
                                   geneSymbol=TRUE, 
                                   fill= "darkblue", col = "darkblue", col.axis="black", 
                                   col.title="black", cex.title=0.5, 
                                   name= "Gencode",
                                   cex=0.1)
    
rna.tracks<- lapply( 1:dim(b.inf)[1] , function(i){
  DataTrack( range =b.inf$bamfile[i], 
             type = 'h', 
             chromosome = chr,
             col = as.character(unname(annot.colors[b.inf$cond1[i]])),
             name = paste(b.inf$cond1[i]),
             window= -1, windowSize = 50, genome = gen,
             col.axis="black", col.title="black", cex.title=0.6,
             ylim=c(0,ymax)) })

ht <- HighlightTrack( trackList = c(rna.tracks) ,
                      start = start(PRange),
                      end   = end(PRange),
                      chromosome = chr,
                      fill = "white",
                      col  = "white" )

track.list<- c(axis, ht,gencode.track)
  }

if(is.na(title)){
  title<-paste0("chr",chr,
                ": ",
               start(PRange)-extend.range,
               " - ", 
               end(PRange)+extend.range
               )
}  
  
  
if(pdf){  
  pdf(figure.name,width=6,height=4)
  plotTracks( track.list, 
              collapseTranscripts = "meta", shape = "arrow", 
              from = start(PRange)-extend.range, 
              to   =   end(PRange)+extend.range,
              title.width = 1,
              col.grid='grey75',
              main = title,
              cex.main = 1)
  dev.off()
}else{
  return(plotTracks( track.list, 
                  collapseTranscripts = "meta", shape = "arrow", 
                  from = start(PRange)-extend.range, 
                  to   =   end(PRange)+extend.range,
                  title.width = 1,
                  col.grid='grey75',
                  main = title,
                  cex.main = 1))
         }
  }


### zUMIs functions to make exon and intron annotation----

makeSAF<-function(gtf, extension_var=FALSE, exon_extension=0, buffer_length=100, scaff_length=0,multi_overlap_var=FALSE){
  print("Loading reference annotation from:")
  print(gtf)
  gtf_file <- plyranges::read_gff(gtf, col_names = c("type","gene_id"))
  
  # Check for gene information in gtf file --------------
  # if "gene" is in type column from GTF, uses those desingations. If not, creates
  # such a field from the max-end and min-start values after grouping by gene_id.
  
  if("gene" %in% levels(gtf_file$type)==TRUE){
    # this represents total length of each gene (unextended)
    gr.gene <- gtf_file %>%
      dplyr::filter(type =="gene")
  }else{
    gr.gtf.filt <- gtf_file %>%
      dplyr::filter(type =="exon")
    
    gr.gene <- gr.gtf.filt %>%
      group_by(gene_id, strand, seqnames) %>%
      summarise(end = max(end),
                start = min(start),
                type = "gene" ) %>%
      plyranges::as_granges()
  }
  
  
  #Calculating total intergenic length for intron-probability
  totalGeneLength<-sum(width(gr.gene))
  #get length of chromosomes from BAM header
  chrL<-system( paste("samtools view -H", abamfile ," | grep '^@SQ' | cut -f3 "), intern = T)
  chrL<-lapply(chrL, function(x){
    x<-as.integer(strsplit(x, ":")[[1]][[2]])
  })
  chrL<-do.call(sum,chrL)
  intergenicBp<- chrL-totalGeneLength
  
  # getting names and length of chromosomes/scaffolds from bam header
  chrom_kept <- .chromLengthFilter(abamfile = abamfile, keep_chr_minLength = scaff_length)
  
  # reading out chromosomes/scaffold names from the GRange gtf object
  chrom_gtf_names <- tibble::as_tibble(x = list("chrom_names" = as.character(seqinfo(gtf_file)@seqnames)),
                                       .name_repair = "universal")
  ## Exon creation ##
  # check to see if any extension is going to occur
  # if so, then do extension for everything. If not still do basic filtering for exons, reducing for
  # only scaffolds chosen to be kept.
  
  # adding check for NULL, which would come from YAML with no extension option
  if(is.null(extension_var)) {
    extension_var <- FALSE
  }
  
  if(extension_var == TRUE) {
    exon <- .exon_extension(gtf_input = gtf_file,
                            exon_extension = exon_extension,
                            buffer_length = buffer_length,
                            chrom_kept = chrom_kept)
  } else {
    filt_gtf <- gtf_file %>%
      plyranges::filter(type == "exon") %>%
      plyranges::group_by(gene_id) %>%
      plyranges::reduce_ranges_directed()
    
    exon <- filt_gtf[as.character((seqnames(filt_gtf))) %in% chrom_kept$name]
  }
  
  
  # this reduces all exons across all genes, and performs gaps to make introns
  gap_obj <- exon %>%
    plyranges::reduce_ranges() %>%
    GenomicRanges::gaps()
  
  # adding code to do disjoin for transcript objects, and then overlapping to get GeneIDs back into introns -CA 03.01.20
  gr.gene.dis <- gr.gene %>%
    plyranges::disjoin_ranges() %>%
    plyranges::find_overlaps_within(.,gr.gene) %>%
    plyranges::mutate(index = 1:length(.))
  
  # additional code to ensure introns are only assigned to one gene annotation
  rm_ind <- plyranges::find_overlaps(gr.gene.dis,               # find overlapping ranges in same object from different genes, assign to "rm_ind"
                                     gr.gene.dis,
                                     minoverlap = 1) %>%
    plyranges::filter(gene_id.x != gene_id.y) %>%               # because this is a self-comparison, filtering out genes that match each other
    data.frame() %>%
    dplyr::pull(index.x) %>%                                    # pull indices, as these represent ranges that are overlapping and belong to
    unique()                                                    # different genes.
  
  gr.gene.dis.filt <- gr.gene.dis %>%
    plyranges::filter(!(index %in% rm_ind)) %>%  # removing based on "rm_ind"
    plyranges::select(-index)                    # getting rid of index after removal, as we don't need anymore
  
  # this finds introns that are only inside of the disjoined gene's transcript object
  intron <- plyranges::find_overlaps_within(gap_obj, gr.gene.dis.filt)
  
  # filtering for introns to be greater than 10 bp and less than 1e5
  intron.filt <- intron %>%
    plyranges::filter(width > 10 & width < 1e5)
  
  # Creation of intron and Exon Saf
  intron.saf <- data.frame(GeneID = intron.filt$gene_id,
                           Chr = seqnames(intron.filt),
                           Start = start(intron.filt),
                           End = end(intron.filt),
                           Len = width(intron.filt),  # needed for intron-probability
                           stringsAsFactors = FALSE)
  
  exon.saf <- data.frame(GeneID = exon$gene_id,
                         Chr    = seqnames(exon),
                         Start  = start(exon),
                         End	  = end(exon),
                         Strand = strand(exon),
                         stringsAsFactors = F)
  
  # Getting strand info back into intron saf
  intron.saf <- dplyr::left_join(intron.saf, unique(exon.saf[,c("GeneID","Strand")]), by = c("GeneID"))
  
  #get intronic bp per gene and intronic reads per gene
  intronGenes <- intron.saf %>% group_by(GeneID) %>% summarize(IntronLengthPerGene = sum(Len), intronsPerGene = length(GeneID), .groups = "drop")
  
  saf <- list(introns = unique(intron.saf), exons = unique(exon.saf), intronsPerGene = intronGenes,intergenicBp = intergenicBp)
  
  
  print("Annotation loaded!")
  
  rm(gr.gene,
     intron,
     exon,
     intron.saf,
     exon.saf,
     intronGenes,
     intergenicBp)
  return(saf)
}

# helper functions

.chromLengthFilter <- function(abamfile, keep_chr_minLength =0){
  bread<- paste("samtools view -H ", abamfile , "| grep '^@SQ' | cut -f2,3 ")
  chr<-data.table::fread(bread,col.names = c("chr","len"), header = F)[
    , chr2 :=gsub("SN:","",chr)][ 
      , len2 :=gsub("LN:","",len)][
        , len2 := as.numeric(len2)][len2>=keep_chr_minLength]
  out <- data.frame("name" = chr$chr2, "length" = chr$len2, stringsAsFactors = FALSE)
  return(out)
}

# Chris Alford's overall function for extending exons. Relies on .extension_funcion, .chromLengthFilter, get_gr,
.exon_extension <- function(gtf_input, exon_extension, buffer_length, chrom_kept){
  suppressWarnings(suppressMessages(require(GenomicRanges)))
  suppressWarnings(suppressMessages(require(GenomicFeatures)))
  suppressWarnings(suppressMessages(require(GenomicAlignments)))
  suppressWarnings(suppressMessages(require(AnnotationDbi)))
  suppressWarnings(suppressMessages(require(plyranges)))
  
  ex_ex <- exon_extension      # amount to extent coding region by in bp
  bl <- buffer_length          # the amount of buffer space in bp between genes
  
  print(paste("starting exon-extension at", Sys.time() ))
  print(paste("extending by", exon_extension, "bp"))
  # reading out chromosomes/scaffold names from the GRange gtf object
  chrom_gtf_names <- tibble::as_tibble(x = list("chrom_names" = as.character(seqinfo(gtf_input)@seqnames)),
                                       .name_repair = "universal")
  
  # doing left join of length data from bam and names from gtf GRange to
  # preserve order of chromosomes/scaffolds from the gtf
  # adding length info back into GRange object from gtf file.
  seqlengths(gtf_input) <- chrom_kept %>%
    dplyr::left_join(chrom_gtf_names, ., by = c( "chrom_names" = "name")) %>%
    dplyr::pull(length)                     # vectorizing length, as seqlengths() requires this
  
  
  # first filtering for exons, making a concensus of transcripts gene wise with reduce, and adding index to keep track of position
  filt_gtf <- gtf_input %>%
    plyranges::filter(type == "exon") %>%
    plyranges::group_by(gene_id) %>%
    plyranges::reduce_ranges_directed()  # adding index and start for manipulations
  
  # Checking whether we will filter out scaffolds based on size here.
  # if the seqnames of our filt_gtf are present in the chrom_kept object from our BAM file, then we keep.
  filt_sub_gtf <- filt_gtf[as.character((seqnames(filt_gtf))) %in% chrom_kept$name]
  
  # was having issue with out of bound indices, so removing scaffolds before adding index.
  filt_sub_gtf <- filt_sub_gtf %>%
    plyranges::mutate(index = 1:length(.),
                      start2 = start(.))
  
  
  # From here we subset out max / min exons gene wise according to strand (+)/(-).
  max_pos <- filt_sub_gtf %>%
    plyranges::filter(strand(.) == "+") %>%
    plyranges::group_by(gene_id) %>%
    plyranges::filter(start2 == max(start2))
  
  max_pos$start2 <- NULL # apparently plyranges no longer likes the option plyranges::select(-blah) , so setting metacolumn to null is more robust?
  
  min_neg <- filt_sub_gtf %>%
    plyranges::filter(strand(.) == "-") %>%
    plyranges::group_by(gene_id) %>%
    plyranges::filter(start2 == min(start2))
  
  min_neg$start2 <-NULL    # apparently plyranges no longer likes the option plyranges::select(-blah) , so setting metacolumn to null is more robust?
  
  # removing start2 column from whole object
  mcols(filt_sub_gtf)$start2 <- NULL
  
  # doing preceed over whole object. This returns vector of the following exon for each index in ref to starting object
  pre_all <- GenomicRanges::precede(filt_sub_gtf,
                                    filt_sub_gtf,
                                    ignore.strand = FALSE) # made a change here.
  
  
  
  # Some entries will be NA in our pre_all vector.
  # Most likely represent exons on scaffolds or last exon in chromosome in either direction
  names(pre_all) <- 1:length(pre_all)                           # naming with position in vector
  na_list <- as.numeric(names(pre_all[is.na(pre_all)]))         # getting locations of NAs out. These correspond to ranges that have no following exon.
  
  # dealing with positive NA's in preceed table
  pos_na_ind <- filt_sub_gtf[na_list] %>%
    plyranges::filter(strand(.) == "+")
  pos_na_ind <- pos_na_ind$index
  
  neg_na_ind <- filt_sub_gtf[na_list] %>%
    plyranges::filter(strand(.) == "-")
  neg_na_ind <- neg_na_ind$index
  
  # actual extension for NAs. I don't have to go piecemeal, because I'm just assigning the same amount each time.
  # also check for chromosome keeping is done above
  end(filt_sub_gtf[pos_na_ind]) <- end(filt_sub_gtf[pos_na_ind]) + ex_ex
  start(filt_sub_gtf[neg_na_ind]) <- start(filt_sub_gtf[neg_na_ind]) - ex_ex
  
  # trimming extension to only go to end of scaffold / chromosome. Also prints which ones are out of bound.
  # internal record of which ones are truncated
  out_of_bound_ind <- GenomicRanges:::get_out_of_bound_index(filt_sub_gtf)
  filt_sub_gtf <- trim(filt_sub_gtf)
  
  # index of the 3' exons in + and - direction
  index_neg <- mcols(min_neg)$index
  index_pos <- mcols(max_pos)$index
  
  
  # Create my list from dataframes
  df <- data.frame(filt_sub_gtf)
  
  # splits data frame into a list based on the gene ID
  df_list<- df %>%
    dplyr::group_split(gene_id)
  
  .extension_function <- function(x, pre_all, df){
    if (x %>%
        dplyr::select(strand) %>%              # check to see if the strand in positive
        unique == "+"){
      strand_var <- "pos"                      # Assign strand variable to "pos"
      ind <- x %>%
        dplyr::filter(end == max(end)) %>%     # finding the exon that has the max "end" value and
        dplyr::pull(index)                     # pulling the index for this entry
    } else {strand_var <- "neg"
    ind <- x %>%                               # likewise, doing the opposite for negative stranded genes
      dplyr::filter(start == min(start)) %>%
      dplyr::pull(index)                       # pulling index for min most exon for negative strand
    }
    if(is.na(pre_all[ind])){                   # pulling following exon of min/max exons in each list entry
      #print("NA value, skipping instance")     # check for NA values in pre_all vector
      return(x)
    } else {
      if(strand_var == "pos"){                 # Check the strand var
        small <- df[ind,]$end                  # If pos, assign preceeding exon as "small" and following exon as "big"
        big <- df[pre_all[ind],]$start
      } else {
        big <- df[ind,]$start                  # if negative assign preceeding exon as big and following exon as small
        small <- df[pre_all[ind],]$end
      }
      # actual extension part
      if(big - small > ex_ex + bl){            # if distance between exons is greater than buffer and extension
        ext <- ex_ex                           # then the full extension is used
      } else {
        if(big - small <= bl){                 # Check if the exons are are closer than buffer already ( ex, if exons are 15 bp away or something)
          ext <- 0                             # if so, setting extension to 0
        } else {
          ext <- big - small - bl              # If all else is false, then set extension to difference between start and end points, minus buffer
        }
      }
      if(strand_var == "pos"){                                       #
        #print(paste("positive gene index",
        #            ind,"extension of",ext))
        x[x$index == ind,]$end <- x[x$index == ind,]$end + ext
        return(x)
      } else{
        #print(paste("negative gene index",
        #            ind,"extension of",ext))
        x[x$index == ind,]$start <- x[x$index == ind,]$start - ext
        return(x)
      }
    }
    return(x)
  }
  
  
  # extension function for the groups
  # each entry in the list will be all exons for the gene
  # index for each range represents position in original Grange.
  df_out <- parallel::mclapply(df_list, .extension_function, pre_all = pre_all, df = df, mc.cores = opt$num_threads)
  
  # unbinding list
  df_output <- dplyr::bind_rows(df_out)
  
  # getting output into a Genomic Range object
  df_out_gr <- get_gr(df_output)
  
  #final return
  print("end extension")
  return(df_out_gr)
}
