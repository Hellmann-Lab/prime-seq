## make intron and exon bed files
# from zUMIs
source("/data/home/wange/ATAC/Coverage_Plot_functions_V1.R")


gtf_file<-"/data/share/htp/prime-seq_Paper/genomes/Hsap/gencode.v35.primary_assembly.annotation.gtf"
abamfile<-bam_files

saf<-makeSAF(gtf = gtf_file)

gr_list<-lapply(saf[1:2],function(le){
       le %>%
         transmute(seqnames=Chr,
                   start=Start-1,
                   end=End,
                   name=GeneID,
                   strand=Strand,
                   length=ifelse(exists("Len"),Len,NA)) 
  })

## split bamfile based on bed files

for (i in c("exons","introns")){
  write.table(gr_list[[i]][,1:3],file=paste(gsub(x = bam_files,pattern=".bam",replacement = ""),i,'granges.bed',sep="."), quote=F, sep="\t", row.names=F, col.names=F)
  
system(
  paste(sep=" ","sbatch --wrap 'samtools view -b -h -L",
        paste(gsub(x = bam_files,pattern=".bam",replacement = ""),i,'granges.bed',sep="."),
       bam_files,
       ">",
       paste(gsub(x = bam_files,pattern=".bam",replacement = ""),i,"bam",sep="."),
       "'")
  )

system(paste(sep=" ","sbatch --wrap 'samtools index",
       paste(gsub(x = bam_files,pattern=".bam",replacement = ""),i,"bam",sep="."),
       "'"))
}
