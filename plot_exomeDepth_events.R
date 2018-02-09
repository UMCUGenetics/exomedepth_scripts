.libPaths("/home/cog/inijman/Rlibs/3.4.1")

library(ExomeDepth)
library(methods)
target.file <-"/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/exomeCopy_crev2_20bp_flat/annotated_ENSEMBL_UCSC_merged_collapsed_sorted_v2_20bpflank_flatnoMT.bed" # set path to BED file
pathToBams <- getwd() 
bam.files <- paste0(pathToBams,"/", dir(pathToBams,"bam$"))
pathToRefBams <- "/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/exomeDepth/refbams/" 
refbam.files <- paste0(pathToRefBams, dir(pathToRefBams, "bam$"))


reference.file<-"/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"


data(exons.hg19.X)
data(Conrad.hg19)

load(file="my.counts")
load(file="my.refcounts")

my.counts.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')
my.counts.dafr$chromosome <- gsub(as.character(my.counts.dafr$space),
                                   pattern = 'chr',
                                   replacement = '') ##remove the annoying chr letters

#print(head(my.counts.dafr))
samplecounts.mat<-as.matrix(my.counts.dafr[,grep(names(my.counts.dafr),pattern='*.bam')])
nsamples<-ncol(samplecounts.mat)

my.refcounts.dafr <- as(my.refcounts[, colnames(my.refcounts)], 'data.frame')
my.refcounts.dafr$chromosome <- gsub(as.character(my.refcounts.dafr$space),
                                      pattern = 'chr',
                                      replacement = '') ##remove the annoying chr letters
my.ref.samples<-colnames(my.refcounts.dafr)[7:131]
#print(head(my.refcounts.dafr))

#loop over samples in my.counts
for (i in 1:nsamples) {
  my.current.samplename <-colnames(my.counts.dafr[6+i])
  print(my.current.samplename)
  load(file = paste(my.current.samplename,"all.exons",sep = "_"))
  
  for (i in 1:nrow(all.exons@CNV.calls)) {
    event.chr=all.exons@CNV.calls[i,]$chromosome
    event.start=all.exons@CNV.calls[i,]$start
    event.end=all.exons@CNV.calls[i,]$end
    event.size=(event.end-event.start)
    event.flank=1000000
    
    if (all.exons@CNV.calls[i,]$BF<10) next
    if (!is.na(all.exons@CNV.calls[i,]$Conrad.hg19)) next
    
    filename<-paste(my.current.samplename,event.chr,event.start,event.end,sep="_")
    
    print(filename)
    png(filename=paste(filename,".png",sep=""),width = 1750, height = 600)
    
    plot (all.exons,
      sequence = event.chr,
      xlim = c(event.start - event.flank, event.end + event.flank),
      count.threshold = 20,
      main = paste(my.current.samplename,' (',event.chr,':',event.start,'-',event.end,')',sep=""),
      cex.lab = 0.8,
      with.gene = TRUE)
    segments(event.start,1,event.end,1,col="red",lwd=8)
    dev.off()
  }
  
}
