.libPaths("/home/cog/inijman/Rlibs/3.4.1")
library(ExomeDepth)
library(methods)

#plot sample specific gene loci from sample all_exons object
my.current.samplename = "U173589CM2016D19618_dedup.realigned.bam"
my.current.gene="SBDS"
event.flank = 10000

#load required files
target.file <-"/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/exomeCopy_crev2_20bp_flat/annotated_ENSEMBL_UCSC_merged_collapsed_sorted_v2_20bpflank_flatnoMT.bed" # set path to BED file
reference.file<-"/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"
genes<-as.data.frame(read.table(file="/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/GENEBODY_LOCATIONS.txt",sep="\t",header=T))
data(exons.hg19.X)
data(Conrad.hg19)


load(file = paste(my.current.samplename,"all.exons",sep = "_"))

#plot specific gene
my.genelocus<-subset(genes, Gene_name==my.current.gene)
plot (all.exons,
      sequence = toString(my.genelocus$Chromosome_name),
      xlim = c( my.genelocus$Gene_start- event.flank, my.genelocus$Gene_end + event.flank),
      count.threshold = 0,
      main = paste(my.current.samplename,my.current.gene,sep="_"),
      cex.lab = 0.8,
      with.gene = TRUE)



#plot all calls
for (i in 1:nrow(all.exons@CNV.calls)) {
    event.chr=all.exons@CNV.calls[i,]$chromosome
    event.start=all.exons@CNV.calls[i,]$start
    event.end=all.exons@CNV.calls[i,]$end
    event.size=(event.end-event.start)
    event.flank=100000
    
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
    dev.off()
}
  

