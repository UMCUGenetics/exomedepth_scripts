.libPaths("/home/cog/inijman/Rlibs/3.4.1")

setwd("/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/wes_set2/exomedepth")
library(ExomeDepth)
target.file <-"/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/exomeCopy_crev2_20bp_flat/annotated_ENSEMBL_UCSC_merged_collapsed_sorted_v2_20bpflank_flatnoMT.bed" # set path to BED file
pathToBams <- "/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/Gijn/test_samples/" 
bam.files <- paste0(pathToBams, dir(pathToBams, "bam$"))
pathToRefBams <- "/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/exomeDepth/refbams/" 
refbam.files <- paste0(pathToRefBams, dir(pathToRefBams, "bam$"))


reference.file<-"/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"

data(exons.hg19)
#data(exons.hg19.X)
data(Conrad.hg19)
my.counts <- getBamCounts(bed.frame = exons.hg19,
                          bam.files = bam.files,
                          include.chr = FALSE,
                          referenceFasta = reference.file)
save(my.counts,file = "my.counts")

#my.refcounts <- getBamCounts(bed.frame = exons.hg19,
#                         bam.files = refbam.files,
#                          include.chr = FALSE,
#                         referenceFasta = reference.file)
#save(my.refcounts,file = "my.refcounts")


load(file="my.counts")
load(file="my.refcounts")

my.counts.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')
my.counts.dafr$chromosome <- gsub(as.character(my.counts.dafr$space),
                                   pattern = 'chr',
                                   replacement = '') ##remove the annoying chr letters


samplecounts.mat<-as.matrix(my.counts.dafr[,grep(names(my.counts.dafr),pattern='*.bam')])
nsamples<-ncol(samplecounts.mat)

my.refcounts.dafr <- as(my.refcounts[, colnames(my.refcounts)], 'data.frame')
my.refcounts.dafr$chromosome <- gsub(as.character(my.refcounts.dafr$space),
                                  pattern = 'chr',
                                  replacement = '') ##remove the annoying chr letters
my.ref.samples<-colnames(my.refcounts.dafr)[7:(ncol(my.refcounts.dafr)-1)]


#loop over samples in my.counts
for (i in 1:nsamples) {
  my.current.samplename <-colnames(my.counts.dafr[6+i])
  print(my.current.samplename)
#  my.current.sample<-samplecounts.mat$my.current.sample
  my.reference.set <- as.matrix(my.refcounts.dafr[,my.ref.samples])
  my.choice<-select.reference.set(test.counts=samplecounts.mat[,i],
                                 reference.counts=(my.reference.set),
                                 bin.length=(my.counts.dafr$end - my.counts.dafr$start)/1000,
                                 n.bins.reduced = 10000)

  print(my.choice[[1]])
  my.matrix <- as.matrix( my.refcounts.dafr[, my.choice$reference.choice, drop = FALSE])
  my.reference.selected <- apply(X = my.matrix,
                               MAR = 1,
                               FUN = sum)

  #CNV calling
  all.exons <- new('ExomeDepth',
                 test = samplecounts.mat[,i],
                 reference = my.reference.selected,
                 formula = 'cbind(test, reference) ~ 1')

  all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 10^-4,
                      chromosome = my.counts.dafr$space,
                      start = my.counts.dafr$start,
                      end = my.counts.dafr$end,
                      name = my.counts.dafr$names)


  all.exons <- AnnotateExtra(x = all.exons,
                           reference.annotation = Conrad.hg19.common.CNVs,
                           min.overlap = 0.5,
                           column.name = 'Conrad.hg19')

  #print(head(all.exons@CNV.calls))
  save(all.exons,file=paste(my.current.samplename,"all.exons",sep = "_"))
  output.file <- paste(my.current.samplename,'exome_calls.csv',sep = "")
  write.csv(file = output.file,
          x = all.exons@CNV.calls,
          row.names = FALSE)
#  event.chr=17
#  event.start=72200074
#  event.end=72667818
#  plot (all.exons,
#      sequence = event.chr,
#      xlim = c(event.start - 1000000, event.end + 100000),
#      count.threshold = 20,
#      main = paste(my.current.samplename,' (',event.chr,':',event.start,'-',event.end,')',sep=""),
#      cex.lab = 0.8,
#      with.gene = TRUE)
}
