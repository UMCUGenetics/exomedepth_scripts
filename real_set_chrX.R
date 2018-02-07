.libPaths("/home/cog/inijman/Rlibs/3.4.1")
getwd()
library(ExomeDepth)
library(methods)
target.file <-"/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/exomeCopy_crev2_20bp_flat/annotated_ENSEMBL_UCSC_merged_collapsed_sorted_v2_20bpflank_flatnoMT.bed" # set path to BED file
pathToBams <- "/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/wes_set2/" 
bam.files <- paste0(pathToBams, dir(pathToBams, "bam$"))
pathToRefBams <- "/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/exomeDepth/refbams/" 
refbam.files <- paste0(pathToRefBams, dir(pathToRefBams, "bam$"))


reference.file<-"/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"


data(exons.hg19.X)
data(Conrad.hg19)
#my.Xcounts <- getBamCounts(bed.frame = exons.hg19.X,
#                          bam.files = bam.files,
#
#                          include.chr = FALSE,
#                          referenceFasta = reference.file)
#save(my.Xcounts,file = "my.Xcounts")

#my.refXcounts <- getBamCounts(bed.frame = exons.hg19.X,
#                         bam.files = refbam.files,
#                          include.chr = FALSE,
#                         referenceFasta = reference.file)
#save(my.refXcounts,file = "my.refXcounts")


load(file="my.Xcounts")
load(file="my.refXcounts")

my.Xcounts.dafr <- as(my.Xcounts[, colnames(my.Xcounts)], 'data.frame')
my.Xcounts.dafr$chromosome <- gsub(as.character(my.Xcounts.dafr$space),
                                   pattern = 'chr',
                                   replacement = '') ##remove the annoying chr letters

print(head(my.Xcounts.dafr))
sampleXcounts.mat<-as.matrix(my.Xcounts.dafr[,grep(names(my.Xcounts.dafr),pattern='*.bam')])
nsamples<-ncol(sampleXcounts.mat)

my.refXcounts.dafr <- as(my.refXcounts[, colnames(my.refXcounts)], 'data.frame')
my.refXcounts.dafr$chromosome <- gsub(as.character(my.refXcounts.dafr$space),
                                  pattern = 'chr',
                                  replacement = '') ##remove the annoying chr letters
my.ref.samples<-colnames(my.refXcounts.dafr)[7:(ncols(my.refXcounts.dafr)-1)]
print(head(my.refXcounts.dafr))

#loop over samples in my.counts
for (i in 1:nsamples) {
  my.current.samplename <-colnames(my.Xcounts.dafr[6+i])
  print(my.current.samplename)
#  my.current.sample<-samplecounts.mat$my.current.sample
  my.reference.set <- as.matrix(my.refXcounts.dafr[,my.ref.samples])
  my.choice<-select.reference.set(test.counts=sampleXcounts.mat[,i],
                                 reference.counts=(my.reference.set),
                                 bin.length=(my.Xcounts.dafr$end - my.Xcounts.dafr$start)/1000,
                                 n.bins.reduced = 10000)

  print(my.choice[[1]])
  my.matrix <- as.matrix( my.refXcounts.dafr[, my.choice$reference.choice, drop = FALSE])
  my.reference.selected <- apply(X = my.matrix,
                               MAR = 1,
                               FUN = sum)

  #CNV calling
  all.exons <- new('ExomeDepth',
                 test = sampleXcounts.mat[,i],
                 reference = my.reference.selected,
                 formula = 'cbind(test, reference) ~ 1')

  all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 10^-4,
                      chromosome = my.Xcounts.dafr$space,
                      start = my.Xcounts.dafr$start,
                      end = my.Xcounts.dafr$end,
                      name = my.Xcounts.dafr$names)


  all.exons <- AnnotateExtra(x = all.exons,
                           reference.annotation = Conrad.hg19.common.CNVs,
                           min.overlap = 0.5,
                           column.name = 'Conrad.hg19')

  #print(head(all.exons@CNV.calls))

  output.file <- paste(my.current.samplename,'exome_Xcalls.csv',sep = "")
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
