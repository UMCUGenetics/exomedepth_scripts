.libPaths("/home/cog/inijman/Rlibs/3.4.1")
library(ExomeDepth)
library(methods)

target.file <-"/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/exomeCopy_crev2_20bp_flat/annotated_ENSEMBL_UCSC_merged_collapsed_sorted_v2_20bpflank_flatnoMT.bed" # set path to BED file
pathToBams <- getwd() 
bam.files <- paste0(pathToBams, "/", dir(pathToBams, "bam$"))
reference.file<-"/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"
data(exons.hg19)
data(Conrad.hg19)


load(file="/home/cog/inijman/scripts/exomedepth_scripts/my.NovaS_refcounts")


my.counts <- getBamCounts(bed.frame = exons.hg19,
                          bam.files = bam.files,
                          include.chr = FALSE,
                          referenceFasta = reference.file)
save(my.counts,file = "my.counts")


my.counts.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')
my.counts.dafr$chromosome <- gsub(as.character(my.counts.dafr$space),
                                   pattern = 'chr',
                                   replacement = '') ##remove the annoying chr letters

print(head(my.counts.dafr))
samplecounts.mat<-as.matrix(my.counts.dafr[,grep(names(my.counts.dafr),pattern='*.bam')])
nsamples<-ncol(samplecounts.mat)

my.refcounts.dafr <- as(my.refcounts[, colnames(my.refcounts)], 'data.frame')
my.refcounts.dafr$chromosome <- gsub(as.character(my.refcounts.dafr$space),
                                  pattern = 'chr',
                                  replacement = '') ##remove the annoying chr letters
my.ref.samples<-colnames(my.refcounts.dafr)[7:(ncol(my.refcounts.dafr)-1)]
print(head(my.refcounts.dafr))

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

  output.file <- paste(my.current.samplename,'exome_calls.csv',sep = "")
  save(all.exons,file=paste(my.current.samplename,"all.exons",sep = "_"))
  write.csv(file = output.file,
          x = all.exons@CNV.calls,
          row.names = FALSE)
}



################## Call X events  ###################################
data(exons.hg19.X)
load(file="/home/cog/inijman/scripts/exomedepth_scripts/my.NovaS_Xrefcounts")
my.Xcounts <- getBamCounts(bed.frame = exons.hg19.X,
                          bam.files = bam.files,

                          include.chr = FALSE,
                          referenceFasta = reference.file)
save(my.Xcounts,file = "my.Xcounts")


my.Xcounts.dafr <- as(my.Xcounts[, colnames(my.Xcounts)], 'data.frame')
my.Xcounts.dafr$chromosome <- gsub(as.character(my.Xcounts.dafr$space),
                                   pattern = 'chr',
                                   replacement = '') ##remove the annoying chr letters

print(head(my.Xcounts.dafr))
sampleXcounts.mat<-as.matrix(my.Xcounts.dafr[,grep(names(my.Xcounts.dafr),pattern='*.bam')])
nsamples<-ncol(sampleXcounts.mat)

my.refXcounts.dafr <- as(my.Xrefcounts[, colnames(my.Xrefcounts)], 'data.frame')
my.refXcounts.dafr$chromosome <- gsub(as.character(my.refXcounts.dafr$space),
                                  pattern = 'chr',
                                  replacement = '') ##remove the annoying chr letters
my.ref.samples<-colnames(my.refXcounts.dafr)[7:(ncol(my.refXcounts.dafr)-1)]
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

  output.file <- paste(my.current.samplename,'exome_Xcalls.csv',sep = "")
  save(all.exons,file=paste(my.current.samplename,"all.Xexons",sep = "_"))
  write.csv(file = output.file,
          x = all.exons@CNV.calls,
          row.names = FALSE)
}
