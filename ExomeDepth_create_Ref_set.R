.libPaths("/home/cog/inijman/Rlibs/3.4.1")

library(ExomeDepth)
library(methods)

target.file <-"/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/exomeCopy_crev2_20bp_flat/annotated_ENSEMBL_UCSC_merged_collapsed_sorted_v2_20bpflank_flatnoMT.bed" # set path to BED file
pathToRefBams <- "/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/NOVASEQ/" 
refbam.files <- paste0(pathToRefBams, dir(pathToRefBams, "bam$"))


reference.file<-"/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"

data(exons.hg19)


my.refcounts <- getBamCounts(bed.frame = exons.hg19,
                         bam.files = refbam.files,
                          include.chr = FALSE,
                         referenceFasta = reference.file)
save(my.refcounts,file = "/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/NOVASEQ/my.NovaS_refcounts")



data(exons.hg19.X)
my.Xrefcounts <- getBamCounts(bed.frame = exons.hg19.X,
                         bam.files = refbam.files,
                          include.chr = FALSE,
                         referenceFasta = reference.file)
save(my.Xrefcounts,file = "/hpc/cog_bioinf/diagnostiek/projects/WES_CNV/NOVASEQ/my.NovaS_Xrefcounts")
