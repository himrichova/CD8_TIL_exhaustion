#! /usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --distribution=block
#SBATCH --mem-per-cpu=64000
#SBATCH --time=14-00:00:00
#SBATCH --partition=longq
#SBATCH --error=output_%j.err
#SBATCH --output=output_%j.out


################################################################################
# (0) Preliminaries
################################################################################
# load the package
library(RnBeads)

# define the directory structure
setwd("/scratch/lab_bsf/projects/BSA_0296_Mouse_Intratumor_CD8/")
dataDir <- file.path(getwd(), "extract")
resultDir <- file.path(getwd(), "RnBeads")

# dataset and file locations
datasetDir = dataDir
bedDir = dataDir
sampleSheet <- file.path(resultDir, "BSA_0296_Mouse_Intratumor_CD8_RnBeads_samples.csv")
reportDir <- file.path(resultDir, "report_new")

options(fftempdir=resultDir)
options(ffcaching="ffeachflush")



# enable parallel computation for critical steps in the analysis
# in the back, RnBeads uses the "foreach" and "doParallel" packages
#parallel.isEnabled()
num.cores <- 4
parallel.setup(num.cores)
#parallel.isEnabled()
#parallel.getNumWorkers()

################################################################################
# (1) Set analysis options
################################################################################

rnb.options(
	analysis.name                                = "RnBeads analysis of DNA methylation project: BSA_0296_Mouse_Intratumor_CD8",
	logging                                      = TRUE,
	email                                        = NULL,
	assembly                                     = "mm10",
	columns.pairing                              = NULL,
	analyze.sites                                = FALSE,
	region.types                                 = c("tiling", "genes", "promoters", "cpgislands"),
	region.aggregation                           = "mean",
	region.subsegments                           = 0,
	region.subsegments.types                     = NULL,
	identifiers.column                           = "sample_name",
	min.group.size                               = 1,
	max.group.count                              = NULL,
	gz.large.files                               = TRUE,
	strand.specific                              = FALSE,
	replicate.id.column                          = NULL,
	import                                       = TRUE,
	import.table.separator                       = ",",
	import.bed.style                             = "Encode",
	import.bed.frame.shift                       = 1,
	import.bed.test                              = TRUE,
	import.bed.test.only                         = FALSE,
	import.skip.object.check                     = FALSE,
	import.idat.chunk.size                       = NULL,
	import.gender.prediction                     = TRUE
)


rnb.options(
	qc                                           = TRUE,
	qc.boxplots                                  = TRUE,
	qc.barplots                                  = TRUE,
	qc.negative.boxplot                          = TRUE,
	qc.snp.heatmap                               = FALSE,
	qc.snp.distances                             = FALSE,
	qc.snp.boxplot                               = FALSE,
	qc.snp.barplot                               = FALSE,
	qc.coverage.plots                            = FALSE,
	qc.coverage.threshold.plot                   = 1:20,
	qc.coverage.histograms                       = FALSE,
	qc.coverage.violins                          = FALSE,
	qc.sample.batch.size                         = 500,
	preprocessing                                = TRUE,
	normalization                                = NULL,
	normalization.method                         = "swan",
	normalization.background.method              = "methylumi.noob",
	normalization.plot.shifts                    = TRUE,
	filtering.whitelist                          = NULL,
	filtering.blacklist                          = NULL,
	filtering.snp                                = "3",
	filtering.cross.reactive                     = FALSE,
	filtering.greedycut                          = FALSE,
	filtering.greedycut.pvalue.threshold         = 0.05,
	filtering.greedycut.rc.ties                  = "row",
	filtering.sex.chromosomes.removal            = FALSE,
	filtering.missing.value.quantile             = 0.5,
	filtering.coverage.threshold                 = 1,
	filtering.low.coverage.masking               = TRUE,
	filtering.high.coverage.outliers             = TRUE,
	filtering.deviation.threshold                = 0,
	imputation.method                            = "none",
	inference                                    = FALSE,
	inference.reference.methylome.column         = NULL,
	inference.max.cell.type.markers              = 10000,
	inference.top.cell.type.markers              = 500,
	inference.sva.num.method                     = "leek",
	exploratory                                  = TRUE,
	exploratory.columns                          = c("sample_name","tissue","group","custom","custom2"),
	exploratory.top.dimensions                   = 1000000,
	exploratory.principal.components             = 8,
	exploratory.correlation.pvalue.threshold     = 0.01,
	exploratory.correlation.permutations         = 10000,
	exploratory.correlation.qc                   = TRUE,
	exploratory.beta.distribution                = TRUE,
	exploratory.intersample                      = FALSE,
	exploratory.region.profiles                  = character(0),
	exploratory.deviation.plots                  = NULL,
	exploratory.clustering                       = "all",
	exploratory.clustering.top.sites             = 1000,
	exploratory.clustering.heatmaps.pdf          = FALSE,
	exploratory.gene.symbols                     = NULL,
	exploratory.custom.loci.bed                  = NULL,
	differential                                 = TRUE,
	differential.site.test.method                = "limma",
	differential.permutations                    = 0,
	differential.comparison.columns              = c("tissue","group","custom","custom2"),
	differential.comparison.columns.all.pairwise = c("tissue","group","custom","custom2"),
	covariate.adjustment.columns                 = NULL,
	differential.adjustment.sva                  = TRUE,
	differential.adjustment.celltype             = TRUE,
	differential.enrichment.go                   = TRUE,
	differential.enrichment.lola                 = TRUE,
	differential.report.sites                    = TRUE,
    export.to.trackhub                           = c("bigBed"),
	export.to.bed                                = TRUE,
	export.to.csv                                = TRUE,
	export.to.ewasher                            = FALSE,
	export.types                                 = c("tiling", "genes", "promoters", "cpgislands")
)

rnb.options(
	logging.memory                               = TRUE,
	logging.disk                                 = TRUE,
	logging.exit.on.error                        = TRUE,
	distribution.subsample                       = 50000,
	disk.dump.big.matrices                       = TRUE,
	disk.dump.bigff                              = TRUE,
	disk.dump.bigff.finalizer                    = "delete",
	enforce.memory.management                    = TRUE,
	enforce.destroy.disk.dumps                   = TRUE
)


################################################################################
# (2) Run the analysis
################################################################################
rnb.run.analysis(
	dir.reports=reportDir,
	sample.sheet=sampleSheet,
	data.dir=bedDir,
	data.type="bs.bed.dir"
)
