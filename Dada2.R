####Preparing the Data####
setwd("Working/Dir")
library("dada2")
library("devtools")
library("dplyr")
library("ggplot2")
library("microbiome")
library("phangorn")
library("phyloseq")
library("Rcpp")
library("reshape2")
library("tidyr")
library("vegan")

##
path = "Reads"

fnFs = sort(list.files(path, pattern="_1.fq.gz"))
fnRs = sort(list.files(path, pattern="_2.fq.gz"))

sample.names = sapply(strsplit(fnFs, "_"), `[`, 1)
show(sample.names)

fnFs = file.path(path, fnFs)
fnRs = file.path(path, fnRs)

filt_path = file.path(path, "filtered")

filtFs = file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs = file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

dev.off()

####Filter and Trimming####
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,160), trimLeft = c(19,20), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)


####Error Correctio Derreplication n####
errF = learnErrors(filtFs, nbases = 1e8, multithread=TRUE, randomize=TRUE, MAX_CONSIST=15)
#100206400 total bases in 622400 reads from 30 samples will be used for learning the error rates.

errR = learnErrors(filtRs, nbases = 1e8, multithread=TRUE, randomize=TRUE, MAX_CONSIST=15)
#101696000 total bases in 726400 reads from 35 samples will be used for learning the error rates.


mergers = vector("list", length(sample.names))
names(mergers) = sample.names
names(filtFs) = sample.names
names(filtRs) = sample.names


plotErrors(errF, nominalQ=TRUE)


for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF = derepFastq(filtFs[[sam]])
  ddF = dada(derepF, err=errF, multithread=TRUE)
  derepR = derepFastq(filtRs[[sam]])
  ddR = dada(derepR, err=errR, multithread=TRUE)
  merger = mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] = merger
}

rm(derepF); rm(derepR)

####Merge Pairs####
seqtab = makeSequenceTable(mergers)
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 9374 bimeras out of 9676 input sequences.


####Exporting ASVs####
getN = function(x) sum(getUniques(x))
track2 = cbind(out, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track2) <- c("input", "filtered", "merged", "tabled", "nonchim")
write.table(track2, "track.txt", sep = "\t", quote = F)
uniquesToFasta(getUniques(seqtab.nochim), "seqtab.nochim.fasta", ids=paste0("sq", seq(length(getUniques(seqtab.nochim)))))
export <- t(seqtab.nochim)
rownames(export) = paste0("sq", seq(length(getUniques(seqtab.nochim))))
export <- cbind('#OTUID' = rownames(export), export)
write.table(export, "otu_table.dada.nochim.txt", sep='\t', row.names=FALSE, quote=FALSE)
write.table(taxa, "taxa_table.dada.nochim.txt", sep='\t', row.names=FALSE, quote=FALSE)


sum(seqtab.nochim)/sum(seqtab)
#[1] 0.8130383


####Assign taxonomy - SILVA####
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa", multithread=TRUE)

taxa <- addSpecies(taxa, allowMultiple = TRUE, "silva_species_assignment_v132.fa")


#inspect the taxonomic assignments#
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#### giving seq headers more manageable names (ASV_1, ASV_2...)####
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}




asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

#Construct phylogenetic tree
seqs <- getSequences(taxa)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

library(phangorn)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

####Phyloseq Object####


sample_info_tab <- read.table("sample_info_shime.tsv", header=T, row.names=1,
                              check.names=F, sep="\t")


Twin_Shime_Longitudinal_Study <- phyloseq(tax_table(taxa), sample_data(sample_info_tab),
                                       otu_table(seqtab, taxa_are_rows = FALSE),phy_tree(fitGTR$tree))

taxa_names(Twin_Shime_Longitudinal_Study) <- paste0("Seq", seq(ntaxa(Twin_Shime_Longitudinal_Study)))


#########################################################################################

summarize_phyloseq(Twin_Shime_Longitudinal_Study)
pseq <- Twin_Shime_Longitudinal_Study

#Pick metadata as data.frame:
meta <- meta(pseq)
#Renomeando Tax
taxa_names(pseq) <- paste0("Seq", seq(ntaxa(pseq)))
#Taxonomy table:
taxonomy <- tax_table(pseq)
#Abundances for taxonomic groups (‘OTU table’) as a TaxaxSamples matrix:
# Absolute abundances
otu.absolute <- abundances(pseq)
# Relative abundances
otu.relative <- abundances(pseq, "compositional")
##Tree
phylogenetic<- phy_tree(pseq)

###Read Filter
##FALTA essa
track2 <-read.table("track.txt", header=T,check.names=F, sep="\t")
ma <- data.matrix(track2)
barplot((ma), main ="Track Reads", font.axis = 1, cex.axis=1 ) 
barplot2(ma)