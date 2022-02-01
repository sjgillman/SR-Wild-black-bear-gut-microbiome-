
### Bioinformatic pipeline for Gillman et al 2020

#### QIIME2 verion: 2019.4
Samples are from Argonne National Laboratory
Pipeline adapted from qiime2 tutorial "Atacama soil microbiome" & "Moving Pictures"
samples are EMP-Paired end multiplexed sequences with new primer set 
w/ barcodes read forward & no longer reversed in demux step

**Import data into QIIME2**
```
## import sequences
qiime tools import \
--type EMPPairedEndSequences \
--input-path Reads \
--output-path paired-end-sequences.qza #you can name this whatever you want
##output artifact: paired-end-sequences.qza
```

**Demultiplexing Sequences**
```
## you will need metadata/mapping file

qiime demux emp-paired \
--m-barcodes-file Meta.tsv \
--m-barcodes-column BarcodeSequence \
--p-no-golay-error-correction \
--i-seqs paired-end-sequences.qza \
--o-per-sample-sequences demuxseq.qza \
--o-error-correction-details demux-detail.qza

##  make a summary visualization 

qiime demux summarize \
--i-data demuxseq.qza \ 
--o-visualization demuxseq.qzv
```
with the .qzv file we will go to qiime2view online and look at the quality of our reads

**Denoising sequences with DADA2 plugin**
prior to denoising, look at demux.qzv to determine if/where to trim sequences
you will need to have r installed in your qiime2 environment
if R is not installed be sure to be in the qiime2 environment
```
## denoising

qiime dada2 denoise-paired \
--i-demultiplexed-seqs demuxseq.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 150 \
--p-trunc-len-r 150 \
--o-table table.qza \
--o-representative-sequencesrep-seqs.qza \
--o-denoising-stats denoising-stats.qza

## output artifacts: table.qza, rep-seqs.qza, denoising-stats.qza
# you will now have artifacts containing the 
# feature table and corresponding feature sequences.
# You can generate summaries of those as follows

##  summary visulization table for determining sample depth for rarifying

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file Meta.tsv
#output visualization: table.qzv
#sampling depth:18257
## remember this is before we remove contaminates

#  make visualization artifacts of rep seq

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv
#output visualization: rep-seq.qzv

## view denoising stats

qiime metadata tabulate \
--m-input-file denoising-stats.qza \
--o-visualization denoising-stats.qzv
#output visualization: denoising-stats.qzv
```
__now we can start using the "moving pictures tutorial" starting at *Taxonomic Analysis sklearn*__
```
## Training classifier
# https://docs.qiime2.org/2019.4/data-resources/
# https://docs.qiime2.org/2019.4/tutorials/feature-classifier/
## we use the SILVA reference database 515/806
# https://www.arb-silva.de/download/archive/qiime
## download SILVA_132_ or whatever the newest version is.
## put SILVA folder in Projects

## Import reference otus

qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path SILVA_132_99_16S.fna \
--output-path SILVA_OTU.qza

## Import reference taxonomy file
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path taxonomy_7_levels.txt \
--output-path ref-taxonomy.qza

##  Extract reference reads
qiime feature-classifier extract-reads \
--i-sequences SILVA_OTU.qza \
--p-f-primer GTGCCAGCMGCCGCGGTAA \
--p-r-primer GGACTACHVGGGTWTCTAAT \
--p-trunc-len 150 \
--p-min-length 100 \
--p-max-length 400 \
--o-reads ref-seqs.qza

## Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy ref-taxonomy.qza \
--o-classifier -classifier.qza

## Test Classifier
qiime feature-classifier classify-sklearn \
--i-classifier -classifier.qza \
--i-reads rep-seqs.qza \
--o-classification taxonomySILVA.qza

#####  fixing white spaces
qiime tools export \
--input-path taxonomySILVA.qza \
--output-path taxonomy-with-spaces

qiime metadata tabulate \
--m-input-file taxonomy-with-spaces/taxonomy.tsv  \
--o-visualization taxonomy-as-metadata.qzv

qiime tools export \
--input-path taxonomy-as-metadata.qzv \
--output-path taxonomy-as-metadata

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-path taxonomy-as-metadata/metadata.tsv \
--output-path -taxonomy-without-spaces.qza


## create visualization
qiime metadata tabulate \
--m-input-file taxonomy-without-spaces.qza \
--o-visualization taxonomySILVA.qzv
```

****Filtering 
```
##  filter out mitochondria and chloroplast
qiime taxa filter-table \
--i-table table.qza \
--i-taxonomy taxonomy-without-spaces.qza \
--p-exclude mitochondria \
--o-filtered-table table-filter.qza

qiime taxa filter-table \
--i-table table-filter.qza \
--i-taxonomy taxonomy-without-spaces.qza \
--p-exclude chloroplast \
--o-filtered-table clean-table.qza

## get rid of unassigned
qiime taxa filter-table \
--i-table clean-table.qza \
--i-taxonomy taxonomy-without-spaces.qza \
--p-exclude Unassigned \
--o-filtered-table clean-table-unassigned-rm.qza


## remove Bacteria only assigned
qiime taxa filter-table \
--i-table clean-table-unassigned-rm.qza \
--i-taxonomy taxonomy-without-spaces.qza \
--p-mode exact \
--p-exclude D_0__Bacteria \
--o-filtered-table clean-table-unassigned_Unknown-rm.qza


## remove Arch only assigned
qiime taxa filter-table \
--i-table clean-table-unassigned_Unknown-rm.qza \
--i-taxonomy taxonomy-without-spaces.qza \
--p-exclude D_0__Archaea \
--o-filtered-table clean-table-unassigned_Unknown_Arch-rm.qza

#  barplot
qiime taxa barplot \
--i-table clean-table-unassigned_Unknown_Arch-rm.qza \
--i-taxonomy taxonomy-without-spaces.qza \
--m-metadata-file Meta.tsv \
--o-visualization taxa-bar-plotsSILVA-clean2.qzv

##  determine depth
qiime feature-table summarize \
--i-table clean--table-unassigned_Unknown_Arch-rm.qza \
--o-visualization clean--table-unassigned_Unknown_Arch-rm.qzv \
--m-sample-metadata-file Meta.tsv
##1050
```

**Generating a tree for phylogenetic diversity analyses with clean data**
```
## Filter
qiime feature-table filter-seqs \
--i-data rep-seqs.qza \
--i-table clean-table-unassigned_Unknown_Arch-rm.qza \
--o-filtered-data filtered-rep-seq.qza

## root
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences filtered-rep-seq.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree -filterd-unrooted-tree.qza \
--o-rooted-tree filtered-rooted-tree.qza
```



