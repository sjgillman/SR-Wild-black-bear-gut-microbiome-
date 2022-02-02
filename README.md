# Wild black bears harbor simple gut microbial communities with little difference between the jejunum and colon (2020)
#### Published Manuscript in Scientific Reports; [DOI: 10.1038](https://doi.org/10.1038/s41598-020-77282-w)
Authors: Sierra J. Gillman, Erin A. McKenney, Diana J.R. Lafferty
All infromation can also be downloaded from [zenodo](https://zenodo.org/record/4060480#.Yfm2x_XMJhE)

**Abstract:**
The gut microbiome (GMB), comprising the commensal microbial communities located in the gastrointestinal tract, has co-evolved in mammals to perform countless micro-ecosystem services to facilitate physiological functions. Because of the complex inter-relationship between mammals and their gut microbes, the number of studies addressing the role of the GMB on mammalian health is almost exclusively limited to human studies and model organisms. Furthermore, much of our knowledge of wildlife–GMB relationships is based on studies of colonic GMB communities derived from the feces of captive specimens, leaving our understanding of the GMB in wildlife limited. To better understand wildlife–GMB relationships, we engaged hunters as citizen scientists to collect biological samples from legally harvested black bears (Ursus americanus) and used 16S rRNA gene amplicon sequencing to characterize wild black bear GMB communities in the colon and jejunum, two functionally distinct regions of the gastrointestinal tract. We determined that the jejunum and colon of black bears do not harbor significantly different GMB communities: both gastrointestinal sites were dominated by Firmicutes and Proteobacteria. However, a number of bacteria were differentially enriched in each site, with the colon harboring twice as many enriched taxa, primarily from closely related lineages.

Directory structure | Description
--- | ---
blackbear-gme-SR/
  README.md
  **data/** | **Description**
  *BearMeta-R.tsv* | data used in R analysis
  *BearMeta.tsv* | data used in QIIME2-- just slightly different format from above
  *BlackBeardemuxsequences.tar.gz* | demultiplexed EMP-paired end sequences demultiplexed on QIIME2
  *physeq.rds* | phyloseq-R object that can be used if wanting to skip rarifying step
  **scripts/** | **Description**
  *QIIME2 Pipeline.md* | bioinformatic pipeline to prepare sequences for analysis in R
  *Statistical Analysis.R* | code required to repeat statistics from manuscript
  **images/**
  *blackbear.png*
  *protocol.pdf* | insctructions on how hunters collected samples

<p align="center">
<img src="images/blackbear.png" width="500" />
  </p>


