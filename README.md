# ArchR_utilies
R functions for a few common purposes:  
- amend UCSC GTF files by adding entries for genes
- split a multi-fasta file for a genome assembly into per-chromosome single fasta files
- automatically generate a seed file for forging a BSgenome package
- forge a BSgenome package from a directory containing single fasta files for individual chromosomes of a genome assembly using the seed file
- check, build and install the BSgenome package
- forge a TxDb database from a GTF file
- get gene IDs and gene symbols from the GTF file or query Biomart databases
- create geneAnnotation and genomeAnnotation for ArchR-based analysis of scATAC-seq data
- forge a OrgDb package from a dataframe containing gene IDs and gene symbols
