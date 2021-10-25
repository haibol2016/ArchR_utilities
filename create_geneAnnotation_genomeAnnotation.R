library("GenomicFeatures")
library("BSgenome.Hsapiens.Ensembl.GRCh38")
library("plyranges")
library("ArchR")
library("tibble")
library("dplyr")
library("biomaRt")
library("collections")
library("BSgenome")

# make TxDb from a GTF file
forgeTxDb <- function(BSgenome, gtf, out_TxDb_dir)
{
    if (missing(BSgenome) | missing(gtf) | missing(out_TxDb_dir))
    {
        stop("All arguments are required: BSgenome, gtf, and out_TxDb_dir!")
    }
    if (!is(BSgenome, "BSgenome")){
        stop(BSgenome, " is not a BSgenome object!")
    }
    if (!file.exists(gtf)){
        stop(gtf, " doesn't exist!")
    }
    if (grepl(".gtf.gz$", gtf))
    {
        in_gtf <- gzfile(gtf, open ="rt")
    } else if (grepl(".gtf$", gtf)) {
        in_gtf <- file(gtf, open = "r")
    } else {
        stop("It seems the GTF file is not a GTF file ",
             "which should with an extension .gtf, or .gtf.gz")
    }
    if (!dir.exists(out_TxDb_dir)){
        dir.create(out_TxDb_dir)
    }
    chrom_len <- seqlengths(BSgenome)
    is_circular <- names(chrom_len) %in% c("chrM", "chrMT", "MT",
                                           "chrPltd", "Pltd")
    chrominfo <- data.frame(chrom = names(chrom_len),
                            length = unname(chrom_len),
                            is_circular = is_circular)
    genome_metadata <- metadata(BSgenome)
    TxDb <- makeTxDbFromGFF(file = in_gtf,
                            format = "gtf",
                            dataSource = genome_metadata$provider,
                            organism = genome_metadata$genome_metadata,
                            taxonomyId = NA,
                            chrominfo = chrominfo,
                            miRBaseBuild = NA)
    close(gtf)
    TxDb_file <- file.path(out_TxDb_dir, 
                           paste0(BSgenome, "TxDb.sqlite"))
    saveDb(TxDb, file = TxDb_file)
    TxDb_file
}


## get gene Symbol from the GTF file if exists, otherwise get the gene symbols
## from biomart
get_geneID_symbol <- function(gtf, species_latin_name)
{
    if (!file.exists(gtf)){
        stop(gtf, " doesn't exist!")
    }
    if (grepl(".gtf.gz$", gtf))
    {
        in_gtf <- gzfile(gtf, open ="r")
    } else if (grepl(".gtf$", gtf)) {
        in_gtf <- file(gtf, open = "r")
    } else {
        stop("It seems the GTF file is not a GTF file ",
             "which should with an extension .gtf, or .gtf.gz")
    }
    require("collections")
    id2symbol_dict <- ordered_dict()
    gtf <- read.delim(in_gtf, header = FALSE, as.is = TRUE,
                      comment.char = "#", quote = "")
    close(in_gtf)
    gtf <- gtf[gtf[, 3] == "gene", 9]
    null <- lapply(gtf, function(.x){
        if (grepl("gene_id", .x)){
            gene_id <- gsub('gene_id\\s+"([^"]+).+', 
                            "\\1", .x, perl = TRUE)
            if (grepl("gene_name", .x))
            {
                gene_symbol <- gsub('.+?gene_name\\s+"([^"]+).+', 
                                    "\\1", .x, perl = TRUE)
            } else {
                gene_symbol <- "NA"
            }
            id2symbol_dict$set(gene_id, gene_symbol)
        }
    })
    
    if (all(unlist(id2symbol_dict$values()) == "NA"))
    {
        message("The GTF file contains no gene symbols. ", 
                "Query Biomart to get gene symbols")
        require("biomaRt")
        species <- tolower(gsub("^(.).+\\W(.+)", "\\1\\2",
                              species_latin_name, perl = TRUE))
        
        ## try different biomart: animal, plant, fungi, metazona
        hosts <-gene_ids <- unlist(id2symbol_dict$keys())
        id_type <- {if (grepl("^ENS.*?G\\d+", gene_ids[1])) {"ensembl_gene_id"} 
            else if (grepl("^\\d+$", gene_ids[1])){"entrezgene_id"} 
            else {"unknown"}}
        if (id_type == "unknow")
        {
            stop("Unkow gene ID type!")
        }
        
        
        hosts <- c("https://www.ensembl.org/", "https://plants.ensembl.org/",
                   "https://fungi.ensembl.org/", "https://metazoa.ensembl.org/")
        
        hosts <- sapply(hosts, function(.x){
            listMarts(host = .x)$biomart[1]
        })
        
        id_symbol <- data.frame()
        for (i in seq_along(hosts))
        {
            ensembl <- useEnsembl(biomart = hosts[i], host = names(hosts)[i])
            datasets <- searchDatasets(ensembl, pattern = species)$dataset
            is_ds <- grepl(paste0(species, "_"), datasets)
            if (any(is_ds))
            {
                dataset <- datasets[is_ds]
                ensembl <- useMart(biomart = hosts[i],
                                   dataset = dataset,
                                   host = names(hosts)[i])
                id_symbol <- select(ensembl, keys = gene_ids,
                       columns = c(id_type,'external_gene_name'),
                       keytype = id_type)
                
                unnamed_geneID_symbol <- 
                    data.frame(id = gene_ids[!gene_ids %in% id_symbol[ ,1]], 
                               external_gene_name = "NA")
                colnames(unnamed_geneID_symbol)[1] <- id_type
                id_symbol <- rbind(id_symbol, unnamed_geneID_symbol)
                id_symbol$external_gene_name <- paste(id_symbol[, 2],
                                                      id_symbol[, 1],
                                                      sep = "_")
                break
            }
        }
        if (nrow(id_symbol) > 0){
            colnames(id_symbol) <- c("gene_id", "symbol")
            return(id_symbol)
        } else {
            message("No gene symbols are get from Biomart!\n",
                    "Check your GTF file!")
            return(NULL)
        }
    } else {
        id2symbol <- data.frame(gene_id = unlist(id2symbol_dict$keys()),
                                     symbol = unlist(id2symbol_dict$values()))
        id2symbol$symbol <- paste(id2symbol$symbol,
                                       id2symbol$gene_id,
                                       sep = "_")
        return(id2symbol)
    }
}

create_ArchR_geneannotation_WO_OrgDb <- function(TxDb = NULL, 
                                                 geneID2Symbol,
                                                 out_dir)
{
    if (is.null(TxDb) || !is(TxDb, "TxDb"))
    {
        stop(TxDb, " must ba a TxDb object!")
    }
    if (missing(geneID2Symbol) || missing(out_dir))
    {
        stop("geneID2Symbol and out are required!")
    }
    if (!is.data.frame(geneID2Symbol) || 
        all(colnames(geneID2Symbol) == c("gene_id", "symbol")))
    {
        stop(geneID2Symbol, 
             "must be a dataframe with colnames:'gene_id', 'symbol'!")
    }
    if (!dir.exists(out_dir))
    {
        dir.create(out_dir)
    }
    ## geneID2Symbol: col1: gene_id; col2: gene_symbol
    symbol <- geneID2Symbol[, 2]
    names(symbol) <- geneID2Symbol[, 1]
    
    ## GRanges for genes
    genes <- GenomicFeatures::genes(TxDb) %>%
        plyranges::remove_names() %>% 
        plyranges::mutate(symbol = symbol[gene_id]) 
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
    
    ## get all transcripts
    tx <- unlist(transcriptsBy(TxDb, by = "gene")) %>%
        plyranges::mutate(gene = names(.)) %>%
        plyranges::remove_names() %>%
        plyranges::select(-c("tx_id")) %>%
        data.frame()
    tx_gene <- tx$gene
    names(tx_gene) <-tx$tx_name
    rm("tx")
    
    ## Create GRanges for exons
    exons <- unlist(exonsBy(TxDb,
                            by = "tx",
                            use.names = TRUE)) %>%
        plyranges::mutate(tx_name = names(.)) %>%
        plyranges::remove_names() %>%
        plyranges::mutate(gene_id = tx_gene[tx_name]) %>%
        plyranges::filter(!is.na(gene_id)) %>%
        plyranges::mutate(symbol = symbol[gene_id]) %>%
        plyranges::select(-c("exon_id", "exon_name", "exon_rank", 
                             "gene_id", "tx_name"))
    exons <- sort(sortSeqlevels(exons), ignore.strand = TRUE)
    
    ## Create GRanges for TSS
    TSS <- unique(resize(GenomicFeatures::transcripts(TxDb), 
                         width = 1, 
                         fix = "start"))
    
    ## Create geneAnnotation
    geneAnnotation <- createGeneAnnotation(
        genome = NULL,
        TxDb = NULL,
        OrgDb = NULL,
        genes = gene,
        exons = exon,
        TSS = tss,
        annoStyle = "ENSEMBL")
    saveRDS(geneAnnotation, file = file.path(out_dir, "geneAnnotation.RDS"))
    geneAnnotation
}

## create genomeAnnotation using a BSgenome and an optional blacklist in the 
## BED format
create_ArchR_genomeannotation <- function(BSgenome, out_dir, 
                                          blacklist_bed = NULL, 
                                          blacklist_hasheader = FALSE)
{
    if (missing(BSgenome) | is.null(BSgenome) | !is(BSgenome, "BSgenome"))
    {
        stop(BSgenome, " must be a BSgenomeobject!")
    }
    if (missing(out_dir))
    {
        stop("out_dir is required!")
    }
    chrom_len <- seqlengths(BSgenome)
    chr_df <- data.frame(seqnames = names(chrom_len),
                         start = 1,
                         end = unname(chrom_len))
    chromSizes <- makeGRangesFromDataFrame(chr_df)
    
    if (!is.null(blacklist_bed))
    {
        if (grepl(".bed.gz$", blacklist_bed))
        {
            blacklist <- gzfile(gtf, open ="rt")
        } else if (grepl(".bed$", gtf)) {
            blacklist <- file(gtf, open = "r")
        } else {
            stop("It seems the Blacklist file is not a bed file ",
                 "which should with an extension .bed, or .bed.gz")
        }
        blacklist_df <- read.delim(blacklist, 
                                   header = blacklist_hasheader)
        if (ncol(blacklist_df) < 3 || 
            !any(blacklist_df[,1] %in% names(chrom_len)) || 
            !any(is.numeric(blacklist_df[,2])) || 
            !any(is.numeric(blacklist_df[,3])))
        {
            stop(blacklist_bed, " is not a valid BED file!")
        }
        
        colnames(blacklist_df)[1:3] <- c("seqnames", "start", "end")
        blacklist <- makeGRangesFromDataFrame(blacklist_df, 
                                              starts.in.df.are.0based = TRUE)
    }
    genomeAnnotation <- createGenomeAnnotation(genome = BSgenome,
                                               chromSizes = chromSizes,
                                               blacklist = blacklist,
                                               filter = TRUE,
                                               filterChr = c("chrM"))
    saveRDS(genomeAnnotation, file = file.path(out_dir, 
                                               "genomeAnnotation.RDS"))
    genomeAnnotation
}


## call the function
# gene_symbol <- 
#    get_geneID_symbol(gtf = "refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz", 
#                      species_latin_name = "Homo sapiens")
