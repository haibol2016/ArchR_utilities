library("rtracklayer")
library("plyranges")
library("GenomicRanges")
library("future.apply")

## Amend UCSC GTF files by adding gene entries
add_gene2gtf <- function(gtf, out_dir = ".", gz_out = TRUE)
{
    if (!file.exists(gtf)){
        stop("gtf doesn't exist!")
    }
    gtf_GR <- rtracklayer::import(gtf)
    if (!dir.exists(out_dir)){
        dir.create(out_dir, recursive = TRUE)
    }
    
    ## split into GRangesList, suppose all gene_ids are unique
    gtf_GRL <- split(gtf_GR, mcols(gtf_GR)$gene_id)
    
    ## using parallel computing
    plan(multisession)
    fixed_gtf <- do.call("rbind", future_lapply(gtf_GRL, function(.x){
        tx_GR <- .x[.x$type == "transcript"]
        gene <- reduce(tx_GR, ignore.strand = FALSE)
        gene$type <- "gene" 
        tx_GR1 <- tx_GR[1]
        gene$source <- tx_GR1$source
        gene$gene_id <- tx_GR1$gene_id
        gene$gene_name <- tx_GR1$gene_name
        data.frame(append(gene, .x))
    })) %>% plyranges::as_granges()
    
    out_filename <- paste0("amended.", Sys.Date(),".", 
                           gsub(".(gz|bz2|xz)$", "", basename(gtf)))
    if (gz_out)
    {
        out_filename <- file.path(out_dir, 
                                  paste0(out_filename, ".gz"))
        gtf_file <- gzfile(out_filename, open = "w")
    } else {
        out_filename <- file.path(out_dir, 
                                  out_filename)
        gtf_file <- file(out_filename, open = "w")
    }

    rtracklayer::export(fixed_gtf, con = gtf_file)
    close(gtf_file)
    out_filename
}
