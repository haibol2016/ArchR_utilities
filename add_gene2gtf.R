library("rtracklayer")

## Amend UCSC GTF files by adding gene entries
add_gene2gtf <- function(gtf, out_dir = ".")
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
    fixed_gtf_GRL <- endoapply(gtf_GRL, function(.x){
        gene <- reduce(.x, ignore.strand=FALSE)
        tx_GR <- .x[.x$type == "transcript"][1]
        mcols(gene) <- mcols(tx_GR)
        gene$type <- "gene" 
        gene$source <- tx_GR$source
        gene$gene_id <- tx_GR$gene_id
        gene$gene_name <- tx_GR$gene_name
        append(gene, this_list[[1]])
    })
    fixed_gtf <- unlist(fixed_gtf_GRL)
    gtf_file <- file.path(out_dir, 
                          paste0("amended.",Sys.Date(),".", basename(gtf)))
    rtracklayer::export(fixed_gtf, 
                        con = gtf_file)
    gtf_file
}
