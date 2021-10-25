library("AnnotationForge")
create_install_OrgDb <- function(geneID_symbom_df, 
                         tax_id = "9606", 
                         genus = "Homo",
                         species = "sapeins",
                         package_dest_dir = "./")
{
    if (!is.data.frame(geneID_symbom_df) || 
        any(colnames(geneID_symbom_df) != c("gene_id", "symbol")))
    {
        stop(geneID_symbom_df, " is not a dataframe or the colnames are not
             correct!")
    }
    id_type <- {if (all(grepl("^ENS.*?G\\d+", geneID_symbom_df[, 1], perl = TRUE)))
    {"ENAEMBL"} else if (all(grepl("^ENS.*?G\\d+", 
                                       geneID_symbom_df[, 1], perl = TRUE)))
    {"GENEID"} else {"unkown"}}
    
    if (id_type == "unkown"){
        stop("gene_id should be either Ensembl gene ID or Entrez Gene ID!")
    }
    
    fSym <- unique(geneID_symbom_df[, c(1,3)])
    colnames(fSym) <- c("GID", "SYMBOL")
    if (all(grepl("_", fSym$SYMBOL, perl= TRUE)))
    {
          fSym$SYMBOL <- gsub("_.+", "", fSym$SYMBOL, perl = TRUE)
    }
    
    ensembl <- unique(geneID_symbom_df[, c(1,1)])
    colnames(ensembl) <- {if (id_type == "ENSEMBL") {c("GID", "ENSEMBL")} else {
        c("GID", "GENEID")
    }}
    if (!dir.exists(package_dest_dir))
    {
        dir.create(package_dest_dir)
    }
    
    makeOrgPackage(gene_info = fSym, 
                   ensembl_trans = ensembl_trans,
                   ensembl = ensembl,
                   version = "0.1",
                   maintainer = "Some One <so@someplace.org>",
                   author = "Some One <so@someplace.org>",
                   outputDir = package_dest_dir,
                   tax_id= tax_id,
                   genus= genus,
                   species= species,
                   goTable=NULL)
    
    package_name <- paste("org", paste0(gsub("^(.).+", "\\1", genus, perl = TRUE),
                           tolower(species)), "eg.db", sep = ".")
    
    install.packages(file.path(tmpdir, package_name), 
                     type = "source", repos=NULL)
}
