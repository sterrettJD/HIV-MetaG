# Converts metaphlan formatted taxonomy names to a taxonomy table for phyloseq
names_to_tax_table <- function(bugs){
    splitted <- strsplit(bugs, split="|", fixed=T)
    # create empty taxonomy matrix
    taxmat <- matrix(NA, 
                     ncol=max(sapply(splitted, length)), 
                     nrow=length(splitted))
    colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
    
    # add split taxonomy to the matrix
    for (i in 1:nrow(taxmat)){
        tax.resolution <- length(splitted[[i]])
        taxmat[i, 1:tax.resolution] <- splitted[[i]]
    }
    # remove the p__, f__, etc to indicate level
    taxmat <- gsub("[a-z]__", "", taxmat)
    
    return(taxmat)
}
