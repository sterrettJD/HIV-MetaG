library(stringi)
library(data.table)
library(tidyverse)
library(r2r)

filter_only_genes <- function(df, with_taxa=F){
    feature.names <- rownames(df)
    # if there's a | in the rowname, it contains a corresponding species
    contains.taxa <- grepl(x=feature.names, pattern="\\|")
    
    if (with_taxa){
        return(gn[contains.taxa,])
    } else {
        return(gn[!contains.taxa,])
    }
}


get_ec_class_db <- function(){
    if(!file.exists("enzclass.txt")){
        download.file("https://ftp.expasy.org/databases/enzyme/enzclass.txt",
                      destfile="enzclass.txt")
        
    } 
    
    enzclass <- read_file("enzclass.txt")
    # split by lines and remove header/footer
    splitted <- str_split(enzclass, pattern="\n")[[1]][12:418]
    # split EC from the corresponding name
    parsed <- sapply(splitted, 
           FUN=function(x) str_split(x, pattern="  ") %>% 
               unlist() %>%
               stri_remove_empty()
           )
    df <- matrix(parsed, nrow=2) %>% t() %>% as.data.frame()
    colnames(df) <- c("EC", "name")
    # trim ECs for hierarchical mapping
    df$EC.trimmed <- df$EC %>% 
        sapply(FUN=function(x) str_split(x, pattern="-")[[1]][1])
    
    # create a hashmap to map keys and values
    mapper <- hashmap()
    mapper[df$EC.trimmed] <- df$name
    return(mapper)
}


