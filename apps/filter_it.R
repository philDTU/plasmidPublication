library(tidyverse)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
#file <- args[1]
#out <- args[2]
#margin <- args[3]


file <- args[1]
out <- args[2]
margin <- 50


check_overlaps <- function(start,end,margin){
    oldEnd <- 0
    res <- c()
    for (row in 1:length(start)) {
        startPoint <- start[row]
        endPoint <- end[row]
        if(oldEnd <= startPoint + margin){
            res <- c(res, FALSE) # str_c(FALSE,startPoint-oldEnd, sep= " "))
        }
        else{
            res <- c(res, TRUE)# str_c(TRUE,startPoint-oldEnd ,sep= " "))
        }
        oldEnd <- end[row]
    }
    return(res)
}


filter_overlaps <- function(file,margin){
    df <- read_tsv(file, col_names = FALSE) %>% dplyr::select(X1,X2,X3,X4)
    colnames(df) <- c("read_id","start","end","contig")
    df$sample <- tools::file_path_sans_ext(basename(file))
    df <- df %>% group_by(read_id) %>% mutate(is_overlaped = check_overlaps(start,end,margin))
    df <-df %>% filter(!any(is_overlaped))
    return(df)
}

df <- filter_overlaps(file,margin)

write_tsv(df,out)
