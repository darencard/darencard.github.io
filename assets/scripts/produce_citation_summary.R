#!/usr/bin/env Rscript --vanilla --default-packages=base,utils,stats

require(rcrossref)
require(scholar)
require(tidyverse)
require(colorblindr)
require(gridExtra)
require(grid)

args<-commandArgs(TRUE)

produce_citation_summary <- function(doi_intable, out_imgpath) {
  # retrieve citation info for Google Scholar
  google_scholar_id <- "umOwsMAAAAAJ"
  
  scholar_profile <- get_profile(google_scholar_id)
  
  # retrieve citation info for CrossRef
  dcc_doi <- read.table(doi_intable, header=TRUE, stringsAsFactors=FALSE)
  
  get_crossref_profile <- function(doi_tbl) {
    cite_vec <- NULL
    for (row in 1:nrow(dcc_doi)) {
      doi <- dcc_doi[row,1]
      cite_num <- cr_citation_count(doi)
      cite_vec <- c(cite_vec, cite_num)
    }
    
    dcc_doi$Cites <- cite_vec
    
    citations <- dcc_doi[order(dcc_doi$Cites,decreasing=TRUE),]
    citations <- cbind(id=rownames(citations),citations)
    citations$id<- as.character(citations$id)
    citations$id<- as.numeric(citations$id)
    
    pub_count <- nrow(citations)
    total_cites <- sum(citations$Cites)
    i10_index <- citations %>% filter(Cites >= 10) %>% nrow()
    h_index <- max(which(citations$id<=citations$Cites))
    
    profile <- list(
      "pub_count"=pub_count,
      "total_cites"=total_cites,
      "i10_index"=i10_index,
      "h_index"=h_index
    )
    
    return(profile)
  }
  
  crossref_profile <- get_crossref_profile(dcc_doi)
  
  # create full citation table
  citation_summary_table <- data.frame(
    "Google Scholar"=c(scholar_profile$total_cites,
                       scholar_profile$i10_index,
                       scholar_profile$h_index),
    "CrossRef"=c(crossref_profile$total_cites,
                 crossref_profile$i10_index,
                 crossref_profile$h_index),
    row.names=c("Total Citations", "i10 Index", "H Index")
  )
  
  citation_table <- tableGrob(citation_summary_table, 
                              cols=c("Google Scholar", "CrossRef"),
                              theme=ttheme_minimal(base_size=10))
  
  # create citation plot
  roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
  }
  
  cite_history <- get_citation_history(google_scholar_id)
  
  upper <- roundUpNice(max(cite_history$cites))
  
  citation_plot <- cite_history %>%
    ggplot(aes(x=year, y=cites)) + 
    geom_bar(stat="identity", fill="#0395DE") +
    geom_text(aes(label=cites), 
              position=position_dodge(width=0.9), vjust=-0.5,
              size=4) +
    scale_y_continuous(limits=c(0, max(cite_history$cites)+10)) +
    labs(x="Year", y="Citations") +
    theme_minimal(base_size=16) +
    theme(panel.grid=element_blank(), axis.text.y=element_blank(),
          axis.title=element_blank(), 
          axis.text.x=element_text(color="black", vjust=5))
  
  # put it all together
  png(out_imgpath, width=6, height=2.5, units="in", res=600)
  grid.arrange(citation_plot, citation_table, nrow=1)
  dev.off()
}

produce_citation_summary(args[1], args[2])