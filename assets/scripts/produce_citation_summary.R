#!/usr/bin/env Rscript --vanilla --default-packages=base,utils,stats

require(rcrossref)
require(scholar)
require(tidyverse)
require(colorblindr)
require(gridExtra)
require(grid)

args<-commandArgs(TRUE)

produce_citation_summary <- function(doi_intable, out_imgpath, out_txtpath) {
  
  get_scholar_profile <- function(google_scholar_id) {
    citations <- get_publications(google_scholar_id, pagesize=500)
    
    pub_count <- citations %>% filter(pubid!="X5YyAB84Iw4C") %>% nrow()
    total_cites <- citations %>% filter(pubid!="X5YyAB84Iw4C") %>% pull(cites) %>% sum()
    i10_index <- citations %>% filter(pubid!="X5YyAB84Iw4C") %>% filter(cites >= 10) %>% nrow()
    h_index <- tail(which(citations$cites >= seq_along(citations$cites)), 1)
    mean <- citations %>% filter(pubid!="X5YyAB84Iw4C") %>% pull(cites) %>% mean(na.rm=TRUE) %>% round(digits=1)
    median <- citations %>% filter(pubid!="X5YyAB84Iw4C") %>% pull(cites) %>% median(na.rm=TRUE)
    
    profile <- list(
      "pub_count"=pub_count,
      "total_cites"=total_cites,
      "i10_index"=i10_index,
      "h_index"=h_index,
      "mean"=mean,
      "median"=median
    )
    
    return(profile)
  }
  
  get_crossref_profile <- function(doi_tbl) {
    cite_vec <- NULL
    for (row in 1:nrow(doi_tbl)) {
      doi <- doi_tbl[row,1]
      cite_num <- cr_citation_count(doi)$count
      cite_vec <- c(cite_vec, cite_num)
    }
    
    doi_tbl$Cites <- cite_vec
    
    citations <- doi_tbl[order(doi_tbl$Cites,decreasing=TRUE),]
    # citations <- cbind(id=rownames(citations),citations)
    # citations$id<- as.character(citations$id)
    # citations$id<- as.numeric(citations$id)
    
    pub_count <- nrow(citations)
    total_cites <- sum(citations$Cites)
    i10_index <- citations %>% filter(Cites >= 10) %>% nrow()
    h_index <- tail(which(citations$Cites >= seq_along(citations$Cites)), 1)
    mean <- round(mean(citations$Cites, na.rm=TRUE), digits=1)
    median <- median(citations$Cites, na.rm=TRUE)
      # max(which(citations$id<=citations$Cites))
    
    profile <- list(
      "pub_count"=pub_count,
      "total_cites"=total_cites,
      "i10_index"=i10_index,
      "h_index"=h_index,
      "mean"=mean,
      "median"=median
    )
    
    return(profile)
  }
  
  # retrieve citation info for Google Scholar
  google_scholar_id <- "umOwsMAAAAAJ"
  
  scholar_profile <- get_scholar_profile(google_scholar_id)
  
  # retrieve citation info for CrossRef
  dcc_doi <- read.table(doi_intable, header=TRUE, stringsAsFactors=FALSE)
  
  crossref_profile <- get_crossref_profile(dcc_doi)
  
  # create full citation table
  citation_summary_table <- data.frame(
    "Google Scholar"=c(scholar_profile$pub_count,
                       scholar_profile$total_cites,
                       scholar_profile$i10_index,
                       scholar_profile$h_index,
                       scholar_profile$mean,
                       scholar_profile$median),
    "CrossRef"=c(crossref_profile$pub_count,
                 crossref_profile$total_cites,
                 crossref_profile$i10_index,
                 crossref_profile$h_index,
                 crossref_profile$mean,
                 crossref_profile$median),
    row.names=c("Publication Count", "Total Citations", "i10 Index", 
                "H Index", "Mean", "Median")
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
    filter(year >= 2015) %>%
    ggplot(aes(x=year, y=cites)) + 
    geom_bar(stat="identity", fill="#0395DE") +
    geom_text(aes(label=cites), 
              position=position_dodge(width=0.9), vjust=-0.5,
              size=4) +
    scale_y_continuous(limits=c(0, max(cite_history$cites)+50)) +
    scale_x_continuous(breaks=cite_history$year) +
    labs(x="Year", y="Citations") +
    theme_minimal(base_size=16) +
    theme(panel.grid=element_blank(), axis.text.y=element_blank(),
          axis.title=element_blank(), 
          axis.text.x=element_text(color="black", vjust=5))
  
  # put it all together
  png(out_imgpath, width=6, height=2.5, units="in", res=600)
  grid.arrange(citation_plot, citation_table, nrow=1)
  dev.off()

  # write out a citation summary text file
  write.table(citation_summary_table, out_txtpath, quote=FALSE,
              sep="\t", row.names=FALSE)

}

produce_citation_summary(args[1], args[2], args[3])

