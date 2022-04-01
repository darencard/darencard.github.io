require(rcrossref)
require(scholar)
require(easyPubMed)
require(tm)
require(tidyverse)
require(gt)
require(ggwordcloud)
require(colorblindr)
require(gridExtra)
require(grid)
require(cowplot)

# inspired by https://twitter.com/alperezqui/status/1365042068435918855
# see also: https://github.com/alperezq/FancyPubFigures

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

# current year
curr_year <- as.numeric(format(Sys.Date(), "%Y"))

# retrieve citation info for Google Scholar
google_scholar_id <- "umOwsMAAAAAJ"

scholar_profile <- get_scholar_profile(google_scholar_id)

# retrieve citation info for CrossRef

### IGNORE LATER
doi_intable <- "~/Cloud_Drives/git/darencard.github.io/assets/scripts/doi.txt"
###

dcc_doi <- read.table(doi_intable, header=TRUE, stringsAsFactors=FALSE)

crossref_profile <- get_crossref_profile(dcc_doi)

# create full citation table
citation_summary_table <- data.frame(
  "GoogleScholar"=c(scholar_profile$pub_count,
                     scholar_profile$total_cites,
                     scholar_profile$mean,
                     scholar_profile$median,
                     scholar_profile$i10_index,
                     scholar_profile$h_index),
  "CrossRef"=c(crossref_profile$pub_count,
               crossref_profile$total_cites,
               crossref_profile$mean,
               crossref_profile$median,
               crossref_profile$i10_index,
               crossref_profile$h_index),
  row.names=c("Publication Count", "Total Citations", "Mean", "Median", "i10 Index",
              "H Index")
)

citation_table <- tableGrob(citation_summary_table,
                            cols=c("Google Scholar", "CrossRef"),
                            theme=ttheme_minimal(base_size=10),
                            labs(title="test"))

# citation_table <- citation_summary_table %>%
#   rownames_to_column("Metric") %>%
#   gt() %>%
#   tab_header(title = md("Citation Metrics"))

# create citation plot
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

cite_history <- bind_rows(data.frame(year=2013, cites=0),
                          get_citation_history(google_scholar_id))

extract_first_author <- function(col) {
  gsub(",.*", "", col)
}

publications <- get_publications(google_scholar_id) %>%
  mutate(first_author=extract_first_author(author)) %>%
  mutate(first=factor(if_else(first_author=="DC Card", 1, 0))) %>%
  filter(pubid!="X5YyAB84Iw4C") %>%
  select(first, journal, year, pubid)

if (2021 %in% publications$year) {
  publications <- publications
} else {
  publications <- bind_rows(data.frame(journal=NA, year=2021, pubid=NA),
                            publications)
}

upper <- roundUpNice(max(cite_history$cites))

citation_plot <- cite_history %>%
  mutate(year=factor(year)) %>%
  # filter(year >= 2015) %>%
  ggplot(aes(x=year, y=cites)) +
    geom_col(aes(fill=cites),  width=0.5) +
    geom_text(aes(label=cites),
              position=position_dodge(width=0.9), vjust=-0.5,
              size=4) +
    scale_y_continuous(limits=c(0, max(cite_history$cites)+50)) +
    scale_fill_viridis_c(direction=-1) +
    labs(title="Citations Per Year",
         x="Year", y="Citations") +
    theme_minimal(base_size=16) +
    theme(panel.grid=element_blank(), axis.text.y=element_blank(),
          axis.title=element_blank(),
          axis.text.x=element_text(color="black", vjust=3),
          legend.position="none")

publish_plot <- publications %>%
  mutate(year=factor(year),
         journal=factor(journal),
         fix=rep(1, nrow(.))) %>%
  ggplot(aes(x=year, y=fix, color=journal)) +
    geom_jitter(aes(shape=first, size=first), width=FALSE) +
    scale_color_viridis_d(direction=-1) +
    scale_size_discrete(range=c(3, 4)) +
    labs(title="Publication Timeline (colored by journal)",
         x="Year", y="Citations") +
    theme_minimal(base_size=16) +
    theme(panel.grid=element_blank(), axis.text.y=element_blank(),
          axis.title=element_blank(),
          axis.text.x=element_text(color="black", vjust=3),
          legend.position="none")

# put it all together
# png(out_imgpath, width=6, height=2.5, units="in", res=600)
# grid.arrange(citation_plot, citation_table, nrow=1)
# dev.off()

top <- plot_grid(publish_plot, citation_plot,
          nrow=2)

bottom <- plot_grid(wordcloud, citation_table, nrow=1)

plot_grid(top, bottom, nrow=2, rel_heights=c(1, 0.5))

# write out a citation summary text file
write.table(citation_summary_table %>% rownames_to_column("metric"),
            "~/Downloads/test_cites.txt", quote=FALSE,
            sep="\t", row.names=FALSE)



pubmed_qry <- "Card DC[Author]"
my_entrez_id <- get_pubmed_ids(pubmed_qry)
my_abstracts_xml <- fetch_pubmed_data(pubmed_id_list = my_entrez_id)
my_PM_list <- articles_to_list(pubmed_data = my_abstracts_xml)

abstracts <- NULL
for (article in 1:length(my_PM_list)) {
  curr_PM_record <- my_PM_list[article]
  my.df <- head(article_to_df(curr_PM_record, max_chars = 10000), 1)
  abstracts <- paste(abstracts, my.df$abstract)
}

abstracts
docs <- Corpus(VectorSource(abstracts))

toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))

# docs <- tm_map(docs, toSpace, "/")
docs <- tm_map(docs, content_transformer(tolower))
docs <- tm_map(docs, removeNumbers)
docs <- tm_map(docs, removeWords, stopwords("english"))
docs <- tm_map(docs, removeWords, c("across",
                                    "also",
                                    "multiple",
                                    "within",
                                    "may",
                                    "provide",
                                    "using",
                                    "can",
                                    "new",
                                    "used",
                                    "among",
                                    "burmese"))
docs <- tm_map(docs, removePunctuation)
docs <- tm_map(docs, stripWhitespace)

dtm <- as.matrix(TermDocumentMatrix(docs))
v <- sort(rowSums(dtm),decreasing=TRUE)
d <- data.frame(word = names(v),freq = v)

wordcloud <- d %>%
  filter(freq >= 15) %>%
  ggplot(aes(label = word, size = freq, color = factor(freq))) +
    geom_text_wordcloud_area(shape = "square") +
    labs(title="Abstract Wordcloud") +
    # scale_size_area(max_size = 5) +
    scale_color_viridis_d(direction=-1) +
    theme_minimal()
