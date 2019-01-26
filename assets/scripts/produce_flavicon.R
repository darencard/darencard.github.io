library(tidyverse)
library(colorblindr)

background <- data.frame(chrom=factor(rep(1:10, each=2000)),
                         pos=1:20000,
                         stat=rnorm(20000, 0.1, 0.1))

outliers <- data.frame(chrom=factor(rep(c(2, 4, 5, 8, 10), each=100)),
                       pos=c(2500:2599, 7300:7399, 9000:9099, 
                             15500:15599, 19400:19499),
                       stat=c(rnorm(100, 0.75, 1),
                              rnorm(100, 0.65, 1),
                              rnorm(100, 0.9, 1),
                              rnorm(100, 0.8, 1),
                              rnorm(100, 0.95, 1)))

plot <- rbind(background, outliers) %>%
  rbind(data.frame(chrom=c(2, 5, 8, 10),
                   pos=c(2500, 11000, 17500, 19400),
                   stat=c(0.75, 0.9, 1, 0.8))) %>%
  ggplot(aes(x=pos, y=stat, color=chrom)) +
    geom_point(size=0.5) +
    scale_y_continuous(limits=c(0, 1)) +
    theme_void() +
    theme(legend.position="none",
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA))

ggsave(plot, filename="~/Desktop/website_flavicon.png", 
       height=1, width=1, units="in", bg="transparent")
