---
layout: posts
title: "Making Axis Text run Downwards in ggplot2"
date: 2022-06-08
excerpt: "How to run text downward instead of rotating it vertically"
---

### Making Axis Text run Downwards in `ggplot2`

I really like `ggplot2` for producing graphics in R and try to use it whenever I can. Over the years I have figured out a few hacks for custom graphics that have helped me and here is another I thought I would share more broadly, since I never did find a clear solution to my problem documented online.

The context here is that I was working with visualizing data based on multi-sequence alignments. There are lots of ways to plot these in R, but I was trying to summarize some data at a subset of alignment sites. Normally, with full alignment data, I would just accept the default `ggplot2` axis text, which has reasonable breaks along the continuous alignment positions allowing the user to understand the location of each site. However, with my subset of sites, it made more sense to set my alignment positions as a factor, which resulted in cluttered axis labels since `ggplot2` will naturally label each categorical factor along the axis. With a larger number of sites, it is logical that this would be messy quickly. Here is an example of what I mean using a quick dummy dataset I created on the fly.

```
# make dataframe of alignment positions (for x-axis) and a second variable for the y-axis

# plot the data as a heatmap using the geom_tile()

```

Now, one obvious solution is to rotate the text 90 degrees so it runs vertically. This setting is well documented and useful but I think it does sacrifice some readability. I have also noticed it is common in published multi-sequence alignments figures to instead run the text vertically while not rotating the letters/words. In other words, instead of writing the letters left-to-right and then rotating everything, the idea is to write the letters up-to-down in a stacked arrangement. Then the reader can read downwards through the axis label. Unfortunately, I did a bunch of Googling in vain and was consistently met with results about rotating the axis labels - it appears the keywords I was using were clashing and not yielding any solutions for what I was looking for. Here is an example plot of what I mean and I will show you how to create this below.



I ended up having the idea that if I can insert a new-line character between each character of the axis text, it may insert those line breaks when it plots and result in the text running downwards like I want. Fortunately, when I looked for solutions for doing this, I found some better options - namely [this StackOverflow thread](https://stackoverflow.com/questions/67428819/create-a-function-that-insert-a-line-break-between-every-letter-of-a-character-v). In particular, I used the nice `stringr` package function for replacing an element of a text string, as outlined in the second thread response.

```
library(dplyr)
library(stringr)
library(ggplot2)

wrap_letters <- function(x) {
      stringr::str_replace_all(x, "(?<=.)(?=.)", "\n")
}
```

This function inserts a new-line character (`\n`) between all characters of a string. Placing it in a function allowed me to call it when going to create my plot, as follows:

```
ggplot(aes(x=wrap_letters(codon), y=reorder(pocket, desc(pocket)), fill=mhc1_Acarolinenesis)) +
    geom_tile(color="black")
```
