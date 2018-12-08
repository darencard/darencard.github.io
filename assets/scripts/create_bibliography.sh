# remove previous PDFs
git rm ../pdfs/*

# create new PDF directory
mkdir ../pdfs

# gather PDFs into directory from my Zotero library
bash ./gather_pdfs.sh dcc_bibliography.csv

# create citation summary graphics
cat dcc_bibliography.csv | awk -F "\"\,\"" '{ print $9 }' | sed '/^$/d' > doi.txt
Rscript produce_citation_summary.R doi.txt ../images/citation_summary.png

# format full publications markdown page and write to proper location
cat \
<(cat yaml_header.txt) \
<(cat citation_summary_markdown.txt) \
<(bash ./format_markdown.sh dcc_bibliography.bib dcc_bibliography.csv) \
<(echo -e "\n\nLast Updated:" `date +%Y-%m-%d`) \
> ../../_pages/publications.md
