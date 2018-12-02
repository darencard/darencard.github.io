# remove previous PDFs
git rm ../pdfs/*

# create new PDF directory
mkdir ../pdfs

# gather PDFs into directory from my Zotero library
bash ./gather_pdfs.sh dcc_bibliography.csv

# format full publications markdown page and write to proper location
cat \
<(cat yaml_header.txt) \
<(bash ./format_markdown.sh dcc_bibliography.bib dcc_bibliography.csv \
> ../../_pages/publications.md
