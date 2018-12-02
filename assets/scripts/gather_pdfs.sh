# script that uses Zotero export to gather PDFs and move them into website repo
# argument 1: CSV file of citation library exported from Zotero

cat ${1} | \
awk -F "\"\,\"" '{ print $38 }' | \
awk -F ";" '{for (i=1;i<=NF;i++){if ($i ~/.pdf$/) {print $i}}}' | \
cut -d "_" -f 1-7 | \
while read path; \
do \
cp ${path}* ../pdfs/; \
done
