# script that takes Zotero outputs and properly formats markdown of citations for website
# argument 1: BIB file of citation library exported from Zotero
# argument 2: CSV file of citation library exported from Zotero

paste -d " " \
<(pandoc tex_guide.tex --columns=100000 --csl=chicago-modified.csl -t markdown_strict --bibliography ${1} | \
sed -e 's/\\\[/\[/g' -e 's/\\\]/\]/g' -e 's/[<>]//g' | \
sed -e 's/Daren\ C\.\ Card/\*\*Daren\ C\.\ Card\*\*/g' -e 's/Card\,\ Daren\ C\./\*\*Card\,\ Daren\ C\.\*\*/g' | \
sed -e 's/DOI/\<i\ class\=\"ai\ ai\-doi\"\>\<\/i\>/g' | sed '/^$/d' | \
sed 's/2018\.\ Perry\,\ Blair\ W\.\,\ \*\*Daren\ C\.\ Card\*\*/2018\.\ Perry\,\ Blair\ W\.\†\/\*\*Daren\ C\.\ Card\†\*\*\ \[\†\ authors\ contributed\ equally\]/g') \
<(pandoc tex_guide.tex --columns=100000 --csl=chicago-modified.csl -t markdown_strict --bibliography ${1} | \
sed -e 's/\\\[/\[/g' -e 's/\\\]/\]/g' -e 's/[<>]//g' | \
sed -e 's/Daren\ C\.\ Card/\*\*Daren\ C\.\ Card\*\*/g' -e 's/Card\,\ Daren\ C\./\*\*Card\,\ Daren\ C\.\*\*/g' | \
sed -e 's/DOI/\<i\ class\=\"ai\ ai\-doi\"\>\<\/i\>/g' | \
awk -F "doi.org/" '{ print $NF }' | cut -d ")" -f 1 | sed '/^$/d' | \
while read doi; \
do \
grep "${doi}" ${2}; done | awk -F "\"\,\"" '{ print $38 }' | \
awk -F ";" '{for (i=1;i<=NF;i++){if ($i ~/.pdf$/) {print $i}}}' | \
awk -F "/" '{ print "[<i class=\"fas fa-file-pdf\"></i>](https://raw.githubusercontent.com/darencard/darencard.github.io/master/assets/pdfs/"$NF")" }') | \
sed -e 'G;' | \
awk '!x{x=sub("^2012. ","## 2012\n\n")}{print $0}' | \
awk '!x{x=sub("^2013. ","## 2013\n\n")}{print $0}' | \
awk '!x{x=sub("^2014. ","## 2014\n\n")}{print $0}' | \
awk '!x{x=sub("^2015. ","## 2015\n\n")}{print $0}' | \
awk '!x{x=sub("^2016. ","## 2016\n\n")}{print $0}' | \
awk '!x{x=sub("^2017. ","## 2017\n\n")}{print $0}' | \
awk '!x{x=sub("^2018. ","## 2018\n\n")}{print $0}' | \
awk '!x{x=sub("^2019. ","## 2019\n\n")}{print $0}' | \
awk '!x{x=sub("^2020. ","## 2020\n\n")}{print $0}' | \
awk '!x{x=sub("^2021. ","## 2021\n\n")}{print $0}' | \
awk '!x{x=sub("^2022. ","## 2022\n\n")}{print $0}' | \
awk '!x{x=sub("^2023. ","## 2023\n\n")}{print $0}' | \
awk '!x{x=sub("^2024. ","## 2024\n\n")}{print $0}' | \
awk '!x{x=sub("^2025. ","## 2025\n\n")}{print $0}' | \
awk '!x{x=sub("^2026. ","## 2026\n\n")}{print $0}' | \
awk '!x{x=sub("^2027. ","## 2027\n\n")}{print $0}' | \
awk '!x{x=sub("^2028. ","## 2028\n\n")}{print $0}' | \
awk '!x{x=sub("^2029. ","## 2029\n\n")}{print $0}' | \
awk '!x{x=sub("^2030. ","## 2030\n\n")}{print $0}' | \
sed -e 's/^20[1-9][0-9]\.\ //g'
