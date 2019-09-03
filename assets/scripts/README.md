## How to have Publications page automatically updated

1. Open Zotero and make sure that it is up to date, such that any new personal publications have been important and properly linked with PDF files using the ZotFile addon.

2. Select all personal citation rows and right click to reveal options. Click `Export Items...`.

3. Export the file first as a `BibTeX` object using default settings and save to the file name `dcc_bibliography.bib` in the path `darencard.github.io/assets/scripts/`. Repeat this step and save also to a CSV file in the same location named `dcc_bibliography.csv`.

4. In the terminal, move into `darencard.github.io/assets/scripts/`

5. Generate the Publications page automatically by executing the `create_bibliography.sh` script: `bash create_bibliography.sh`. Ensure the script runs with no errors and troubleshoot as necessary.

6. Return to the `darencard.github.io` base directory.

7. Stage, commit, and push the updates. After a few minutes, updates should be visible online.

```
git add *
git commit -m "updated publications"
git push origin master
```
