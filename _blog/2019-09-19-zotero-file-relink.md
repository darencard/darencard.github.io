---
layout: posts
title: "Relinking ZotFile Attachments"
date: 2019-09-19
excerpt: "How to re-establish ZotFile links to PDF attachments"
---

## Background

For the past couple of years, I have used [Zotero](https://www.zotero.org/) to manage my bibliography. I highly recommend it because it does a great job with storing bibliography, has plugins to automatically download bibliographic info and PDFs in Chrome and to insert citations in MS Word, and because it is open source! Bibliography entries can by synced with Zotero for free (to a limit), making them accessible anywhere through the internet. The wonderful plugin [ZotFile](http://zotfile.com/) can be used to automatically rename and organize PDF files in a custom way and to link bibliography entries with PDFs stored in arbitrary locations, making it useful for storing your PDF library in a custom location. I currently keep my PDFs in a local [Box](https://www.box.com/) directory that is synced with the cloud, so I can easily view my files through Zotero locally and also retrieve them from anywhere later through the internet.

I am currently migrating to a new computer and am able to sync down my Zotero bibliography into the Zotero app and my PDF library to my local folder. My bibliography entries in Zotero have the links still associated with them and my PDFs are in the same absolute path as my last computer. However, all of the links are broken, which is a major issue. Rather than manually re-establishing the links, I decided to look into a more automated solution and I was fortunately able to find one! I document it below in case I, or others, need to do this again.

## The Problem

1. As described, my full bibliography is downloaded and viewable in the Zotero app. So this is a prerequisite.
2. My Box directory containing my full PDF library is also synced using the Box Sync (not Box Drive!) app - a second prerequisite. I have it stored in the same path as my last computer, but it turns out you could store it anywhere. All of my PDFs are in a single directory and I currently do not have them organized into subdirectories (Zotero does allow this), so I cannot speak to how subdirectories may complicate this procedure.
3. I had also re-installed the ZotFile plugin on my fresh computer and had reconfigured the settings so that the file renaming matched what my renaming scheme is and the file location path was pointing to where my PDF library was stored.
4. The big issue is that whenever I double click on a bibliography entry, which should automatically open the linked PDF like it did on my old computer, I get the following error message, indicating that my links are all broken.

```
File Not Found

The attached file could not be found.

It may have been moved or deleted outside of Zotero, or, if the file was added on another computer, it may not yet have been synced to zotero.org.
```

## The Solution

1. In looking for a solution, I came across several threads in the Zotero Forums that were useful in pointing me in the right direction. See [thread #1](https://forums.zotero.org/discussion/64933/batching-editing-selected-file-links-with-zutilo) and [thread #2](https://forums.zotero.org/discussion/71618/scan-for-multiple-files-after-renaming-folder). Turns out the quick solution is another Zotero plugin called [Zutilo](https://github.com/willsALMANJ/Zutilo). This needs to be download and installed, which is pretty easy.
2. The documentation for link re-establishment in Zutilo could be better but I eventually figured it out. In essence, all you need to do is provide the "old" link name and the "new" link name for a rename "old" --> "new" type command.
3. Some settings first need to be configured with the Zutilo plugin. In Zotero, click on the `Tools` menu and on `Zutilo Preferences...`. Then, under the `User Interface` tab, make sure you adjust the settings so that `Show attachments` and `Modify attachments` are set to `Zutilo context menu`.
4. Next, we need to find the "old" link name. To do so, click the arrow button and then right-click on the linked PDF attachment and go to `Zutilo` and `Show attachment paths`. A prompt will appear with this path. Mine had this message:

```
Zutilo: showing path of attachment 1 of 1

attachments:<filename>.pdf
```

5. Based on this message, the "old" path is `attachments:`. You can copy that out.
6. Next you need to know the full path to where your associated PDF is stored. I did this in the terminal by navigating to the folder and running `pwd`. On my system, this ended up being `/Users/darencard/Cloud_Drives/Box Sync/Zotero_library/` but it will vary for everyone else.
7. Now, you should first test the relinking with a bibliography entry or two before doing everything. So click on one or a few bibliography entries (the collapsed parent line, not the linked file lines) to select them. Then right-click and go to `Zutilo` and `Modify attachment paths`.
8. You will now get a prompt that displays the following text.

```
Zutilo: Modify file attachments

Old partial path to attachments' directory to be replaced:
```

9. This is basically asking for the "old" path. In the blank input line, you should paste in the "old" path, which will probably be `attachments:`. Also check the option `Replace all instances of partial path string`. Click OK.
10. The next prompt has the following text.

```
Zutilo: Modify file attachments

New partial path to be used for attachments with paths matching the old partial path:
```

11. This is basically asking for the "new" path. In the blank input line, you should paste in the "new" path, which will probably something like `/Users/darencard/Cloud_Drives/Box Sync/Zotero_library/` (don't forget the last `/`). Click OK.
12. You should now be able to double click on those bibliographic entries and have the linked PDFs properly open. You can repeat the above steps on the rest of your files to re-establish those links as well. Note that this only relinks the PDFs and other file links (e.g., to HTML files/PubMed entries downloaded by the Zotero Chrome plugin) will still not function. Given I don't need those, this does not bother me.

That's it! Hopefully this is helpful to others!
