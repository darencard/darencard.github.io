---
layout: posts
title: "Removing Space in Box Sync Folder"
date: 2019-09-27
excerpt: "How to remove the annoying space in the default Box Sync folder"
---

Currently, I am taking advantage of all the free cloud storage I can by scattering my files across different services. I use [Box](https://www.box.com/) to store my [Zotero](https://www.zotero.org/) bibliography library and PDFs. I also have some scripts that automatically build my Publications page for my academic website. Upon doing a fresh install of [Box Sync](), I was finding that the default name for the directory where the files were stored is `Box Sync`. The space is very annoying for anyone who does a lot of command-line work and the scripts I use for my website rely on the directory name `Box_Sync`. There is no way to easily change this in the Box GUI, so apparently I came up with something before and forgot what it was.

After looking around, I came across some [excellent instructions](https://community.box.com/t5/Using-Box-Sync/For-Box-Admins-Changing-the-Default-Location-for-the-Box-Sync/ta-p/84) on the [Box Community Forums](https://community.box.com/t5/Box-Community/ct-p/English). I'll resummarize them briefly here. Note that these instructions are for a Mac computer only and probably won't work on other systems.

1. Uninstall Box Sync completely: On a Mac, this includes deleting the app but also deleting some other files. See information linked in the instructions above for more info.

```bash
sudo rm -rf /Applications/Box\ Sync.app
sudo rm -rf /Library/PrivilegedHelperTools/com.box.sync*
rm -rf ~/Library/Logs/Box
rm -rf ~/Library/Application\ Support/Box
rm -rf ~/Library/Caches/com.box.sync
```

2. Create new plist file: We must create a new `.plist` file that Box will use to set the default Box Sync storage path and directory name. This involves creating a file `/Library/Preferences/com.box.sync.plist` that contains some important text. You can do so with `sudo nano`. Here is what the output should look like.

```bash
cat /Library/Preferences/com.box.sync.plist
# <!--?xml version="1.0" encoding="UTF-8"?--><plist version="1.0"><dict><key>SyncRootFolder</key><string>/Users/<USERNAME>/<OTHER_DIRS>/Box_Sync</string></dict></plist>
```

Your path will vary based on your username (replace `<USERNAME>` with your username) and any other directories you want to use to store the synced folder (replace `<OTHER_DIRS>` with path or nothing). Also note that the sync directory is now `Box_Sync` instead of `Box Sync`, which is what I want (you could change it to anything you want).

3. Reinstall Box Sync as you normally would by copying it to the Applications directory.

4. When Box Sync starts up, you will first log in. If all went well, the next page should give you the option "Start Syncing". Down below, there is a small link to "Customize Folder Settings". If you click on this, you should see the default path that matches what you set in the .plist file above. Importantly, you could change the parent directory if you would like and Box will create a directory called `Box_Sync` (or your desired name) for you, which is what we are really after.
