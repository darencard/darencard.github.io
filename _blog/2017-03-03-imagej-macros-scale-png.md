---
layout: posts
title: "ImageJ macro for adding scale bars automatically"
date: 2017-03-03
excerpt: "Useful ImageJ macro that automatically adds scale bars to images"
---

Below is a ImageJ macro that will read a user-provided image file, add a scale bar, and then output a PNG image. The user must specify two arguments, the input file and the output file, as one quote-enclosed argument (see example below). The scale bar characteristics have been hard coded and can be changed by hand, or the script can be modified to allow argument specifications.

This macro is designed to be called from the command line using the ImageJ executable. With my Mac OSX computer running Fiji, the path is `/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx`. This has not been tested elsewhere and may not work without some effort. It relies on the Bio-Formats plugin to read the file and was written to convert from Zeiss's .czi files, so no guarantee that it works with others as desired. It is especially important to note that this does not set the scale, but infers it based on the metadata stored in the .czi files. Therefore, it will probably not work well with other file types.

If you store this macro file in the designated `macros` directory for your ImageJ/Fiji installation, you do not need to specify the full path to the macro file. In the example below I have named the macro `czi_scale_to_png.ijm` but you could use anything. Also note the quotes around the two arguments, which are necessary for this to run correctly.

#### Usage:

`/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --headless -macro czi_scale_to_png.ijm 'input output'`

#### Macro:

```
// read arguments from command-line
args = getArgument();

// load Bio-Formats macro extensions
run("Bio-Formats Macro Extensions");

// split arguments by space and assign to input/output file variable
arglist = split(args, " ");

infile = arglist[0];
outfile = arglist[1];

// open the input file
Ext.openImagePlus(infile);

// add a scale bar of unit 10, in black 12pt font in the bottom right corner
run("Scale Bar...", "width=10 height=4 font=12 color=Black background=None location=[Lower Right] bold overlay");

// save as a png
saveAs("Png", outfile);

// close input file
close();

// evaluate
eval("script", "System.exit(0);");
```

#### Batching:

It is easy to batch process a set of files using simple Shell scripting. Here is an example that will process all .czi files in a directory. I find this easier to understand and more powerful that hard-coding a loop into the macro.

```bash
for i in *.czi; do echo $i; cmd="/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --headless -macro czi_scale_to_png.ijm '$i $(echo $i | sed 's/.czi/.png/g')'"; echo $cmd; eval $cmd; done
```
