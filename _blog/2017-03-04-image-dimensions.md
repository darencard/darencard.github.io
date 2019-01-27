---
layout: posts
title: "Image dimensions from ImageJ"
date: 2017-03-04
excerpt: "Automatically detecting and outputing image dimensions from ImageJ."
---

Below is a ImageJ Python script that will read a user-provided image file and output the width and height, with units, to stdout. This macro is designed to be called from the command line using the ImageJ executable. With my Mac OSX computer running Fiji, the path is `/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx`. This has not been tested elsewhere and may not work without some effort. It relies on the Bio-Formats plugin to read the file and was written to convert from Zeiss's .czi files, so no guarantee that it works with others as desired. It is especially important to note that this does not set the scale, but infers it based on the metadata stored in the .czi files. Therefore, it will probably not work well with other file types.

#### Usage:

`/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --headless get_image_dims.py input`

#### Python script:

```python
#!/usr/bin/env python

import sys
from loci.plugins import BF

infile = sys.argv[1]
imps = BF.openImagePlus(infile)
imp = imps[0]

cal = imp.getCalibration()

out_list = []
out_list.append(infile)
out_list.append(cal.pixelWidth * imp.getDimensions()[0])
out_list.append(cal.pixelHeight * imp.getDimensions()[1])
out_list.append(cal.unit)

sys.stdout.write('\t'.join(map(str, out_list)) + '\n')
```

#### Batching:

It is easy to batch process a set of files using simple Shell scripting. Here is an example that will process all .czi files in a directory. Given my installation always throws up a bunch of Java errors to stderr, I decided just to direct that to `/dev/null`.

`for i in *.czi; do /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --headless get_image_dims.py $i 2> /dev/null; done`

#### Additonal info:

Here are some links to pages where I extracted information useful for creating this script.
- [Extract calibration](http://imagej.1557.x6.nabble.com/Getting-the-scale-td3686642.html)
- [Reading data with Bio-Format](https://ilovesymposia.com/2014/02/26/fiji-jython/)
- [General scripting with Fiji](https://www.ini.uzh.ch/~acardona/fiji-tutorial/)
- [Calling script from command-line](http://stackoverflow.com/questions/41189104/running-jython-script-from-terminal-with-parameter)
