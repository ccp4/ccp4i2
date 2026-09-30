from __future__ import print_function
import sys
import os
import glob
import subprocess
htmlFiles=glob.glob(os.path.join(os.path.abspath('*.html')))

for htmlFile in htmlFiles:
    root, ext = os.path.splitext(htmlFile)
    dirName, fileName = os.path.split(root)
    rstFile = root+".rst"
    subprocess.call(['pandoc', htmlFile, '-o', rstFile])
    print(fileName+'.rst')
