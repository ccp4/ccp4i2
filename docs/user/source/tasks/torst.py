from __future__ import print_function
import sys
import os
import glob
import subprocess

hashes = '################################################################################'

directories=glob.glob(os.path.abspath('*'))

for directory in directories:
    dirRoot = os.path.split(directory)[1]
    htmlFiles=glob.glob(os.path.join(directory,'*.html'))
    for htmlFile in htmlFiles:
        root, ext = os.path.splitext(htmlFile)
        dirName, fileName = os.path.split(root)
        rstTempFile = root+"_tmp.rst"
        rstFile = root+".rst"
        subprocess.call(['pandoc', htmlFile, '-o', rstTempFile])
        with open(rstTempFile,'r') as inputFile:
            text = inputFile.read()
            with open(rstFile, 'w') as outputFile:
                outputFile.write('{}\n'.format(hashes[:len(dirRoot)]))
                outputFile.write('{}\n'.format(dirRoot))
                outputFile.write('{}\n'.format(hashes[:len(dirRoot)]))
                outputFile.write('\n{}'.format(text))
        os.unlink(rstTempFile)
        print(os.path.join(os.path.split(os.path.dirname(rstFile))[1],fileName+'.rst'))
