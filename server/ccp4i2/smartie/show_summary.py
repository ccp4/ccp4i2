#     baubles.py: a smarter CCP4 logfile browser
#     Copyright (C) 2006-2007 Peter Briggs, Wanjuan Yang, CCLRC 
#
#     This code is distributed under the terms and conditions of the
#     CCP4 licence agreement as `Part 1' (Annex 2) software.
#     A copy of the CCP4 licence can be obtained by writing to the
#     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
#
########################################################################
#
# show_summary.py
#
########################################################################

"""show_summary

Given the name of a CCP4 log file, print the text enclosed in summary
tags to standard out."""

__cvs_id__ = "$Id: show_summary.py,v 1.1 2008/08/26 14:23:39 pjx Exp $"

from ccp4i2.smartie import smartie
import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python "+str(sys.argv[0])+" <file>")
        sys.exit(1)
    logfile = smartie.parselog(str(sys.argv[1]))
    for i in range(0,logfile.nsummaries()):
        print(smartie.strip_logfile_html(logfile.summary(i).retrieve()))
