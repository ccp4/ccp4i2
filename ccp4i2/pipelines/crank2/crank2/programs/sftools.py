#!/usr/bin/python
import os,sys
from program import program
import common

class sftools(program):
  name="SFtools"
  binary="sftools"
  never_merge=True
  # setting rundir immediately so that relative path can be always used.  sftools has a limit of 200 characters
  always_rundir=True
  # in case of cell parameters, the first regexp is for global cell, the second for particular crystal cell
  stat={}
  stat['a'] = common.stats(regexp=[r" a\s+=\s+(\S+)", r" ({0})\s+{0}\s+{0}\s+{0}\s+{0}\s+{0}\s+".format('\d+\.\d+')+r"{0}"])
  stat['b'] = common.stats(regexp=[r" b\s+=\s+(\S+)", r" {0}\s+({0})\s+{0}\s+{0}\s+{0}\s+{0}\s+".format('\d+\.\d+')+r"{0}"])
  stat['c'] = common.stats(regexp=[r" c\s+=\s+(\S+)", r" {0}\s+{0}\s+({0})\s+{0}\s+{0}\s+{0}\s+".format('\d+\.\d+')+r"{0}"])
  stat['alpha'] = common.stats(regexp=[r" alpha\s+=\s+(\S+)", r" {0}\s+{0}\s+{0}\s+({0})\s+{0}\s+{0}\s+".format('\d+\.\d+')+r"{0}"])
  stat['beta'] = common.stats(regexp= [r" beta\s+=\s+(\S+)",  r" {0}\s+{0}\s+{0}\s+{0}\s+({0})\s+{0}\s+".format('\d+\.\d+')+r"{0}"])
  stat['gamma'] = common.stats(regexp=[r" gamma\s+=\s+(\S+)", r" {0}\s+{0}\s+{0}\s+{0}\s+{0}\s+({0})\s+".format('\d+\.\d+')+r"{0}"])
  stat['spacegroup'] = common.stats(regexp=r"Space group name\s+:\s+(.+)")
  stat['spacegroup_num'] = common.stats(regexp=r"Space group number\s+:\s+(\S+)")
  stat['xname'] = common.stats(regexp=r"\s+\d+\s+\S\s+{0}\s*\ncrystal name:\s+(\S+)\s*\ndataset name:\s+\S+",convert=False)
  stat['dname'] = common.stats(regexp=r"\s+\d+\s+\S\s+{0}\s*\ncrystal name:\s+\S+\s*\ndataset name:\s+(\S+)",convert=False)
  stat['resolution'] = common.stats(regexp=r"\s+\d+\s+\S\s+{0}\s+\-?\d+\.\d+\s*(?:(?:\-?\d+\.\d+)|(?:\*+))\s+\-?\d+\.\d+\s+(\d+\.\d+)")
  stat['low_resol'] = common.stats(regexp=r"\s+\d+\s+\S\s+{0}\s+\-?\d+\.\d+\s*(?:(?:\-?\d+\.\d+)|(?:\*+))\s+\-?\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+)")
  stat['max'] = common.stats(regexp=r"\s+\d+\s+\S\s+{0}\s+-?\d+\.\d\d\s*(\d+\.\d+)")
  stat['min'] = common.stats(regexp=r"\s+\d+\s+\S\s+{0}\s+-?(\d+\.\d\d)\s*\d+\.\d+")
  stat['aver'] = common.stats(regexp=r"\s+\d+\s+\S\s+{0}\s+-?\d+\.\d\d\s*\d+\.\d+\s*(-?\d+\.\d+)")
  stat['wavelength'] = common.stats(regexp=r"\s+(\d+\.\d+)\s+{0}\s+{1}")  # xname,dname
  stat['illegal_cell'] = common.stats(regexp=r"!!! ERROR: illegal cell parameters !!!")
  stat['zero_sigma'] = common.stats(regexp=r"!!! WARNING, zero sigma values !!!")
  stat['correl'] = common.stats(regexp=r"OVERALL STATISTICS\s+\S+\s+(-?\d+\.\d+)")
  stat['all_absent'] = common.stats(regexp=r"{0}\s+<  A  L  L    A  B  S  E  N  T  >")  # label
  stat['all_zero'] = common.stats(regexp=r"{0}\s+0.0+\s+0.0+\s+0.0+")  # label
  stat['num_refl'] = common.stats(regexp=r"(\d+) reflections now stored in memory")
