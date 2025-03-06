from __future__ import print_function


import os, sys, time
import smartie as SM
from xml.etree import ElementTree as ET
import shutil,glob

# ------------------------------------------------------------------------------

class LogConverter:

   def convert(self, logfile, job_info):
      cell_tag = "D"
      job_tag = "Job"
      job_path = "/Job"
      xrtns = "{http://www.ccp4.ac.uk/xrt}%s"

      self.xrttree = ET.ElementTree(ET.Element("report"))
      e0 = self.xrttree.getroot()
      e0.text = "\n"
      e0.tail = "\n"
      e0 = ET.SubElement(e0, xrtns %"results")
      e0.text = "\n\n"
      e0.tail = "\n"

      self.xmltree = ET.ElementTree(ET.Element(job_tag))
      f0 = self.xmltree.getroot()
      f0.text = "\n"
      f0.tail = "\n"

      e1 = ET.SubElement(e0, xrtns %"title", select=job_path + "/Title")
      e1.tail = "\n"
      e1.tail = "\n"

      f1 = ET.SubElement(f0, "Title")
      f1.text = job_info.title
      f1.tail = "\n"
      f1.tail = "\n"

      for ip in range(0, logfile.nprograms()):
         program = logfile.program(ip)
         program_tag = "SubJob_%03d" %ip
         program_path = job_path + "/" + program_tag

         all_attributes = program.attributes()
#        print all_attributes
#        for attribute in all_attributes:
#           print attribute, program.get_attribute(attribute)

         progname = ""
         if "name" in all_attributes:
            progname = program.get_attribute("name")

         elif "termination_name" in all_attributes:
            progname = program.get_attribute("termination_name")

         progtitle = "Run"
         if progname:
            progtitle += " of %s" %(progname)

         rundate = ""
         if "rundate" in all_attributes:
            rundate = program.get_attribute("rundate")
            progtitle += " on %s" %(rundate)

         runtime = ""
         if "runtime" in all_attributes:
            runtime = program.get_attribute("runtime")
            progtitle += " at %s" %(runtime)

         e1 = ET.SubElement(e0, xrtns %"section", title=progtitle)
         e1.text = "\n\n"
         e1.tail = "\n\n"

         f1 = ET.SubElement(f0, program_tag)
         f1.text = "\n"
         f1.tail = "\n"

         attributes = program.attributes()
         if "name" in attributes:
            e1.attrib["progname"] = program.get_attribute("name")

         if "rundate" in attributes:
            info_tag = "RunDate"
            e1.attrib["rundate"] = program_path + "/" + info_tag
            f2 = ET.SubElement(f1, info_tag)
            f2.text = program.get_attribute("rundate")
            f2.tail = "\n"

         if "runtime" in attributes:
            info_tag = "RunTime"
            e1.attrib["runtime"] = program_path + "/" + info_tag
            f2 = ET.SubElement(f1, info_tag)
            f2.text = program.get_attribute("runtime")
            f2.tail = "\n"

         e2 = None
         for ik in range(0, program.nkeytexts()):
            keytext = program.keytext(ik)
            keytext_tag = "KeyText_%03d" %(ik + 1)
            keytext_path = program_path + "/" + keytext_tag

            e2 = ET.SubElement(e1, xrtns %"keytext")
            e2.attrib["folded"] = "false"
            e2.attrib["name"] = keytext.name()
            e2.attrib["select"] = keytext_path
            e2.tail = "\n"

            f2 = ET.SubElement(f1, keytext_tag)
            f2.text = "\n" + keytext.message().strip() + "\n"
            f2.tail = "\n"

         if e2 != None:
            e2.tail = "\n\n"

         for xrt_table_tag in ("graph", "table"):
            jt = 0
            e1a = e1

            if xrt_table_tag == "graph" and program.ntables() > 0:
               e1a = ET.SubElement(e1, xrtns %"graph")
               e1a.text = "\n\n"
               e1a.tail = "\n\n\n"

            for it in range(0, program.ntables()):
               jt += 1
               table = program.table(it)
               table_tag = "Table_%03d" %jt
               table_path = program_path + "/" + table_tag

               e2 = ET.SubElement(e1a, xrtns %"table", select=table_path)
               e2.attrib["type"] = "plain"
               e2.attrib["title"] = table.title().strip()
               e2.text = "\n\n"
               e2.tail = "\n\n"
               if xrt_table_tag == "table":
                  e2.attrib["folded"] = "true"

               for ic in range(0, table.ncolumns()):
                  column = table.table_column(ic)
                  column_tag = "Column_%03d" %ic
                  cell_path = column_tag + "/" + cell_tag

                  e3 = ET.SubElement(e2, xrtns %"data", title=column.title())
                  e3.tail = "\n"

               e3.tail = "\n\n"

               if xrt_table_tag == "graph":
#                 cells = iter(table.data())
#                 print table.data()
#                 plaintable = [""]
#                 for ir in range(0, table.nrows()):
#                    plainline = list()
#                    for ic in range(0, table.ncolumns()):
#                       plainline.append(str(cells.next()).strip())

#                    plaintable.append(" ".join(plainline))

#                 plaintable.append("")
                  f2 = ET.SubElement(f1, table_tag)
#                 f2.text = "\n".join(plaintable)
                  f2.text = table.data().rstrip() + "\n"
                  f2.tail = "\n"

                  for ig in range(0, table.ngraphs()):
                     graph = table.table_graph(ig)
                     columns = graph.columns()

                     e3 = ET.SubElement(e2, xrtns %"plot")
                     e3.attrib["title"] = graph.title().strip()
                     e3.text = "\n"
                     e3.tail = "\n\n"

                     for column in columns[1:]:
                        e4 = ET.SubElement(e3, "plotline", xcol=str(columns[0]), ycol=str(column), colour="auto")
                        e4.tail = "\n"

                     scaling = graph.scaling()
                     if scaling == "N":
                        e3.attrib["ymin"] = "0"

                     elif scaling and scaling != "A":
                        xrange, sep, yrange = scaling.partition("x")
                        if sep == "x":
                           xmin, sep, xmax = xrange.partition("|")
                           if sep == "|":
                              e3.attrib["xmin"] = xmin.strip()
                              e3.attrib["xmax"] = xmax.strip()

                           ymin, sep, ymax = yrange.partition("|")
                           if sep == "|":
                              e3.attrib["ymin"] = ymin.strip()
                              e3.attrib["ymax"] = ymax.strip()


      f1 = ET.SubElement(f0, "Files")
      f1.text = "\n"
      f1.tail = "\n"

      rfiles = list()
      ifiles = list()
      for cou in range(0, len(job_info.ifiles)):
         file = job_info.ifiles[cou]
         if file[0].startswith("RESTRAINT"):
            rfiles.append(file)

         else:
            ifiles.append(file)

      ofiles = list()
      for cou in range(0, len(job_info.ofiles)):
         file = job_info.ofiles[cou]
         if file[0].startswith("RESTRAINT"):
            rfiles.append(file)

         else:
            ofiles.append(file)

      cou = 0
      for title, files in (("Input Files", ifiles), ("External Restraints", rfiles), ("Output Files", ofiles)):
         if files:
            e1 = ET.SubElement(e0, xrtns %"section", title=title)
            e1.text = "\n\n"
            e1.tail = "\n\n"
            e2 = ET.SubElement(e1, xrtns %"files")
            e2.text = "\n"
            e2.tail = "\n\n"
            for key, type, path, title in files:
               exists = False
               if key.startswith("STRUCTURE"):
                  exists =  os.path.isfile(path + ".pdb") and os.path.isfile(path + ".mtz")

               else:
                  exists =  os.path.isfile(path)

               if exists:
                  tag = key
                  cou += 1
                  tag = "File_%03d" %cou
                  xml_path = job_path + "/Files/" + tag
                  e3 = ET.SubElement(e2, xrtns %"file", key=key, type=type, select=xml_path)
                  e3.tail = "\n"
                  if title:
                     e3.attrib["title"] = title

                  always_folded = ["RESTRAINT_HTML"]
                  if key in always_folded:
                     e3.attrib["folded"] = "always"

                  f2 = ET.SubElement(f1, tag, path=path)
                  f2.tail = "\n"

# ------------------------------------------------------------------------------
# a patch for QuickScale and QuickSymm in imosflm and for a log-file generated without using the interface

class CLReader:

   def __init__(self, imap):
      self.ifiles = list()
      if "HKLIN" in imap:
         self.ifiles.append(("HKLIN", "hkl:hkl", os.path.realpath(imap.get("HKLIN")), ""))

      self.ofiles = list()
      if "HKLOUT" in imap:
         self.ofiles.append(("HKLOUT", "hkl:hkl", os.path.realpath(imap.get("HKLOUT")), ""))

      self.xrt = imap.get("REP_XRT")
      self.xml = imap.get("REP_XML")
      self.log = imap.get("LOGFILE")
      self.title = imap.get("TITLE")
#     print "\n".join(("", "LOG-input:", self.log, "", "XRT-output:", self.xrt, "", "XML-output:", self.xml))

# ------------------------------------------------------------------------------

class DBReader:

   def read_def(self, in_def_path, typed):
      def_path = os.path.realpath(in_def_path)
      def_map = dict()
      def_file = open(def_path)
      for line in def_file:
         if not line.startswith("#"):
            key, sep, value = line.strip().replace("\t", " ").partition(" ")
            if sep == " ":
               if typed:
                  type, sep, value = value.lstrip().partition(" ")

            if sep == " ":
               value = value.lstrip()
               lquote = ""
               if value.startswith('"'):
                  lquote = '"'

               elif value.startswith("'"):
                  lquote = "'"

               rquote = ""
               if value.endswith('"'):
                  rquote = '"'

               elif value.endswith("'"):
                  rquote = "'"

               if lquote == rquote:
                  value = value.strip(lquote)
                  if value == ".":
                     value = os.path.dirname(def_path)

                  elif value == "..":
                     value = os.path.dirname(os.path.dirname(def_path))

                  def_map[key] = value


      def_file.close()
      return def_map


   def __init__(self, jobid, projid, dir_def_path):
      dir_def_map = self.read_def(dir_def_path, True)
      proj_path_map = dict()
      db_path_map = dict()
#     print dir_def_map
      for rec_id, key in list(dir_def_map.items()):
         if rec_id.startswith("PROJECT_ALIAS,"):
            rec, sep, id = rec_id.partition(",")
            proj_path_map[key] = dir_def_map.get("PROJECT_PATH," + id)
            db_path_map[key] = dir_def_map.get("PROJECT_DB," + id)

         elif rec_id.startswith("DEF_DIR_ALIAS,"):
            rec, sep, id = rec_id.partition(",")
            proj_path_map[key] = dir_def_map.get("DEF_DIR_PATH," + id)

      proj_path = proj_path_map.get(projid)
      db_path = db_path_map.get(projid)
#     print db_path_map

#     print "AAAAAA", db_path
      proj_def_path = os.path.join(db_path, "database.def")
      proj_def_map = self.read_def(proj_def_path, False)

      self.log = os.path.join(proj_path, proj_def_map.get("LOGFILE," + jobid))
      self.title = proj_def_map.get("TITLE," + jobid)
      taskname = proj_def_map.get("TASKNAME," + jobid)
      jobname = jobid + "_" + taskname

      job_def_path = os.path.join(db_path, jobname + ".def")
      job_def_map = self.read_def(job_def_path, False)
      self.ifiles = list()
      self.ofiles = list()
      for record, files in (("INPUT_FILES", self.ifiles), ("OUTPUT_FILES", self.ofiles)):
         for key in job_def_map.get(record, "").split():
            name = job_def_map.get(key)
            dir = job_def_map.get("DIR_" + key)
            if name and dir:
               dir = proj_path_map.get(dir, dir)
               path = os.path.join(dir, name)
               title = ""
#              if os.path.exists(path):
               if True:
                  type = "other"
                  if key.startswith("XYZ"):
                     type = "xyz"

                  elif key.startswith("HKL"):
                     type = "hkl:hkl"
                     exceptions = ["refmac5", "refmac5_review", "mrbump", "dimple"]
                     if record == "OUTPUT_FILES" and taskname in exceptions:
                        type = "hkl:map"

                  elif key == "RESTRAINTFILE":
                     type = "text"
                     title = "Restraints"

                  elif key == "OUTPUTS":
                     key = "Analysis"
                     type = "html"
                     title = "Prosmart analysis"

                  elif path.endswith(".seq"):
                     type = "text"

                  elif path.endswith(".pir"):
                     type = "text"

                  elif path.endswith(".fasta"):
                     type = "text"

                  elif path.endswith(".log") or path.endswith(".loggraph"):
                     type = "log"

                  files.append((key, type, path, title))

      if taskname == "refmac5":
         htmlfile = os.path.join(proj_path, jobid + "_prosmart", "ProSMART_Results.html")
         self.ofiles.append(("RESTRAINT_HTML", "html", htmlfile, "Prosmart analysis"))

      if taskname == "phaser_MR":
         fllist = proj_def_map.get("OUTPUT_FILES," + jobid, "").split()
         fllist.sort()
         pdblist = list()
         mtzlist = list()
         for flname in fllist:
            if flname.endswith(".pdb"):
               pdblist.append(os.path.join(proj_path, flname[:-4]))

            elif flname.endswith(".mtz"):
               mtzlist.append(os.path.join(proj_path, flname[:-4]))

            elif flname.endswith(".sol"):
               self.ofiles.append(("SOLUTIONS", "text", os.path.join(proj_path, flname), "Solutions"))

            elif flname.endswith(".sum"):
               self.ofiles.append(("SUMMARY", "summary", os.path.join(proj_path, flname), "Summary"))

         cou = 0
         for pdbname in pdblist:
            if pdbname in mtzlist:
               cou += 1
               self.ofiles.append(("STRUCTURE%d" %(cou), "xyz:map", pdbname, "Solution %d and maps" %(cou)))

            else:
               self.ofiles.append(("XYZOUT_PART", "xyz", pdbname + ".pdb", "Partial Solution Ensemble"))


#        print self.ofiles

      self.xrt = imap.get("REP_XRT")
      self.xml = imap.get("REP_XML")
#     print "\n".join(("", "LOG-input:", self.log, "", "XRT-output:", self.xrt, "", "XML-output:", self.xml))

# ------------------------------------------------------------------------------

if __name__ == '__main__':

      imap = dict()
      ifile = open(sys.argv[1])
      for line in ifile:
         key, sep, value = line.strip().partition(" ")
         if sep == " ":
            imap[key] = value.lstrip()

      ifile.close()
      job_info = None
      if "JOB_ID" in imap:
         jobid = imap.get("JOB_ID")
         projid = imap.get("PROJECT")
         dir_def_path = imap.get("DIR_DEF")
         job_info = DBReader(jobid, projid, dir_def_path)
         if os.path.realpath(job_info.log) != os.path.realpath(imap.get("LOGFILE")):
            raise Exception()
         logs = [os.path.realpath(p) for k,t,p,i in job_info.ofiles if t=='log']

      else:
         job_info = CLReader(imap)
         logs = sorted(glob.glob(os.path.join(os.path.dirname(job_info.log), '*.loggraph')))

      if os.path.realpath(job_info.log) not in logs:
         logs.insert(0, job_info.log)

      use_log = job_info.log
      # in case of multiple logs, merge all into one file
      if len(logs) > 1:
         use_log = os.path.join(os.path.dirname(job_info.xrt), 'merged.log')
         with open(use_log,'w') as g:
            for log in logs:
               if os.path.isfile(log):
                  with open(log) as f:
                     shutil.copyfileobj(f,g)

      log_file = open(use_log)
      cou = 0
      report = None
      for line in log_file:
         cou += 1
         if cou > 100:
            break

         if line.startswith("_JOB_DIRECTORY:"):
            key, sep, value = line.partition(":")
            report = os.path.join(os.path.realpath(value.strip()), "report")
            break

      if report:
         if os.path.isfile(report + ".xrt") and os.path.isfile(report + ".xml"):
            shutil.copy(report + ".xrt", job_info.xrt)
            shutil.copy(report + ".xml", job_info.xml)

      else:
         converter = LogConverter()
         converter.convert(SM.parselog(use_log), job_info )
         converter.xmltree.write(job_info.xml)
         ofile = open(job_info.xrt, "w")
         ofile.write(ET.tostring(converter.xrttree.getroot()).decode().replace("ns0", "xrt"))
         ofile.close()

# ------------------------------------------------------------------------------

