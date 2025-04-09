import os
import shutil
import subprocess
import sys
import tempfile


# ------------------------------------------------------------------------------

class DefMap(dict):

   def __init__(self, def_path, typed):
      if os.path.isfile(def_path):
         if os.access(def_path, os.R_OK):
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
                        self[key] = value.strip(lquote)

            def_file.close()

   def get_key(self, query):
      for key, value in list(self.items()):
         if query == value:
            return key

      for key, value in list(self.items()):
         if query == value.strip().replace("/", "\\"):
            return key

   def get_viewers(self):
      if "EXPORT_TO_QTRVIEW" in self and "RUN_COOT" in self and "RUN_CCP4MG" in self:
         if self["EXPORT_TO_QTRVIEW"]:
            view_coot = self["RUN_COOT"]
            if os.path.isabs(view_coot):
               view_coot = os.path.abspath(view_coot)

            view_qtmg = self["RUN_CCP4MG"]
            if os.path.isabs(view_qtmg):
               view_qtmg = os.path.abspath(view_qtmg)

            return(view_coot, view_qtmg)

# ------------------------------------------------------------------------------

class TmpDir:

   def __init__(self):
      self.path = os.path.abspath("report_dir")
      if not os.path.isdir(self.path):
         ccp4_scr = os.environ["CCP4_SCR"]
         if ccp4_scr:
            if not os.access(ccp4_scr, os.W_OK):
               ccp4_scr = None

         try:
            self.path = tempfile.mkdtemp(suffix="ccp4qtrv", dir=ccp4_scr)

         except:
            print("Cannot create or write into temporary directory. Exiting.")
            sys.exit(1)

   def clear(self):
      try:
         shutil.rmtree(self.path)

      except:
         print("Cannot remove temporary directory", self.path)
         sys.exit(1)

# ------------------------------------------------------------------------------

def from_command_line(log_file_in):
   log_file_path = os.path.abspath(log_file_in)
   if not os.path.isfile(log_file_path):
      print(log_file_path, "is not a file. Exiting.")
      sys.exit(1)

   if not log_file_path.endswith(".log") and not log_file_path.endswith(".txt") :
      print(log_file_path.endswith, "dos not have extension log. Exiting.")
      sys.exit(1)

   log_dir_path, log_file_name = os.path.split(log_file_path)
   log_root_name = os.path.splitext(log_file_name)[0]
   log_root_parts = log_root_name.split("_")

   ccp4path = os.getenv("CCP4")
   if not ccp4path:
      print("Environement variable CCP4 is not set. Exiting.")
      sys.exit(1)

   ccp4path = os.path.abspath(ccp4path)

   viewer = os.path.join(ccp4path, "qtrview.app", "Contents", "MacOS", "qtrview")
   if not os.path.isfile(viewer):
      viewer = os.path.join(ccp4path, "bin", "qtrview")
      if not os.path.isfile(viewer):
         viewer = viewer + ".exe"
         if not os.path.isfile(viewer):
            print("qtrview is not found at expected location. Exiting.")
            sys.exit(1)

   if not os.access(viewer, os.X_OK):
      print(viewer, "is not executable. Exiting.")
      sys.exit(1)

   generator = os.path.join(os.path.dirname(__file__),"qtrgeneric.py")
   if not os.path.isfile(generator):
      print(generator, "does not exist. Exiting.")
      sys.exit(1)

   home_path = os.path.abspath(os.path.expanduser("~"))
   ddef_file_path = os.path.abspath("directories.def")
   if not os.path.isfile(ddef_file_path):
      ddef_file_path = os.path.join(log_dir_path, "CCP4_DATABASE", "directories.def")
      if not os.path.isfile(ddef_file_path):
         ddef_file_path = os.path.join(home_path, ".CCP4", "unix", "directories.def")
         if not os.path.isfile(ddef_file_path):
            ddef_file_path = os.path.join(home_path, "CCP4", "windows", "directories.def")
            if not os.path.isfile(ddef_file_path):
               ddef_file_path = None

   rep_dir = TmpDir()

   try:
      rep_rootname = os.path.join(rep_dir.path, "report")
      rep_inp_data = list()
      rep_inp_data.append("PYTHON " + sys.executable)
      rep_inp_data.append("XRT_GEN " + generator)
      rep_inp_data.append("LOGFILE " + log_file_path)
      rep_inp_data.append("REP_DIR " + rep_dir.path)
      rep_inp_data.append("REP_XRT " + rep_rootname + ".xrt")
      rep_inp_data.append("REP_XML " + rep_rootname + ".xml")
      rep_inp_data.append("CLEANUP NO")
      job_status = None
      if ddef_file_path:
         dir_def_map = DefMap(ddef_file_path, True)
         log_dir_key = dir_def_map.get_key(log_dir_path)
         if log_dir_key:
            var_name, var_sep, var_index = log_dir_key.partition(",")
            if var_sep == "," and var_name == "PROJECT_PATH":
               project = dir_def_map["PROJECT_ALIAS," + var_index]
               proj_db_dir = os.path.abspath(dir_def_map["PROJECT_DB," + var_index])
               if proj_db_dir:
                  if os.path.isfile(os.path.join(proj_db_dir, log_root_name + ".def")):
                     proj_db_map = DefMap(os.path.join(proj_db_dir, "database.def"), False)
                     if len(log_root_parts) > 1:
                        job_id = log_root_parts[0]
                        if log_file_name == proj_db_map["LOGFILE," + job_id]:
                           job_status = proj_db_map["STATUS," + job_id]
                           if job_status:
                              rep_inp_data.append("JOB_ID " + job_id)
                              rep_inp_data.append("PROJECT " + project)
                              rep_inp_data.append("DIR_DEF " + ddef_file_path)
                              rep_inp_data.append("STATUS " + job_status)

      if not job_status:
         rep_inp_data.append("STATUS RUNNING")
         job_title = None
         hkl_out_prefix = None
         if log_root_name == "pointandscale_sca":
            job_title = "Quick Scale"
            hkl_out_prefix = "ctruncate_"

         elif log_root_name == "pointandscale_sym":
            job_title = "Quick Symmetry"
            hkl_out_prefix = "pointless_"

         if hkl_out_prefix:
            log_file_unit = open(log_file_path)
            log_file_data = log_file_unit.read().split()
            log_file_unit.close()
            hkl_file_list = list()
            for word in log_file_data:
               if word.endswith(".mtz"):
                  hkl_file_path = word
                  if not os.path.isabs(word):
                     hkl_file_path = os.path.join(log_dir_path, word)

                  hkl_file_path = os.path.abspath(hkl_file_path)

                  if hkl_file_path not in hkl_file_list:
                     hkl_file_list.append(hkl_file_path)

            hkl_out_path = None
            hkl_in_name = None
            for hkl_file_path in hkl_file_list:
               hkl_file_name = os.path.split(hkl_file_path)[1]
               if hkl_file_name.startswith(hkl_out_prefix):
                  hkl_out_path = hkl_file_path
                  hkl_in_name = hkl_file_name.partition("_")[2]

            hkl_in_path = None
            for hkl_file_path in hkl_file_list:
               hkl_file_name = os.path.split(hkl_file_path)[1]
               if hkl_file_name == hkl_in_name:
                  hkl_in_path = hkl_file_path

            if hkl_in_path and hkl_out_path:
               rep_inp_data.append("HKLIN " + hkl_in_path)
               rep_inp_data.append("HKLOUT " + hkl_out_path)

            rep_inp_data.append("TITLE " + job_title)

         else:
            rep_inp_data.append("TITLE " + log_file_name)


      conf_viewers = DefMap(os.path.join(log_dir_path, "CCP4_DATABASE", "configure.def"), True).get_viewers()
      if not conf_viewers:
         conf_viewers = DefMap(os.path.join(home_path, ".CCP4", "unix", "configure.def"), True).get_viewers()
         if not conf_viewers:
            conf_viewers = DefMap(os.path.join(home_path, ".CCP4", "windows", "configure.def"), True).get_viewers()
            if not conf_viewers:
               conf_viewers = DefMap(os.path.join(ccp4path, "share", "ccp4i", "etc", "UNIX", "configure.def"), True).get_viewers()
               if not conf_viewers:
                  conf_viewers = DefMap(os.path.join(ccp4path, "share", "ccp4i", "etc", "WINDOWS", "configure.def"), True).get_viewers()

      if conf_viewers:
         rep_inp_data.append("RUN_COOT " + conf_viewers[0])
         rep_inp_data.append("RUN_CCP4MG " + conf_viewers[1])

      if os.path.exists("nograph"):
         print("\n".join(rep_inp_data))

   except:
      print("Cannot generate input data for QtRView. Exiting.")
      rep_dir.clear()
      sys.exit(1)

   try:
      rep_inp_data.append("")
      rep_inp_path = rep_rootname + ".inp"
      rep_inp_unit = open(rep_inp_path, "w")
      rep_inp_unit.write("\n".join(rep_inp_data))
      rep_inp_unit.close()

   except:
      print("Cannot write input file for QtRView. Exiting.")
      rep_dir.clear()
      sys.exit(1)

   if os.path.exists("nograph"):
      pass

   else:
      try :
         subprocess.Popen([viewer, "--inp-file", rep_inp_path]).wait()

      except :
         print("Cannot launch QtRView. Exiting.")
         rep_dir.clear()
         sys.exit(1)

   rep_dir.clear()
   sys.exit(0)

# ------------------------------------------------------------------------------

if __name__ == '__main__':
   if len(sys.argv) != 2:
      print("Usage:", sys.argv[0], "log-file")
      sys.exit(1)

   from_command_line(sys.argv[1])

# ------------------------------------------------------------------------------

