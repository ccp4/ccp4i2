import fs from "fs";
import path from "path";
import os from "os";

/**
 * Sets up the CCP4 environment variables.
 */
export function ccp4_setup_windows(CCP4Dir) {
  // Set CCP4 to the current directory
  let CCP4 = CCP4Dir;
  process.env.CCP4 = CCP4;
  // Remove trailing backslash
  CCP4 = CCP4.endsWith(path.sep) ? CCP4.slice(0, -1) : CCP4;

  // Check if CCP4 has spaces (optional part, as per original batch script)
  if (CCP4.includes(" ")) {
    console.error(
      "CCP4 path shall not have spaces. Move it to a different directory."
    );
    process.exit(1);
  }

  // Set CCP4_MASTER based on CCP4 directory
  let CCP4_MASTER = path.dirname(CCP4);

  // Define other environment variables
  process.env.CBIN = path.join(CCP4, "bin");
  process.env.CLIB = path.join(CCP4, "lib");
  process.env.CLIBD = path.join(CCP4, "lib", "data");
  process.env.CEXAM = path.join(CCP4, "examples");
  process.env.CHTML = path.join(CCP4, "html");
  process.env.CINCL = path.join(CCP4, "include");
  process.env.CCP4I_TOP = path.join(CCP4, "share", "ccp4i");
  process.env.CLIBD_MON = path.join(CCP4, "lib", "data", "monomers", path.sep);
  process.env.MMCIFDIC = path.join(CCP4, "lib", "ccp4", "cif_mmdic.lib");
  process.env.CRANK = path.join(CCP4, "share", "ccp4i", "crank");
  process.env.CCP4_OPEN = "unknown";
  process.env.GFORTRAN_UNBUFFERED_PRECONNECTED = "Y";
  process.env.HDF5_PLUGIN_PATH = path.join(
    CCP4,
    "lib",
    "base",
    "lib",
    "plugins"
  );
  process.env.SSL_CERT_FILE = path.join(CCP4, "etc", "ssl", "cacert.pem");
  process.env.QT_QPA_PLATFORM_PLUGIN_PATH = path.join(CCP4, "QtPlugins");
  process.env.CCP4I_TCLTK = path.join(CCP4, "bin");
  process.env.PDB_EXTRACT = path.join(CCP4, "share");

  // Set CCP4_BIN
  process.env.CCP4_BIN = path.join(CCP4, "bin");

  // Modify PATH variable
  const newPath = [
    path.join(CCP4, "nodejs"),
    path.join(CCP4, "bin"),
    path.join(CCP4, "etc"),
    path.join(CCP4, "Scripts"),
    //process.env.PATH.replace("CCP4", "C__4"),
  ].join(path.delimiter);
  process.env.PATH = newPath;

  // Clear old CCP4 version variables
  delete process.env.CCP4_LIB;
  delete process.env.CCP4_BIN;
  delete process.env.IMOSFLM_WISH;
  delete process.env.CCP4_BROWSER;
  delete process.env.CDOC;
  delete process.env.CETC;
  delete process.env.CLIBS;
  delete process.env.CPROG;
  delete process.env.BINSORT_SCR;
  delete process.env.IMOSFLM_VERSION;
  delete process.env.CCP4_HELPDIR;
  delete process.env.CCP4I_HELP;
  delete process.env.DBCCP4I_TOP;
  delete process.env.MOLREPLIB;
  delete process.env.MOSFLM_WISH;
  delete process.env.PISA_CONF_FILE;
  delete process.env.PUBLIC_FONT84;
  delete process.env.GFORTRAN_UNBUFFERED_ALL;

  // Set CCP4_SCR if it's not defined
  if (!process.env.CCP4_SCR) {
    process.env.CCP4_SCR = path.join(os.homedir(), "AppData", "Local", "Temp");
  }

  // Check if CCP4_SCR folder exists and create it if necessary
  const TW = path.join(
    process.env.CCP4_SCR,
    "._CCP4SCR_" + process.env.USERNAME
  );
  try {
    if (!fs.existsSync(TW)) {
      fs.mkdirSync(TW);
    }
    fs.rmdirSync(TW);
  } catch (err) {
    console.error("Error: The folder cannot be created or is not writable");
    process.exit(1);
  }

  // Check if CCP4_SCR has unexpected characters
  const regex = /^[A-Za-z0-9_\-.:\\]*$/;
  if (!regex.test(process.env.CCP4_SCR)) {
    console.error(
      `The scr path ${process.env.CCP4_SCR} contains unexpected characters.`
    );
    process.exit(1);
  }

  // Output environment variables
  console.log(`CCP4: ${CCP4}`);
  console.log(`CCP4_SCR: ${process.env.CCP4_SCR}`);
}
