import fs from "fs";
import path from "path";
import os from "os";

// Setting CCP4 directories and environment variables
export function ccp4_setup_sh(CCP4Dir) {
  const CCP4_MASTER = path.dirname(CCP4Dir);
  const CCP4 = CCP4Dir;
  const CCP4_SCR = path.join(
    "/tmp",
    os.userInfo().username.replace(/[^a-zA-Z0-9]/g, "_")
  );

  // Setting environment variables
  process.env.CCP4_MASTER = CCP4_MASTER;
  process.env.CCP4 = CCP4;
  process.env.CCP4_SCR = CCP4_SCR;
  process.env.CCP4I_TCLTK = path.join(CCP4, "bin");

  // Optional: set HTTP proxy (uncomment if needed)
  process.env.HTTP_PROXY = "wwwblah.blah.ac.uk:xxxx/blah.blah";

  // Gfortran-related variable
  process.env.GFORTRAN_UNBUFFERED_PRECONNECTED = "Y";

  // Setting paths to executables and libraries
  process.env.CBIN = path.join(CCP4, "bin");
  process.env.CLIB = path.join(CCP4, "lib");
  process.env.CLIBD = path.join(CCP4, "lib", "data");
  process.env.CETC = path.join(CCP4, "etc");
  process.env.CINCL = path.join(CCP4, "include");
  process.env.CHTML = path.join(CCP4, "html");
  process.env.CEXAM = path.join(CCP4, "examples");
  process.env.CCP4I_TOP = path.join(CCP4, "share", "ccp4i");
  process.env.MMCIFDIC = path.join(CCP4, "lib", "ccp4", "cif_mmdic.lib");
  process.env.CRANK = path.join(process.env.CCP4I_TOP, "crank");
  process.env.CLIBD_MON = path.join(CCP4, "lib", "data", "monomers", path.sep);
  process.env.CCP4_HELPDIR = path.join(CCP4, "help");

  // Check if CCP4 exists
  if (!fs.existsSync(CCP4)) {
    console.log(`WARNING: The directory ${CCP4} does not exist.`);
    console.log("WARNING: The CCP4 programs will not run correctly.");
  }

  // Ensure CCP4_SCR directory exists, create it if not
  if (!fs.existsSync(CCP4_SCR)) {
    fs.mkdirSync(CCP4_SCR, { recursive: true });
  }

  console.log(`CCP4_SCR directory set to: ${CCP4_SCR}`);

  // HARVESTHOME can be set if needed (uncomment if required)
  // process.env.HARVESTHOME = os.homedir();

  // Set MOSFLM_WISH if needed (uncomment and modify if required)
  // process.env.MOSFLM_WISH = path.join(process.env.CCP4I_TCLTK, 'wish');

  // Set PDB_EXTRACT path
  process.env.PDB_EXTRACT = path.join(CCP4, "share");

  // CCP4_OPEN is set to 'UNKNOWN' by default
  process.env.CCP4_OPEN = "UNKNOWN";

  // Add directories to PATH
  process.env.PATH = `${path.join(CCP4, "etc")}:${path.join(CCP4, "bin")}:${
    process.env.PATH
  }`;

  // Add to MANPATH if defined
  if (process.env.MANPATH) {
    process.env.MANPATH = `${path.join(CCP4, "share", "man")}:${
      process.env.MANPATH
    }`;
  }

  // Set SSL_CERT_FILE if it's not defined
  if (!process.env.SSL_CERT_FILE) {
    process.env.SSL_CERT_FILE = path.join(CCP4, "etc", "ssl", "cacert.pem");
  }

  // Traditional aliases (this is just for display; they can't be used directly in Node.js)
  console.log(
    "The following aliases are set in the shell environment (not usable in Node.js):"
  );
  console.log("ccp4, xtal, cbin, cetc, cincl, clib, clibd, cexam, chtml");

  // Clean up (in case older versions of CCP4 are sourced)
  delete process.env.DBCCP4I_TOP;

  // Source external setup file (if needed)
  // This would require invoking a shell script from within Node.js if needed
  const arpwarpSetupPath = `${CCP4_MASTER}/arpwarp_setup.bash`;
  if (fs.existsSync(arpwarpSetupPath)) {
    console.log(`Sourcing external setup file: ${arpwarpSetupPath}`);
    // You can execute this script by invoking a shell command using child_process if needed
  } else {
    console.log(`External setup file not found: ${arpwarpSetupPath}`);
  }

  console.log("CCP4 environment setup complete!");
}
