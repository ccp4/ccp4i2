import sys

# Analysis of XDS header
class Xdstype():
    #  Items and flags:
    #  program    XDS program which wrote the file
    #  xtype      file type  (and program)
    #    'XDS-INTEGRATE'   INTEGRATE
    #    'XDS-ASCII'       CORRECT
    #    'nXDS'            nXDS
    #    'XDS-XSCALE'      XSCALE
    #  merged      True if file is merged (usually False)
    #  scaled      True if file is scaled (ie not INTEGRATE)
    #  scaleable   True if file contains all information for scaling
    #                 (INtEGRATE or XDS_ASCII)
    #  friedel     True if Friedel symmetry has been imposed (usually False)
    #  wavelength  list of wavelengths, maybe multiple (XSCALE) or empty (nXDS)
    
    def __init__(self, xdsfile):
        # Determine type of XDS file
        
        xdsin = open(xdsfile, 'r')

        # Get all header lines
        header = []
        while(9):
            header.append(xdsin.readline())
            if 'END_OF_HEADER' in header[-1]:
                break
        xdsin.close()

        # Is it merged?
        self.merged = False   #  default for INTEGRATE
        if 'MERGE=TRUE' in header[0]:
            self.merged = True

        #  Friedel
        self.friedel = False   # usual case
        if "FRIEDEL'S_LAW=TRUE" in header[0]:
            self.friedel = True

        generate = [s for s in header if 'Generated' in s]
        #  Name of program which wrote the file
        self.program = generate[0].split()[2]

        self.xtype = None
        self.scaled = True  # except for INTEGRATE
        self.scaleable = True  # does file contain information for (re)scaling?
        if self.program == 'INTEGRATE':
            self.xtype = 'XDS-INTEGRATE'
            self.scaled = False
        elif self.program == 'CORRECT':
            self.xtype = 'XDS-ASCII'
        elif self.program == 'nXDS':
            self.xtype = 'nXDS'
            self.scaleable = False
        elif self.program  == 'XSCALE':
            self.xtype = 'XDS-XSCALE'
            self.scaleable = False

        # Wavelength(s)
        wavelengthlines = [s for s in header if 'X-RAY_WAVELENGTH' in s]
        self.wavelengths = []
        for l in wavelengthlines:
            ls = l.split()
            idx = [idx for idx,s in enumerate(ls) if 'WAVELENGTH' in s]
            if len(idx) > 0:
                wavelength = ls[idx[0]+1]
                self.wavelengths.append(wavelength)

    def dump(self):
        print(self.program, self.xtype,
              '\n   Merged:', self.merged,
              ', Scaled:', self.scaled,
              ', Scaleable:', self.scaleable,
              ', Friedel:', self.friedel,
              ', Wavelengths:', self.wavelengths)
            
            

if __name__ == "__main__":
    filename = sys.argv[1]
    Xt = Xdstype(filename)
    Xt.print()
    
