import math
import gemmi


class ColumnSet():
    # a set of CIF columns
    def __init__(self, columnlabels, flagpresent, text, typecode):
        # setlabels     the column labels
        # flagpresent   bool flags, False if missing
        # text          brief text describing the column set content
        # typecode      int code for data type

        self.columnlabels = columnlabels
        self.flagpresent = flagpresent
        self.text = text
        self.typecode = typecode
        self.OK = False not in flagpresent

    def missingColumns(self):
        if self.OK:
            return ""
        return ", ".join(
            label
            for label, present in zip(self.columnlabels, self.flagpresent)
            if not present
        )

    def formatcolumnlabels(self):
        return ", ".join(self.columnlabels)


class CIFLabelSets():
    # class to manage appropriate sets of reflection labels in an mmcif file
    # Each set should be all present or all absent

    # Note that Gemmi doesn't directly distinguish between columns
    #   from _refln and _diffrn_refln loops
    # For each set in list:
    #   1) list of CIF column names
    #   2) data type string
    #   3)  data type code
    #       = 0 freeR,
    #       > 0 suitable for input as observed data
    #             +1 intensities Imean, +2 I+-, +3 Fmean, +4 F+-
    #       < 0 other data,
    #             -1 DelAnom,  -2 just phases (HL), -3 map coefficients
    ACCEPTED_SETS = [[['F_meas_au','F_meas_sigma_au'],'Fmean', +3],
                     [['intensity_meas','intensity_sigma'],'Imean', +1],
                     [['pdbx_F_plus','pdbx_F_plus_sigma', 
                      'pdbx_F_minus','pdbx_F_minus_sigma'], 'F+- anomalous', +4],
                     [['pdbx_I_plus','pdbx_I_plus_sigma',
                      'pdbx_I_minus', 'pdbx_I_minus_sigma'], 'I+- anomalous', +2],
                     [['status'], 'FreeR status', 0],
                     [['pdbx_r_free_flag'], 'FreeR flag', 0],
                     [['pdbx_anom_difference', 'pdbx_anom_difference_sigma'],
                      'DelF anomalous', -1],
                     [['pdbx_HL_A_iso', 'pdbx_HL_B_iso',
                       'pdbx_HL_C_iso', 'pdbx_HL_D_iso'],'HL phases', -2],
                     [['pdbx_FWT', 'pdbx_PHWT',
                       'pdbx_DELFWT', 'pdbx_DELPHWT',
                       'phase_meas',  'fom'],'Map coefficients', -3]
                     ]
    
    ACCEPTED_SETS_UNMERGED = \
             [[['intensity_net','intensity_sigma'],'unmerged I', +1]]

    # Codes <= this are not monitored for missing  items
    TYPECODE_IGNORE_MISSING = -3

    def __init__(self, clabels, mergedtype):
        #print("CIF column labels", clabels, mergedtype)

        # Which set to use
        if mergedtype:
            accsets = CIFLabelSets.ACCEPTED_SETS
        else:
            accsets = CIFLabelSets.ACCEPTED_SETS_UNMERGED

        self.OK  = True  #  OK if all sets are complete
        self.columnsets = []  # list of column sets found

        # I'm sure there is a clever way of doing this, but ...
        for accset in accsets:  # loop possible sets
            ncols = len(accset[0])
            foundcolumns = []
            for cl in clabels:  # loop actual column labels
                if cl in accset[0]:
                    foundcolumns.append(cl)

            if len(foundcolumns) == 0:
                continue

            # At least some found, flag missing ones
            flagpresent = []
            for i in range(ncols):
                present = False
                if accset[0][i] in foundcolumns:
                    present = True
                flagpresent.append(present)

            #print("flagpresent", flagpresent, " accset", accset)
            colset = ColumnSet(accset[0], flagpresent, accset[1], accset[2])
            self.columnsets.append(colset)
            if not colset.OK:
                self.OK = False

    # - - - - - - - - - - - - - - -
    def getTypeCodes(self):
        #   3)  data type code
        #       = 0 freeR,
        #       > 0 suitable for input as observed data
        #             +1 intensities Imean, +2 I+-, +3 Fmean, +4 F+-
        #       < 0 other data,
        #             -1 DelAnom,  -2 just phases (HL), -3 map coefficients
        codes = []
        for cset in self.columnsets:
            codes.append(int(cset.typecode))
        return codes
    # - - - - - - - - - - - - - - -
    def getcolumnsets(self):
        return self.columnsets

    # - - - - - - - - - - - - - - -
    def columnsetstext(self):
        # Return text labels as list of strings
        s = []
        for cset in self.columnsets:
            s.append(cset.text)
        return s
    
    # - - - - - - - - - - - - - - -
    def numcolumnsets(self):
        return len(self.columnsets)

    # - - - - - - - - - - - - - - -
    def columnnames(self):
        names = {}
        if len(self.columnsets) == 0: return names
        for colset in self.columnsets:
            names[colset.text] = colset.columnlabels
        return names
    # - - - - - - - - - - - - - - -
    def formatContent(self):
        # Column content information
        s = ''
        if len(self.columnsets) == 0: return s
        for colset in self.columnsets:
            s += colset.text  
            if not colset.OK and \
                   colset.typecode > CIFLabelSets.TYPECODE_IGNORE_MISSING:
                s += ' (missing ' + colset.missingColumns() + ')'
            s += ', '
        return s[:-2]


class CifBlockInfo:
    # Extract and store information from a Gemmi rblock
    def __init__(self, rblock):
        if not rblock:
            # not a reflection block
            self.OK = False
            return
        
        self.OK = True
        self.rblock = rblock

        self.bname = self.rblock.block.name
        self.entry_id = rblock.entry_id
        self.unmerged = self.rblock.is_unmerged()
        # cell as a list
        self.cell = [rblock.cell.a,rblock.cell.b,rblock.cell.c,
                     rblock.cell.alpha,rblock.cell.beta,rblock.cell.gamma]
        # the name in the file, in case we need it
        self.spacegroup_name_file = rblock.spacegroup.hm
        #  maybe convert eg "R 3" to "R 3 :R"
        self.spacegroup_name = self.fixRhombohedralSpacegroupName()

        wvl = rblock.wavelength
        if math.isnan(wvl):
            wvl = 0.0
        self.wavelength = wvl
        self.havecolumnlabels = True
        try:
            self.columnlabels = rblock.column_labels()
        except:
            self.havecolumnlabels = False

        # Resolution range, tuple of (dmax, dmin)
        da = rblock.make_d_array()
        self.resolutionrange = (max(da), min(da))
        # Highest resolution
        self.highres = self.resolutionrange[1]

        self.categories = rblock.block.get_mmcif_category_names()
        # Reflection loops contain either items
        #   1) _refln.*          usually merged data, True, or 
        #   2) _diffrn_refln.*   unmerged data, False
        self.reflnlooptype = True
        if '_diffrn_refln.' in self.categories:
            self.reflnlooptype = False
        
        self.details = None
        difinfo = self.rblock.block.get_mmcif_category('_diffrn')
        if difinfo:
            self.details = ''
            if 'details' in difinfo:
                self.details = difinfo['details'][0]

        #  Get column label groups
        self.labelsets = CIFLabelSets(self.columnlabels, self.reflnlooptype)
        # formatted column information
        self.info = self.formatInfo()
        # dictionary of column names
        self.columnnames = self.labelsets.columnnames()

        anomalous = False
        if 'anomalous' in self.labelsets.formatContent():
            anomalous = True
        hklcheck = HKLcheck(rblock, anomalous)
        self.hklstatus = hklcheck.hklstatus
        self.hklcheckformat = self.hklstatus.formathklstatus()
        self.nunique = hklcheck.nrefunique() # Number of unique reflections

        # inspect FreeR
        self.freerStatus = None
        self.freerMissing = None
        if self.haveFreeR():
            self.freerStatus = FreerStatus(rblock, self.columnlabels)
            self.freerMissing = self.freerStatus.freermissing
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fixRhombohedralSpacegroupName(self):
        #  Add ":H" or "R" to R if needed
        name = self.spacegroup_name_file
        if name[0] == 'R':
            if ':' not in name:
                # No ':' character, check cell
                #  H lattice has alpha = beta = 90, gamma = 120
                lat = 'R'
                if math.isclose(self.cell[5], 120.0):
                    if math.isclose(self.cell[4], 90.0):
                        if math.isclose(self.cell[3], 90.0):
                            lat = 'H'
            name = name+' :'+lat
        return name
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def rfreeFlagmissing(self):
        return self.freerMissing
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def haveFreeR(self):
        # True if the file has a FreeR column (status or pdbx_r_free_flag)
        if ('status' in self.columnlabels) or \
           ('pdbx_r_free_flag' in self.columnlabels):
            return True
        return False
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def validFreeR(self):
        # Returns: None if no FreeR; True if valid; False if invalid (eg all the same)
        if self.freerStatus is None: return None
        valid = self.freerStatus.valid()
        return valid
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def freerWarning(self):
        # Return warning messages if there are problems with the FreeR data
        # Returns None if there is no FreeR, empty string if no warnings
        if not self.haveFreeR(): return None
        return self.freerStatus.warning()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def standardasu(self):
        return self.hklstatus.standardasu
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def ismerged(self):
        return self.hklstatus.merged
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def nrefunique(self):
        # Number of unique reflections
        return self.nunique

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def formatInfo(self):
        # brief formatted column information
        s = self.labelsets.formatContent()
        return s

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def formatColumnInfo(self):
        # formatted column information in more detail
        s = ''
        # list of ColumnSet objects
        columnsets = self.labelsets.getcolumnsets()
        for colset in columnsets:
            # List of data types and columns
            s += colset.text + ": [" +\
                 colset.formatcolumnlabels() + "], "

        return s[:-2]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def allowsimplewrite(self, reducehkl=True):
        # Is it valid to use a simple CifToMtz procedure?
        # Should be a single asymmetric unit,
        #   and CCP4 standard asu unless reducehkl == False
        if self.freerMissing:
            return False
        if self.hklstatus.oneasu:
            if self.hklstatus.standardasu or (not reducehkl):
                return True
        return False
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def allowanomalouswrite(self):
        # Is it valid to write MTZ generating anomalous pairs from
        if self.hklstatus.anomlinepairs:
            return True
        return False

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def merged_diffn_data(self):
        # return True if block contains merged diffraction data,
        #  ie Fs or Is
        if not self.ismerged():
            return False  # unmerged data

        codes = self.labelsets.getTypeCodes()
        for code in codes:
            if int(code) > 0:
                return True
        return False

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def columnsOK(self, formerged=True):
        # Check column sets for validity for merged or unmerged data
        # Returns status, message
        # status  = 0 OK, -1 not
        #  message list of data types and columns

        # Column sets may be
        #     importable     - typecode >= 0
        #  or not importable - typecode < 0
        #           
        # Cases:
        #  1) all or sole column sets are OK
        #  2) at least one importable column set is complete

        message = ''
        if self.labelsets.numcolumnsets() == 0:
            message += " appears to contain no suitable reflection data items"
            status = -1
        elif self.labelsets.OK:
            # No problems
            status = 0
        else:
            # list of ColumnSet objects
            columnsets = self.labelsets.getcolumnsets()
            status = -1
            for colset in columnsets:
                typecode = colset.typecode
                if typecode >= 0:
                    if colset.flagpresent:
                        # at least one set OK
                        status = +1
            if status < 0:
                message += " has an incomplete set of data items, missing items: "
        return status, message

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
class HKLstatus():
    # Properties from HKLcheck
    def __init__(self, merged, oneasu, standardasu, anomlinepairs,
                 maxmultiplicity):
        self.merged = merged
        self.oneasu = oneasu
        self.standardasu = standardasu       # Not standard CCP4 asu
        # True if we have anomalous pairs on separate lines
        self.anomlinepairs = anomlinepairs
        self.maxmultiplicity = maxmultiplicity

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def formathklstatus(self):
        s = ""
        if self.oneasu:
            if self.standardasu:
                s = "One asymmetric unit, standard CCP4 setting"
            else:
                s = "One asymmetric unit, not standard CCP4 setting"
        else:
            if self.anomlinepairs:
                s = "Anomalous pairs I/F +- on different lines"
            else:
                s = "Unmerged data"
        return s


class HKLcheck():
    # class to check the hkl list from an mmcif block, to diagnose any issues

    def __init__(self, rblock, anomalous=False):
        # Construct from a cif block
        # anomalous = True if columns have explicit I/F+-
        # hkl_list is of type numpy.ndarray
        self.hkl_list = rblock.make_miller_array()

        sg = rblock.spacegroup   # gemmi.SpaceGroup
        asu = gemmi.ReciprocalAsu(sg)
        ops = sg.operations()
        
        # Put reduced indices into a new list
        nrefs = len(self.hkl_list)
        rhkl = [None]*nrefs
        isymcounts = {}
        nplus = 0
        nminus = 0

        for i, hkl in enumerate(self.hkl_list):
            hklisym = asu.to_asu(hkl, ops)
            rhkl[i] = tuple(hklisym[0])       # for sort
            isym = hklisym[1]
            if str(isym) in isymcounts:
                isymcounts[str(isym)] += 1
            else:
                isymcounts[str(isym)] = 0

            # isym odd are I/F+, isym even are I/F-, count acentric only
            if not ops.is_reflection_centric(hkl):
                if isym%2 == 0:
                    nminus += 1
                else:
                    nplus += 1
            
        #print("NPM", nplus, nminus)

        rhkl.sort()
        previoushkl = rhkl[0]
        nsame = 0
        hist = {}  # histogram multiplicity
        ncount = 1
        nunique = 1
        for i in range(1, nrefs):
            if rhkl[i] == previoushkl:
                nsame += 1
                ncount += 1
            else:
                nunique += 1
                previoushkl = rhkl[i]
                if ncount > 0:
                    hist[ncount] = hist.get(ncount, 0) + 1
                    ncount = 1

        self.maxhist =  max(hist);
        self.unique = nunique

        # Assess reflection list
        # Cases:
        #  1. just one asymmetric unit present (nsame == 0)
        #   a) hkl in "standard" CCP4 asu
        #   b) non-standard asu
        #  2. more than one asu present (nsame > 0)
        #   a) I/F+ and - on different lines
        #      no explicit anomalous, but should be generated
        #      Maximum in histogram of multiplicities == 2
        #   b) unmerged data

        merged = True
        oneasu = True
        standardasu = False      # Not standard CCP4 asu
        # True if we have anomalous pairs on separate lines
        anomlinepairs = False

        if nsame == 0:
            oneasu = True
            if len(isymcounts) == 1:
                if list(isymcounts.keys())[0] == '1':
                    # Standard CCP4 asu
                    standardasu = True
        else:
            oneasu = False
            #  Case 2
            #  Do we have explicit anomalous?
            if anomalous:
                #print("explicit anomalous")
                merged = False
            else:
                # treat as unmerged
                merged = False
                if self.maxhist == 2:
                    # Check that the two Isym values are oven and odd
                    isymsum = 0
                    for key in isymcounts:
                        isymsum += int(key)
                    #print("isyms", isymcounts, isymsum)
                    if isymsum%2 == 1:
                        # Sum(Isym) should be odd
                        anomlinepairs = True
                        merged = True
                        if isymsum == 3:
                            # Standard CCP4 asu
                            standardasu = True

        # Store hklstatus
        self.hklstatus = HKLstatus(merged, oneasu, standardasu,
                                   anomlinepairs, self.maxhist)
        #print(self.hklstatus .formathklstatus())

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def nrefunique(self):
        return self.unique  # Number of unique reflections
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def merged(self):
        # False if unmerged
        return self.hklstatus.merged
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def oneasu(self):
        # True if just one asymmetric unit present
        return self.hklstatus.oneasu
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def print(self):
        if 'anomalous' in self.labeltypes:
            print("explicit anomalous")
        else:
            print("No explicit anomalous")
        print(self.hklcheck.formathklstatus())


def getcolumn(rb, tag):
    # extract a column with label including tag
    # return list of strings
    
    s = rb.block.as_string()
    lines = s.split('\n')

    j = 0  # line number
    jloop = -1

    for line in lines:
        # Find tag and record line number of start of loop
        if 'loop_' in line: jloop = j
        if tag in line:
            break
        j += 1

    # Count columns, ie labels in loop_
    ncol = 0
    jtag = -1
    for line in lines[jloop+1:]:
        if line[0] != '_':
            break
        if tag in line:
            jtag = ncol
        ncol += 1

    jlinestart = jloop+ncol+1 # first line of data

    coldata = []   # required data
    for line in lines[jlinestart:]:
        d = line.split()
        if len(d) != ncol:
            break
        coldata.append(d[jtag])
    return coldata


class FreerStatus():
    def __init__(self, rblock, columnlabels):
        # columnlabels contain either 'status' or 'pdbx_r_free_flag'
        self.OK = True
        self.rblock = rblock
        # Do we have both? (must have at least one if we have got here)
        freercolumns = 'rfreeflag'
        if ('status' in columnlabels) and \
               ('pdbx_r_free_flag' in columnlabels):
            freercolumns = 'both'
        elif ('status' in columnlabels):
            freercolumns = 'status'

        self.statussame = None
        self.rfreeflagsame = None
        self.flagssame = None
        # True if some FreeR flags are missing, eg NaN
        self.freermissing = False
                
        #print('freercolumns', freercolumns)
        if freercolumns != 'rfreeflag':
            self.statussame = self.checkStatus()
        if freercolumns != 'status':
            self.rfreeflagsame = self.readRfreeFlag()
        if freercolumns == 'both':
            # print(self.statuslist,  self.freerlist)
            if self.statuslist is not None and self.freerlist is not None:
                if len(self.statuslist) != len(self.freerlist):
                    # Shouldn't happen
                    self.OK = False
                    print("!!! FreerStatus list lengths don't match,",
                          len(self.statuslist), len(self.freerlist))
                    return
                self.flagssame = self.sameFlags()

    def warning(self):
        # Returns list of warning messages if the is a problem, else ''
        s = []
        if self.statussame is not None:
            if self.statussame:
                s.append("WARNING: all reflection status flags are the same")
        if self.rfreeflagsame is not None:
            if self.rfreeflagsame:
                s.append("WARNING: all FreeR flags are the same")
        if self.flagssame is not None:
            if not self.flagssame:
                s.append("WARNING: status flags and FreeR flags do not match")
        if self.freermissing:
            s.append("WARNING: some FreeR flags are flagged as missing")

        return s

    def checkStatus(self):
        # check to see if status flags are all the same
        # returns list of ints, 'f' -> 0, or 'o' -> 1
        self.statuslist = self.make_freer_array()
        # return True if all values are the same
        return self.allTheSame(self.statuslist)

    def make_freer_array(self):
        #  status is a single character string, convert to int
        scol = getcolumn(self.rblock, "status")
        # Check for valid status flag
        s = self.isStatusValid(scol[0])
        self.freerStatusType = s
        if s == 'valid':
            freer = []
            for sc in scol:
                freer.append(self.status_to_freeflag(sc))
            return freer
        elif s == 'integer':
            self.statussame = self.readRfreeFlag('status')
            return self.freerlist

    def isStatusValid(self, status):
        # Check for valid status flag
        vf = self.status_to_freeflag(status)
        s = 'valid'
        if math.isnan(vf):
            # invalid, maybe integer
            try:
                int(status)
                s = 'integer'
            except:
                s = 'invalid'
        return s

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def status_to_freeflag(self, status):
        # Convert status to freerflag, cf gemmi/include/cif2mtz.hpp
        # 'o' -> 1, 'f' -> 0, 'x' -> -1 else NaN
        c = status[0]
        #print("status_to_freeflag status", status)
        if c == '\'' or c == '"':  # remove potential leading junk
            c = status[1]
        if (c == 'o'):
            return 1;
        if (c == 'f'):
            return 0;
        if (c == 'x'):
            return -1;
        return math.nan;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def readRfreeFlag(self, tag='pdbx_r_free_flag'):
        # tag = 'status' or 'rfreeflag'
        # returns list of floats
        self.freerlist = self.rblock.make_float_array(tag)
        for i, x in enumerate(self.freerlist):
            if math.isnan(x):
                self.freerlist[i] = -1
                self.freermissing = True
            else:
                # Convert to int
                self.freerlist[i] = int(x+0.01)
                
        # check to see if RfreeFlags are all the same
        # return True if all values are the same
        return self.allTheSame(self.freerlist)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def allTheSame(self, frlist):
        # Return True if all values are the same, ignore negatives
        if frlist is None: return False
        first = frlist[0]
        for rf in frlist[1:]:
            if rf >= 0 and rf != first:
                return False
        return True
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def sameFlags(self):
        # Return True if values are equivalent in both lists
        #  self.statuslist self.freerlist
        # statuslist is 0 or 1, freerlist is 0 or >0
        #  statuslist may also = NaN if status flag was eg 'x'
        for i in range(len(self.statuslist)):
            if not math.isnan(self.statuslist[i]):
                if self.freerlist[i] > 0:
                    if self.statuslist[i] == 0:
                        return False  # Flag > 0, status = 0
                else:
                    if self.statuslist[i] != 0:
                        return False  # Flag = 0, status = 1
        return True
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def format(self):
        # Return formatted assessment
        s = []
        header = ""
        if self.statussame is not None:
            header = "status flag"
            if self.statussame:
                s.append("WARNING: all reflection status flags are the same")
            else:
                s.append("Reflection status flags are not all the same")
        if self.rfreeflagsame is not None:
            header = "Rfree flag number"
            if self.rfreeflagsame:
                s.append("WARNING: all FreeR flags are the same")
            else:
                s.append("FreeR flags are not all the same")
        if self.flagssame is not None:
            header = "both status flag and Rfree flag number"
            if self.flagssame:
                s.append("Status flags and FreeR flags match")
            else:
                s.append("WARNING: status flags and FreeR flags do not match")

        s.insert(0,"Block contains "+header)
        return s
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def valid(self):
        # Return True if valid ie not all the same
        #   else False
        #   unmatched is treated as valid
        if self.statussame is not None:
            if self.statussame:
                return False
        if self.rfreeflagsame is not None:
            if self.rfreeflagsame:
                return False
        return True

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
def printLines(slist):
    # Print a list of strings sensibly
    WIDTH = 76  # characters
    sl = " "
    n = 1
    for i, string in enumerate(slist):
        n += len(string) + 2
        if n > WIDTH:
            sl += '\n  '
            n = 1
        sl += string
        if i < len(slist)-1:
            sl += ', '
    print(sl)

# - - - - - - - - - - - - - - -
def printBlockInfo(cifblockinfo):
    #  Print block info
    print("\n* * * * *")
    print("Information for entryID",cifblockinfo.entry_id,
          "from reflection block ", cifblockinfo.bname)

    if cifblockinfo.unmerged:
        print("Unmerged data")
    else:
        print("Merged data")

    print("Unit cell:", *cifblockinfo.cell)
    print("Space group:", cifblockinfo.spacegroup_name)
    print("Wavelength:", cifblockinfo.wavelength)
    print("Resolution range:",
          "{:.3f} to {:.3f}".format(cifblockinfo.resolutionrange[0],
                                  cifblockinfo.resolutionrange[1]))
            
    if cifblockinfo.reflnlooptype:
        looptype = '"_refln"'
    else:
        looptype = '"_diffrn_refln"'
    print("Column labels, ", len(cifblockinfo.columnlabels),
          " columns, reflection loop type", looptype)
    printLines(cifblockinfo.columnlabels)

    if cifblockinfo.details is None:
        print("No _diffrn block")
    else:
        if cifblockinfo.details != '':
            print("Details:", cifblockinfo.details)
        else:
            print("No details")

    print("Data contents:", cifblockinfo.labelsets.formatContent())

    print("Assessment of hkl list:",cifblockinfo.hklcheckformat)

    frs = cifblockinfo.freerStatus
    if frs is not None:
        print("Assessment of FreeR flags:") 
        s = frs.format()
        for l in s:
            print(l)
        if len(frs.warning()) > 0:
            for s in frs.warning():
                print(s)
