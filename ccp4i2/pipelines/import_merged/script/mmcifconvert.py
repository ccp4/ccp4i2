import math
import sys

from lxml import etree
import gemmi
import numpy

from .importutils import addXMLelement, ReflectionDataTypes
from .mmcifutils import CifBlockInfo, printBlockInfo


class GetColumn():
    def __init__(self):
        pass

    def getcolumn(self, rb, tag):
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

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
class Taggedhkl():
    # reduced tuple(hkl), index into original list, isym value
    def __init__(self, rhkl, idx, isym):
        self.rhkl = rhkl
        self.idx = idx
        self.isym = isym

    def format(self):
        return self.rhkl, self.idx, self.isym
        #return str(self.rhkl), str(self.idx), str(self.isym)
    
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
class CIFReflectionData:
    # class to contain reflection data from an mmcif file

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__(self, rblock):
        # Construct from a cif block
        self.rblock = rblock
        # hkl_list is of type numpy.ndarray
        self.hkl_list = rblock.make_miller_array()

        self.cifblockinfo = CifBlockInfo(rblock)

        self.labelsets = self.cifblockinfo.labelsets  # CIFLabelSets
        self.labeltypes = self.labelsets.columnsetstext()

        self.outfile = 'job_obsout.mtz'
        self.freerfile = 'job_freer.mtz'
        self.dataout = True
        # Resolution limits are in self.cifblockinfo.resolutionrange
        # Start XML report
        self.convertXML = etree.Element('MMCIF_CONVERT')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def setMTZnames(self, outfile, freerfile):
        self.outfile = outfile
        if outfile is None: self.dataout = False
        self.freerfile = freerfile

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def gettype(self, typelist):
        # Loop possible types in order of priority
        # Return first one found in file
        for labtype in typelist:
            if labtype in self.labeltypes:
                foundtype = labtype
                break            
        return foundtype

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def gettypes(self, typelist):
        # Loop possible types in order of priority
        # Return list of found ones
        foundtypes = []
        for labtype in typelist:
            if labtype in self.labeltypes:
                foundtypes.append(labtype)
        return foundtypes

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def gettypesinfile(self, typelist):
        # Return all relevant types found in file
        foundtypes = []
        for labtype in typelist:
            if labtype in self.labeltypes:
                foundtypes.append(labtype)
        return foundtypes

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def simpleMTZwrite(self, cifdatatype, haveFreeR=True, resorange=None):
        # use CifToMtz to just create MTZ file of data
        #  either just the data itself (Is or Fs) if haveFreeR is False
        #    else both data and FreeR
        # hkl will NOT be reduced to the standard CCP4 asymmetric unit
        #  other data are unchanged

        print(">> simpleMTZwrite, mtztype", cifdatatype)
        #print(ReflectionDataTypes.REFLECTION_DATA[cifdatatype])

        self.cifdatatype = cifdatatype
        freerspecs = None

        if haveFreeR:
            # FreeR
            freertype = self.gettype(ReflectionDataTypes.FREER_TYPES)
            freerspecs = ReflectionDataTypes.REFLECTION_DATA[freertype]
            if freerspecs is None:
                print("No valid FreeR set found")
                return False
            
        # Main data
        spec_lines = None
        if self.dataout:
            # I or F data
            #  Check for all types
            if cifdatatype not in ReflectionDataTypes.REFLECTION_DATA:
                cifdatatype = self.gettype(ReflectionDataTypes.DATA_PRIORITY)
            self.cifdatatype = cifdatatype
            spec_lines = ReflectionDataTypes.REFLECTION_DATA[cifdatatype]
            if spec_lines is None: return False
            #print("cifdatatype", cifdatatype,",  spec_lines", spec_lines)
            conv = gemmi.CifToMtz()
            conv.spec_lines = spec_lines
            mtz = conv.convert_block_to_mtz(self.rblock)

            # Set dataset info
            dts = mtz.datasets[1]
            dts.project_name = self.cifblockinfo.entry_id
            dts.crystal_name = self.cifblockinfo.bname
            dts.dataset_name = self.cifblockinfo.bname
            mtz.datasets[1] = dts
            if resorange is not None:
                mtz = self.setMtzResolutionLimits(mtz, resorange)
            highres = min(mtz.make_d_array())
            mtz.write_to_file(self.outfile)
            print("File ", self.outfile, " written")
            
        # Now FreeR
        if haveFreeR:
            conv = gemmi.CifToMtz()
            conv.spec_lines = freerspecs
            mtz = conv.convert_block_to_mtz(self.rblock)
            mtz.datasets.pop()
            if resorange is not None:
                mtz = self.setMtzResolutionLimits(mtz, resorange)
            mtz.write_to_file(self.freerfile)
            print("File ", self.freerfile, " written")

        self.storeColumnNames(spec_lines, freerspecs)
        reducemessage = self.makereducemessage(False)
        nreduced = 0
        self.addtoXML(cifdatatype, len(self.hkl_list),
                      'ciftomtz convert', reducemessage,
                      haveFreeR, nreduced, resorange)

        return True
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def setMtzResolutionLimits(self, mtz, resorange):
        dmin = resorange[1]
        if resorange[0] > 0.0:
            dmax = resorange[0]
            mtz.set_data(mtz.array[mtz.make_d_array() >= dmin])
            mtz.set_data(mtz.array[mtz.make_d_array() <= dmax])
        else:
            mtz.set_data(mtz.array[mtz.make_d_array() >= dmin])
        return mtz
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def outputMTZ(self, mtz, data, outfile, resorange):
        # Set data into mtz, trim if necessary, write to outfile
        mtz.set_data(data)
        if resorange is not None:
            mtz = self.setMtzResolutionLimits(mtz,resorange)
        mtz.write_to_file(outfile)
        
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def makereducemessage(self, reducehkl=True):
        reducemessage = ''
        if not self.cifblockinfo.standardasu():
            # Not standard asu
            if  reducehkl:
                reducemessage = 'hkl reduced to standard asymmetric unit'
            else:
                reducemessage = 'no hkl change from non-standard asymmetric unit'
        return reducemessage

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def samehkl(self, hkl1, hkl2):
        if hkl1 == tuple(hkl2):
            return True
        return False
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def anomMTZwritedata(self, cifdatatype,
                         reducehkl=True, haveFreeR=True, resorange=None):
        # disentangle anomalous data with I/F +- on different lines
        # This type of file seems to have I/F +- on successive lines,
        #  but we can't rely on that
        #
        #  cifdatatype 
        #  reducehkl=True to reduce to standard CCP4 asu
        #  haveFreeR=True to write FreeR file as well as data
        #
        # We have on entry:
        # self.rblock, self.hkl_list, self.cifblockinfo,
        # self.labelsets, self.labeltypes
        # print("anomMTZwritedata")

        if self.dataout:
            # Normal data
            if cifdatatype not in ReflectionDataTypes.DATA_PRIORITY:
                cifdatatype = self.gettype(ReflectionDataTypes.DATA_PRIORITY)
            i = ReflectionDataTypes.DATA_PRIORITY.index(cifdatatype)
            #  Equivalent anomalous type
            anomcifdatatype = ReflectionDataTypes.DATA_PRIORITY[i-1]
            self.cifdatatype = anomcifdatatype

            if 'anomalous' in cifdatatype:
                print("Explicit anomalous not allowed for anomMTZwrite")
                return False

            # cifdatatype will be Imean or Fmean, or FreeR
            # each spec contains "ciflabel MTZlabel MTZtype 1"
            cifspecs = ReflectionDataTypes.REFLECTION_DATA[cifdatatype]
            if len(cifspecs) != 2:
                print("Wrong specs", cifspecs)
                return False
            # Data arrays, I or F, sigma
            cifdata  = self.rblock.make_float_array(cifspecs[0].split()[0])
            cifsigma = self.rblock.make_float_array(cifspecs[1].split()[0])

            # Check sizes
            nref = len(self.hkl_list)
            if len(cifdata) != nref or len(cifsigma) != nref:
                print("Wrong array lengths", nref, len(cifdata), len(cifsigma))
                return False


        freerspecs = None
        if haveFreeR:
            # Loop possible FreeR types in order of priority
            freertype = self.gettype(ReflectionDataTypes.FREER_TYPES)
            if freertype is None:
                print("No valid set found")
                return False
            # freertype will be 'FreeR flag' or 'FreeR status'
            # each spec contains "ciflabel MTZlable MTZtype 1"
            freerspecs = ReflectionDataTypes.REFLECTION_DATA[freertype]
            # FreeR
            tag = freerspecs[0].split()[0]
            if tag == 'status':
                #  status is a single character string, convert to float
                ciffreer = self.make_freer_array(tag)
            else:
                ciffreer  = self.rblock.make_float_array(tag)

        taggedlist = self.maketaggedlist()

        mtzspecs = None
        nrefunique = self.cifblockinfo.nrefunique()
        if self.dataout:
            # Now for the main cif data
            mtzspecs = ReflectionDataTypes.REFLECTION_DATA[self.cifdatatype]
            ncolumns = len(mtzspecs) + 3
            data = numpy.zeros((nrefunique, ncolumns), dtype=numpy.float)
            dataline0 = [math.nan]*ncolumns  # Missing items set to NaN
            dataline = list(dataline0)
        if haveFreeR:
            freerdata = numpy.zeros((nrefunique, 4), dtype=numpy.float)
            frdataline0 = [math.nan]*4
            frdataline = list(frdataline0)

        self.storeColumnNames(mtzspecs, freerspecs)

        juniq = 0
        previoushkl = None
        nreduced = 0
        ops = self.rblock.spacegroup.operations()

        for i, taggedhkl in enumerate(taggedlist):
            # Loop reflections
            idx = taggedhkl.idx    # index into original list
            isym = taggedhkl.isym  # isym, odd for I/F+
            if reducehkl:
                rhkl = taggedhkl.rhkl  # Unique hkl
                if not self.samehkl(rhkl, self.hkl_list[idx]):
                    nreduced += 1            
            else:
                rhkl = self.hkl_list[idx]
                            
            if previoushkl is None:
                previoushkl = rhkl
            elif rhkl != previoushkl:
                # New hkl, output old one
                if haveFreeR:
                    frdataline[3] = ciffreer[idx]
                    freerdata[juniq] = frdataline
                    frdataline = list(frdataline0)

                if self.dataout:
                    # print("*", juniq, idx,rhkl,  isym, dataline)
                    data[juniq] = dataline
                    dataline = list(dataline0)

                juniq += 1
                previoushkl = rhkl

            if haveFreeR:
                frdataline[0] = rhkl[0]
                frdataline[1] = rhkl[1]
                frdataline[2] = rhkl[2]
            if self.dataout:
                dataline[0] = rhkl[0]
                dataline[1] = rhkl[1]
                dataline[2] = rhkl[2]
                IFdata = cifdata[idx]
                IFsigma = cifsigma[idx]
                # Insert data
                if ops.is_reflection_centric(rhkl):
                    # Centric + = -
                    dataline[3] = IFdata
                    dataline[4] = IFsigma
                    dataline[5] = IFdata
                    dataline[6] = IFsigma
                else:
                    j = 3+2*((isym+1)%2)   # j = 3 or 5
                    dataline[j] = IFdata
                    dataline[j+1] = IFsigma

        # End loop reflections

        if juniq < nrefunique:
            # Add last one if not done
            if self.dataout: data[juniq] = dataline
            if haveFreeR: freerdata[juniq] = frdataline

        # print("Unique reflections", juniq)
        #print("nreduced", nreduced)

        if self.dataout:
            # Start to construct the MTZ file for main data
            mtz = gemmi.Mtz(with_base=True)  # with h,k,l
            self.start_mtz(mtz, mtzspecs)
            self.outputMTZ(mtz, data, self.outfile, resorange)
            print("File ", self.outfile, " written, number of reflections",
                  len(data))

        # and the FreeR file
        if haveFreeR:
            # Start to construct the MTZ file
            mtz = gemmi.Mtz(with_base=True)  # with h,k,l
            mtz.spacegroup = self.rblock.spacegroup
            mtz.set_cell_for_all(self.rblock.cell)
            mtz.history = ['From cifconvert, blockname: '+ self.cifblockinfo.bname]
            for spec in freerspecs:
                #  MTZ column label, column type
                self.addMTZcolumnLabelType(mtz, spec)

            self.outputMTZ(mtz, freerdata, self.freerfile, resorange)
            print("File ", self.freerfile, " written, number of reflections",
                  len(freerdata))

        reducemessage = self.makereducemessage(reducehkl)
        self.addtoXML(self.cifdatatype, juniq,
                      'combine anomalous lines',
                      reducemessage, haveFreeR, nreduced, resorange)

        return True

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def nonstandardMTZwritedata(self, cifdatatype,
                                reducehkl=True, haveFreeR=True,
                                freeRmissing=False, resorange=None):
        # Write MTZ from cif,
        #  reducing hkl to standard CCP4 asu if reducehkl True
        # Loop possible types in order of priority
        # 
        # We have on entry:
        # self.rblock, self.hkl_list, self.cifblockinfo,
        # self.labelsets, self.labeltypes

        print(">>> nonstandardMTZwritedata", cifdatatype)
        self.cifdatatype = cifdatatype
        mtzspecs = None
        nref = len(self.hkl_list)

        if self.dataout:
            if cifdatatype not in ReflectionDataTypes.REFLECTION_DATA:
                cifdatatype = self.gettype(ReflectionDataTypes.DATA_PRIORITY)
                self.cifdatatype = cifdatatype

            anomalous = False
            if 'anomalous' in cifdatatype:
                anomalous = True

            mtzspecs = ReflectionDataTypes.REFLECTION_DATA[cifdatatype]

            ciftags = [cs.split()[0] for cs in mtzspecs]
            ndatacols = len(ciftags)

            cifdata = [self.rblock.make_float_array(tag) for tag in ciftags]

            for i in range(ndatacols):
                if len(cifdata[i]) != nref:
                    print("Wrong array length", nref, i, len(cifdata[i]))
                    return False

        freerspecs = None
        if haveFreeR:
            # Loop possible FreeR types in order of priority
            freertype = self.gettype(ReflectionDataTypes.FREER_TYPES)
            if freertype is None:
                print("No valid set found")
                return False
            # freertype will be 'FreeR flag' or 'FreeR status'
            # each spec contains "ciflabel MTZlable MTZtype 1"
            freerspecs = ReflectionDataTypes.REFLECTION_DATA[freertype]
            # FreeR
            tag = freerspecs[0].split()[0]
            if tag == 'status':
                #  status is a single character string, convert to float
                ciffreer = self.make_freer_array(tag)
            else:
                ciffreer  = self.rblock.make_float_array(tag)

        self.storeColumnNames(mtzspecs, freerspecs)

        taggedlist = self.maketaggedlist()

        # Now for the main cif data
        nrefunique = self.cifblockinfo.nrefunique()
        if self.dataout:
            ncolumns = len(mtzspecs) + 3
            data = numpy.zeros((nrefunique, ncolumns), dtype=numpy.float)
            dataline = [math.nan]*ncolumns
            anomswap = [2,3, 0,1]  # indices to swap + and - if needed
        if haveFreeR:
            freerdata = numpy.zeros((nrefunique, 4), dtype=numpy.float)
            frdataline = [math.nan]*4

        ops = self.rblock.spacegroup.operations()

        juniq = 0

        nswaplist = 10
        nswapped = 0
        nreduced = 0
        nallnan = 0
        nobsnanfrvalid = 0
        
        for i, taggedhkl in enumerate(taggedlist):
            # Loop reflections
            idx = taggedhkl.idx    # index into original list
            isym = taggedhkl.isym  # isym, odd for I/F+
            if reducehkl:
                rhkl = taggedhkl.rhkl  # Unique hkl
                if self.samehkl(rhkl, self.hkl_list[idx]):
                    nreduced += 1            
            else:
                rhkl = self.hkl_list[idx]  # original hkl

            if self.dataout:
                dataline[0] = rhkl[0]
                dataline[1] = rhkl[1]
                dataline[2] = rhkl[2]

                d = [math.nan]*ndatacols
                allnan = True
                for i in range(ndatacols):
                    d[i] = cifdata[i][idx]
                    if not math.isnan(cifdata[i][idx]):
                        allnan = False
                #print("*",rhkl, idx, d)
                if allnan:
                    nallnan += 1
                    if haveFreeR and not math.isnan(ciffreer[idx]):
                        nobsnanfrvalid += 1
                        #    print("*allNan",rhkl, idx, d, ciffreer[idx])
                    continue  # no output
                
                if anomalous:
                    if not ops.is_reflection_centric(rhkl):
                        doswap = (isym%2 == 0)
                        if doswap and nswaplist > 0:
                            print("Swapping plus and minus, reflection",
                                  rhkl,  " changed from", self.hkl_list[idx])
                            nswaplist -= 1
                            if nswaplist == 0:
                                print("... additional swaps omitted")
                        for i in range(ndatacols):
                            # Swap plus and minus if original hkl has changed hand
                            if doswap:
                                dataline[3+i] = d[anomswap[i]]
                            else:
                                dataline[3+i] = d[i]
                        if doswap:
                            nswapped += 1
                    else:
                        for i in range(ndatacols):
                            dataline[3+i] = d[i]
                        
                    # two pairs of columns 3,4 and 5,6, set Nan if zero
                    self.settonan(dataline, 5, 2)
                else:
                    for i in range(ndatacols):
                        dataline[3+i] = d[i]

                self.settonan(dataline, 3, 2)
                
                #print("*", juniq, idx,rhkl, juniq,  dataline)
                data[juniq] = dataline

            if haveFreeR:
                frdataline[0] = rhkl[0]
                frdataline[1] = rhkl[1]
                frdataline[2] = rhkl[2]
                frdataline[3] = ciffreer[idx]
                freerdata[juniq] = frdataline

            juniq += 1
    
        if nswapped > 0:
            print(nswapped,
                  "reflections have swapped anomalous plus and minus")
            addXMLelement(self.convertXML, 'nanomswapped', str(nswapped))

        if nallnan > 0:
            addXMLelement(self.convertXML, 'OmittedMissing', str(nallnan))
            if nobsnanfrvalid > 0:
                addXMLelement(self.convertXML, 'OmittedMissingValidFreeR',
                              str(nobsnanfrvalid))
            print("Omitted reflections with all values missing:", nallnan)
            print("  Number of these with valid FreeRflag:", nobsnanfrvalid)

        #print("nreduced", nreduced)
            
        if self.dataout:
            # Start to construct the MTZ file
            mtz = gemmi.Mtz(with_base=True)  # with h,k,l
            self.start_mtz(mtz, mtzspecs)
            self.outputMTZ(mtz, data, self.outfile, resorange)
            print("File ", self.outfile, " written, number of reflections",
                  len(data))

        # and the FreeR file
        if haveFreeR:
            # Start to construct the MTZ file
            mtz = gemmi.Mtz(with_base=True)  # with h,k,l
            mtz.spacegroup = self.rblock.spacegroup
            mtz.set_cell_for_all(self.rblock.cell)
            mtz.history = ['From cifconvert, blockname: '+ self.cifblockinfo.bname]
            for spec in freerspecs:
                #  MTZ column label, column type
                self.addMTZcolumnLabelType(mtz, spec)

            self.outputMTZ(mtz, freerdata, self.freerfile, resorange)
            print("File ", self.freerfile, " written,  number of reflections",
                  len(freerdata))

        reducemessage = self.makereducemessage(reducehkl)
        message = 'non-standard asymmetric unit'
        if freeRmissing:
            message = "missing FreeR values"
        self.addtoXML(cifdatatype, juniq, message,
                      reducemessage, haveFreeR, nreduced, resorange)
        return True
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def settonan(self, dataline, i1, n):
        if abs(dataline[i1+1]) < 1.0e-10:
            for i in range(i1, i1+n):
                dataline[i] = math.nan
        return dataline
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def maketaggedlist(self):
        # Make hkl list for sorting, to give access to I/F+- mates
        sg = self.rblock.spacegroup   # gemmi.SpaceGroup
        ops = sg.operations()
        asu = gemmi.ReciprocalAsu(sg)

        #  Taggedhkl(rhkl, idx, isym)
        taggedlist = []
        for i, hkl in enumerate(self.hkl_list):
            hklisym = asu.to_asu(hkl, ops)
            #            if not self.samehkl(tuple(hklisym[0]), hkl):
            #                print("maketaggedlist", hklisym[0], hkl)
        
            taggedlist.append(Taggedhkl(tuple(hklisym[0]), i, hklisym[1]))

        taggedlist.sort(key=lambda x: x.rhkl)
        return taggedlist

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def start_mtz(self, mtz, specs):
        mtz.spacegroup = self.rblock.spacegroup
        mtz.set_cell_for_all(self.rblock.cell)
        mtz.add_dataset(self.cifblockinfo.bname)
        dts = mtz.datasets[1]
        dts.project_name = self.cifblockinfo.entry_id
        dts.crystal_name = self.cifblockinfo.bname
        dts.dataset_name = self.cifblockinfo.bname
        dts.wavelength =  self.cifblockinfo.wavelength
        mtz.datasets[1] = dts

        mtz.history = ['From cifconvert, blockname: '+ self.cifblockinfo.bname]
        
        # Add columns
        for spec in specs:
            #  MTZ column label, column type
            self.addMTZcolumnLabelType(mtz, spec)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def make_freer_array(self, tag):
        #  status is a single character string, convert to int
        gc = GetColumn()
        scol = gc.getcolumn(self.rblock, tag)
        freer = []
        for sc in scol:
            freer.append(self.status_to_freeflag(sc))
        return freer

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def  status_to_freeflag(self, status):
        # Convert status to freerflag, cf gemmi/include/cif2mtz.hpp
        # Anything other than 'o' or 'f' is set as NaN
        c = status[0]
        if c == '\'' or c == '"':  # remove potential leading junk
            c = status[1]
        if (c == 'o'):
            return 1.0;
        if (c == 'f'):
            return 0.0;
        return math.nan;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getMTZcol(self, spec):
        if spec is None: return spec
        s = ''
        for sp in spec:
            s += sp.split()[1] + ","
        return s[:-1]  # remove trailing comma

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def storeColumnNames(self, mtzspecs, freerspecs):
        # Send output column names to XML
        mtzcols = self.getMTZcol(mtzspecs)      # may be None
        freercols = self.getMTZcol(freerspecs)  # may be None
        
        # Output columns to XML
        if mtzcols is not None:
            addXMLelement(self.convertXML, 'outputcolumnnames',  mtzcols)
        if freercols is not None:
            addXMLelement(self.convertXML, 'freercolumnname',  freercols)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def addMTZcolumnLabelType(self, mtz, spec):
        #  MTZ column label, column type, fix Freer type
        clabel = spec.split()[1]
        ctype = spec.split()[2]
        if ctype == 's':
            ctype = 'I'
        mtz.add_column(clabel, ctype)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def addtoXML(self, datatype, nunique, typemessage,
                 reducehklmessage='', haveFreeR=True, nreduced=0,
                 resorange=None):
        #  Record stuff in XML
        if datatype is not None:
            addXMLelement(self.convertXML, 'datatype', datatype)
            addXMLelement(self.convertXML, 'nrefoutput', str(nunique))
        addXMLelement(self.convertXML, 'message', typemessage)
        if reducehklmessage != '':
            addXMLelement(self.convertXML,
                          'reducemessage', reducehklmessage)
        # Resolution
        rrngxml = self.resorangeXML(self.cifblockinfo.resolutionrange)
        self.convertXML.append(rrngxml)
        if resorange is not None:
            if resorange[0] > 0.0 or resorange[1] > 0.0:
                rrngxml = self.resorangeXML(resorange, 'cutresolution')
                self.convertXML.append(rrngxml)

        if haveFreeR:
            addXMLelement(self.convertXML,
                          'freeRwritten', str(nunique))
            # List of all FreeRtypes types found in file,
            #  ie pdbx_r_free_flag and/or status
            freerinfile = self.gettypesinfile(ReflectionDataTypes.FREER_TYPES)
            s = ''
            for ftype in freerinfile:
                freerspecs = ReflectionDataTypes.REFLECTION_DATA[ftype]
                s += freerspecs[0].split()[0] + ' '
            addXMLelement(self.convertXML,
                          'FreeRinFile', s[:-1])
            # Type used for MTZ freeR
            ftype = self.gettype(ReflectionDataTypes.FREER_TYPES)
            freerspecs = ReflectionDataTypes.REFLECTION_DATA[ftype]
            addXMLelement(self.convertXML,
                          'FreeRused', freerspecs[0].split()[0])

        if nreduced != 0:
            addXMLelement(self.convertXML,
                          'numreducedhkl', str(nreduced))
                      
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def resorangeXML(self, resorange, id=None):
        # XML version of resolution range
        resoxml = etree.Element('ResolutionRange')
        if id is not None:
            resoxml.set('id', id)
        if resorange[0] > 0.0:
            addXMLelement(resoxml, 'min', "{:.3f}".format(resorange[0]))
        if resorange[1] > 0.0:
            addXMLelement(resoxml, 'max', "{:.3f}".format(resorange[1]))
        return resoxml

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getXML(self):
        # Return XML report
        return self.convertXML

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def contentFlag(self):
        contentflag = ReflectionDataTypes.CONTENT_FLAGS[self.cifdatatype]
        return contentflag

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
class ConvertCIF():
    def __init__(self, inputfilename, blockname, cifdatatype,
                 outfile, freerfile,
                 reducehkl=True, resorange=None):

        print("ConvertCIF", cifdatatype)
        mmcif = gemmi.cif.read(inputfilename)

        rblocks = gemmi.as_refln_blocks(mmcif)
        nblocks = len(rblocks)

        if nblocks == 0:
            print("File ", inputfilename, " contains no reflection blocks")
            sys.exit()

        print("\nEntry ", rblocks[0].entry_id,
              "contains", nblocks, "reflection blocks")

        idxblock = 0
        rblock = rblocks[0]
        for rb in rblocks:
            blkinfo = CifBlockInfo(rb)

            bname = blkinfo.bname
            if (blockname is not None) and (bname == blockname):
                rblock = rb
                break

        blkinfo = CifBlockInfo(rblock)
        printBlockInfo(blkinfo)

        refdata = CIFReflectionData(rblock)
        refdata.setMTZnames(outfile, freerfile)

        # Do we have a FreeR column?  (status or pdbx_r_free_flag)
        haveFreeR = blkinfo.haveFreeR()

        # Data types for cif and mtz names
            
        # Cases:
        #  1. just one asymmetric unit present
        #   a) hkl in "standard" CCP4 asu
        #      Use simple cif2mtz call, for data and FreeR
        #   b) non-standard asu
        #      Reduce hkl to standard CCP4 asu unless reducehkl == False
        #  2. more than one asu present (nsame > 0)
        #   a) I/F+ and - on different lines
        #      no explicit anomalous, but should be generated
        #   b) unmerged data
        #      Can't do conversion, bail out

        # Can we do simple write (using CifToMtz)?
        self.status = False
        if blkinfo.allowsimplewrite(reducehkl):
            print("simple write")
            # Now write the main data
            self.status = refdata.simpleMTZwrite(cifdatatype,
                                                 haveFreeR, resorange)
        else:
            # Other cases need special treatment
            if not blkinfo.ismerged():
                # Unmerged data
                print("Can't convert unmerged data")
            elif blkinfo.allowanomalouswrite():
                print("anom pairs")
                # Data with anomalous pairs on different lines
                self.status = refdata.anomMTZwritedata(cifdatatype,
                                                       reducehkl,
                                                       haveFreeR, resorange)
            else:
                if (not blkinfo.standardasu()) or blkinfo.rfreeFlagmissing():
                    print("Non-standard asu or rfreeFlagmissing")
                    # Data in non-standard asu, reduce hkl
                    self.status = \
                       refdata.nonstandardMTZwritedata(cifdatatype,
                                        reducehkl, haveFreeR,
                                        blkinfo.rfreeFlagmissing(), resorange)
                else:
                    # Shouldn't happen
                    print("Shouldn't happen: standard asu and more than one asu")

        self.contentflag = refdata.contentFlag()
        self.XMLreport = refdata.getXML()
        #print(etree.tostring(self.XMLreport, pretty_print=True))

        if not self.status:
            print("Data write failed")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def contentFlag(self):
        return self.contentflag

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getstatus(self):
        return self.status

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getXML(self):
            return self.XMLreport
