import sys
import os
import argparse

from ccp4i2.core.mgimports import mmdb2
import ccp4srs

from gemmi import cif

from .dictFileToMonomer import dictFileToMonomer

def changeDictionaryAtomNames(doc, changes):
    for block in doc:
        for item in block:
            if item.loop is None:
                continue
            for tag in item.loop.tags:
                if "atom_id" in tag:
                    col = block.find_values(tag)
                    for i, atomId in enumerate(col):
                        if atomId in changes:
                            newId = changes[atomId]
                            if "'" in newId or '"' in newId:
                                col[i] = cif.quote(newId)
                            else:
                                col[i] = newId

def replaceMatchesInDict(matches,theDict,outfile):
    #Gemmi returns the vaules in the cif document with no leading/trailing whitespace.
    matches = {k.strip(): v.strip() for k, v in matches.items()}
    try:
        doc = cif.read_file(theDict)
        changeDictionaryAtomNames(doc,matches)
        doc.write_file(outfile)
    except Exception as e:
        print("Failed to change CIF dictionary names %s" % e)

def initSRS():
    SRS = ccp4srs.Manager()
    SRS.loadIndex(os.path.join(os.environ["CCP4"],"share","ccp4srs"))

    for i in range(mmdb2.nAminoacidNames):
        SRS.loadStructure(mmdb2.getAAProperty(mmdb2.AAProperties,i).name)
    for i in range(mmdb2.nNucleotideNames):
        SRS.loadStructure(mmdb2.getPstr(mmdb2.NucleotideName,i))

    return SRS

def getMatches(g1,g2,reverse=False,tryAllSizes=False):

    # This seems to be necessary or else things get broken .... Hmm.
    g1.RemoveChirality()
    g2.RemoveChirality()

    U = ccp4srs.GraphMatch();
    U.SetTimeLimit(2)

    g1.Build(False)
    g2.Build(False)
    print("g1.GetNofVertices()",g1.GetNofVertices())
    print("g2.GetNofVertices()",g2.GetNofVertices())
    g1.Print();
    print()
    g2.Print();
    print()

    if reverse:
        maxAtoms = g2.GetNofVertices()
    else:
        maxAtoms = g1.GetNofVertices()

    if tryAllSizes:
        while maxAtoms > 3:
            U.MatchGraphs(g1,g2,maxAtoms,True,ccp4srs.EXTTYPE_Ignore)
            if U.GetNofMatches() > 0:
                break
            maxAtoms -= 1
    else:
        U.MatchGraphs(g1,g2,maxAtoms,True,ccp4srs.EXTTYPE_Ignore)

    return U

def getBestMatch(U):
    mindist = sys.maxsize
    minMatch = -1
    for i in range(U.GetNofMatches()):
        match = U.GetMatch(i);
        dist = 0
        for j in range(len(match[0])):
            dist += abs(match[0][j] - match[1][j])
        if dist < mindist:
            mindist = dist
            minMatch = i
    if minMatch > -1:
        return U.GetMatch(minMatch)

def matchAtoms(ifname,ofname=None,dictMatchName=None,selection="",dictFileName=None):
    if dictFileName is not None:
        SRS = ccp4srs.Manager()
    else:
        SRS = initSRS()
    
    mmdb2.InitMatType()
    molHnd = mmdb2.Manager()
    RC = molHnd.ReadCoorFile(ifname)

    if selection == "" or selection == "all" or selection is None:
        # Default chain/res/atom
        selection = "*/*/*"

    selHnd = molHnd.NewSelection()
    selindexp = mmdb2.intp()
    print("Selection",selection)
    molHnd.Select(selHnd,mmdb2.STYPE_RESIDUE,selection,mmdb2.SKEY_NEW)
    selRes = mmdb2.GetResidueSelIndex(molHnd,selHnd,selindexp)
    print("Selected",selindexp.value(),"residues")
    if selindexp.value()>1:
        print("Selected more that 1 residue, will use the first.")
    if selindexp.value()==0:
        print("Selected no residues, giving up.")
        return {}

    retMatches = {}

    res = mmdb2.getPCResidue(selRes,0)
    gPDB = ccp4srs.Graph()
    gPDB.MakeGraph(res)
    
    excludeH = False
    
    testFullDB = True
    
    if dictMatchName is not None:
        testFullDB = False
    
    if dictFileName is not None:
        testFullDB = False

    if not testFullDB:
        if dictFileName is not None:
            monomer = dictFileToMonomer(dictFileName)
            if monomer is None:
                return {}

            rcSRS = mmdb2.intp()
            gSRS = monomer.getGraph(rcSRS)

            if gSRS is None:
                return {}
            
        else:
            monomer = SRS.getMonomer(dictMatchName)
            rcSRS = mmdb2.intp()
            gSRS = monomer.getGraph(rcSRS)
        
        if gSRS and gPDB:
            if excludeH:
                gPDB.ExcludeType(mmdb2.getElementNo(("H")))
                gSRS.ExcludeType(mmdb2.getElementNo(("H")))

        if gPDB.GetNofVertices() <= gSRS.GetNofVertices():
            U = getMatches(gPDB,gSRS,tryAllSizes=True)
        else:
            print("Reversing ....")
            U = getMatches(gPDB,gSRS,reverse=True,tryAllSizes=True)

        if U.GetNofMatches() >0:

            bestMatch = getBestMatch(U)
        
            print("Best Match")
            print(bestMatch,len(bestMatch[0]),"(",gPDB.GetNofVertices(),",",gSRS.GetNofVertices(),")")
            allMatchNames = []

            for i in range(len(bestMatch[0])):
                ipdb,isrs = bestMatch[0][i], bestMatch[1][i]
                if monomer.atom(isrs-1).element().strip() == res.GetAtom(ipdb-1).element.strip():
                    matchName = monomer.atom(isrs-1).name()
                    if len(monomer.atom(isrs-1).element().strip()) == 1 and len(matchName)<4 and matchName[0] != " ":
                        matchName = " " + matchName
                    if len(matchName)<4:
                        matchName = matchName.ljust(4," ")
                    retMatches[res.GetAtom(ipdb-1).GetAtomName()] = matchName
                    allMatchNames.append(matchName)
                    #res.GetAtom(ipdb-1).SetAtomName(matchName)
            print(retMatches, len(retMatches),res.GetNumberOfAtoms())

            for i in range(res.GetNumberOfAtoms()):
                matched = False
                for k,v in list(retMatches.items()):
                    if res.GetAtom(i).name == k:
                        matched = True
                        break
                if not matched:
                    print(res.GetAtom(i).name, "not matched")
                    if res.GetAtom(i).name in allMatchNames:
                        print("Problem this name is already used",res.GetAtom(i).name)
                        element = res.GetAtom(i).element.strip()
                        if len(element) == 1:
                            for iguess in range(1, 1000):
                                if iguess < 100:
                                    guessName = f" {element}{iguess:<2}"
                                else:
                                    guessName = f"{element}{iguess}"
                                if guessName not in allMatchNames:
                                    for j in range(res.GetNumberOfAtoms()):
                                        if res.GetAtom(j).name == guessName:
                                            continue
                                    print(guessName,"is plausible")
                                    retMatches[res.GetAtom(i).GetAtomName()] = guessName
                                    allMatchNames.append(guessName)
                                    break
                        elif len(element) > 2:
                            print("bad element",res.GetAtom(i).element)

            print(allMatchNames)
            print(retMatches, len(retMatches),res.GetNumberOfAtoms())

            atoms = []
            for k,v in list(retMatches.items()):
                atoms.append(res.GetAtom(k))
            i = 0
            for k,v in list(retMatches.items()):
                atoms[i].SetAtomName(v)
                i += 1

        else:
            print("No reverse matches, cannot rename atoms.")
            return {}


        print("Writing",ofname, "after matching with",dictMatchName,".")
        if ofname is not None:
            molHnd.FinishStructEdit()
            molHnd.WritePDBASCII(ofname);
    
    srsStructFile = SRS.getStructFile()
    
    if testFullDB:
        # Now we need some intelligence of when to not try to match
        rcSRS = mmdb2.intp()
        matches = []
        for i in range(0,SRS.n_entries()):
            monomer = SRS.getMonomer(i,srsStructFile)
            #print"Checking with", monomer.ID()
            if monomer:
                gSRS = monomer.getGraph(rcSRS)
                # FIXME ? -3 or less? More?
                if gSRS and gPDB:
                    if excludeH:
                        gPDB.ExcludeType(mmdb2.getElementNo(("H")))
                        gSRS.ExcludeType(mmdb2.getElementNo(("H")))
                if gSRS and (gPDB.GetNofVertices() - gSRS.GetNofVertices()) > -3:
                    U = getMatches(gPDB,gSRS)
                    if U.GetNofMatches() > 0:
                        #print gSRS.GetNofVertices(), gPDB.GetNofVertices(), U.GetNofMatches(), monomer.chem_name(), monomer.ID()
                        matches.append((gSRS.GetNofVertices(),gSRS,U,monomer))
                        #sys.exit()
    
        matches = sorted(matches)
        matches.reverse()
    
        atomMatches = {}
        imatch = 0
        for m in matches:
            if abs(m[1].GetNofVertices()-gPDB.GetNofVertices()) < 3:
                nVertices,gSRS,U,monomer = m[0], m[1], m[2], m[3]
                bestMatch = getBestMatch(U)
                imatch += 1
                nmatched = 0
                for i in range(len(bestMatch[0])):
                    ipdb,isrs = bestMatch[0][i], bestMatch[1][i]
                    if monomer.atom(isrs-1).element().strip() == res.GetAtom(ipdb-1).element.strip():
                        matchName = str(monomer.atom(isrs-1).name())
                        pdbResName = str(res.GetAtom(ipdb-1).name)
                        if len(monomer.atom(isrs-1).element().strip()) == 1 and len(matchName)<4 and matchName[0] != " ":
                            matchName = " " + matchName
                        if len(matchName)<4:
                            matchName = matchName.ljust(4," ")
                        if pdbResName in atomMatches:
                            atomMatches[pdbResName].append(matchName)
                        else:
                            atomMatches[pdbResName] = [matchName]
    
        # Work out consensus.
        atomMatches_new = {}
        for orig,new in list(atomMatches.items()):
            theSet = set(new)
            maxC = 0
            maxN = ""
            for name in theSet:
                if new.count(name) > maxC:
                    maxC = new.count(name)
                    maxN = name
            atomMatches_new[orig] = maxN
            print("Possible consensus match ",orig, maxN)
    
        print("There were",imatch,"matches")
    
        # Check which match best matches the consensus
        maxnhit = 0
        overallBestMatch = None
        for m in matches:
            if abs(m[1].GetNofVertices()-gPDB.GetNofVertices()) < 3:
                nVertices,gSRS,U,monomer = m[0], m[1], m[2], m[3]
                bestMatch = getBestMatch(U)
                nhit = 0
                for i in range(len(bestMatch[0])):
                    ipdb,isrs = bestMatch[0][i], bestMatch[1][i]
                    if monomer.atom(isrs-1).element().strip() == res.GetAtom(ipdb-1).element.strip():
                        matchName = str(monomer.atom(isrs-1).name())
                        pdbResName = str(res.GetAtom(ipdb-1).name)
                        if len(monomer.atom(isrs-1).element().strip()) == 1 and len(matchName)<4 and matchName[0] != " ":
                            matchName = " " + matchName
                        if len(matchName)<4:
                            matchName = matchName.ljust(4," ")
                        if atomMatches_new[pdbResName] == matchName:
                            nhit += 1
                if nhit > maxnhit:
                    print("new max nhit", monomer.ID(),nhit)
                    maxnhit = nhit
                    overallBestMatch = m
    
        # Finally define the consensus
        if overallBestMatch is not None:
            nVertices,gSRS,U,monomer = overallBestMatch[0], overallBestMatch[1], overallBestMatch[2], overallBestMatch[3]
            bestMatch = getBestMatch(U)
            for i in range(len(bestMatch[0])):
                ipdb,isrs = bestMatch[0][i], bestMatch[1][i]
                if monomer.atom(isrs-1).element().strip() == res.GetAtom(ipdb-1).element.strip():
                    matchName = str(monomer.atom(isrs-1).name())
                    pdbResName = str(res.GetAtom(ipdb-1).name)
                    if len(monomer.atom(isrs-1).element().strip()) == 1 and len(matchName)<4 and matchName[0] != " ":
                        matchName = " " + matchName
                    if len(matchName)<4:
                        matchName = matchName.ljust(4," ")
                    print("Probable best match",pdbResName,matchName)
                    retMatches[res.GetAtom(ipdb-1).GetAtomName()] = matchName
                    res.GetAtom(ipdb-1).SetAtomName(matchName)
    
        if ofname is not None:
            print("Writing",ofname, "after whole ccp4srs search.")
            molHnd.FinishStructEdit()
            molHnd.WritePDBASCII(ofname);

    return retMatches

if __name__ == "__main__":
    ofname = None
    ifname = None
    dictMatchName = None
    dictFileName = None
    selection = None
    parser = argparse.ArgumentParser(description='Specify atom matching input file, optional output file, dictionary to match and atom selection from input file.')
    parser.add_argument('-i', help='input PDB or mmcif file',metavar="pdb/cif filename",required=True)
    parser.add_argument('-o', help='output PDB or mmcif file',metavar="pdb/cif filename")
    parser.add_argument('-s', help='atoms selection as CID, e.g. A/22',metavar="A/1")
    parser.add_argument('-d', help='Optional specific dictionary entry to match. Without this argument or -f the whole of ccp4srs is searched.',metavar="UNK")
    parser.add_argument('-f', help='Optional dictionary file. Without this argument or -d the whole of ccp4srs is searched.',metavar="UNK")
    args = parser.parse_args()
    ifname = args.i
    ofname = args.o
    selection = args.s
    dictMatchName = args.d
    dictFileName = args.f
    retMatches = matchAtoms(ifname,ofname,dictMatchName,selection,dictFileName)
    print(retMatches)
