from __future__ import print_function

import sys
from lxml import etree

class MDLAtom(object):
    def __init__(self, molLine=None):
        super(MDLAtom,self).__init__()
        self.x = -1e30
        self.y = -1e30
        self.z = -1e30
        self.element = 'C'
        self.massDifference = 0
        self.charge = 0
        self.nbonds = 0
        if molLine is not None: self.fromMolLine(molLine)
    def fromMolLine(self, molLine):
        self.x = float(molLine[ 0:10])
        self.y = float(molLine[10:20])
        self.z = float(molLine[20:30])
        self.element = molLine[30:34].strip()
        self.massDifference =int(molLine[34:36])
        self.charge =int(molLine[36:39])
    def __str__(self):
        result = ''
        result += '%10.5f%10.5f%10.5f %3s%2d%3d\n'%(self.x, self.y, self.z, self.element, self.massDifference, self.charge)
        return result

class MDLBond(object):
    def __init__(self, molLine=None):
        super(MDLBond,self).__init__()
        self.at1 = -1
        self.at2 = -1
        self.order = -1
        self.stereo = -1
        if molLine is not None: self.fromMolLine(molLine)
    def fromMolLine(self, molLine):
        self.at1 = int(molLine[ 0:3])
        self.at2 = int(molLine[3:6])
        self.order = int(molLine[6:9])
        self.stereo = int(molLine[9:12])
    def __str__(self):
        result = ''
        result += '%3d%3d%3d%3d\n'%(self.at1, self.at2, self.order, self.stereo)
        return result

class MDLMolecule(object):
    def __init__(self, filePath = None, molBlock=None):
        super(MDLMolecule,self).__init__()
        self.nAtoms = 0
        self.nBonds = 0
        self.atoms = []
        self.bonds  = []
        if filePath is not None: self.fromMolFile(filePath)
        elif molBlock is not None: self.fromMolBlock(molBlock)
    
    def fromMolFile(self, filePath):
        with open(filePath,'r') as file:
            molBlock = file.read()
        return self.fromMolBlock(molBlock)

    def fromMolBlock(self, molBlock):
        lines = molBlock.split('\n')
        countsRead = False
        nAtomsRead = 0
        for line in lines:
            if line.startswith('#'): pass
            elif not countsRead and '999' in line:
                self.nAtoms = int(line[0:3])
                self.nBonds = int(line[4:6])
                countsRead = True
            elif len(self.atoms) < self.nAtoms:
                atom = MDLAtom(line)
                self.atoms.append(atom)
            elif len(self.bonds) < self.nBonds:
                bond=MDLBond(line)
                self.bonds.append(bond)
                self.atoms[bond.at1-1].nbonds += bond.order
                self.atoms[bond.at2-1].nbonds += bond.order
        print(self.atoms, self.bonds)
        return

    def __str__(self):
        a  = 'Natoms: '+str( self.nAtoms) + '\n'
        a += 'Nbonds: '+str( self.nBonds) + '\n'
        for atom in self.atoms:
            a+=str(atom)
        for bond in self.bonds:
            a+=str(bond)
        return a
            
    def normalize(self, size = None):
        xs = [atom.x for atom in self.atoms]
        xmin = min(xs)
        xmax = max(xs)
        ys = [-atom.y for atom in self.atoms]
        ymin = min(ys)
        ymax = max(ys)
        self.factor = 12.
        if size is not None:
            xFactor = size[0]/((xmax-xmin)+24)
            yFactor = size[1]/((ymax-ymin)+24)
            self.factor = min(xFactor,yFactor)
        zs = [atom.z for atom in self.atoms]
        zmin = min(zs)
        zmax = max(zs)
        for atom in self.atoms:
            atom.x = 15.+self.factor*(atom.x-xmin)
            atom.y = 15.+self.factor*(-atom.y-ymin)
            atom.z = 15.+self.factor*(atom.z-zmin)

    def svgXML(self, size = None):
        import math
        if size is None: size = (240,180)
        self.normalize(size)
        svgNode = etree.fromstring('''
<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="WIDTHpx" height="HEIGHTpx">
</svg>'''.replace('WIDTH',str(size[0])).replace('HEIGHT',str(size[1])))
        styleNode = etree.SubElement(svgNode,'style',type='text/css')
        fontSize = 8.*(self.factor/8.)
        bondStrokeWidth = 1.*self.factor/8.
        textOutlineStrokeWidth = fontSize/2.
        styleNode.text = etree.CDATA('''
            polygon.forwards { fill:black; stroke:black; stroke-width:0; stroke-linecap:round;stroke-linejoin:round;}
            polygon.backwards { fill:black; stroke:black; stroke-width:0; stroke-linecap:round;stroke-linejoin:round;}
            line {  stroke:black;stroke-width:BONDSTROKEWIDTH;stroke-linecap:round;}
            text.outline {font-family:Arial, Helvetica, sans-serif; font-size:FONTSIZEpx;stroke-width:TEXTOUTLINESTROKEWIDTH; stroke:white; fill:white; }
            text.foreground {font-family:Arial, Helvetica, sans-serif; font-size:FONTSIZEpx;stroke-width:0; stroke:black; fill:black;}
'''.replace('BONDSTROKEWIDTH',str(bondStrokeWidth)).replace('FONTSIZE',str(fontSize)).replace('TEXTOUTLINESTROKEWIDTH',str(textOutlineStrokeWidth)) )
        svgNode.append(styleNode)
        for bond in self.bonds:
            at1 = self.atoms[bond.at1-1]
            at2 = self.atoms[bond.at2-1]
            norm = [at2.y-at1.y, at1.x-at2.x]
            normLength = math.sqrt(norm[0]*norm[0]+norm[1]*norm[1])
            #Try to preempt zeroLength bonds with kludge
            try:
                norm = [(self.factor/8.)*comp/normLength for comp in norm]
            except:
                norm = [(self.factor/8.)*(at2.y-at1.y), (self.factor/8.)*(at1.x-at2.x)]
            if bond.order == 1 and bond.stereo == 1:
                points = [[at1.x,at1.y],[at2.x+1.5*norm[0],at2.y+1.5*norm[1]], [at2.x-1.5*norm[0],at2.y-1.5*norm[1]]]
                pointsString = " ".join([str(point[0])+','+str(point[1]) for point in points])
                polygonNode = etree.SubElement(svgNode,'polygon',points=pointsString)
                polygonNode.set('class','forwards')
            elif bond.order == 1 and bond.stereo == 6:
                offsetMax = 1.5
                nSegs = 6
                vec12 = [at2.x-at1.x, at2.y-at1.y]
                for iSeg in range(nSegs):
                    factorStart = float((2*iSeg)  ) / float((2*nSegs)-1)
                    factorEnd   = float((2*iSeg)+1) / float((2*nSegs)-1)
                    at1Crd = [at1.x, at1.y]
                    centreStart = [at1Crd[i] + (factorStart  * vec12[i]) for i in range(2)]
                    centreEnd =   [at1Crd[i] + (factorEnd    * vec12[i]) for i in range(2)]
                    offsetStart = factorStart * offsetMax
                    offsetEnd   = factorEnd   * offsetMax
                    points = [[centreStart[0] + offsetStart*norm[0],  centreStart[1] + offsetStart*norm[1]],
                              [centreEnd[0]   + offsetEnd  *norm[0],  centreEnd[1]   + offsetEnd  *norm[1]],
                              [centreEnd[0]   - offsetEnd  *norm[0],  centreEnd[1]   - offsetEnd  *norm[1]],
                              [centreStart[0] - offsetStart*norm[0],  centreStart[1] - offsetStart*norm[1]],
                              ]
                    pointsString = " ".join([str(point[0])+','+str(point[1]) for point in points])
                    polygonNode = etree.SubElement(svgNode,'polygon',points=pointsString)
                    polygonNode.set('class','backwards')
            elif bond.order == 2:
                lineNode = etree.SubElement(svgNode,'line',x1=str(at1.x-norm[0]), y1=str(at1.y-norm[1]), x2=str(at2.x-norm[0]), y2=str(at2.y-norm[1]))
                lineNode = etree.SubElement(svgNode,'line',x1=str(at1.x+norm[0]), y1=str(at1.y+norm[1]), x2=str(at2.x+norm[0]), y2=str(at2.y+norm[1]))
            
            else:
                lineNode = etree.SubElement(svgNode,'line',x1=str(at1.x), y1=str(at1.y), x2=str(at2.x), y2=str(at2.y))
        for atom in self.atoms:
            if atom.element != 'C' or atom.charge != 0:
                element = atom.element
                valence = 0
                if element == 'C':
                    valence = 4
                elif element.upper() == 'SI':
                    valence = 4
                elif element == 'N':
                    valence = 3
                elif element == 'B':
                    valence = 3
                elif element == 'O':
                    valence = 2
                elif element == 'S':
                    valence = 2
                if valence>0:
                    Hes = (valence+atom.charge)-atom.nbonds
                    if   Hes == 4: element += 'H4'
                    elif Hes == 3: element += 'H3'
                    elif Hes == 2: element += 'H2'
                    elif Hes == 1: element += 'H'
                if atom.charge == 1: element += '+'
                if atom.charge == -1: element += '-'
                elementOutlineNode = etree.SubElement(svgNode, 'text', x=str(atom.x), y=str(atom.y), dy='0.5em', dx='-0.5em')
                elementOutlineNode.set('class','outline')
                elementOutlineNode.text = element
                elementForegroundNode = etree.SubElement(svgNode, 'text', x=str(atom.x), y=str(atom.y), dy='0.5em', dx='-0.5em')
                elementForegroundNode.set('class','foreground')
                elementForegroundNode.text = element

        return svgNode

if __name__ == '__main__':
    from ccp4i2.core import CCP4Utils
    newLigand = MDLMolecule('job_1/prodrg-in.mdl')
    newLigand.normalize()
    aNode = newLigand.svgXML()
    a = etree.tostring(aNode, pretty_print=True)
    with open('new.svg','w') as outputSVG:
        CCP4Utils.writeXML(outputSVG,a)

