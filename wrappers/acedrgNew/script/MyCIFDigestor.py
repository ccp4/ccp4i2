from __future__ import print_function

import sys, os

class MyCIFException(BaseException):
    pass

class MyCIFBlock(object):
    def __init__(self, *args,**kws):
        pass

class MyCIFLoopLine(object):
    pass

from collections.abc import Sequence
class MyCIFLoop(Sequence):
    TAG = '_loop'
    def __init__(self, text):
        super(MyCIFLoop,self).__init__()
        self.keys = []
        self.data = []
        self.category = None
        haveDataCompList = False
        for line in text.split('\n'):
            if (not line.startswith('#')) and len(line.strip()) > 0 and (not line.startswith('loop_')):
                if line.lstrip().startswith("data_comp_list"): haveDataCompList = True
                if haveDataCompList:
                    if self.category is None: self.category = line.split('.')[0]
                    if line.startswith(self.category):
                        self.keys.append(line.split('.')[1])
                    elif len(self.keys) > 0:
                        tokens = self.parseLine(line)
                        self.data.append(tokens)

    def parseLine(self, line):
        tokens = []
        inToken = False
        delimiters = [' ','\t']
        token = ''
        for i in range(len(line)):
            if inToken:
                if line[i] in delimiters:
                    tokens.append(token)
                    inToken = False
                    delimiters = [' ','\t']
                else:
                    token += line[i]
            else:
                if line[i] == ' ' or line[i] == '\t':
                    pass
                elif line[i] == "'" or line[i] == '"':
                    inToken = True
                    delimiters = [line[i]]
                    token = ''
                else:
                    inToken = True
                    token = line[i:i+1]
        if inToken:
            tokens.append(token)
        return tokens

    def __getitem__(self, key):
        dataline = self.data[key]
        return self.dataAsDict(dataline)

    def __len__(self):
        return len(self.data)
    
    def setLinePropertyValue(self,linedict,property,value):
        theLineIndex = self.index(linedict)
        theKeyIndex = self.keys.index(property)
        self.data[theLineIndex][theKeyIndex] = str(value)

    def dataAsDict(self, dataline):
        dataAsDict = {}
        for i in range(len(self.keys)):
            #dataAsDict['.'.join([self.category, self.keys[i]])] = dataline[i]
            dataAsDict[self.keys[i]] = dataline[i]
        return dataAsDict

    def __str__(self):
        selfAsString = 'loop_\n'
        for key in self.keys:
            selfAsString += (self.category+'.'+key+'\n')
        for dataLine in self.data:
            quotedText = []
            for value in dataLine:
                if ' ' in value:
                    quotedText.append("'"+value+"'")
                else:
                    quotedText.append(value)
            selfAsString += ('\t'.join(quotedText) + '\n')
        return selfAsString

class MyCIFDataBlock(MyCIFBlock):
    category='Banana'
    def __init__(self, *args,**kws):
        super(MyCIFDataBlock,self).__init__(*args, **kws)
        self.elements = []
        if kws.get('text',None) is not None:
            text=kws.get('text')
            inLoop = False
            loopBlock = ''
            for line in text.split('\n'):
                if line.strip().lower().startswith('data_'):
                    self.category = line.split()[0]
                else:
                    #Identify start of loop
                    if not inLoop:
                        if line.startswith('loop_'):
                            inLoop = True
                            loopBlock = ''
                    else:
                    #Identify end of loop
                        if line.startswith('loop_'):
                            loop = MyCIFLoop(loopBlock)
                            self.elements.append(loop)
                            #print len(loop)
                            loopBlock = ''
                        elif len(line.strip()) == 0:
                            loop = MyCIFLoop(loopBlock)
                            self.elements.append(loop)
                            inLoop = False
                            loopBlock = ''
                        loopBlock += (line + '\n')
            if inLoop and len(loopBlock) > 0:
                loop = MyCIFLoop(loopBlock)
                self.elements.append(loop)
                loopBlock = ''
        if kws.get('category', None) is not None:
            self.category = kws.get('category')
                
    def loops(self):
        return [element for element in self.elements if type(element) == MyCIFLoop]
    def addElement(self, **kws):
        newElement = MyCIFElement(**kws)
        self.elements.append(newElement)
        return newElement
    def __str__(self):
        selfAsString = self.category + '\n'
        for element in self.elements:
            selfAsString += str(element)
        return selfAsString

class MyCIFGlobalBlock(MyCIFBlock):
    category='global_'
    def __init__(self, text):
        super(MyCIFGlobalBlock,self).__init__()
        self._lib_name = '?'
        self._lib_version = '?'
        self._lib_update = '?'
        for line in text.split('\n'):
            tokens = line.split()
            if len(tokens) > 1 and tokens[0].lower() == '_lib_name': self._lib_name = tokens[1]
            if len(tokens) > 1 and tokens[0].lower() == '_lib_version': self._lib_version = tokens[1]
            if len(tokens) > 1 and tokens[0].lower() == '_lib_update': self._lib_update = tokens[1]
    def __str__(self):
        selfAsString = ''
        selfAsString += 'global_\n'
        selfAsString += '_lib_name %s\n'%(self._lib_name)
        selfAsString += '_lib_version %s\n'%(self._lib_version)
        selfAsString += '_lib_update %s\n'%(self._lib_update)
        return selfAsString

class MyCIFElement(MyCIFBlock):
    types = ['MyCIFString','MyCIFTextBlock']
    def __init__(self, *args, **kws):
        super(MyCIFElement,self).__init__(*args, **kws)
        self.myCIFType = kws.get('myCIFType', 'MYCIFString')
        #Provide default formatters
        if self.myCIFType == 'MyCIFString':
            self.formatter = kws.get('formatter',"'{:s}'")
        elif self.myCIFType == 'MyCIFTextBlock':
            self.formatter = kws.get('formatter',"\n;{:s}\n;")
        self.name = kws.get('name', '_undefined')
        self.value = kws.get('value', None)
    def __str__(self):
        selfAsString = self.name+' '
        if self.myCIFType == 'MyCIFString':
            selfAsString += self.formatter.format(str(self.value))
        elif self.myCIFType == 'MyCIFTextBlock':
            selfAsString += self.formatter.format(str(self.value))
        return selfAsString+'\n'

class MyCIFFile(object):
    def __init__(self,*args,**kw):
        self.blocks = []
        if kw.get('filePath',None) is not None:
            self.initWithFilePath(kw.get('filePath'))
    def initWithFilePath(self,filePath=None):
        #print 'Opening file ', filePath
        if filePath is None or not os.path.isfile(filePath):
            raise BaseException('Oops')
        with open(filePath,'r') as file:
            lines = file.readlines()
            inGlobal = False
            inData = False
            inChod = True
            blockText=''
            for line in [longLine.strip() for longLine in lines]:
                if line.lower().startswith('#'):
                    pass
                    #print 'Comment', line
                elif inChod:
                    if line.lower().startswith('global_'):
                        inChod = False
                        inData = False
                        inGlobal = True
                    elif line.lower().startswith('data_'):
                        inGlobal = False
                        inChod = False
                        inData = True
                elif inGlobal:
                    if line.lower().startswith('data_'):
                        globalBlock = MyCIFGlobalBlock(blockText)
                        self.blocks.append(globalBlock)
                        blockText = ''
                        inGlobal = False
                        inChod = False
                        inData = True
                elif inData:
                    if line.lower().startswith('data_'):
                        dataBlock = MyCIFDataBlock(text=blockText)
                        self.blocks.append(dataBlock)
                        blockText = ''
                        inGlobal = False
                        inChod = False
                        inData = True
                blockText += (line+'\n')
            if inGlobal:
                globalBlock = MyCIFGlobalBlock(blockText)
                self.blocks.append(globalBlock)
            if inData:
                dataBlock = MyCIFDataBlock(text=blockText)
                self.blocks.append(dataBlock)
    def dataBlocks(self):
        return [aBlock for aBlock in self.blocks if type(aBlock) is MyCIFDataBlock]
    
    def addDataBlock(self, **kws):
        newBlock = MyCIFDataBlock(**kws)
        self.blocks.append(newBlock)
        return newBlock

    def addElement(self, **kws):
        newElement = MyCIFElement(**kws)
        self.blocks.append(newElement)
        return newElement
    
    def __str__(self):
        return '#\n'.join([str(block) for block in self.blocks])


if __name__ == '__main__':
    print(sys.argv)
    a = MyCIFFile(filePath=sys.argv[1])
    if len(sys.argv) > 2:
        with open(sys.argv[2],'w') as outputFile:
            outputFile.write(str(a))
    for block in a.dataBlocks():
        for loop in [aLoop for aLoop in block.loops if aLoop.category == '_chem_comp_atom']:
            for loopline in loop:
                if 'x' in loopline:
                    loop.setLinePropertyValue(loopline,'x',0.0)
    print(str(a))


