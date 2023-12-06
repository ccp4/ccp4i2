from __future__ import print_function

"""
     CCP4ComTemplate.py: CCP4 GUI Project
     Copyright (C) 2010 University of York

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

"""
   Liz Potterton Nov 2010 - Handle template files that generate program command files
"""
import os
import re
from core.CCP4ErrorHandling import *
from core.CCP4Config import QT,XMLPARSER
if QT():
    from core.CCP4QtObject import CObject
else:
    from core.CCP4Object import CObject

class CComTemplateElement(CObject):

    ERROR_CODES = {1 : {'description' : 'Unknow data specified in command template'},
                   2 : {'description' : 'Error attempting to get text representation of data object'},
                   3 : {'description' : 'Unknown error attempting to find value for word'},
                   4 : {'description' : 'Error attempting to evaluate code in command template'},
                   5 : {'description' : 'Evaluating code in command template did not give Boolean result'},
                   6 : {'description' : 'Variable can not be converted to a string for output to com file'}}

    def __init__(self, parent=None, diagnostic=False):
        CObject.__init__(self, parent)
        self.contents = []
        self.diagnostic = diagnostic
        # Hold data that needs to be visible to child classes - presently
        # just applies to the loop var in CComTemplateLoop
        self.publicData = {}

    def splitLines(self,text):
        textLinesIn = text.split('\n')
        textLines = []
        for lineIn in textLinesIn:
            line = lineIn.strip()
            if len(line) > 0 and line[0] != '#':
                textLines.append(line)
        #print 'splitLines',textLines
        return textLines

    def parseTemplate(self, template):
        templateLines = self.splitLines(template)
        nextLine, err = self.interpretTemplateLines(templateLines)
        return err

    def interpretTemplateLines(self, templateLines=[], nextLine=0):
        errorReport = CErrorReport()
        while self.testToContinue(templateLines, nextLine):
            err = None
            if templateLines[nextLine][0:2] == 'IF':
                self.contents.append(CComTemplateIf(self, diagnostic=self.diagnostic))
                errorReport.extend( self.contents[-1].setParameters(line=templateLines[nextLine]))
                newNextLine, err = self.contents[-1].interpretTemplateLines(templateLines, nextLine+1)
            elif templateLines[nextLine][0:4] == 'LOOP':
                self.contents.append(CComTemplateLoop(self, diagnostic=self.diagnostic))
                errorReport.extend(self.contents[-1].setParameters(line=templateLines[nextLine]))
                newNextLine, err = self.contents[-1].interpretTemplateLines(templateLines, nextLine+1)
            elif templateLines[nextLine][0:5] == 'CASE ':
                self.contents.append(CComTemplateCase(self, diagnostic=self.diagnostic))
                errorReport.extend(self.contents[-1].setParameters(line=templateLines[nextLine]))
                newNextLine, err = self.contents[-1].interpretTemplateLines(templateLines, nextLine+1)
            elif templateLines[nextLine][0:3] == 'AT ':
                self.contents.append(CComTemplateAt(self, diagnostic=self.diagnostic))
                errorReport.extend(self.contents[-1].setParameters(line=templateLines[nextLine]))
                newNextLine = nextLine+1
            elif templateLines[nextLine][0:9] == 'LABELLINE':
                self.contents.append(CComTemplateMtzLabel(self, diagnostic=self.diagnostic))
                newNextLine, err=self.contents[-1].interpretTemplateLines(templateLines, nextLine)
            elif self.testForStatement(templateLines,nextLine):
                self.contents.append(templateLines[nextLine])
                newNextLine = nextLine+1
            else:
                self.contents.append(CComTemplateLine(self, diagnostic=self.diagnostic))
                newNextLine, err = self.contents[-1].interpretTemplateLines(templateLines, nextLine)
            if err is not None:
                errorReport.extend(err)
            if newNextLine > 0:
                nextLine = newNextLine
            else:
                nextLine = nextLine + 1
            if self.diagnostic:
                print('At line:', nextLine, end=' ')
                if nextLine < len(templateLines):
                    print(templateLines[nextLine])
        return nextLine + 1, errorReport


    def testToContinue(self, templateLines, nextLine):
        if nextLine < len(templateLines):
            return True
        else:
            return False

    def testForStatement(self, templateLines, nextLine):
        return False

    def __str__(self):
        text = ''
        for item in self.contents:
            text = text + str(item)
        return text

    def __len__(self):
        return len(self.contents)

    def getPublicData(self, name=''):
        tmplObj = self
        while isinstance(tmplObj, CComTemplateElement):
            #print 'getPublicData', name, tmplObj.__class__, tmplObj.publicData
            if name in tmplObj.publicData:
                return tmplObj.publicData[name]
            else:
                tmplObj = tmplObj.parent()
        return NotImplemented

    def substituteValues(self, line=None, container=None):
        errorReport = CErrorReport()
        wordList = line.split()
        for word in wordList:
            dataValue = None
            if word[0] == '$':
                try:
                    dataValue = self.getValue(word, container=container)
                    #print 'substituteValues', word, dataValue
                except CException as e:
                    errorReport.extend(e)
                except:
                    errorReport.append(self.__class__, 3, word)
                #print errorReport
                if dataValue is None or dataValue == 'None':
                    dataValue = ''
                #print 'substituteValues', word, dataValue, line
                line = re.sub('\\' + word, dataValue, line, 1)
                #print 'substituteValues after sub', line
        return line,errorReport
  
    def getValue(self, word, container=None):
        # Test if a parameter maintained by CComTemplateElement (eg loop var)
        dataValue = self.getPublicData(word[1:])
        #print 'getValue',word,dataValue,type(dataValue)
        if dataValue is not NotImplemented:
            try:
                dataValue = str(dataValue)
            except:
                raise CException(self.__class__, 6, word)
            return dataValue
        dataValue = None
        dataName = word[1:].split('.')
        dataObj = container.find(dataName[0])
        #print 'getValue',word,dataName,dataObj
        if dataObj is None:
            raise CException(self.__class__,1,word)
        if dataName[-1] == 'quote':
            quote = True
            dataName = dataName[0:-1]
        else:
            quote = False
        if len(dataName) > 1:
            try:
                for item in dataName[1:]:
                    dataObj = dataObj.get(item)
            except:
                raise CException(self.__class__, 1, word)
        try:
            dataValue = str(dataObj)
        except:
            raise CException(self.__class__, 2, word)
        if quote:
            return '"' + dataValue + '"'
        else:
            return dataValue

    def evaluateCode(self, code='', container=None):
        localVars,errorReport = self.getLocalVars(code, container=container)
        #print 'evaluateCode',code,localVars
        testCode = re.sub('\$','',code)
        #print 'CComTemplateElement.evaluateCode localVars', localVars, testCode
        result = NotImplemented
        #try:
        result = eval(testCode, {}, localVars)
        #except:
        #  errorReport.append(self.__class__, 4, code)
        #print 'CComTemplateElement.evaluateCode result', result
        return result,errorReport

    def evaluateCodeToBoolean(self, code='', container=None):
        result, errorReport = self.evaluateCode(code=code, container=container)
        boolResult = False
        try:
            boolResult = bool(result)
        except:
            errorReport.append(self.__class__, 5, code)
        return boolResult, errorReport

    def evaluateCodeFragments(self, line, container):
        # Search the input line for code fragments enclosed in curly brace
        # and replace these with an evaluation of the code
        errorReport = CErrorReport()
        m = re.search(r'(.*?)\{(.*?)\}(.*)', line)
        if m is None:
            text = line
        else:
            g = m.groups()
            #print 'CComTemplateLine.evaluateCodeFragments', g
            text = g[0]
            result,err = self.evaluateCode(g[1], container)
            if len(err) > 0:
                text = text + '{' + g[1] + '}'
                errorReport.extend(err)
            else:
                try:
                    resultText = str(result)
                    text = text + resultText
                except:
                    text = text + '{' + g[1] + '}'
                    errorReport.append(self.__class__, 103, g[1])
                if len(g[2]) > 0:
                    moreText, err = self.evaluateCodeFragments(g[2], container)
                    text = text + moreText
                    errorReport.extend(err)
        return text, errorReport

    def getLocalVars(self, line=None, container=None):
        errorReport = CErrorReport()
        localVars = {}
        if container is not None and line is not None:
            r = re.search(r'(.*?)(\$[a-zA-Z\-_]*)(.*)', line)
            while r is not None:
                #print 'getLocalVars', r.groups()
                dataName = r.groups()[1][1:].split('.')[0]
                dataObj = container.find(dataName)
                #print 'getLocalVars', dataName, dataObj
                if dataObj is not None:
                    localVars[dataName] = dataObj
                else:
                    errorReport.append(self.__class__, 1, dataName)
                r = re.search(r'(.*?)(\$[a-zA-Z\-_]*)(.*)', r.groups()[2])
        tmplObj = self
        while isinstance(tmplObj, CComTemplateElement):
            localVars.update(tmplObj.publicData)
            tmplObj = tmplObj.parent()
        return localVars,errorReport


class CComTemplate(CComTemplateElement):
    ERROR_CODES = {}
    ERROR_CODES.update(CComTemplateElement.ERROR_CODES)
    ERROR_CODES.update({101 : {'description' : 'Command template file does not exist'},
                        102 : {'description' : 'Failed to read command template file'},})

    def __init__(self, parent=None, fileName=None, template=None, diagnostic=False):
        CComTemplateElement.__init__(self, parent=None, diagnostic=diagnostic)
        if fileName is not None:
            self.loadTemplateFromFile(fileName)
        elif template is not None:
            self.parseTemplate(template)

    def loadTemplateFromFile(self, fileName):
        from core import CCP4Utils
        template = ''
        if not os.path.exists(fileName):
            raise CException(self.__class__, 101, fileName)
        if os.path.splitext(fileName)[1] == '.xml':
            from core import CCP4File
            xFile = CCP4File.CI2XmlDataFile(fileName)
            header = xFile.loadFile()
            # !!  Should we be checking plugin/version etc.
            body = xFile.getBodyEtree()
            template = str(body.text)
        else: 
            template = CCP4Utils.readFile(fileName)
        #print 'CComTemplate.loadTemplateFromFile', template
        if len(template) == 0:
            raise CException(self.__class__, 102, fileName)
        return self.parseTemplate(template)

    def makeComScript(self, container=None):
        text = ''
        errorReport = CErrorReport()
        for item in self.contents:
            itemText, err = item.makeComScript(container=container)
            text = text + itemText
            errorReport.extend(err)
        return text,errorReport

    def makeComLine(self, container=None):
        wordList = []
        errorReport = CErrorReport()
        for item in self.contents:
            itemText, err = item.makeComScript(container=container)
            wordList.append(itemText.split())
            errorReport.extend(err)
        return wordList, errorReport


class CComTemplateIf(CComTemplateElement):

    ERROR_CODES = {}
    ERROR_CODES.update(CComTemplateElement.ERROR_CODES)

    def __init__(self, parent=None, diagnostic=False, line=''):
        CComTemplateElement.__init__(self, parent, diagnostic)

    def setParameters(self, line=''):
        self.testCode = re.sub('IF', '', line).strip()
        if self.testCode[0] == '{' and self.testCode[-1] == '}':
            self.testCode = self.testCode[1:-2]
        return CErrorReport()

    def testToContinue(self, templateLines, nextLine):
        if templateLines[nextLine] == 'ENDIF':
            return False
        else:
            return True

    def testForStatement(self, templateLines, nextLine):
        if templateLines[nextLine] == 'ELSE':
            return True
        else:
            return False

    def __str__(self):
        text = 'CComTemplateIf: IF\n'
        for item in self.contents:
            if item == 'ELSE':
                text = text + 'CComTemplateIf: ELSE\n'
            else:
                text = text + str(item)
        text = text + 'CComTemplateIf: ENDIF\n'
        return text

    def makeComScript(self, container=None):
        errorReport = CErrorReport()
        text = ''
        writeStatus,err = self.evaluateCodeToBoolean(self.testCode, container)
        errorReport.extend(err)
        if self.contents.count('ELSE') > 0:
            elseIndex = self.contents.index('ELSE')
        else:
            elseIndex = len(self.contents)
        if writeStatus:
            for item in self.contents[0:elseIndex]:
                itemText,err = item.makeComScript(container)
                errorReport.extend(err)
                text = text + itemText
        elif elseIndex < len(self.contents):
            for item in self.contents[elseIndex+1:]:
                itemText,err = item.makeComScript(container)
                errorReport.extend(err)
                text = text + itemText
        return text,errorReport


class CComTemplateLoop(CComTemplateElement):
                      
    ERROR_CODES = {}
    ERROR_CODES.update(CComTemplateElement.ERROR_CODES)
    ERROR_CODES.update({101 : {'description' : 'LOOP statement has wrong number of arguments'},
                        102 : {'description' : "LOOP variable name should not have '$'"},
                        103 : {'description' : 'LOOP start value not interpretable as an integer'},
                        104 : {'description' : 'LOOP end value not interpretable as an integer'},
                        105 : {'description' : 'LOOP increment not interpretable as an integer'},
                        106 : {'description' : 'LOOP increment is zero'}})

    def __init__(self,parent=None, diagnostic=False):
        CComTemplateElement.__init__(self, parent, diagnostic)
        self.loopVarName = None
        self.startVarCode = None
        self.endVarCode = None
        self.loopIncrCode = None
        self.startValue = 1
        self.endValue = 1
        self.loopIncr = 1
        self.loopVarValue = 0

    def setParameters(self, line):
        errorReport = CErrorReport()
        text = re.sub('LOOP','',line).strip()
        wordList = text.split()
        if len(wordList) < 3 or len(wordList) > 4 :
            errorReport.append(self.__class__, 101, line)
            return errorReport
        self.loopVarName = wordList[0]
        # Assume the following are 'code' that need evaluating when we
        # have dat from a container
        self.startVarCode = wordList[1]
        self.endVarCode = wordList[2]
        if len(wordList) > 3:
            self.loopIncrCode = wordList[3]
        return errorReport

    def testToContinue(self, templateLines, nextLine):
        if templateLines[nextLine] == 'ENDLOOP':
            return False
        else:
            return True

    def makeComScript(self, container=None):
        #print 'makeComScript',self.contents[0]
        errorReport = CErrorReport()
        text = ''
        v, err = self.evaluateCode(self.startVarCode, container)
        if len(err) > 0:
            errorReport.append(self.__class__, 103, self.startVarCode)
        else:
            try:
                self.startValue = int(v)
            except:
                errorReport.append(self.__class__, 103, self.startVarCode)
        v,err = self.evaluateCode(self.endVarCode, container)
        if len(err) > 0:
            errorReport.append(self.__class__, 104 ,self.endVarCode)
        else:
            try:
                self.endValue = int(v)
            except:
                errorReport.append(self.__class__, 104, self.endVarCode)
        if self.loopIncrCode is not None:
            v,err = self.evaluateCode(self.loopIncrCode, container)
            if len(err) > 0:
                errorReport.append(self.__class__, 105, self.loopIncrCode)
            else:
                try:
                    self.loopIncr = int(v)
                except:
                    errorReport.append(self.__class__, 105, self.loopIncrCode)
        if self.loopIncr == 0:
            errorReport.append(self.__class__, 106, self.loopIncrCode)
        if len(errorReport) > 0:
            return '', errorReport
        self.loopVarValue = self.startValue - self.loopIncr
        while self.incrementLoopVar():
            for item in self.contents:
                itemText,err = item.makeComScript(container)
                # 'makeComScript',self.loopVarValue,itemText
                errorReport.extend(err)
                text = text + itemText
        return text,errorReport

    def incrementLoopVar(self):
        self.loopVarValue = self.loopVarValue + self.loopIncr
        if (self.loopIncr > 0 and self.loopVarValue > self.endValue) or \
           (self.loopIncr < 0 and self.loopVarValue < self.endValue):
            return False
        else:
            self.publicData[self.loopVarName] = self.loopVarValue
            #print 'incrementLoopVar',self.publicData
            return True

class CComTemplateCase(CComTemplateElement):
    ERROR_CODES = {}
    ERROR_CODES.update(CComTemplateElement.ERROR_CODES)
    ERROR_CODES.update ({101 : {'description' : 'Error attempting to compare CASE values'}})

    def __init__(self, parent=None, diagnostic=False):
        CComTemplateElement.__init__(self, parent, diagnostic)
        self.caseCode = None

    def setParameters(self, line):
        self.caseCode  = line[4:].strip()
        return CErrorReport()

    def testToContinue(self,templateLines,nextLine):
        if templateLines[nextLine] == 'ENDCASE':
            return False
        else:
            return True

    def testForStatement(self, templateLines, nextLine):
        if len(templateLines[nextLine]) > 8 and templateLines[nextLine][0:9] == 'CASEMATCH':
            return True
        else:
            return False

    def makeComScript(self, container=None):
        errorReport = CErrorReport()
        text = ''
        value,err = self.substituteValues(self.caseCode, container)
        errorReport.extend(err)
        if value is NotImplemented:
            return '', errorReport
        import types
        i = -1
        for item in self.contents:
            i = i + 1
            if isinstance(item, str):
                testValue,err = self.substituteValues(item[9:].strip(), container)
                errorReport.extend(err)
                if testValue is not NotImplemented:
                    try:
                        equal = value.__eq__(testValue)
                    except:
                        errorReport.append(self.__class__, 101, item)
                        equal = False
                    #print 'CComTemplateCase',testValue,type(testValue),value,type(testValue),equal
                    if equal:
                        i = i + 1
                        while i < len(self.contents) and not isinstance(self.contents[i], str):
                            itemText,err = self.contents[i].makeComScript(container)
                            errorReport.extend(err)
                            text = text + itemText
                            i = i + 1
                        return text, errorReport
        return text, errorReport

class CComTemplateAt(CComTemplateElement):
    ERROR_CODES = {}
    ERROR_CODES.update(CComTemplateElement.ERROR_CODES)
    ERROR_CODES.update({101 : {'description' : 'Can not interpret file name'},
                        102 : {'description' : 'File not found'}})

    def __init__(self, parent=None, diagnostic=False):
        CComTemplateElement.__init__(self,parent,diagnostic)
        self.fileName = None
    
    def setParameters(self, line):
        from core.CCP4Utils import interpretPath
        self.fileName = interpretPath(line[2:].strip())
        return CErrorReport()
    
    def makeComScript(self, container=None):
        errorReport = CErrorReport()
        text = ''
        if not os.path.exists(self.fileName):
            errorReport.append(self.__class__, 102, self.fileName)
            return text, errorReport
        self.comTemplate = CComTemplate(parent=self, fileName=self.fileName)
        fileText,fileErr = self.comTemplate.makeComScript(container)
        errorReport.extend(fileErr)
        text = text + fileText
        return text, errorReport
  
    
class CComTemplateLine(CComTemplateElement):
    # This is the leaf on the tree - the contents holds
    # template lines rather than other CComTemplateElement objects
    ERROR_CODES = {}
    ERROR_CODES.update(CComTemplateElement.ERROR_CODES)
    ERROR_CODES.update({101 : {'description' : 'Failed to parse line with {..code..}'},
                        102 : {'description' : 'Error evaluating {..code..}'},
                        103 : {'description' : 'Error converting result of {..code..} to a string'}})

    def __init__(self, parent=None, diagnostic=False, continuation=0):
        CComTemplateElement.__init__(self, parent, diagnostic=diagnostic)
        self.continuationDepth = continuation
        p = '-'
        for i in range(self.continuationDepth):
            p = p + '-'
        self.testPattern = p
  
    def testToContinue(self, templateLines, nextLine):
        if re.match(self.testPattern, templateLines[nextLine]) is None:
            return False
        else:
            return True

    def stripContinuation(self, text):
        return text.lstrip('-').strip()

    def interpretTemplateLines(self, templateLines, nextLine):
        errorReport = CErrorReport()
        self.contents.append(self.stripContinuation(templateLines[nextLine]))
        nextLine = nextLine + 1
        while nextLine < len(templateLines):
            if re.match(self.testPattern, templateLines[nextLine]) is not None:
                l = CComTemplateLine(self, continuation=self.continuationDepth+1)
                nextLine, err = l.interpretTemplateLines(templateLines, nextLine)
                #print 'from interpretTemplateLines', self, nextLine
                self.contents.append(l)
                errorReport.extend(err)
            else:
                return nextLine, errorReport
        # We only get here if have reached end of template
        return nextLine + 1, errorReport

    def makeComScript(self,container=None):
        #print 'CComTemplateLine.makeComScript', self.contents[0]
        errorReport = CErrorReport()
        text = ''
        if self.contents[0][0] == '{':
            terminator = self.contents[0].index('}')
            if terminator < 0:
                errorReport.append(self.__class__, 101, self.contents[0])
                return '', errorReport
            testCode = self.contents[0][1:terminator]
        else:
            if self.contents[0].count(' ') > 0:
                terminator = self.contents[0].index(' ')
                testCode = self.contents[0][0:terminator]
            else:
                testCode = self.contents[0]
        #print 'makeComScript',testCode
        writeStatus,err = self.evaluateCodeToBoolean(testCode, container)
        errorReport.extend(err)
        templText = self.contents[0][terminator + 1:]
        if templText.count('|'):
            idx = templText.index('|')
            altText = templText[idx + 1:]
            templText = templText[0:idx - 1]
        else:
            altText = None
        #print 'makeComScript', templText, altText
        if writeStatus or altText is not None:
            if writeStatus:
                subText,err = self.substituteValues(templText.strip(), container)
            else:
                subText,err = self.substituteValues(altText.strip(), container)
            #print 'CComTemplateLine.makeComScript subText',subText
            errorReport.extend(err)
            if len(subText) > 0:
                text,err = self.evaluateCodeFragments(subText, container)
                errorReport.extend(err)
                if len(text) > 0:
                    text = self.continuationMarker() + text
                    if text[-1] != '\n':
                        text = text + '\n'
        if writeStatus and len(self.contents) > 1:
            for line in self.contents[1:]:
                lineText, err = line.makeComScript(container=container)
                errorReport.extend(err)
                if len(lineText) > 0:
                    text = text + lineText
        #print 'CComTemplateLine.makeComScript',text
        return text, errorReport

    def continuationMarker(self):
        if self.continuationDepth > 0:
            return '  - '
        else:
            return ''

    def __str__(self):
        text = 'CComTemplateLine ' + str(self.continuationDepth) + ': ' + self.contents[0] + '\n'
        if len(self.contents) > 1:
            for line in self.contents[1:]:
                text = text + str(line)
        return text

class CComTemplateMtzLabel(CComTemplateElement):

    ERROR_CODES = {101 : {'description' : 'Pair of words in LABELLINE are not label and value'},
                   102 : {'severity' : SEVERITY_WARNING, 'description' : 'No value for LABELLINE parameter'}}
    ERROR_CODES.update(CComTemplateElement.ERROR_CODES)

    def __init__(self, parent=None, diagnostic=False, line=''):
        CComTemplateElement.__init__(self, parent, diagnostic)
        self.testPattern = '-'

    def testToContinue(self,templateLines,nextLine):
        if nextLine >= len(templateLines):
            return False
        if re.match(self.testPattern, templateLines[nextLine]) is None:
            return False
        else:
            return True

    def __str__(self):
        text = 'CComTemplateMtzLabel: LABELLINE\n'
        for item in self.contents:
            text = text + str(item)
        return text
  
    def interpretTemplateLines(self, templateLines=[], nextLine=0):
        errorReport = self.appendContent(re.sub('LABELLINE', '', templateLines[nextLine]))
        while self.testToContinue(templateLines, nextLine + 1):
            nextLine = nextLine + 1
            #print 'CComTemplateMtzLabel.interpretTemplateLines', nextLine, templateLines[nextLine]
            errorReport.append(self.appendContent(re.sub('-', '', templateLines[nextLine], 1)))
        return nextLine + 1, errorReport

    def appendContent(self, line):
        items = line.split()
        errorReport = CErrorReport()
        for ii in range(len(items)/2):
            if items[ii*2][0] != '$' and items[ii*2 + 1][0] == '$':
                self.contents.append([items[ii*2], items[ii*2 + 1]])
            else:
                errorReport.append(self.__class__, 101, items[ii*2] + ' ' + items[ii*2 + 1])
        return errorReport

    def makeComScript(self, container=None):
        errorReport = CErrorReport()
        text = ''
        #print 'CComTemplateMtzLabel.makeComScript contents',self.contents
        for name, param in self.contents:
            try:
                value = self.getValue(param, container=container)
            except CException as e:
                errorReport.extend(e)
            else:
                if value is None or value == 'None':
                    errorReport.append(self.__class__, 102, name)
                    if len(text) > 0:
                        text = 'LABIN ' + text.strip() + '\n'
                    return text, errorReport
                else:
                    text = text + name + '=' + value + ' '
        if len(text) > 0:
            text = 'LABIN ' + text.strip() + '\n'
        return text,errorReport

#===========================================================================================
import unittest
def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testComTemplate)
    #suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testQObject))
    return suite

def testModule():
    suite = TESTSUITE()
    #print 'suite',suite
    unittest.TextTestRunner(verbosity=2).run(suite)

class testComTemplate(unittest.TestCase):

    def makeContainer(self):    
        from core import CCP4Container
        from core.CCP4Utils import getCCP4I2Dir
        c = CCP4Container.CContainer()
        c.loadContentsFromXml(os.path.join(getCCP4I2Dir(),'test','data','test_com_template_1.def.xml'))
        c.loadDataFromXml(os.path.join(getCCP4I2Dir(),'test','data','test_com_template_1.params.xml'))
        #print 'makeContainer',c.dataOrder()
        return c

    def test_1(self):
        template = '''
1 whatever $WHATEVER
$TEST_PAR test $TEST_VALUE
'''
        t = CComTemplate(template=template)
        #print t
        self.assertEqual(len(t),2,'Failed to load template')

    def test_2(self):
        template = '''
1 title $TITLE
IF $HOWEVER
  1 whatever $WHATEVER
ELSE
  $TEST_PAR test $TEST_VALUE
ENDIF
1 trailer $TRAILER
'''
        t = CComTemplate(template=template)
        #print t
        self.assertEqual(len(t),3,'Failed to load template')
        self.assertTrue(isinstance(t.contents[1],CComTemplateIf),'Failed to load template CComTemplateIf')

    def test_3(self):
        c = self.makeContainer()
        t = CComTemplate(template='')
        v = t.getValue('$WHATEVER',container=c)
        self.assertEqual(v,'12','Failed to convert $WHATEVER to data value')
        v = t.getValue('$TITLE.text',container=c)
        self.assertEqual(v,'Try this title','Failed to convert $TITLE.text to data value')
        try: 
            v = t.getValue('$DONTEXIST',container=c)
        except CException as e:
            errcode = e[0]['code']
        self.assertEqual(errcode,1,'Failed to convert $DONTEXIST to correct error code')


    def test_4(self):
        c = self.makeContainer()
        t = CComTemplate(template='')
        v,err = t.substituteValues('1 whatever $WHATEVER',container=c)
        self.assertEqual(v,'1 whatever 12','Failed to substituteValues in WHATEVER line')    
        v,err = t.substituteValues('$TEST_PAR test $TEST_VALUE',container=c)
        self.assertEqual(v,'True test 24.0','Failed to substituteValues in TEST_VALUE line')

    def test_5(self):
        c = self.makeContainer()
        t = CComTemplate(template='')
        v,err = t.evaluateCodeToBoolean('$WHATEVER.isSet()',container=c)
        self.assertTrue(v,'Failed to evaluate $WHATEVER.isSet()')
        v,err = t.evaluateCodeToBoolean('10 < $WHATEVER < 15',container=c)
        self.assertTrue(v,'Failed to evaluate 10 < $WHATEVER < 15')

    def test_6(self):
        c = self.makeContainer()
        template = '''1 whatever $WHATEVER\n$TEST_PAR test $TEST_VALUE\n{$WHATEVER >20} greater $WHATEVER\n{$WHATEVER <20} less than $WHATEVER\n'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_6',text,err
        self.assertEqual(len(err),0,'Unexpected errors in report')
        self.assertEqual(text,'''whatever 12\ntest 24.0\nless than 12\n''','Failed to create correct script')

    def test_7(self):
        c = self.makeContainer()
        template = '''1 whatever $WHATEVER\n- $TEST_PAR test $TEST_VALUE\n-- {$WHATEVER >20} greater $WHATEVER\n-- {$WHATEVER <20} less than $WHATEVER\n'''
        t = CComTemplate(template=template)
        #print 'test_7'
        #print t
        text,err = t.makeComScript(container=c)
        #print 'test_7'
        #print text
        #print err
        self.assertEqual(len(err),0,'Unexpected errors in report')
        self.assertEqual(text,'''whatever 12\n  - test 24.0\n  - less than 12\n''','Failed to create correct script')

    def test_8(self):
        c = self.makeContainer()
        template = '''
1 whatever $WHATEVER
IF not $HOWEVER
  $TEST_PAR test $TEST_VALUE
  {$WHATEVER >20} greater $WHATEVER
ELSE
  {$WHATEVER <20} less than $WHATEVER\n
ENDIF
{len( $TRAILER )>2} $TRAILER
'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_8'
        #print text
        #print err
        self.assertEqual(len(err),0,'Unexpected errors in report')
        self.assertEqual(text,'''whatever 12\nless than 12\nfoo\n''','Failed to create correct script')

    def test_9(self):
        c = self.makeContainer()
        template = '''
1 whatever $WHATEVER
IF $HOWEVER
  IF $TEST_PAR
    1 test $TEST_VALUE
  ELSE
    1 greater $WHATEVER
  ENDIF
ELSE
  {$WHATEVER <20} less than $WHATEVER\n
ENDIF
{len( $TRAILER )>2} $TRAILER
'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_9'
        #print text
        #print err
        self.assertEqual(len(err),0,'Unexpected errors in report')
        self.assertEqual(text,'''whatever 12\ntest 24.0\nfoo\n''','Failed to create correct script')


    def test_10(self):
        c = self.makeContainer()
        template = '''
1 HOWEVER $HOWEVER
LOOP N 1 $WHATEVER
  1 RUN $N
ENDLOOP
{len( $TRAILER )>2} $TRAILER
'''
        t = CComTemplate(template=template)
        text, err = t.makeComScript(container=c)
        #print 'test_10'
        #print text
        #print err
        self.assertEqual(len(err),0,'Unexpected errors in report')

    def test_11(self):
        l = CComTemplateLine()
        t,err = l.evaluateCodeFragments('test { 3 * 4 } this',None)
        self.assertEqual(t,'test 12 this','Error in CComTemplateLine.evaluateCodeFragments')
        self.assertEqual(len(err),0,'Unexpected errors in test 11')
 
    def test_12(self):
        c = self.makeContainer()
        template = '''
1 HOWEVER $HOWEVER
LOOP N 1 $WHATEVER
  1 RUN $N
  - 1 TEST { $N * $TEST_VALUE }
ENDLOOP
{len( $TRAILER )>2} $TRAILER
'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_12'
        #print text
        #print err
        self.assertEqual(len(err),0,'Unexpected errors in report')
 
    def test_13(self):
        c = self.makeContainer()
        template = '''
1 HOWEVER $HOWEVER
CASE $TRAILER
CASEMATCH bar
  1 matched bar
CASEMATCH foo
  1 matched foo
ENDCASE
1 WHATEVER $WHATEVER
'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_13'
        #print text
        #print err
        self.assertEqual(text,'HOWEVER True\nmatched foo\nWHATEVER 12\n')
        self.assertEqual(len(err),0,'Unexpected errors in report')

    def test_14(self):
        c = self.makeContainer()
        template = '''
1 HOWEVER $HOWEVER
CASE $WHATEVER
CASEMATCH 12
  1 matched 12
CASEMATCH  24
  1 matched 24
ENDCASE
1 WHATEVER $WHATEVER
'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_14'
        #print text
        #print err
        self.assertEqual(text,'HOWEVER True\nmatched 12\nWHATEVER 12\n')
        self.assertEqual(len(err),0,'Unexpected errors in report')

    def test_15(self):
        c = self.makeContainer()
        template = '''
1 HOWEVER $HOWEVER
CASE $TRAILER
CASEMATCH bar
  1 matched bar
CASEMATCH foo
  1 matched foo
  AT $CCP4I2_TOP/test/data/test_com_template_1.com
ENDCASE
'''
        t = CComTemplate(template=template)
        text,err = t.makeComScript(container=c)
        #print 'test_15'
        #print text
        #print err
        self.assertEqual(text,'HOWEVER True\nmatched foo\nwhatever 12\nless than 12\n')
        self.assertEqual(len(err),0,'Unexpected errors in report')

    def test_16(self):
        c = self.makeContainer()
        template = '''
LABELLINE F $F_SIGF.F SIGF $F_SIGF.SIGF
LABELLINE PHI $PHI_W.PHI W $PHI_W.W
'''
        t = CComTemplate(template=template,diagnostic=False)
        text,err = t.makeComScript(container=c)
        #print 'test_16'
        #print text,text.split('\n')
        #print err
        self.assertEqual(len(err),1,'Wrong number of errors in LABELLINE test')
        self.assertEqual(text,'LABIN F=ffoo SIGF=ffoosig\nLABIN PHI=foophi\n','LABELLINE output incorrect')
