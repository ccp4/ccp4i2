"""
    ZZPipelineNameZZ.py: CCP4 GUI Project
    
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

import os
from ccp4i2.core.CCP4PluginScript import CPluginScript

class ZZPipelineNameZZ(CPluginScript):
    TASKNAME = 'ZZPipelineNameZZ'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'ZZPluginMaintainerZZ'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                        ['log_mtzjoin.txt', 0]
                       ]
    RUNEXTERNALPROCESS = False

    #Uncomment the following if this pipeline will run plugins asynchronously
    #ASYNCHRONOUS=True

    def __init__(self, *args, **kws):
        super(ZZPipelineNameZZ, self).__init__(*args, **kws)

    #The startProcess method is where you build in the pipeline logic
    def startProcess(self, command, **kws):
        self.ZZFirstPluginNameZZPlugins = []
        self.completedPlugins = []
        for xyzin in self.container.inputData.XYZIN_LIST:
            ZZFirstPluginNameZZPlugin = self.makePluginObject("ZZFirstPluginNameZZ")
            ZZFirstPluginNameZZPlugin.container.inputData.XYZIN = xyzin
            ZZFirstPluginNameZZPlugin.container.inputData.F_SIGF = self.container.inputData.F_SIGF
            
            #Uncomment the following lines if you want to run plugins asynchronously
            '''
            ZZFirstPluginNameZZPlugin.doAsync = True
            import functools
            self.connectSignal(ZZFirstPluginNameZZPlugin,'finished',
                               functools.partial(self.pluginFinished,ZZFirstPluginNameZZPlugin))
            self.ZZFirstPluginNameZZPlugins.append(ZZFirstPluginNameZZPlugin)
            '''
            
            ZZFirstPluginNameZZResult = ZZFirstPluginNameZZPlugin.process()
        return CPluginScript.SUCCEEDED

    #This method will be called as each plugin completes if the pipeline is run asynchronously
    def pluginFinished(self, whichPlugin):
        self.completedPlugins.append(whichPlugin)
        if len(self.ZZFirstPluginNameZZPlugins) == len(self.completedPlugins):
            postProcessStaus = super(ZZPipelineNameZZ, self).postProcess(processId=self._runningProcessId)
            self.reportStatus(postProcessStatus)
            
    def processOutputFiles(self):
        #Create (dummy) PROGRAMXML
        from lxml import etree
        from ccp4i2.core import CCP4Utils
        import sys, os, shutil
        pipelineXMLStructure = etree.Element("ZZPipelineNameZZ")
        
        for iPlugin, ZZFirstPluginNameZZPlugin in enumerate(self.ZZFirstPluginNameZZPlugins):
            out = self.container.outputData
            
            #Copy output XYZs
            out.XYZOUT_LIST.append(out.XYZOUT_LIST.makeItem())
            out.XYZOUT_LIST[-1].setFullPath(os.path.join(self.getWorkDirectory(),"XYZOUT_"+str(iPlugin)+".pdb"))
            shutil.copyfile(ZZFirstPluginNameZZPlugin.container.outputData.XYZOUT.fullPath.__str__(),
                            out.XYZOUT_LIST[-1].fullPath.__str__())
            
            #Copy output 2FoFcs
            out.FPHIOUT_LIST.append(out.FPHIOUT_LIST.makeItem())
            out.FPHIOUT_LIST[-1].setFullPath(os.path.join(self.getWorkDirectory(),"FPHIOUT_"+str(iPlugin)+".mtz"))
            shutil.copyfile(ZZFirstPluginNameZZPlugin.container.outputData.FPHIOUT.fullPath.__str__(),
                            out.FPHIOUT_LIST[-1].fullPath.__str__())
            
            #Copy output FoFcs
            out.DIFFPHIOUT_LIST.append(out.DIFFPHIOUT_LIST.makeItem())
            out.DIFFPHIOUT_LIST[-1].setFullPath(os.path.join(self.getWorkDirectory(),"DIFFPHIOUT_"+str(iPlugin)+".mtz"))
            shutil.copyfile(ZZFirstPluginNameZZPlugin.container.outputData.DIFFPHIOUT.fullPath.__str__(),
                            out.DIFFPHIOUT_LIST[-1].fullPath.__str__())
            
            #Catenate output XMLs
            pluginXMLStructure = CCP4Utils.openFileToEtree(ZZFirstPluginNameZZPlugin.makeFileName("PROGRAMXML"))
            cycleElement = etree.SubElement(pluginXMLStructure,"Cycle")
            cycleElement.text = str(iPlugin)
            pipelineXMLStructure.append(pluginXMLStructure)
        
        with open(self.makeFileName("PROGRAMXML"),"w") as pipelineXMLFile:
            CCP4Utils.writeXML(pipelineXMLFile,etree.tostring(pipelineXMLStructure))
        
        return CPluginScript.SUCCEEDED
