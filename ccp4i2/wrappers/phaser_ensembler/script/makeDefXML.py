import os
import re

from lxml import etree
from phaser.ensembler import PHIL_MASTER


# define mappings for .type to <className>
phil_type_as_class_name = {
    "str" : "CString",
        "int" : "CInt",
            "bool" : "CBoolean",
                "float" : "CFloat",
                    "choice" : "CString",
}

def parse_choice_options (phil_def) :
    return ",".join([ re.sub(r"\*", "", word.value) for word in phil_def.words ])

def parse_captions (phil_def) :
    return ",".join([re.sub("_"," ",item) for item in phil_def.caption.split()])

def definitionToElement(keyword):
    value = keyword.extract()
    elem = etree.SubElement(root, "content")
    elem.set("id", keyword.full_path().replace('.','__'))
    elem_class = etree.SubElement(elem, "className")
    elem_class.text = phil_type_as_class_name.get(keyword.type.phil_type,
      "CString")
    qualifiers = etree.SubElement(elem, "qualifiers")
    if (keyword.short_caption is not None) : # control label
      gui_label = etree.SubElement(qualifiers, "guiLabel")
      gui_label.text = keyword.short_caption
    if (keyword.help is not None) : # tooltip
      toolTip = etree.SubElement(qualifiers, "toolTip")
      toolTip.text = keyword.help
    if (keyword.expert_level is not None) :
      expert_level = etree.SubElement(qualifiers, "expertLevel")
      expert_level.text = str(keyword.expert_level)
    if (keyword.type.phil_type == "bool") and (value is not None) :
      default = etree.SubElement(qualifiers, "default")
      default.text = str(value)
    elif (keyword.type.phil_type == "choice") :
      enum = etree.SubElement(qualifiers, "enumerators")
      enum.text = parse_choice_options(keyword)
      default = etree.SubElement(qualifiers, "default")
      default.text = str(value)
      if (keyword.caption is not None) :
        menu = etree.SubElement(qualifiers, "menuText")
        menu.text = parse_captions(keyword)
      if (not keyword.type.multi) :
        only = etree.SubElement(qualifiers, "onlyEnumerators")
        only.text = "1"
    elif (keyword.type.phil_type in ["int", "float"]) :
      default = etree.SubElement(qualifiers, "default")
      default.text = str(value)
      if (keyword.type.value_min is not None) :
        min = etree.SubElement(qualifiers, "min")
        min.text = str(keyword.type.value_min)
      if (keyword.type.value_max is not None) :
        max = etree.SubElement(qualifiers, "max")
        max.text = str(keyword.type.value_max)
    return elem

def flattenScope(scope, root):
    for object in scope.objects:
        if object.is_definition:
            keywordElement = definitionToElement(object)
            root.append(keywordElement)
        elif object.is_scope:
            print('Object is scope', object.full_path())
            flattenScope(object, root)


if __name__ == "__main__":
    CCP4=os.environ['CCP4']
    print(CCP4)
    from ....core import CCP4Container

    paramsContainer = CCP4Container.CContainer()
    header = paramsContainer.addHeader()
    header.setCurrent()
    header.function='DEF'
    header.pluginName.set('phaser_ensembler')
    element, errors = paramsContainer.saveContentsToEtree()
    
    root = etree.Element('container',id='phaser_ensembler')

    #Here introduce the CCP4i2 classes needed to handle input and output files
    inputDataXML = etree.fromstring(
'''
<dummyRoot>
    <container id="inputData">
        <content id="XYZIN_LIST">
            <className>CList</className>
            <qualifiers>
                <listMinLength>2</listMinLength>
            </qualifiers>
            <subItem>
                <className>CPdbDataFile</className>
                <qualifiers>
                    <mustExist>True</mustExist>
                    <allowUndefined>False</allowUndefined>
                    <fromPreviousJob>False</fromPreviousJob>
                    <ifAtomSelection>True</ifAtomSelection>
                </qualifiers>
            </subItem>
        </content>
        <content id="ALIGNIN">
            <className>CSeqAlignDataFile</className>
        </content>
      <content id="OVERRIDEID">
          <className>CFloat</className>
          <qualifiers>
              <default>90.0</default>
          </qualifiers>
      </content>
    </container>
    <container id="outputData">
        <content id="XYZOUT">
            <className>CPdbDataFile</className>
            <qualifiers>
                <default><subType>3</subType></default>
                <saveToDb>True</saveToDb>
            </qualifiers>
        </content>
    </container>
</dummyRoot>
        ''')
    for child in inputDataXML: root.append(child)

    #Get the "phil" information  as a scope to translate into def.xml
    #And flatten this into a container called keywords
    keywordElement = etree.SubElement(root,'container',id='keywords')
    flattenScope(PHIL_MASTER, keywordElement)

    paramsContainer.loadContentsFromEtree(root, overwrite=True)
    paramsContainer.saveContentsToXml(fileName='phaser_ensembler.def.xml')
