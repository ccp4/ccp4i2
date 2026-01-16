# Based on create_def_xml.py for xia2 written by David Waterman

"""create phaser_MR_AUTO.def.xml from PHIL parameters"""
import os
import re
import io

import phaser
from lxml import etree

from ccp4i2.utils.phil_handlers import Phil2Etree, PhilTaskCreator


class PhaserPhil2Etree(Phil2Etree):
  
  def __init__(self, phil_scope, modes):
    self.phil_scope = phil_scope
    self.modes = modes
  
  def __call__(self, root_id):
    phil_params = etree.Element('container', id=root_id)
    self.convertScope(self.phil_scope, phil_params, self.modes)
    return phil_params
    
  def is_selected_mode(self, obj, modes):
    if 'composition' not in obj.full_path():
      style = obj.style
      parent = obj.primary_parent_scope.name
      if parent == 'keywords' and style is not None:
        return len(set(re.split('[:, ]', style)) & set(modes)) > 0
      else: # some scopes have 2 styles
        if style is not None:
          if 'phaser:mode' in style: 
            return len(set(re.split('[:, ]', style)) & set(modes)) > 0
        return True
    else:
      return False
  
  def make_keyword(self, phil_path):
    import re
    if 'general' in phil_path:
      return '_'.join(kw.upper()[:4] for kw in re.split('[._]',phil_path.replace('phaser.keywords.general.','')))
    else:
      return '_'.join(kw.upper()[:4] for kw in re.split('[._]',phil_path.replace('phaser.keywords.','')))
    
  def make_menu_text(self, words, caption):
    if len(words) == len(caption.split()):
      return ','.join(caption.split()).replace('_',' ')
    else:
      return None
    
  def attributesToQualifiers(self, qualifiers, phil_obj, modes, force_label=True):
    # guiLabel
    if phil_obj.short_caption is not None:
      gui_label = etree.SubElement(qualifiers, "guiLabel")
      gui_label.text = self.sanitize_text(phil_obj.short_caption)
    elif force_label:
      gui_label = etree.SubElement(qualifiers, "guiLabel")
      gui_label.text = self.sanitize_text(phil_obj.name)
    # toolTip
    if phil_obj.help is not None:
      toolTip = etree.SubElement(qualifiers, "toolTip")
      toolTip.text = self.sanitize_text(phil_obj.help)
    # guiDefinition
    guiDefinition = etree.SubElement(qualifiers, "guiDefinition")
    if phil_obj.style is not None and 'phaser:mode' in phil_obj.style:
      phaser_mode = etree.SubElement(guiDefinition, "phaserMode")
      phaser_mode.text = str(','.join([mode for mode in re.split('[:, ]', phil_obj.style) if mode in modes]))
    if phil_obj.expert_level is not None:
      expert_level = etree.SubElement(guiDefinition, "expertLevel")
      expert_level.text = str(phil_obj.expert_level)
    if phil_obj.multiple is not None:
      multiple = etree.SubElement(guiDefinition, "multiple")
      multiple.text = self.sanitize_text(str(phil_obj.multiple))
    return

  def definitionToElement(self, keyword, phil_params):
    value = keyword.extract()
    elem = etree.SubElement(phil_params, "content")
    elem.set("id", self.make_keyword(keyword.full_path()))
    elem_class = etree.SubElement(elem, "className")

    # Map phil type to class and qualifiers
    phil_type = keyword.type.phil_type
    if phil_type == "bool" and str(value) not in ["True", "False"]:
      phil_type = "ternary"
    elem_class.text = self.phil_type_as_class_name.get(phil_type, "CString")
    qualifiers = etree.SubElement(elem, "qualifiers")
    self.attributesToQualifiers(qualifiers, keyword, modes=modes)

    # Set defaults for strings and bools
    if (phil_type in ["bool", "str"]) and (value is not None):
      default = etree.SubElement(qualifiers, "default")
      default.text = self.sanitize_text(str(value))

    # Set default for ternary logic, treated like a choice
    elif (phil_type == "ternary"):
      if value is None:
        value = "Auto"
      enum = etree.SubElement(qualifiers, "enumerators")
      enum.text = "True,False," + self.sanitize_text(str(value))
      menu_text = etree.SubElement(qualifiers, "menuText")
      menu_text.text = "Yes,No," + value if value ==  "Auto" else self.sanitize_text(str(value))
      default = etree.SubElement(qualifiers, "default")
      default.text = self.sanitize_text(str(value))
      only = etree.SubElement(qualifiers, "onlyEnumerators")
      only.text = "True"

    # Set attributes for choices
    elif (phil_type == "choice"):
      if (keyword.type.multi):
        # enumerators do not map to PHIL's multi choice well. In that case
        # just use a string.
        default = etree.SubElement(qualifiers, "default")
        default.text = keyword.as_str().split('=')[1].strip()
      else:
        enum = etree.SubElement(qualifiers, "enumerators")
        enum.text = self.parse_choice_options(keyword)
        if keyword.caption is not None:
          menuText = self.make_menu_text(keyword.words, keyword.caption)
          if menuText is not None:
            menu_text = etree.SubElement(qualifiers, "menuText")
            menu_text.text = self.sanitize_text(str(menuText))
        default = etree.SubElement(qualifiers, "default")
        if isinstance(value, list):
          value = " ".join(["*" + self.sanitize_text(str(v)) for v in value])
        default.text = self.sanitize_text(str(value))
        only = etree.SubElement(qualifiers, "onlyEnumerators")
        only.text = "True"

    # Set defaults and limits for numerics
    elif (phil_type in ["int", "float"]) :
      default = etree.SubElement(qualifiers, "default")
      default.text = self.sanitize_text(str(value))
      if (keyword.type.value_min is not None) :
        min = etree.SubElement(qualifiers, "min")
        min.text = str(keyword.type.value_min)
      if (keyword.type.value_max is not None) :
        max = etree.SubElement(qualifiers, "max")
        max.text = str(keyword.type.value_max)
    return elem
    
  def convertScope(self, scope, container, modes=None):
    for obj in scope.objects:
      if obj.is_definition and self.is_selected_mode(obj, modes) and 'keywords' in  obj.full_path():
        keywordElement = self.definitionToElement(obj, container)
        if self.valid_keyword('set'+self.make_keyword(obj.full_path())):
          container.append(keywordElement)
          print('adding: ', obj.full_path())
        else:
          print('not adding: ', obj.full_path())
      elif obj.is_scope and self.is_selected_mode(obj, modes):
        if obj.primary_parent_scope.name == 'keywords':
          sub_container = etree.Element('container', id=obj.full_path().split('.')[-1])
          qualifiers = etree.SubElement(sub_container, "qualifiers")
          self.attributesToQualifiers(qualifiers, obj, modes=modes, force_label=False)
          container.append(sub_container)
          self.convertScope(obj, sub_container, modes=modes)
        else:
          self.convertScope(obj, container, modes=modes)
    return
  
  def valid_keyword(self, kw): 
    inp_obj = phaser.InputMR_AUTO()
    return hasattr(inp_obj, kw) and callable(getattr(inp_obj, kw))
  
class PhaserKeywordsCreator(PhilTaskCreator):
  
  def __init__(self, phaser_phil, defaults_file, debug=False):

    # import PHIL scope
    from iotbx.phil import parse
    self.defaults_file =  defaults_file
    self.phaser_phil = phaser_phil
    defaults_phil = parse(file_name=self.defaults_file)
    master_phil = parse(file_name=self.phaser_phil).fetch(source=defaults_phil)
    PhilTaskCreator.__init__(self, master_phil, debug)
    self.fmt_dic['PLUGINNAME'] = "phaser_MR"
    # Translate the PHIL scope into an etree
    p2e = PhaserPhil2Etree(master_phil, modes)
    self.phil_tree = p2e(root_id='keywords')

    # Minimal boilerplate XML string
    self.boilerPlateXML = '''<ccp4i2>
    <ccp4i2_header>
        <function>DEF</function>
        <pluginName>{PLUGINNAME}</pluginName>
    </ccp4i2_header>
    <ccp4i2_body id="{PLUGINNAME}">
        <container id="inputData">
        </container>
        <container id="guiParameters">
          <content id="EXPERT_LEVEL">
            <className>CInt</className>
            <qualifiers>
              <enumerators>0,1,2,3</enumerators>
              <menuText>Basic,Intermediate,Expert,Developer</menuText>
              <default>1</default>
              <onlyEnumerators>True</onlyEnumerators>
            </qualifiers>
          </content>
        </container>
        <container id="keywords">
        </container>
        <container id="outputData">
        </container>
    </ccp4i2_body>
    </ccp4i2>
    '''
    
    self.inputDataXML = etree.fromstring('''
    <container id="inputData">
        <content id="F_OR_I">
            <className>CString</className>
            <qualifiers>
                <onlyEnumerators>True</onlyEnumerators>
                <enumerators>F,I</enumerators>
                <guiLabel>Use Fs or Is (depending on availability)</guiLabel>
                <default>I</default>
            </qualifiers>
        </content>
        <content id="F_SIGF">
            <className>CObsDataFile</className>
            <qualifiers>
                <mustExist>True</mustExist>
                <allowUndefined>False</allowUndefined>
                <fromPreviousJob>True</fromPreviousJob>
                <requiredContentFlag>1,2,3,4</requiredContentFlag>
            </qualifiers>
        </content>
        <content id="I_SIGI">
            <className>CObsDataFile</className>
            <qualifiers>
                <mustExist>True</mustExist>
                <allowUndefined>True</allowUndefined>
                <fromPreviousJob>True</fromPreviousJob>
                <requiredContentFlag>1,3</requiredContentFlag>
            </qualifiers>
        </content>
        <content id="FPHI">
            <className>CMapCoeffsDataFile</className>
            <qualifiers>
                <mustExist>False</mustExist>
                <allowUndefined>True</allowUndefined>
                <fromPreviousJob>True</fromPreviousJob>
            </qualifiers>
        </content>
        <content id="ABCD">
            <className>CPhsDataFile</className>
            <qualifiers>
                <mustExist>False</mustExist>
                <allowUndefined>True</allowUndefined>
                <fromPreviousJob>False</fromPreviousJob>
                <sameCrystalAs>F_SIGF</sameCrystalAs>
            </qualifiers>
        </content>
        <content id="KILLFILEPATH">
            <className>CFilePath</className>
            <qualifiers>
                <mustExist>False</mustExist>
                <allowUndefined>True</allowUndefined>
            </qualifiers>
        </content>
        <content id="ASU_PROTEIN_MW">
            <className>CFloat</className>
            <qualifiers>
                <min>1</min>
            </qualifiers>
        </content>
        <content id="ASU_NUCLEICACID_MW">
            <className>CFloat</className>
            <qualifiers>
                <min>1</min>
            </qualifiers>
        </content>
        <content id="ASUFILE">
            <className>CAsuDataFile</className>
            <qualifiers>
                <mustExist>True</mustExist>
                <fromPreviousJob>True</fromPreviousJob>
                <allowUndefined>True</allowUndefined>
                <selectionMode>2</selectionMode>
            </qualifiers>
        </content>
       <content id="ENSEMBLES">
            <className>CEnsembleList</className>
            <qualifiers>
                <listMinLength>1</listMinLength>
                <saveToDb>True</saveToDb>
            </qualifiers>
        </content>
        <content id="FIXENSEMBLES">
            <className>CList</className>
            <subItem>
                <className>CString</className>
                <qualifiers>
                    <allowUndefined>False</allowUndefined>
                </qualifiers>
            </subItem>
            <qualifiers>
                <listMinLength>0</listMinLength>
            </qualifiers>
        </content>
    </container>
    ''')

  def __call__(self):

    task_xml = etree.fromstring(self.boilerPlateXML.format(**self.fmt_dic))

    # Insert inputData
    if self.inputDataXML is not None:
      self._replace_element(task_xml, 'inputData', self.inputDataXML)

    # Replace keywords in task_xml with the phil_tree
    self._replace_element(task_xml, 'keywords', self.phil_tree)

    # Write out prettified version
    out_file = self.fmt_dic['PLUGINNAME'] + '.def.xml'
    parser = etree.XMLParser(remove_blank_text=True)
    tree = etree.parse(io.StringIO(etree.tostring(task_xml).decode("utf-8")), parser)
    try:
      with open(out_file, 'wb') as f:
        print('Writing def.xml to %s' % out_file)
        f.write(etree.tostring(tree, pretty_print=True, xml_declaration=True))
    except OSError as exception:
      if exception.errno == errno.EACCES:
        raise RuntimeError('No write permission to this directory')
    return

if __name__ == '__main__':
  modes = ['MR_AUTO', 'MR*', 'MR_FRF', 'MR_FTF', 'MR_PAK', 'MR_RNP', '*']
  phaser_phil = os.path.join(os.environ['CCP4'],'lib','py2','site-packages','phaser','phaser','phenix_interface','__init__.params')
  defaults_file = os.path.join(os.environ['CCP4'],'lib','py2','cctbx','include','phaser_defaults.params')
  pkc = PhaserKeywordsCreator(phaser_phil=phaser_phil, defaults_file=defaults_file)
  pkc()
