#!/usr/bin/env python3
"""
CCP4i2 CData Container System Architecture Diagram
Shows the class hierarchy, value state tracking, qualifiers, and XML serialization.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
import numpy as np

# Set up the figure
fig, ax = plt.subplots(1, 1, figsize=(24, 18))
ax.set_xlim(0, 24)
ax.set_ylim(0, 18)
ax.set_aspect('equal')
ax.axis('off')

# Color palette
colors = {
    'base_purple': '#8B5CF6',
    'container_blue': '#3B82F6',
    'fundamental_green': '#10B981',
    'file_orange': '#F59E0B',
    'xtal_cyan': '#06B6D4',
    'state_pink': '#EC4899',
    'qualifier_red': '#EF4444',
    'xml_indigo': '#6366F1',
    'text_dark': '#1F2937',
    'white': '#FFFFFF',
    'arrow': '#6B7280',
    'section_purple': '#EDE9FE',
    'section_blue': '#DBEAFE',
    'section_green': '#D1FAE5',
    'section_orange': '#FEF3C7',
}

def draw_class_box(ax, x, y, width, height, color, class_name, attributes=None,
                   methods=None, is_abstract=False):
    """Draw a UML-style class box"""
    # Main box
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.01,rounding_size=0.05",
                         facecolor='white', edgecolor=color,
                         linewidth=2)
    ax.add_patch(box)

    # Title bar
    title_height = 0.4
    title_box = FancyBboxPatch((x, y + height - title_height), width, title_height,
                                boxstyle="round,pad=0.01,rounding_size=0.05",
                                facecolor=color, edgecolor='none')
    ax.add_patch(title_box)

    # Class name
    style = 'italic' if is_abstract else 'normal'
    ax.text(x + width/2, y + height - title_height/2, class_name,
            fontsize=9, ha='center', va='center', color='white',
            fontweight='bold', style=style)

    # Attributes
    if attributes:
        attr_y = y + height - title_height - 0.25
        for attr in attributes:
            ax.text(x + 0.1, attr_y, attr, fontsize=6.5, va='center',
                    color=colors['text_dark'])
            attr_y -= 0.22

    # Methods
    if methods:
        # Separator line
        sep_y = y + height - title_height - 0.25 - (len(attributes or []) * 0.22) - 0.1
        ax.plot([x + 0.05, x + width - 0.05], [sep_y, sep_y],
                color='#D1D5DB', linewidth=1)

        method_y = sep_y - 0.2
        for method in methods:
            ax.text(x + 0.1, method_y, method, fontsize=6, va='center',
                    color=colors['container_blue'])
            method_y -= 0.2

def draw_rounded_box(ax, x, y, width, height, color, label, sublabel=None,
                     text_color='white', fontsize=9):
    """Draw a rounded rectangle with label"""
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.02,rounding_size=0.1",
                         facecolor=color, edgecolor='none', alpha=0.95)
    ax.add_patch(box)

    label_y = y + height*0.65 if sublabel else y + height/2
    ax.text(x + width/2, label_y, label, fontsize=fontsize,
            ha='center', va='center', color=text_color, fontweight='bold')

    if sublabel:
        ax.text(x + width/2, y + height*0.3, sublabel, fontsize=fontsize-2,
                ha='center', va='center', color=text_color, alpha=0.9)

def draw_section_box(ax, x, y, width, height, color, title):
    """Draw a section box with title"""
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.01,rounding_size=0.1",
                         facecolor=color, edgecolor='#9CA3AF',
                         linewidth=1, alpha=0.4)
    ax.add_patch(box)
    ax.text(x + 0.15, y + height - 0.2, title, fontsize=11,
            ha='left', va='top', color=colors['text_dark'], fontweight='bold')

def draw_inheritance_arrow(ax, start, end):
    """Draw inheritance arrow (hollow triangle)"""
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle='-|>', color=colors['arrow'],
                               lw=1.5, fc='white'))

def draw_composition_arrow(ax, start, end, label=None):
    """Draw composition arrow (filled diamond)"""
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle='->', color=colors['arrow'], lw=1.5))
    if label:
        mid_x = (start[0] + end[0]) / 2
        mid_y = (start[1] + end[1]) / 2
        ax.text(mid_x + 0.1, mid_y, label, fontsize=6, color=colors['arrow'])

# Title
ax.text(12, 17.6, 'CCP4i2 CData Container System', fontsize=20,
        ha='center', va='center', fontweight='bold', color=colors['text_dark'])
ax.text(12, 17.2, 'Hierarchical Data Model with XML Serialization', fontsize=12,
        ha='center', va='center', color=colors['arrow'])

# ==================== CLASS HIERARCHY ====================

draw_section_box(ax, 0.2, 11, 11.5, 6, colors['section_purple'], 'Class Hierarchy')

# HierarchicalObject (base)
draw_class_box(ax, 0.5, 14.8, 3.2, 1.8, colors['base_purple'],
               'HierarchicalObject',
               ['_parent: weakref', '_children: list'],
               ['set_parent()', 'add_child()'],
               is_abstract=True)

# CData
draw_class_box(ax, 4.2, 14.3, 3.5, 2.3, colors['base_purple'],
               'CData',
               ['_value', '_value_states: dict', '_qualifiers: dict', 'name: str'],
               ['getEtree()', 'setEtree()', 'get_qualifier()', 'validity()'],
               is_abstract=True)

# CContainer
draw_class_box(ax, 8.2, 14.5, 3, 2, colors['container_blue'],
               'CContainer',
               ['_data_order: list', 'CONTENT_ORDER'],
               ['addContent()', 'addObject()', 'dataOrder()', 'copyData()'])

# Fundamental types
draw_class_box(ax, 0.5, 11.5, 2.2, 2.3, colors['fundamental_green'],
               'CInt',
               ['value: int', 'min', 'max'],
               ['validity()'])

draw_class_box(ax, 2.9, 11.5, 2.2, 2.3, colors['fundamental_green'],
               'CFloat',
               ['value: float', 'min', 'max'],
               ['validity()'])

draw_class_box(ax, 5.3, 11.5, 2.4, 2.3, colors['fundamental_green'],
               'CString',
               ['value: str', 'maxLength', 'enumerators'],
               ['validity()'])

draw_class_box(ax, 7.9, 11.5, 2, 2.3, colors['fundamental_green'],
               'CBoolean',
               ['value: bool'],
               ['validity()'])

draw_class_box(ax, 10.1, 11.5, 1.4, 2.3, colors['fundamental_green'],
               'CList',
               ['items[]'],
               ['append()'])

# Inheritance arrows
draw_inheritance_arrow(ax, (5.95, 14.3), (5.05, 16.6))  # CData -> HierarchicalObject
draw_inheritance_arrow(ax, (9.7, 14.5), (6.95, 14.5))    # CContainer -> CData

# From CData to fundamental types
for x in [1.6, 4.0, 6.5, 8.9, 10.8]:
    draw_inheritance_arrow(ax, (x, 13.8), (x, 14.3))

# ==================== FILE TYPES ====================

draw_section_box(ax, 12, 11, 5.5, 6, colors['section_orange'], 'File Types')

draw_class_box(ax, 12.3, 14.8, 2.5, 1.8, colors['file_orange'],
               'CDataFile',
               ['path', 'projectId', 'annotation'],
               ['exists()', 'fullPath()'])

# Specialized file types
file_types = [
    ('CMtzDataFile', 'MTZ reflection'),
    ('CGenericReflFile', 'Generic refl'),
    ('CMapDataFile', 'CCP4 map'),
    ('CMmcifDataFile', 'mmCIF'),
]
for i, (name, desc) in enumerate(file_types):
    x = 12.3 + (i % 2) * 2.6
    y = 12.8 - (i // 2) * 1.2
    draw_rounded_box(ax, x, y, 2.4, 0.9, colors['file_orange'],
                     name, desc, fontsize=7, text_color='#1F2937')

draw_inheritance_arrow(ax, (13.55, 14.8), (5.95, 14.5))  # CDataFile -> CData

# ==================== XTAL DATA ====================

draw_section_box(ax, 17.8, 11, 6, 6, colors['section_blue'], 'Crystallography Types')

xtal_types = [
    ('CFloatRange', 'start, end: float'),
    ('CIntRange', 'start, end: int'),
    ('CCell', 'a,b,c,alpha,beta,gamma'),
    ('CSpaceGroup', 'name, number'),
    ('CResolutionRange', 'resLow, resHigh'),
    ('CSeqData', 'sequence: str'),
]

for i, (name, attrs) in enumerate(xtal_types):
    x = 18 + (i % 2) * 3
    y = 15.8 - (i // 2) * 1.3
    draw_class_box(ax, x, y, 2.8, 1.1, colors['xtal_cyan'],
                   name, [attrs])

# ==================== VALUE STATE SYSTEM ====================

draw_section_box(ax, 0.2, 5.5, 7.5, 5, colors['section_pink'] if 'section_pink' in colors else '#FCE7F3',
                 'Value State Tracking')

# ValueState enum
draw_rounded_box(ax, 0.5, 8.5, 3, 1.5, colors['state_pink'],
                 'ValueState (Enum)', 'Tracks field assignment', fontsize=8)

states = [
    ('NOT_SET', 'Never assigned - excluded from XML'),
    ('DEFAULT', 'Using qualifier default'),
    ('EXPLICITLY_SET', 'User assigned - always serializes'),
]
for i, (state, desc) in enumerate(states):
    y = 7.9 - i * 0.7
    draw_rounded_box(ax, 0.6, y, 1.8, 0.5, '#DB2777',
                     state, fontsize=6.5)
    ax.text(2.5, y + 0.25, desc, fontsize=6, va='center', color=colors['text_dark'])

# State tracking in CData
ax.text(4, 10, 'State Tracking:', fontsize=9, fontweight='bold', color=colors['text_dark'])
code = '''obj._value_states = {
  'start': ValueState.NOT_SET,
  'end': ValueState.EXPLICITLY_SET,
}

# Check if field is set
if obj.isSet('start'):
  process(obj.start)'''
ax.text(4, 9.6, code, fontsize=6, family='monospace', va='top', color=colors['text_dark'])

# ==================== QUALIFIER SYSTEM ====================

draw_section_box(ax, 8, 5.5, 8, 5, colors['section_orange'], 'Qualifier System (Metadata)')

# Qualifier categories
qual_cats = [
    ('Numeric', ['min', 'max']),
    ('String', ['minLength', 'maxLength', 'enumerators', 'onlyEnumerators']),
    ('File', ['mustExist', 'fromPreviousJob', 'mimeTypeName', 'saveToDb']),
    ('UI', ['guiLabel', 'toolTip', 'charWidth']),
    ('List', ['listMinLength', 'listMaxLength']),
]

y_start = 9.8
for cat, quals in qual_cats:
    ax.text(8.2, y_start, f'{cat}:', fontsize=8, fontweight='bold', color=colors['text_dark'])
    ax.text(9.5, y_start, ', '.join(quals), fontsize=7, color=colors['text_dark'])
    y_start -= 0.5

# Decorator example
ax.text(12.5, 10, '@cdata_class Decorator:', fontsize=9, fontweight='bold', color=colors['text_dark'])
decorator_code = '''@cdata_class(
  qualifiers={
    'min': 0,
    'max': 100,
    'guiLabel': 'Count',
  }
)
class CCount(CInt):
  pass'''
ax.text(12.5, 9.6, decorator_code, fontsize=6, family='monospace', va='top', color=colors['text_dark'])

# ==================== XML SERIALIZATION ====================

draw_section_box(ax, 16.3, 5.5, 7.5, 5, colors['section_purple'], 'XML Serialization')

ax.text(16.5, 10, 'getEtree() / setEtree():', fontsize=9, fontweight='bold', color=colors['text_dark'])

xml_example = '''<container id="inputData">
  <content id="HKLIN">
    <className>CMtzDataFile</className>
    <qualifiers>
      <guiLabel>Input MTZ</guiLabel>
      <mustExist>True</mustExist>
    </qualifiers>
    <value>/path/file.mtz</value>
  </content>
  <content id="RESOLUTION">
    <className>CFloatRange</className>
    <start>2.0</start>
    <end>50.0</end>
  </content>
</container>'''
ax.text(16.5, 9.5, xml_example, fontsize=5.5, family='monospace', va='top', color=colors['text_dark'])

# ==================== DATA FLOW ====================

draw_section_box(ax, 0.2, 0.5, 23.6, 4.5, colors['section_blue'], 'Data Flow: Definition to Execution')

# Flow boxes
flow_steps = [
    ('1', '.def.xml', 'Task Definition', colors['xml_indigo']),
    ('2', 'DefXmlParser', 'Load Definition', colors['base_purple']),
    ('3', 'CContainer', 'inputData\nparameters\noutputData', colors['container_blue']),
    ('4', 'User Sets\nValues', 'API / GUI', colors['fundamental_green']),
    ('5', 'getEtree()', 'Serialize to XML', colors['xml_indigo']),
    ('6', 'Storage', 'DB / File / Queue', colors['file_orange']),
    ('7', 'setEtree()', 'Deserialize', colors['xml_indigo']),
    ('8', 'Execute', 'Task Wrapper', colors['xtal_cyan']),
]

for i, (num, title, desc, color) in enumerate(flow_steps):
    x = 0.5 + i * 2.9
    draw_rounded_box(ax, x, 1.5, 2.6, 2.5, color, title, desc, fontsize=8)
    # Number badge
    circle = plt.Circle((x + 0.3, 4.2), 0.25, color=color)
    ax.add_patch(circle)
    ax.text(x + 0.3, 4.2, num, fontsize=8, ha='center', va='center',
            color='white', fontweight='bold')

# Arrows between flow steps
for i in range(len(flow_steps) - 1):
    x = 0.5 + i * 2.9 + 2.6
    ax.annotate('', xy=(x + 0.2, 2.75), xytext=(x, 2.75),
                arrowprops=dict(arrowstyle='->', color=colors['arrow'], lw=1.5))

# ==================== LEGEND ====================

ax.text(0.5, 4.7, 'Color Legend:', fontsize=9, fontweight='bold', color=colors['text_dark'])
legend_items = [
    (colors['base_purple'], 'Base Classes'),
    (colors['container_blue'], 'Containers'),
    (colors['fundamental_green'], 'Fundamental Types'),
    (colors['file_orange'], 'File Types'),
    (colors['xtal_cyan'], 'Xtal Data'),
    (colors['state_pink'], 'Value States'),
]

for i, (color, label) in enumerate(legend_items):
    x = 3 + i * 2.8
    rect = Rectangle((x, 4.55), 0.3, 0.3, facecolor=color, edgecolor='none')
    ax.add_patch(rect)
    ax.text(x + 0.4, 4.7, label, fontsize=7, va='center', color=colors['text_dark'])

# Key patterns
ax.text(20, 4.7, 'Key Patterns:', fontsize=8, fontweight='bold', color=colors['text_dark'])
ax.text(20, 4.3, '- Hierarchical parent-child (weak refs)', fontsize=6.5, color=colors['text_dark'])
ax.text(20, 4.0, '- Metadata-driven via qualifiers', fontsize=6.5, color=colors['text_dark'])
ax.text(20, 3.7, '- XML round-trip preservation', fontsize=6.5, color=colors['text_dark'])

plt.tight_layout()
plt.savefig('/Users/nmemn/Developer/ccp4i2/docs/diagrams/cdata_architecture.png',
            dpi=200, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig('/Users/nmemn/Developer/ccp4i2/docs/diagrams/cdata_architecture.svg',
            format='svg', bbox_inches='tight', facecolor='white', edgecolor='none')
print("Saved: cdata_architecture.png and .svg")
