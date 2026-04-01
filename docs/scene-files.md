# Scene Files: Molecular Graphics in Reports

This document describes the CCP4i2 scene file mechanism for embedding molecular
graphics views in job reports. It covers the legacy CCP4mg-based system, the
current state of the code, and the path toward a Moorhen-based replacement.

## Overview

Scene files (`.scene.xml`) are XML descriptions of a molecular graphics view:
which data to load, how to draw it, and where to point the camera. In legacy
CCP4i2 they were consumed by CCP4mg to produce both static images and
interactive 3D views inside reports.

The mechanism has three layers:

1. **Scene templates** stored alongside each wrapper (`wrappers/<task>/script/<task>.scene.xml`)
2. **Report-time generation** where templates are cloned, customised with actual
   job data, and saved into the job directory
3. **A viewer** that reads the scene file, loads the data, and renders the view

Layers 1 and 2 are largely intact in the codebase. Layer 3 (CCP4mg) is gone in
the Django/React architecture and needs to be replaced by Moorhen.

---

## Scene XML Format

### File structure

A scene file has the standard CCP4i2 XML envelope when saved to disk:

```xml
<ccp4:ccp4i2>
  <ccp4i2_header>
    <function>MGSCENE</function>
    <projectName>...</projectName>
    <jobId>...</jobId>
    <!-- etc. -->
  </ccp4i2_header>
  <ccp4i2_body>
    <scene>
      <!-- scene content -->
    </scene>
  </ccp4i2_body>
</ccp4:ccp4i2>
```

Templates stored in the source tree omit the envelope and contain only the
`<scene>` element. The `Picture` class (`report/pictures.py`) wraps the scene
content in the envelope when saving.

A file may contain **multiple `<scene>` elements** to define a presentation
with several views.

### The `<scene>` element

Each `<scene>` contains:

| Section | Purpose |
|---------|---------|
| `<data>` | Data objects to load (`<MolData>`, `<MapData>`) |
| `<View>` | Camera: centre, orientation, scale, slab, fog |
| `<wizard>` | *(Optional)* Named CCP4mg display template |

If a `<wizard>` is present, it controls how data is drawn. Otherwise,
display is controlled by `<MolDisp>` / `<MapDisp>` children of the data
elements.

### `<data>` section

#### `<MolData>` -- coordinate data

```xml
<MolData id="id1">
  <filename database="filenames/XYZOUT" />

  <!-- Optional: custom monomer dictionaries for non-standard residues -->
  <customResCIFFiles>
    <cifmonomer>
      <name>LIG</name>
      <filename>/path/to/monomer.cif</filename>
    </cifmonomer>
  </customResCIFFiles>

  <!-- Zero or more display objects -->
  <MolDisp>...</MolDisp>
  <MolDisp>...</MolDisp>

  <!-- Optional: surface display -->
  <SurfaceDispobj>...</SurfaceDispobj>
</MolData>
```

The `id` attribute is a unique identifier (e.g. `id1`, `id2`) used to
cross-reference from `<View>`, `<wizard>`, and `<MapData model="...">`.

#### `<MapData>` -- electron density

```xml
<MapData id="id3">
  <filename database="filenames/FPHIOUT" />
  <columns>
    <F>F</F>
    <PHI>PHI</PHI>
  </columns>
  <model>id1</model>            <!-- MolData id to mask around -->
  <gridSize>0.5</gridSize>
  <contourUnits>sigma</contourUnits>
  <isDifferenceMap>False</isDifferenceMap>

  <MapDisp>
    <contourLevel>1.0</contourLevel>
    <difference>0</difference>
    <colour>ice blue</colour>
  </MapDisp>
</MapData>
```

For difference maps, use `<difColumns>` instead of `<columns>` and set
`<isDifferenceMap>True</isDifferenceMap>`:

```xml
<MapData id="id4">
  <filename database="filenames/DIFFPHIOUT" />
  <isDifferenceMap>True</isDifferenceMap>
  <difColumns>
    <F>F</F>
    <PHI>PHI</PHI>
  </difColumns>
  <model>id1</model>
  <contourUnits>sigma</contourUnits>
  <MapDisp>
    <contourLevel>3.0</contourLevel>
    <difference>1</difference>
    <colour>green</colour>
    <second_colour>red</second_colour>
  </MapDisp>
</MapData>
```

#### Data source binding

The `<filename>` element supports two mechanisms for binding to job data:

| Mechanism | Syntax | Resolution |
|-----------|--------|------------|
| **Database attribute** | `<filename database="filenames/XYZOUT" />` | Resolved at view-time by looking up the named parameter in the job container |
| **String substitution** | `<filename>SUBSTITUTEME</filename>` | Replaced in Python code during report generation |

The `database=` approach is preferred -- it keeps templates generic and
delegates resolution to the report generator.

### `<MolDisp>` -- molecular display object

Each `<MolDisp>` defines one representation layer for a `<MolData>` object:

```xml
<MolDisp>
  <select>all</select>                    <!-- Atom selection -->
  <style>CYLINDERS</style>                <!-- Rendering style -->
  <colour_parameters>
    <colour_mode>one_colour</colour_mode> <!-- or bychain, bymodel, etc. -->
    <one_colour>yellow</one_colour>
    <non_C_atomtype>1</non_C_atomtype>    <!-- Colour non-C atoms by element -->
  </colour_parameters>

  <!-- Optional: bond drawing -->
  <drawing_style>
    <show_multiple_bonds>1</show_multiple_bonds>
    <deloc_ring>1</deloc_ring>
  </drawing_style>

  <!-- Optional: atom labels -->
  <label_params>
    <label_select>label_all</label_select>
    <label_colour>complement</label_colour>
    <label_user_oneperres>1</label_user_oneperres>
    <label_style>
      <!-- 15 entries controlling which label fields to show -->
      <label_style_entry>0</label_style_entry>
      <!-- ... -->
    </label_style>
  </label_params>
</MolDisp>
```

**Rendering styles** used in the codebase:

| Style | Description |
|-------|-------------|
| `SPLINE` | Smooth ribbon/spline trace |
| `CYLINDERS` | Bond cylinders (stick model) |
| `BALLSTICK` | Ball-and-stick |
| `FATBONDS` | Thick bonds |

**Colour modes:**

| Mode | Description |
|------|-------------|
| `one_colour` | Single colour (set via `<one_colour>`) with optional element colouring for non-carbon (`<non_C_atomtype>1</non_C_atomtype>`) |
| `bychain` | Colour by chain ID |
| `bymodel` | Colour by model number |

**Atom selections** use CCP4mg selection syntax:

| Selection | Meaning |
|-----------|---------|
| `all` | All atoms |
| `peptide` | Peptide residues |
| `{not /*/*/*/*[H]}` | Exclude hydrogens |
| `neighb cid="A:ATP:501" maxd=10.0 group=all excl=central` | Neighbours within 10A of a residue, excluding the central residue |
| `A/1170(DUP)` | Specific residue by chain/number/name |

### `<SurfaceDispobj>` -- molecular surface

```xml
<SurfaceDispobj>
  <visible>1</visible>
  <transparent>1</transparent>
  <opacity>0.7</opacity>
  <select>neighb cid="A:LIG:1" maxd=8.0 group=all excl=central,solvent</select>
  <style_parameters>
    <select>cid</select>
    <cid>not {A:LIG:1 or solvent}</cid>
  </style_parameters>
</SurfaceDispobj>
```

### `<View>` -- camera setup

Two approaches: **automatic** (derived from data) or **explicit** (fixed values).

#### Automatic view

```xml
<View>
  <centre_MolData>id1</centre_MolData>
  <centre_selection>A:ATP:501</centre_selection>  <!-- Optional: centre on specific atoms -->
  <orientation_auto>
    <molData>id1</molData>
    <selection>A:ATP:501</selection>               <!-- Best orientation for this selection -->
  </orientation_auto>
  <scale_auto>1</scale_auto>                       <!-- Zoom to fit the selection -->
  <slab_enabled>0</slab_enabled>
</View>
```

If `<centre_selection>` is omitted or set to `all`, the view centres on the
entire molecule.

#### Explicit view

```xml
<View>
  <scale>58.0125</scale>
  <orientation>
    <q0>0.0134139428983</q0>
    <q1>0.336496398576</q1>
    <q2>-0.507403789037</q2>
    <q3>-0.793178186004</q3>
  </orientation>
  <centre_xyz>
    <x>-23.676</x>
    <y>8.538</y>
    <z>5.262</z>
  </centre_xyz>
  <fog_enabled>1</fog_enabled>
  <fog_near>-1.1</fog_near>
  <fog_far>26.01</fog_far>
  <slab_offset>27.5</slab_offset>
  <slab_enabled>false</slab_enabled>
  <slab_width>55.0</slab_width>
</View>
```

### `<wizard>` -- display template

The wizard mechanism delegates display setup to a named CCP4mg preset:

```xml
<wizard>
  <template>ribbons:colour chains</template>
  <parameters>
    <MolData1>id1</MolData1>           <!-- Cross-ref to <MolData id="id1"> -->
    <MapData1>id3</MapData1>           <!-- Cross-ref to <MapData id="id3"> -->
    <MapData2>id4</MapData2>
  </parameters>
</wizard>
```

The `<template>` value is in the form `folder:template name`. Templates used
in the codebase:

| Template | Used by |
|----------|---------|
| `ribbons:colour chains` | molrep_mr, molrep_den |
| `Ribbons:Colour by chain` | modelcraft |
| `map and model:cylinders map symmetry sticks` | coot_rebuild |

#### Wizard template internals (from CCP4mg documentation)

Wizard templates are Python scripts within CCP4mg stored in
`ccp4mg/data/picture_wizard_files/`, organised into style folders (e.g.
`ribbons/`, `ligand_binding_site/`, `nucleic_acid/`, `interfaces/`).

Each template file (`style_name.mgpic.py`) has four sections:

```
#SECTION TITLE
Short description shown in the Picture Wizard UI

#SECTION NOTES
User-facing documentation with hints on customising

#SECTION CHOICES
Python list defining the parameters the wizard accepts:
CHOICES = [
    { 'name': 'MolData1', 'type': 'MolData', 'label': 'Select model' },
    { 'name': 'ligand1', 'type': 'monomer', 'label': 'Ligand 1',
      'parent': 'MolData1' },
    { 'name': 'label_site_residues', 'type': 'logical',
      'label': 'Label binding site residues', 'initial': '0' },
]

#SECTION SCRIPT
Python script using picture definition file syntax (MolDisp, SurfaceDispObj,
etc.) with variables from CHOICES substituted at runtime.
```

The SCRIPT section creates display objects using picture definition file
syntax. For example a ligand binding site wizard might create:

- A `MolDisp` with `neighb` selection for residues within 3A (main/side chain)
- A `MolDisp` with `neighb` selection for H-bond capable groups within 4A
- A `MolDisp` showing context (whole protein as ribbons, or local main chain
  as worms within 6-7A)
- A `SurfaceDispObj` for the binding pocket (initially hidden)

The wizard applies intelligence beyond static templates -- it examines the
loaded data to determine chains, monomers, secondary structure, etc. and
adapts the display accordingly. Templates that need user input (e.g. interface
definitions requiring two sets of atoms) declare this via CHOICES.

The `<wizard>` element in scene XML references these templates. The
`<template>` value is in the form `folder:template name`. The `<parameters>`
element provides values for the template's CHOICES, with data object
cross-references (`MolData1`, `MapData1`, etc.) linking to the `id` attributes
of `<MolData>` and `<MapData>` elements.

Templates used in CCP4i2 scene files:

| Template | Used by |
|----------|---------|
| `ribbons:colour chains` | molrep_mr, molrep_den |
| `Ribbons:Colour by chain` | modelcraft |
| `map and model:cylinders map symmetry sticks` | coot_rebuild |

The wizard is an alternative to explicit `<MolDisp>` elements -- use one or
the other, not both.

---

## CCP4mg Picture Definition File Format

In addition to the scene XML format used by CCP4i2, CCP4mg has a human-readable
Python-based serialization format: **picture definition files** (`.mgpic.py`). This
format is important for the transliterator because:

- It is the non-binary persistence format (as opposed to `.pkl` status files)
- It exposes the full CCP4mg data model in a scriptable form
- Wizard templates are written in this format
- It can be read and written by external programs

### Data model overview

A picture definition file defines a scene as a collection of Python objects:

| Object | Purpose |
|--------|---------|
| `MolData` | Coordinate data (filename, transformation) |
| `MapData` | Electron density (filename, F/PHI columns, grid size) |
| `MolDisp` | How to draw a molecule (selection, colour, style) |
| `MapDisp` | How to draw a map contour (radius, level, style) |
| `SurfaceDispObj` | Molecular surface (selection, context, colour) |
| `HBonds` / `Contacts` | Hydrogen bonds / close contacts |
| `Annotation` | 3D text labels |
| `Crystal` | Symmetry display |
| `View` | Camera (centre, orientation, zoom) |
| `Wizard` | Reference to a wizard template |
| `ColourScheme` | Named colour rules |
| `SelectionScheme` | Named selection rules |
| `ParamsManager` | Drawing style parameters |
| `Colours` | Custom colour definitions (RGB) |

### Example

```python
MolData (
    filename = ['demo', '2ins.pdb', '/home/lizp/demo_data/2ins.pdb']
)

MolDisp (
    selection = 'all',
    colour = 'bychain',
    style = 'SPLINE'
)
```

### MolDisp attributes (key subset)

| Attribute | Type | Description |
|-----------|------|-------------|
| `selection` | string | Selection command (e.g. `'all'`, `'catrace'`, `'neighb cid="A/27" maxd=5.0'`) |
| `selection_parameters` | dict | GUI state: `select`, `ranges`, `monomers`, `neighb_*`, etc. |
| `colour` | string | Scheme name: `'atomtype'`, `'bychain'`, `'secstr'`, `'one_colour'`, `'bvalue'`, etc. |
| `colour_parameters` | dict | Details: `colour_mode`, `one_colour`, `non_C_atomtype`, `colour_rules`, `colour_blend` |
| `style` | string | `'BONDS'`, `'SPLINE'`, `'CYLINDERS'`, `'BALLSTICK'`, `'SPHERES'`, `'Surface'` |
| `drawing_style` | string | Name of associated `ParamsManager` for ribbon/cylinder widths etc. |

### Colour schemes

Colours can be specified as:
- Named schemes: `atomtype`, `bychain`, `secstr`, `bvalue`, `occupancy`, `mainside`, `thru_chain`
- Single colour: `one_colour` with `colour_parameters = { 'one_colour': 'yellow' }`
- Custom rules: `colour_rules` with list of `[colour, selection]` pairs
- Blend through model: `colour_blend` with per-chain start/end colours
- Non-carbon by element: `non_C_atomtype = 1` overlays element colours on non-C atoms

### Selection commands (full syntax from CCP4mg docs)

Selections are text strings composed of coordinate IDs (CIDs), aliases, commands,
and operators:

**CID format:** `/mdl/chn/seq(res).ic/atm[elm]:aloc`

**Aliases:**

| Alias | Meaning |
|-------|---------|
| `amino_acid` | All amino acid residues |
| `nucleic` | All nucleic acid residues |
| `solvent` | All solvent |
| `catrace` | CA atoms in amino acids |
| `main` | Main chain atoms |
| `side` | Side chain atoms |
| `peptide_link` | N, C, O, H atoms |

**Commands:**

| Command | Arguments | Description |
|---------|-----------|-------------|
| `neighb` | `cid`, `maxd`, `group`, `hbonded`, `excl` | Neighbourhood selection |
| `sphere` | `x`, `y`, `z`, `radius`, `type` | Sphere around a point |
| `range` | `start`, `end` | Residue range |

**Operators:** `or`, `and`, `not`, `xor`, `excl` (= and not)

Curly braces `{...}` group sub-expressions.

### SurfaceDispObj

```python
SurfaceDispObj (
    selection = 'neighb cid="//500(NDP)" maxd=5.0 group=residue',
    colour = 'electrostatic',
    style = 'amino_acid',          # context atoms
    surface_style = 'solid',       # or 'dots', 'mesh'
)
```

The `style` attribute is the context selection (atoms that shape the surface
but are not themselves surfaced). The `colour` can be `'electrostatic'` for
Poisson-Boltzmann potential colouring.

### View

```python
View (
    centre_xyz = [-6.78, -17.32, -0.62],
    zoom = 0.51,
    orientation = [0.572, 0.209, -0.517, 0.599]  # quaternion
)
# Or centre on a model/selection:
View (
    centre_MolData = '1abc',
    centre_selection = '//A/27',
    zoom = 0.51,
    orientation = [0.572, 0.209, -0.517, 0.599]
)
```

### Relevance to the transliterator

The picture definition file format is the primary input for the transliterator.
The scene XML used by CCP4i2 is a subset of the full CCP4mg data model --
it covers `MolData`, `MapData`, `MolDisp`, `MapDisp`, `SurfaceDispobj`, `View`,
and `Wizard`. The transliterator needs to handle this subset, but understanding
the full `.mgpic.py` format is important because:

1. **Wizard templates are written in it** -- expanding a wizard means executing
   the SCRIPT section which produces `MolDisp`, `SurfaceDispObj`, etc.
2. **The full colour/selection vocabulary** is richer than what appears in scene
   XML -- the transliterator should support the common cases
3. **It is the format CCP4mg uses for non-binary persistence** -- any tool
   migrating CCP4mg scenes will encounter `.mgpic.py` files

---

## Report Integration

### Python classes

| Class | File | Role |
|-------|------|------|
| `Picture` | `report/pictures.py` | Wraps a scene for inclusion in a report |
| `Launch` | `report/actions.py` | Button to open scene in external viewer |
| `Container.addPicture()` | `report/core.py` | Convenience method on all report containers |
| `Container.getPictures()` | `report/core.py` | Recursively collects scene file paths |

### How reports add pictures

**Simple: reference a static template**

```python
# modelcraft_report.py
baseScenePath = os.path.join(ccp4i2_root, 'wrappers', 'modelcraft',
                             'script', 'modelcraft_1.scene.xml')
pictureFold.addPicture(
    label="Autobuilt structure",
    sceneFile=str(baseScenePath),
    id="autobuild_1",
)
```

The `database=` attributes in the template will be resolved at view-time.

**Medium: clone template, substitute filenames**

```python
# qtpisa_report.py
baseSceneXML = CCP4Utils.openFileToEtree(baseScenePath)
et = etree.ElementTree(baseSceneXML)
filename_element = et.findall(".//scene/data/MolData/filename")[0]
del filename_element.attrib["database"]      # Remove lazy binding
filename_element.text = actual_path           # Set concrete path
sceneFilePath = os.path.join(jobDirectory, 'qtpisa_scene0.scene.xml')
et.write(sceneFilePath)
pic = pictureGallery.addPicture(sceneFile=sceneFilePath, label='...')
```

**Complex: build scenes programmatically**

The refmac report (`wrappers/refmac/script/refmac_report.py`) builds per-ligand
scenes from scratch using helper methods:

```python
# For each ligand in the structure:
molDataNode = etree.SubElement(dataNode, 'MolData', id='id1')
etree.SubElement(molDataNode, 'filename', database='filenames/XYZOUT')

# Draw neighbours as sticks
self.appendMolDisp(molDataNode, selectionText='neighb cid="A:LIG:1" maxd=10.0 ...',
                   carbonColour='green', style='CYLINDERS')

# Draw transparent surface of surroundings
surfDispObj = etree.SubElement(molDataNode, "SurfaceDispobj")
# ... configure opacity, selection ...

# Draw the ligand itself as ball-and-stick
self.appendMolDisp(molDataNode, selectionText='A:LIG:1',
                   carbonColour='yellow', style='BALLSTICK', bondOrder=True)

# Centre view on the ligand
viewNode = etree.SubElement(sceneNode, 'View')
etree.SubElement(viewNode, 'centre_MolData').text = 'id1'
etree.SubElement(viewNode, 'centre_selection').text = 'A:LIG:1'
```

Helper methods in `refmac_report.py`:

- **`appendMolDisp()`** -- builds a `<MolDisp>` with selection, colour, style
- **`addMap()`** -- builds a `<MapData>` with columns, contour, colour
- **`subElementWithText()`** -- utility for creating elements with text content

### Object gallery pattern

Multiple scenes are typically presented in an `ObjectGallery`, which renders
as a tabbed or selectable list of views:

```python
pictureFold = parent.addFold(label='Picture', initiallyOpen=False)
pictureGallery = pictureFold.addObjectGallery(
    style='float:left;', height='550px',
    tableWidth='260px', contentWidth='450px',
)
# Add multiple scenes to the gallery
for ligand in ligands:
    scene = build_scene_for(ligand)
    pictureGallery.addPicture(sceneFile=scene, label=f'Ligand {ligand}')
```

---

## Current State

### What exists

- **14 scene template files** across `wrappers/` and `pipelines/`
- **`Picture` and `Launch` classes** in the report framework
- **Report code** in ~10 wrappers that creates scenes
- **Moorhen integration** in `client/renderer/components/moorhen/` with full
  molecule/map loading via the REST API

### What is broken or disabled

| Wrapper | Status | Notes |
|---------|--------|-------|
| refmac | `return` at line 450 | "disabled while reconsidering their role" |
| qtpisa | `drawPictures()` commented out | |
| mrparse | `return` before picture code | "FIXME - XML PICTURE" |
| modelcraft | Active | References static template |
| molrep_mr, molrep_den | Active | Reference static templates |
| ccp4mg_*, pdbview_edit | Active | But depend on CCP4mg viewer |

### What is missing

1. **`Picture.as_data_etree()`** -- the `Picture` class has no serialization
   method, so scenes are not included in the report XML sent to the frontend
2. **Frontend scene component** -- no React component consumes scene data from
   report XML
3. **Scene-to-Moorhen interpreter** -- nothing translates scene XML concepts
   into Moorhen API calls

---

## Existing Scene Templates

Quick reference of all templates and what they display:

| Template | Data | Display | View |
|----------|------|---------|------|
| `modelcraft_1.scene.xml` | XYZOUT | Wizard: Ribbons colour by chain | Auto, centred on molecule |
| `molrep_mr_1.scene.xml` | XYZOUT + FPHIOUT + DIFFPHIOUT | Wizard: Ribbons colour chains + 2 maps | Auto, centred on molecule |
| `coot_1.scene.xml` | XYZOUT + FPHIIN + DELFPHIIN | Wizard: map and model, cylinders | Auto, centred on molecule |
| `acedrg.scene.xml` | PDB (SUBSTITUTEME) | Cylinders, yellow, no-H, labels, bond orders | Auto |
| `acedrgNew.scene.xml` | PDB (SUBSTITUTEME) | Same as acedrg | Auto |
| `qtpisa.scene.xml` | PDB (SUBSTITUTEME) | Same as acedrg | Auto |
| `ccp4mg_general.scene.xml` | XYZOUT | Spline | Auto |
| `ccp4mg_edit_model.scene.xml` | XYZOUT | Spline, colour by model | Auto |
| `ccp4mg_edit_nomrbump.scene.xml` | XYZOUT | Spline | Auto |
| `pdbview_edit.scene.xml` | XYZOUT | Spline | Auto |
| `prosmart_refmac_1.scene.xml` | *(empty)* | Used as base by refmac report | -- |
| `MakeLink.scene.xml` | *(empty)* | Placeholder | -- |
| `dr_mr_modelbuild_pipeline_1.scene.xml` | *(empty)* | Placeholder | -- |

---

## Architecture for Moorhen Migration

### Project boundaries

The work splits across three projects with distinct responsibilities:

```
CCP4i2 (this project)          CCP4mg scene transliterators        Moorhen
─────────────────────          ────────────────────────────         ───────
Scene templates &              Python library that converts         Scene interpreter:
report-time generation         CCP4mg scene XML (including          declarative JSON
                               wizard templates) into a             → Moorhen draw
Resolves:                      renderer-neutral intermediate        commands
 • database= refs → file IDs  format
 • neighb selections → CIDs                                        Offered as a
 • SUBSTITUTEME → paths        Offered to CCP4mg project as a      Moorhen-native
                               migration path away from Qt          capability
Emits resolved scene XML
       │
       ▼
Transliterator (import or call)
       │
       ▼
JSON scene description ──────────────────────────────────────────► Moorhen
```

**CCP4i2** owns scene generation and data resolution. It produces fully-resolved
scene XML where all filenames are concrete (file IDs or URLs), all selections
are standard CID strings, and all wizard references are intact.

**CCP4mg scene transliterators** are a Python library (potentially contributed
to or maintained alongside the CCP4mg project) that converts CCP4mg scene XML
-- including wizard template expansion -- into a renderer-neutral intermediate
format. This serves two audiences:

- CCP4i2 uses it at report generation time to produce Moorhen-ready JSON
- The CCP4mg project can use it as part of their own migration away from Qt,
  since the intermediate format can target any renderer built on a compatible
  data model

The transliterators need access to the CCP4mg wizard template definitions
(Python scripts within CCP4mg) to expand wizard references into explicit
display objects. The underlying data model of CCP4mg's scene graph has
significant overlap with MoleculesToTriangles (which underpins Moorhen's
rendering via Coot), making the mapping relatively natural.

**Moorhen** owns the scene interpreter: a client-side (TypeScript) module that
takes a declarative JSON scene description and translates it into Moorhen API
calls. This is a general Moorhen capability -- any application embedding Moorhen
can use it, not just CCP4i2.

### Conceptual mapping: CCP4mg → intermediate → Moorhen

| CCP4mg Scene Concept | Intermediate (JSON) | Moorhen API |
|---------------------|--------------------|----|
| `<MolData>` + `<filename>` | `{ type: "molecule", url: "..." }` | `MoorhenMolecule.loadToCootFromURL()` |
| `<MapData>` + columns | `{ type: "map", url: "...", columns: {...} }` | `MoorhenMap.loadToCootFromMtzURL()` |
| `<MolDisp style="SPLINE">` | `{ representation: "ribbon", ... }` | Ribbon representation |
| `<MolDisp style="CYLINDERS">` | `{ representation: "bonds", ... }` | Bonds representation |
| `<MolDisp style="BALLSTICK">` | `{ representation: "ballAndStick", ... }` | Ball-and-stick representation |
| `<MapDisp>` + contour/colour | `{ contourLevel: 1.0, colour: "..." }` | `MoorhenMap.setContourLevel()` |
| `<SurfaceDispobj>` | `{ representation: "surface", ... }` | Surface representation |
| `<View>` (auto) | `{ centre: { cid: "//A/501" }, autoOrient: true }` | `setViewOn()` |
| `<View>` (explicit) | `{ centre: {x,y,z}, orientation: {q0..q3} }` | Direct camera setup |
| `<customResCIFFiles>` | `{ dictionaries: [{ url: "..." }] }` | `MoorhenMolecule.addDict()` |
| `<wizard>` | `{ wizard: "ribbons_by_chain", params: {...} }` | Moorhen-native recipe |
| `database="filenames/XYZOUT"` | Resolved by CCP4i2 to file ID / URL | -- |

### Design decisions

**Selection syntax: server-side expansion**

CCP4mg selections like `neighb cid="A:LIG:1" maxd=10.0 group=all excl=central`
have no equivalent in Moorhen. Rather than implementing a custom selection
parser on the client, these are **expanded server-side at report generation
time** into concrete CID strings that Moorhen can consume directly.

At report generation the coordinate file is available via gemmi, so the
expansion is straightforward:

1. Parse the `neighb` directive (central CID, distance cutoff, exclusions)
2. Use `gemmi.NeighborSearch` to find residues within the cutoff
3. Emit a union of per-residue CID selectors (e.g. `//A/501-503,510,515/`)

This keeps the Moorhen scene interpreter simple -- it passes CID strings
through without needing to understand crystallographic geometry. The refmac
report, which is the main consumer of `neighb` selections, already loads
the structure to enumerate ligands, so adding a neighbour search there is
natural.

The trade-off is that expanded selections are verbose and frozen at generation
time, but since scenes are produced once per completed job this is fine.

**Wizards: a Moorhen-native concept, not just a compatibility shim**

CCP4mg wizard templates are Python scripts that apply display intelligence --
"colour each chain differently", "show map at 1 sigma with model as cylinders",
"highlight this ligand in its environment". They avoid requiring the user (or
the embedding application) to write extensive imperative code to produce
context-appropriate views.

If Moorhen aspires to be an illustrator of molecular scenes as well as an
editor, then a wizard-like abstraction belongs in Moorhen itself. Moorhen
wizards would be named recipes ("ligand environment", "ribbons by chain",
"map and model") that interrogate the loaded molecule and set up appropriate
representations. Because they have access to the live data, they can be
smarter than static expansions -- adapting contour levels to map statistics,
choosing colours that contrast with the chain count, deciding whether to show
a surface based on ligand size, etc.

This changes the transliterator's role: rather than expanding CCP4mg wizards
into primitive draw commands, it **maps CCP4mg wizard names to Moorhen wizard
names** -- a vocabulary translation. The actual display logic lives in Moorhen
and can evolve independently.

The CCP4mg wizard template source code (Python scripts within CCP4mg) is the
reference for what each recipe should do. Only 3 templates are used in current
CCP4i2 scene files, but Moorhen's recipe library could grow to cover the full
CCP4mg wizard vocabulary over time, supporting broader migration.

**Intermediate format: JSON, not XML**

The intermediate format between the transliterator and Moorhen should be JSON.
It maps naturally to Moorhen's TypeScript API, is easy to embed in REST
responses, and avoids coupling to CCP4mg's XML schema. CCP4i2 scene XML
remains the authoring format on the server side.

### The transliterator as a CCP4mg migration tool

The transliterator library has value beyond CCP4i2. CCP4mg is deeply Qt-based,
and CCP4 has a strategic priority to move away from the cross-platform support
burden that Qt carries. The transliterator provides a path:

1. **Parse** CCP4mg scene XML (including wizard template expansion) into a
   Python data model
2. **Emit** a renderer-neutral intermediate representation

### MoleculesToTriangles: the rendering target

MoleculesToTriangles (M2T) underpins Moorhen's rendering layer (via Coot's
management of the API). Its core class `MolecularRepresentation` is exactly
three things composed together:

```
MyMolecule + CompoundSelection + ColorScheme + renderStyle
```

This is the declarative model that CCP4mg's `.mgpic.py` format should have
been targeting. A `MolecularRepresentation` says: "for this molecule, select
these atoms, colour them with these rules, draw them in this style."

**Key M2T classes** (in `MoleculesToTriangles/CXXClasses/`):

| Class | Role |
|-------|------|
| `MolecularRepresentation` | Selection + colour + style → display primitives |
| `CompoundSelection` | Atom selection (MMDB CID strings with `&`, `\|`, `!`) |
| `ColorScheme` | Ordered list of `ColorRule`s, last match wins |
| `ColorRule` | `CompoundSelection` → colour (solid or property ramp) |
| `SolidColorRule` | Fixed colour for a selection |
| `AtomPropertyRampColorRule` | Colour ramp by B-factor, etc. |
| `Representation` | Base: holds `DisplayPrimitive`s and render params |
| `DisplayPrimitive` | Triangles, cylinders, spheres, lines, etc. |

**Render styles** in M2T (`MolecularRepresentation::redraw()`):

| M2T renderStyle | CCP4mg equivalent | Method |
|----------------|-------------------|--------|
| `"Ribbon"` | SPLINE | `drawRibbon()` |
| `"Cylinders"` | CYLINDERS | `drawBondsAsCylinders()` |
| `"Sticks"` | BONDS / FATBONDS | `drawBondsAsNewSticks()` |
| `"Spheres"` | SPHERES | `drawSpheres()` |
| `"MolecularSurface"` | Surface | `drawMolecularSurface()` |
| `"VdWSurface"` | -- | `drawVdWSurface()` |
| `"AccessibleSurface"` | -- | `drawAccessibleSurface()` |
| `"DishyBases"` | NUCLEICBASEBLOCKS | `drawDishyBases()` |
| `"StickBases"` | -- | `drawStickBases()` |
| `"HydrogenBonds"` | HBonds | `drawHydrogenBonds()` |
| `"Calpha"` | CA trace | `drawCalphas()` |

The mapping from CCP4mg to M2T is almost 1:1:

| CCP4mg MolDisp | M2T MolecularRepresentation |
|---|---|
| `selection` | `CompoundSelection` |
| `colour` / `colour_parameters` | `ColorScheme` (list of `ColorRule`s) |
| `style` | `renderStyle` string |
| Drawing params (ribbon width, etc.) | `floatParameters` / `intParameters` |

**The intermediate JSON format from the transliterator is essentially serialised
`MolecularRepresentation` descriptors.** Each element in the JSON scene
corresponds to one `MolecularRepresentation` that Moorhen will instantiate.

### Coot API: what's exposed today

M2T is not exposed directly to Moorhen. It sits behind Coot's
`molecules_container_t` API, which Moorhen calls via WASM. The relevant
methods (from `molecules-container.hh`):

**Bonds / sticks / ball-and-stick:**

| Method | Notes |
|--------|-------|
| `get_bonds_mesh(imol, mode, ...)` | Modes: `"COLOUR-BY-CHAIN-AND-DICTIONARY"`, `"CA+LIGANDS"`, `"VDW-BALLS"` |
| `get_bonds_mesh_instanced(imol, mode, ...)` | As above, instanced. Params: bond_width, atom_radius_to_bond_width_ratio, smoothness |
| `get_bonds_mesh_for_selection_instanced(imol, cid, mode, ...)` | Bonds for a CID selection only |
| `get_goodsell_style_mesh_instanced(imol, ...)` | Goodsell rendering |

**Ribbons / surfaces (via M2T):**

| Method | Notes |
|--------|-------|
| `get_molecular_representation_mesh(imol, cid, colour_scheme, style, ss_flag)` | Styles: `"Ribbon"`, `"MolecularSurface"`. Colour schemes: `"colorRampChainsScheme"`, `"colorBySecondaryScheme"`, `"Chain"` |
| `get_gaussian_surface(imol, ...)` | Gaussian surface |
| `get_map_contours_mesh(imol, x, y, z, radius, level)` | Map contours |

**Colour rules:**

| Method | Notes |
|--------|-------|
| `add_colour_rule(imol, selection_cid, colour)` | CID + hex colour. For M2T representations |
| `add_colour_rules_multi(imol, combo_string)` | Bulk rules: `"//A/1^#cc0000\|//A/2^#cb0002"` |
| `delete_colour_rules(imol)` / `get_colour_rules(imol)` | Manage rule set |
| `set_user_defined_atom_colour_by_selection(imol, ...)` | CID → colour index mapping (for bonds) |
| `set_user_defined_bond_colours(imol, colour_map)` | Colour index → RGBA |
| `set_use_bespoke_carbon_atom_colour(imol, state)` | Override carbon colour |
| `set_base_colour_for_bonds(imol, r, g, b)` | Base for colour wheel |

**M2T parameters:**

| Method | Notes |
|--------|-------|
| `M2T_updateFloatParameter(imol, name, value)` | e.g. ribbon widths, cylinder radius |
| `M2T_updateIntParameter(imol, name, value)` | e.g. angular sampling |

### Gaps for scene rendering

1. **Two separate rendering paths** -- bonds go through `get_bonds_mesh`,
   ribbons/surfaces through `get_molecular_representation_mesh`. A scene
   with "ribbon for protein + sticks for ligand" requires calling both APIs
   with different selections. There is no composite scene API.

2. **Limited colour scheme vocabulary for M2T** --
   `get_molecular_representation_mesh` accepts only 3 named colour schemes.
   The `add_colour_rule` mechanism is more flexible (arbitrary CID→colour)
   but it is unclear whether it feeds into `get_molecular_representation_mesh`
   or only into bond meshes.

3. **Limited style vocabulary for M2T** -- only `"Ribbon"` and
   `"MolecularSurface"` via `get_molecular_representation_mesh`. Cylinders,
   ball-and-stick, sticks go through `get_bonds_mesh` (a different code path
   with different colour handling).

4. **No `"bychain"` / `"atomtype"` / `"one_colour"` for M2T path** -- the
   named colour schemes for `get_molecular_representation_mesh` don't include
   the common CCP4mg colouring modes. These would need to be constructed from
   `add_colour_rule` calls.

5. Paul's own comment in the code: *"This API will change - we want to
   specify surfaces and ribbons too"* (on `export_model_molecule_as_gltf`).

### Extension point: `molecules_container_js`

Moorhen extends Coot's API via a subclass `molecules_container_js` (in
`Moorhen/wasm_src/moorhen-wrappers.cc`) that adds Moorhen-specific methods.
This class has direct access to Coot internals (`get_mol()` returns
`mmdb::Manager*`) and can use gemmi, MMDB, and any Coot data structures.

Existing Moorhen extensions in this file include `DrawGlycoBlocks`,
`DrawMoorhenMetaBalls`, `GetSecondaryStructure`, and various export methods
(`export_molecular_representation_as_obj`, `export_molecular_representation_as_3mf_xml`).

**This is the natural place to add a scene setup API.** A method like
`setup_scene_from_json(const std::string &scene_json)` could:

1. Parse the JSON scene description
2. For each representation element:
   - Set up colour rules via `add_colour_rule`
   - Set M2T parameters via `M2T_updateFloatParameter`
   - Generate the appropriate mesh (bonds or molecular representation)
3. Return the collection of meshes for Moorhen to render

Alternatively, individual methods could be added to fill the specific gaps:

- `get_molecular_representation_mesh` with expanded colour scheme support
  (accepting `add_colour_rule`-defined schemes, not just 3 named ones)
- A unified representation API that handles bonds and ribbons through one
  path with consistent colour handling
- Support for the M2T render styles that aren't currently reachable
  (`"Cylinders"`, `"Sticks"`, `"DishyBases"`, `"HydrogenBonds"`)

The choice between a monolithic scene API and incremental extensions depends
on whether the Moorhen TypeScript layer or the WASM layer should own the
scene interpretation logic. Putting it in WASM (C++) means the scene
interpreter runs close to the data and can use gemmi/MMDB directly for
things like neighbour searches. Putting it in TypeScript means more
flexibility for the Moorhen UI to participate in scene setup.

### Coot API: architectural constraints

**Per-molecule vs per-representation colour schemes**

M2T's design gives each `MolecularRepresentation` its own `ColorScheme`.
A scene with "ribbon by chain for protein" and "ball-and-stick by element
for ligand" needs two representations with different colour schemes. But
the Coot API exposes `add_colour_rule` per molecule (`imol`), not per
representation. There is no API to create named colour schemes, track them,
or associate a specific scheme with a specific representation request.

This means you can't currently set up two differently-coloured views of the
same molecule through the Coot API. Extending Coot to support per-
representation colour schemes -- or having the Moorhen WASM layer construct
M2T `ColorScheme` objects directly -- would be needed for full scene support.

**Coot's colour intelligence**

Bypassing Coot's colour management entirely would lose non-trivial logic:

- Adjusting saturation/brightness for dark vs light backgrounds
- Rotating the colour wheel based on how many models are loaded (so chains
  across different molecules get distinct colours)
- Element colour tables and carbon colour overrides

This logic applies when Coot *constructs* colour schemes. For scene
rendering, the ideal separation is:

1. **Coot builds the `ColorScheme`s** (fast, applies background/model-count
   intelligence) -- this is the "deciding what to draw" phase
2. **M2T triangulates** (expensive, produces geometry) -- this is the
   "computing the triangles" phase

The scene JSON feeds into phase 1; phase 2 is the existing M2T pipeline.

**Worker thread contention**

In Moorhen, the Coot API runs in a single WebWorker thread. This is not a
UI-blocking problem (the main thread stays responsive), but it means
representations cannot be created or updated while Coot is performing
computationally expensive operations like real-space refinement. The worker
thread serialises all operations -- refinement, mesh generation, colour
setup -- so scene rendering competes for time with interactive editing.

For the scene rendering use case (setting up a view of a completed job) this
is manageable: there is no concurrent refinement. But it means scene setup
should be designed to minimise worker thread occupancy:

- **Batch all Coot calls** -- set up all colour schemes and M2T parameters
  in one burst, then request all meshes, rather than interleaving scene
  setup with other operations
- **Pre-compute what you can outside Coot** -- selection expansion (e.g.
  `neighb` → concrete CIDs via gemmi) should happen in the transliterator
  or scene evaluator, not in the worker thread
- **Avoid redundant mesh regeneration** -- if multiple representations
  share the same colour scheme, set it up once

For Moorhen as an interactive illustrator (where representations need live
updating during editing), the deeper issue is that M2T triangulation is
invoked through Coot's API and therefore queued behind whatever Coot is
doing. A longer-term architectural change -- allowing M2T to work from a
snapshot of coordinate data independently of Coot's live molecule -- would
decouple representation updates from editing operations. This is a broader
Moorhen architecture question beyond the scope of the scene format design,
but the scene format should not introduce assumptions that make it harder.

### Future direction: M2T on gemmi with thin adapters

M2T is currently coupled to MMDB throughout (`mmdb::Atom*`, selection
handles, `mmdb::Manager`). A migration path under consideration:

1. **Thin adapter layer**: present the `mmdb::Atom` API but with a
   `gemmi::Atom` behind it. M2T rendering code (ribbons, cylinders,
   surfaces, etc.) continues unchanged
2. **AST-based selection**: replace MMDB handle-based `CompoundSelection`
   with an AST evaluator operating on adapter atoms. Prior art exists in
   CCP4i2's `CPdbDataFile` which already implements boolean-combineable
   selection over atoms without MMDB handles
3. **Snapshot decoupling**: M2T works from a `gemmi::Structure` copy,
   independent of Coot's live MMDB molecule. Enables parallel
   triangulation during refinement

This would let M2T be used with or without Coot, and would make the scene
rendering pipeline direct: JSON → AST selections + colour rules → gemmi-
backed M2T → triangles. Coot's colour intelligence (background adaptation,
colour wheel rotation) could be consulted during scene setup without
M2T depending on Coot at triangulation time.

#### Selection model: toward an AST

M2T's current `CompoundSelection` parses selection strings and composes them
via MMDB's handle system (`SKEY_AND`, `SKEY_OR`, `SKEY_NEW`). This ties
selection semantics to MMDB's runtime and makes selections opaque strings
that can only be evaluated, not inspected or transformed.

A better approach -- and the direction for the new scene format -- is an
**AST-based selection model** where:

- **Parse once, evaluate anywhere** -- the same selection tree could be
  evaluated against MMDB, gemmi, or any atom iterator
- **Transform selections** -- expanding `neighb` into concrete CIDs becomes
  a tree rewrite, not string manipulation
- **Compose cleanly** -- `and` / `or` / `not` are tree nodes, not string
  concatenation with `{` `}` brackets
- **Serialise naturally** -- the AST maps directly to JSON, no string parsing
  on the receiving end
- **Validate statically** -- check well-formedness without running against a
  molecule

The CCP4mg selection language (CIDs, aliases, `neighb`, `sphere`, `range`,
operators) is a reasonable vocabulary. The problem was never the
expressiveness -- it was that selections lived as opaque strings parsed at
the last moment. An AST representation of the same concepts is strictly
better and would serve as the selection representation in the intermediate
JSON format.

**Prior art:** CCP4i2's `CPdbDataFile` already implements logical-operator
combineable selection syntax over atoms without MMDB's handle system. This
is not hypothetical -- the design for composing selections with boolean
operators over atom properties has been done in this project and could be
adapted for a general-purpose AST selection evaluator.

This also connects to the "declarative not imperative" principle. A wizard
recipe expressed as an AST of selections, styles, and colours is data all
the way down -- no Python execution needed. The recipe *is* the AST.

### CCP4mg wizard templates: local reference

The 45 CCP4mg wizard templates are distributed with CCP4 and can be found at:

```
/Applications/ccp4-9/.../site-packages/ccp4mg/data/picture_wizard_files/
```

Organised into 10 style folders:

| Folder | Templates | Key patterns |
|--------|-----------|-------------|
| `Bonds/` | 5 | All atoms, CA trace + ligands, colour by file/chain |
| `Ribbons/` | 9 | By chain, by secondary structure, blend N→C, broken ribbons |
| `Ligand_binding_site/` | 10 | Site ± ribbons ± surface, map around site, electrostatics |
| `map_and_model/` | 2 | Cylinders + map + symmetry sticks |
| `nucleic_acid/` | 6 | Base blocks, ribbons by strand/restype |
| `interfaces/` | 5 | Buried surface, contact surfaces, residues |
| `surfaces/` | 1 | Dot region |
| `PISA_data/` | 4 | Assembly, interface, surface |
| `NMR_or_symmetry_models/` | 1 | CA trace across models |
| `utils/` | 3 | Small molecules, metals (shared helpers) |

**All 45 templates are fundamentally declarative recipes wrapped in
unnecessary imperative Python.** The Python is conditional logic interrogating
the loaded structure ("has peptide?", "has monomers?", "how many chains?")
to decide which display objects to create. Each template produces a fixed
pattern of `MolDisp`, `SurfaceDispObj`, `HBonds` objects with parameterised
selections and colours.

For example, `Ribbons/colour_chains.mgpic.py` (used by CCP4i2's
`modelcraft_1.scene.xml`) reduces to:

1. Has protein? → `MolDisp(selection=peptide, colour=bychain, style=SPLINE)`
2. Has disulphides? → `MolDisp(selection=ssbond, style=CYLINDERS, colour=bychain)`
3. Has nucleic acid? → ribbon + base cylinders + base blocks
4. Show small molecules? → standard water/metal/solute template

`Ligand_binding_site/Site_and_ribbons_by_chain.mgpic.py` (the most complex)
reduces to:

1. Protein ribbon for context
2. One `MolDisp` per ligand (ball-and-stick, atomtype colours)
3. Side chains within 4A of ligands (cylinders, bychain, non-C by element)
4. Main chain within 4A that H-bonds to ligands (cylinders)
5. Outer neighbourhood at 7A (thin bonds, hidden by default)
6. H-bonds: ligand↔protein, sidechain↔sidechain, mainchain, water
7. Electrostatic surface around ligand pocket
8. Whole-protein surface (hidden by default)

Every one of these can be expressed as a declarative JSON recipe. The Python
is compensating for the lack of a declarative way to say "neighbourhood of
parameter X" -- which an AST-based selection model handles directly.

### Summary: the rendering pipeline

```
Wizard recipe + coordinate file
        │
        ▼
    Evaluate recipe against structure
    (inspect: chains, ligands, SS, neighbours via gemmi)
        │
        ▼
    Expanded drawing primitives (JSON)
    [{ selection: AST, colorScheme: [...], renderStyle: "Ribbon", params: {...} }, ...]
        │
        ▼
    Moorhen receives JSON
        │
        ▼
    Coot API: instantiate MolecularRepresentations
    (CompoundSelection + ColorScheme + renderStyle)
        │
        ▼
    MoleculesToTriangles
    (selection → atoms → geometric primitives → coloured triangles)
        │
        ▼
    WebGL render
```

This means the transliterator could serve as infrastructure for:

- CCP4i2 report scenes → Moorhen (the immediate use case)
- Standalone CCP4mg scene files → Moorhen (for users migrating workflows)
- Any future renderer built on M2T or a compatible data model

### Open questions

- **Transliterator scope**: should start as a separate package from the outset,
  with no dependency on CCP4i2's codebase.
- **Wizard template access**: the 45 templates ship with CCP4. The
  transliterator can reverse-engineer declarative equivalents from the
  `.mgpic.py` files -- the imperative Python is boilerplate structure
  interrogation, not complex logic.
- **Selection AST design**: needs a concrete schema. The vocabulary (CIDs,
  aliases, `neighb`, `sphere`, `range`, operators) is well-defined; the
  representation as JSON AST nodes needs specification.
- **Surfaces**: M2T supports molecular, VdW, and accessible surfaces.
  Moorhen's exposure of these via Coot needs verification.
- **Static image generation**: CCP4mg could render to PNG for thumbnails.
  Moorhen is browser-only. Do we need server-side rendering for thumbnails,
  or are live-only views acceptable?
