# `def.xml` Reference

A task's **`def.xml`** declares its data model: the input files it consumes, the
output files it produces, and its control parameters. The server uses it to build
the task's parameter container, drive validation, and — when no bespoke UI is
registered — **auto-render the task interface** (see
[Authoring a Task](authoring-a-task.md)).

`def.xml` is one of **two** ways to declare parameters. For tools that already
expose a PHIL `master_phil` (Phenix, PhaserTNG, DIALS, …) you can skip the
control-parameter boilerplate entirely and ingest the PHIL instead — see
[PHIL_TASK_GUIDE.md](../server/ccp4i2/wrappers/PHIL_TASK_GUIDE.md). Even PHIL
tasks keep a small `def.xml` for their `inputData`/`outputData`.

> This page documents the common, load-bearing cases. For the full set of types
> and qualifiers, read a few real wrappers under `server/ccp4i2/wrappers/*/script/*.def.xml`
> and the `CData` subclasses in `core/`.

---

## Skeleton

Every `def.xml` has the same shape: a header, then a body with **three
containers** — `inputData`, `outputData`, `controlParameters`.

```xml
<?xml version='1.0' encoding='ASCII'?>
<ccp4:ccp4i2 xmlns:ccp4="http://www.ccp4.ac.uk/ccp4ns">
  <ccp4i2_header>
    <function>DEF</function>
    <pluginName>mytask</pluginName>        <!-- must match the task name -->
  </ccp4i2_header>
  <ccp4i2_body id="mytask">
    <container id="inputData">        <!-- files/data the task consumes -->
      ...
    </container>
    <container id="outputData">       <!-- files the task produces -->
      ...
    </container>
    <container id="controlParameters"><!-- scalar options / keywords -->
      ...
    </container>
  </ccp4i2_body>
</ccp4:ccp4i2>
```

Each parameter is a `<content id="NAME">` with a `<className>` (its `CData` type)
and a `<qualifiers>` block. The `id` is the name you address the parameter by
everywhere else (wrapper code, i2run, the frontend, validation `name` fields).

---

## Worked example (`freerflag`)

A complete, real task — good to copy from
([source](../server/ccp4i2/wrappers/freerflag/script/freerflag.def.xml)):

```xml
<container id="inputData">
  <content id="F_SIGF">
    <className>CObsDataFile</className>
    <qualifiers>
      <mustExist>True</mustExist>
      <allowUndefined>False</allowUndefined>
      <fromPreviousJob>True</fromPreviousJob>
    </qualifiers>
  </content>
  <content id="FREERFLAG">
    <className>CFreeRDataFile</className>
    <qualifiers>
      <mustExist>True</mustExist>
      <sameCrystalAs>F_SIGF</sameCrystalAs>   <!-- cross-parameter constraint -->
    </qualifiers>
  </content>
</container>

<container id="outputData">
  <content id="FREEROUT">
    <className>CFreeRDataFile</className>
  </content>
</container>

<container id="controlParameters">
  <content id="COMPLETE">
    <className>CBoolean</className>
    <qualifiers><default>True</default></qualifiers>
  </content>
  <content id="FRAC">
    <className>CFloat</className>
    <qualifiers>
      <toolTip>Fraction in freeR set</toolTip>
      <allowUndefined>True</allowUndefined>
      <default>0.05</default>
      <min>0.0</min>
      <max>1.0</max>
    </qualifiers>
  </content>
  <content id="GEN_MODE">
    <className>CString</className>
    <qualifiers>
      <onlyEnumerators>True</onlyEnumerators>
      <menuText>Generate a new set,complete an existing set</menuText>
      <enumerators>GEN_NEW,COMPLETE</enumerators>
      <default>GEN_NEW</default>
    </qualifiers>
  </content>
</container>
```

---

## Common `CData` types (`className`)

| `className` | Use for | Notable qualifiers |
|---|---|---|
| `CBoolean` | on/off flag → checkbox | `default` |
| `CInt` | integer → number field | `default`, `min`, `max` |
| `CFloat` | real number → number field | `default`, `min`, `max` |
| `CString` | free text, or a **menu** when `onlyEnumerators` set | `enumerators`, `menuText`, `onlyEnumerators`, `default` |
| `CMtzDataFile` | reflection MTZ | file qualifiers below |
| `CObsDataFile` | observed intensities/amplitudes | `contentFlag`, `subType` |
| `CFreeRDataFile` | free-R flags | `sameCrystalAs` |
| `CMapCoeffsDataFile` / `CPhsDataFile` | map coefficients / phases | `contentFlag` |
| `CPdbDataFile` | coordinates (PDB/mmCIF) | file qualifiers below |
| `CSeqDataFile` / `CAsuDataFile` | sequence / ASU contents | |
| `CDictDataFile` | monomer/restraint dictionary (CIF) | |
| `CList` | a repeatable list of any of the above | `listMinLength`, `listMaxLength` |

> For repeatable inputs (a tool that takes *N* coordinate files), use a `CList`
> whose sub-item is the file type. PHIL tasks handle multiplicity automatically
> via `.multiple = True` + a list shim — see the PHIL guide.

---

## Common qualifiers

**Validation / presence**
- `mustExist` — the referenced file must exist on disk (checked at run time).
- `allowUndefined` — if `False`, the parameter is **required** (blocks the run
  when unset). This is the main lever for "required input" errors.
- `default` — initial value.
- `min` / `max` — numeric bounds (also enforced in the UI).

**Menus (enumerated `CString`)**
- `onlyEnumerators` — restrict to the listed values (renders as a dropdown).
- `enumerators` — comma-separated stored values (e.g. `GEN_NEW,COMPLETE`).
- `menuText` — comma-separated human labels, positionally matched to
  `enumerators`.

**UI hints**
- `guiLabel` — label shown in the interface (falls back to the `id`).
- `toolTip` — hover help.

**File data**
- `fromPreviousJob` — offer outputs of earlier jobs in this project as inputs.
- `contentFlag` — which representation of a multi-form data object is required
  (e.g. intensities vs amplitudes); `<min>0</min>` means "any".
- `subType` — sub-classification of the data (see per-type `enumerators`).

**Cross-parameter constraints**
- `sameCrystalAs` / `sameCrystalLevel` — tie two data objects to the same
  crystal (used to keep `F_SIGF` and `FREERFLAG` consistent).

---

## How `def.xml` maps to the UI

- Each `content` becomes a widget chosen by its `className` (checkbox, number
  field, text field, dropdown, file browser…).
- `controlParameters` group under a **Parameters/Keywords** area;
  `inputData`/`outputData` under an **Input/Output** area.
- With **no registered interface**, this is rendered automatically by
  `GenericInterface`. A bespoke React interface is only needed for conditional
  visibility, custom layout, or reactive behaviour — see the
  [Task Interface Implementation Guide](../client/renderer/components/task/task-elements/TASK_INTERFACE_IMPLEMENTATION_GUIDE.md).

## Validation

`allowUndefined`, `mustExist`, and `requiredContentFlag` are checked
automatically by the base `validity()`. For task-specific or cross-parameter
rules (and expensive pre-flight checks), override `validity()` /
`runTimeValidity()` in the wrapper — see
[Validity Patterns](../mddocs/pipeline/VALIDITY_PATTERNS.md) and the validation
section of the project `CLAUDE.md`.

---

**Next:** [Authoring a Task](authoring-a-task.md) walks the whole journey —
def.xml → wrapper → registration → UI.
