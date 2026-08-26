# Coordinate format fidelity

How CCP4i2 decides whether a coordinate file is PDB or mmCIF, what it writes
when a wrapper asks for "the selected atoms", and why those two questions are
currently answered by the file's *name* rather than its *content*.

**Status:** audit complete, nothing implemented. Written 2026-08-25 from the
diagnosis of `test_substitute_ligand`, one of the ten failures in the
[pre-C1 i2run baseline](error-handling-remediation.md#the-pre-c1-baseline).

**Do not apply this before the C1 work has been measured.** The pre-C1 baseline
is the comparator for C1's impact, and every change proposed here moves i2run
results too. Landing them together would make both unreadable. Sequencing is set
out under [When to apply this](#when-to-apply-this).

---

## The finding in one paragraph

`CPdbDataFile` offers two methods that read as a matched pair —
`getSelectedAtomsPdbFile` and `getSelectedAtomsFile` — and the names promise
exactly what a caller would want: one gives you a PDB whatever you started with,
the other preserves the input format. Neither promise is kept. The
"format-preserving" one is a thin wrapper that computes a filename and calls the
"PDB" one; the "PDB" one decides what to write by looking at the extension of
the path it was handed, and skips writing altogether — copying the input
verbatim — when no atom selection is set. So the format of the output is decided
by a filename, and 28 call sites have each worked out their own way of coping.

---

## What the two methods actually do

`CPdbDataFile.getSelectedAtomsPdbFile(fileName)` —
[CCP4ModelData.py:2095](../server/ccp4i2/core/CCP4ModelData.py):

| Selection set? | Output name | What is written |
|---|---|---|
| no | anything | **verbatim copy of the input** — the content keeps the *input* format |
| yes | `*.cif`, `*.mmcif` | mmCIF |
| yes | anything else | PDB |

The format branch lives in `CPdbData.writeSelection`
([CCP4ModelData.py:1742](../server/ccp4i2/core/CCP4ModelData.py)), which ends:

```python
suffix = Path(file_path).suffix.lower()
if suffix in ['.cif', '.mmcif']:
    output_structure.make_mmcif_document().write_file(str(file_path))
else:
    output_structure.write_pdb(str(file_path))
```

`CPdbDataFile.getSelectedAtomsFile(baseName, workDirectory)` —
[CCP4ModelData.py:2151](../server/ccp4i2/core/CCP4ModelData.py) — picks
`.cif` or `.pdb` from `isMMCIF()`, then calls `getSelectedAtomsPdbFile` with
that name. It is format-preserving *because* the other method is not.

So the real contract of the pair is:

- `getSelectedAtomsPdbFile` — "write in whatever format this filename implies,
  unless there is no selection, in which case copy the bytes".
- `getSelectedAtomsFile` — "compute a filename that implies the input format,
  then call the above".

The intended contract — and the one the names, the docstrings and every caller
assume — is:

- `getSelectedAtomsPdbFile` — **a PDB file, always**, for programs that read
  nothing else.
- `getSelectedAtomsFile` — **the input's format, always**, for programs that
  read both and lose information when down-converted.

## Why nobody noticed

Only **three** i2run tests set an atom selection at all
(`test_coordinate_selector`, `test_i2run`, `test_phaser_rnp_pipeline`). Every
other test of these 28 call sites exercises the no-selection row of the table
above — the row that copies bytes and therefore cannot disagree with itself.
The selection paths are close to untested.

---

## How format is inferred today

Every read in `CCP4ModelData.py` (lines 1621, 1960, 2008, 2284) and in the
wrappers calls `gemmi.read_structure(path)` with the default arguments. The
default is **not** what the surrounding comments assume ("gemmi auto-detects the
input format", "handles PDB and mmCIF automatically"):

```
read_structure(path, merge_chain_parts=True, format=CoorFormat.Unknown, ...)
```

`CoorFormat.Unknown` means *detect from the extension*. `CoorFormat.Detect`,
passed explicitly, means *detect from the content*. Measured on gemmi 0.7.5,
the version in `ccp4-20260702`:

| file | actual content | `read_structure(p)` | `read_structure(p, format=Detect)` |
|---|---|---|---|
| `mis.pdb` | mmCIF | **FAIL** — "perhaps it is cif not pdb?" | OK → `Mmcif` |
| `mis.cif` | PDB | **FAIL** — ValueError | OK → `Pdb` |
| `x.tmp` | PDB | **FAIL** — "Unknown format" | OK → `Pdb` |
| `x.xyz` | PDB | **FAIL** — "Unknown format" | OK → `Pdb` |
| `a.ent` | PDB | OK → `Pdb` | OK → `Pdb` |
| `a.mmcif` | mmCIF | OK → `Mmcif` | OK → `Mmcif` |

Content-based inference — the thing this whole area needs — is one keyword
argument. Exactly one call site in the codebase already passes it:
`density_calculator` ([density_calculator.py:21](../server/ccp4i2/wrappers/density_calculator/script/density_calculator.py)),
which is why that wrapper can write its input to `xyzin.xyz` and still work.

`CPdbDataFile.isMMCIF()` ([CCP4ModelData.py:2035](../server/ccp4i2/core/CCP4ModelData.py))
does sniff content, but by hand: it reads the first two lines looking for
`data_` or `loop_`, and gives up on anything else — an mmCIF file that opens
with comments or blank lines is misclassified. It should ask gemmi.

---

## Consumer audit — 28 call sites

What each caller does, and what its program actually needs. **Verify the
"needs" column against the installed program before changing that caller** —
several wrappers carry a comment about a program's limitations that was true
when it was written and is not true now (see [dimple](#worked-example-dimple)).

### Group A — needs a real PDB

| Caller | Output name | Called when |
|---|---|---|
| [chainsaw.py:16](../server/ccp4i2/wrappers/chainsaw/script/chainsaw.py) | `XYZIN_sel.pdb` | selection only |
| [acorn.py:259](../server/ccp4i2/wrappers/acorn/script/acorn.py) | `XYZIN_selected_atoms.pdb` | selection only |
| [molrep_mr.py:53](../server/ccp4i2/wrappers/molrep_mr/script/molrep_mr.py) | `XYZIN_selected_atoms.pdb` | selection only |
| [molrep_den.py:37](../server/ccp4i2/wrappers/molrep_den/script/molrep_den.py) | `XYZIN_selected_atoms.pdb` | selection only |
| [mrbump_basic.py:172](../server/ccp4i2/wrappers/mrbump_basic/script/mrbump_basic.py) | `selected_N.pdb` | selection only |
| [phaser_ensembler.py:19](../server/ccp4i2/wrappers/phaser_ensembler/script/phaser_ensembler.py) | `selected_N.pdb` | selection only |
| [phaser_MR.py:167](../server/ccp4i2/pipelines/phaser_pipeline/wrappers/phaser_MR/script/phaser_MR.py) | `selectedAtomModel_i_j.pdb` | selection only |
| [phasertng_picard.py:93](../server/ccp4i2/wrappers/phasertng_picard/script/phasertng_picard.py) | `<tag>.pdb` | selection only |
| [phasertng_riker.py:99](../server/ccp4i2/wrappers/phasertng_riker/script/phasertng_riker.py) | `<tag>.pdb` | selection only |
| [zanuda.py:24](../server/ccp4i2/wrappers/zanuda/script/zanuda.py) | `model_selected.tmp` | selection only |
| [findmyseq.py:34](../server/ccp4i2/wrappers/findmyseq/script/findmyseq.py) | `XYZIN_selected_atoms.pdb` | selection only |
| [SubstituteLigand.py:138](../server/ccp4i2/pipelines/SubstituteLigand/script/SubstituteLigand.py) | `selected_atoms.pdb` | **always** |

Two things are wrong across this group.

**The no-selection branch bypasses the API entirely.** `chainsaw`, `acorn`,
`findmyseq`, `phaser_ensembler` and `mrbump_basic` all pass
`XYZIN.fullPath` *unmodified* to the program when no selection is set. An mmCIF
input therefore reaches a PDB-only program as mmCIF. `molrep_mr`/`molrep_den`
guard this (`elif getExt() != '.pdb': convertFormat('pdb', ...)`), and the
`phasertng` pair guard it with an extension test — so four wrappers handle the
case, five do not.

**`zanuda` is broken whenever a selection is set.** It writes to
`model_selected.tmp`, with the comment "writes in source format" — wrong twice:
`writeSelection` sees `.tmp`, so it writes **PDB**; then
`gemmi.read_structure(inputPath)` sees `.tmp` and raises **"Unknown format"**.
No test sets a selection, so it has never been seen.

### Group B — wants the input's format preserved

| Caller | How it gets there |
|---|---|
| [refmac.py:124](../server/ccp4i2/wrappers/refmac/script/refmac.py) | hand-rolled suffix swap |
| [servalcat.py:78](../server/ccp4i2/wrappers/servalcat/script/servalcat.py) | hand-rolled suffix swap |
| [sheetbend.py:33](../server/ccp4i2/wrappers/sheetbend/script/sheetbend.py) | hand-rolled suffix swap |
| [csymmatch.py:25,36](../server/ccp4i2/wrappers/csymmatch/script/csymmatch.py) | hand-rolled suffix swap (×2) |
| [coordinate_selector.py:20](../server/ccp4i2/wrappers/coordinate_selector/script/coordinate_selector.py) | hand-rolled suffix swap |
| [gesamt.py:41,51](../server/ccp4i2/wrappers/gesamt/script/gesamt.py) | `getSelectedAtomsFile` ✔ |
| [modelcraft.py:32](../server/ccp4i2/wrappers/modelcraft/script/modelcraft.py) | `getSelectedAtomsFile` ✔ |
| [lorestr_i2.py:42](../server/ccp4i2/wrappers/lorestr_i2/script/lorestr_i2.py) | `getSelectedAtomsFile` ✔ |

Five of the eight write the same four lines by hand:

```python
if self.container.inputData.XYZIN.isMMCIF():
    self.inputCoordPath = str(pathlib.Path(self.inputCoordPath).with_suffix('.cif'))
self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.inputCoordPath)
```

That is `getSelectedAtomsFile` spelled out longhand, and it only works because
`getSelectedAtomsPdbFile` is not a PDB method. Any fix that makes the PDB method
honest **must** convert these five in the same commit, or they will start
writing PDB content into `.cif` files — the same bug, mirrored.

### Group C — the consumer reads whatever it is given

| Caller | Consumer |
|---|---|
| [density_calculator.py:20](../server/ccp4i2/wrappers/density_calculator/script/density_calculator.py) | gemmi, with `format=Detect` — the only content-based read in the tree |
| [add_fractional_coords.py:18](../server/ccp4i2/wrappers/add_fractional_coords/script/add_fractional_coords.py) | its own `script.py`, which tries `gemmi.cif.read` and falls back to PDB |
| [i2Dimple.py:42,45](../server/ccp4i2/wrappers/i2Dimple/script/i2Dimple.py) | dimple ≥ 2.7 — see below |

## Worked example: dimple

`i2Dimple.processInputFiles` casts mmCIF to PDB by hand, explaining:

> dimple needs PDB input. If XYZIN is mmCIF, convert it to PDB first … gemmi
> auto-detects the input format, replacing the previous mmdb2
> ReadCIFASCII/WritePDBASCII cast

Both sentences are now false, and between them they produced the failure that
started this audit. Checked against `ccp4-20260702`, dimple 2.7:

- its CLI whitelist accepts `.cif`, `.cif.gz`, `.mmcif`, `.mmcif.gz`
  (`dimple/main.py:763`);
- `wf.copy_uncompressed(opt.pdbs[0], 'ini.pdb')` (`main.py:71`) is not a copy —
  `dimple/workflow.py:954` is `st = gemmi.read_structure(src); st.write_pdb(dst)`.

dimple reads the model itself and writes its own `ini.pdb`; everything after
that point is PDB-internal, but the conversion happens **inside dimple, once**.
The wrapper's cast is redundant. It is also harmful: it writes mmCIF content to
`selected_atoms.pdb` (via the copy branch of `getSelectedAtomsPdbFile`), and
dimple's own `gemmi.read_structure` — extension-based, like everyone's — then
fails with *"Incorrect file format (perhaps it is cif not pdb?)"*.

The correct wrapper is one call to `getSelectedAtomsFile`, no gemmi, no cast:
i2Dimple belongs in Group B. The one condition is that the file we hand dimple
is **named honestly**, which is precisely the invariant this document
establishes.

The general lesson: check the program, not the comment. Several Group A entries
may have been re-plumbed upstream the same way, while the wrapper still carries
a note written for a version that is years gone.

---

## What to change

1. **`detect_coordinate_format(path)`** — a module-level helper in
   `CCP4ModelData` returning `'mmcif'` / `'pdb'` / `None`, decided by content
   via `gemmi.read_structure(path, format=gemmi.CoorFormat.Detect)`, with the
   extension used only when the content is unreadable. `isMMCIF()` and the
   `contentFlag` inference both delegate to it.
2. **Every `gemmi.read_structure` gains `format=gemmi.CoorFormat.Detect`** — in
   `CCP4ModelData` and in the wrappers. A file is then trusted for what it
   contains, not what it is called. This alone removes the entire class of
   mis-named-file read failures, including i2Dimple's.
3. **`writeSelection(atoms, path, fmt)` takes an explicit format.** No function
   in this area may infer an output format from an output filename.
4. **`getSelectedAtomsPdbFile` guarantees PDB**, on both paths — the
   no-selection path converts instead of copying — and fails loudly when the
   structure cannot be represented in PDB (>62 chains, chain IDs longer than
   two characters, atom-serial or residue-number overflow) rather than writing
   something truncated.
5. **`getSelectedAtomsFile` guarantees the source format** with a matching
   extension, on both paths.
6. **Callers.** The five Group B hand-rolled suffix swaps become
   `getSelectedAtomsFile` (**required** — see above). The five Group A wrappers
   that pass the original file through when no selection is set call the
   PDB-guaranteeing method **unconditionally**. `i2Dimple` moves to Group B and
   loses its cast. `zanuda` stops writing `.tmp`.
7. **Report failures.** These methods return an int that almost every caller
   discards; `phasertng_picard`/`riker` are the only ones that check. Whatever
   lands should raise or report through the plugin error report — the same
   argument as C1 in [error-handling-remediation.md](error-handling-remediation.md).

## Test plan

A CCP4-free unit matrix (gemmi only), asserting the **content** of every output
by re-reading it with an explicit `format=`, never its name:

- inputs: PDB, mmCIF, **PDB named `.cif`**, **mmCIF named `.pdb`**, and a file
  with an extension gemmi does not know (`.tmp`);
- × selection set / not set;
- × `getSelectedAtomsPdbFile` / `getSelectedAtomsFile`;
- assert: the PDB method always yields PDB content; the format-preserving method
  always yields the input's format and an extension that matches it; a
  mis-named input is read correctly in every case.

Plus: `detect_coordinate_format` against the six-row gemmi table above, and a
regression test for each of the specific defects named here (the SubstituteLigand
copy, the zanuda `.tmp`, the five unguarded Group A no-selection paths).

The i2run evidence should be a targeted run of the affected wrappers —
`substitute_ligand`, `dimple`, `chainsaw`, `csymmatch`, `gesamt`, `modelcraft`,
`refmac`, `servalcat`, `sheetbend`, `zanuda`, `coordinate_selector` — **with a
selection set**, since that is the untested half and the whole point.

## What was implemented, and what was not

**Landed 2026-08-25** on `coordinate-format-fidelity`, against the post-C1
baseline.

Done:

1. `detect_coordinate_format()` — content first (a cheap scan of the leading
   lines), gemmi's `CoorFormat.Detect` second, the extension last. `isMMCIF` and
   the `contentFlag` introspection both go through it.
2. Every `gemmi.read_structure` in `CCP4ModelData` passes
   `format=gemmi.CoorFormat.Detect`.
3. `writeSelection` takes the format as an argument.
4. `getSelectedAtomsPdbFile` guarantees PDB on both paths; the no-selection path
   converts instead of copying. Where PDB cannot express the structure, gemmi
   raises and the error is returned.
5. `getSelectedAtomsFile` guarantees the source format on both paths, with a
   matching name.
6. Callers: the five hand-rolled suffix swaps (`refmac`, `servalcat`,
   `sheetbend`, `csymmatch`, `coordinate_selector`) now call
   `getSelectedAtomsFile`; `i2Dimple` and `SubstituteLigand` do the same and
   drop the cast; `zanuda` makes one call and no longer writes `.tmp`.

**Not done: the five Group A wrappers that pass the input through unconverted
when no selection is set** — `chainsaw`, `acorn`, `findmyseq`,
`phaser_ensembler`, `mrbump_basic`. Routing those through the PDB-guaranteeing
call would fix an mmCIF reaching a PDB-only program, but it is a behaviour
change on the most-used path, and it rests on the assumption that each program
needs PDB. That assumption is exactly what proved false for dimple. Changing
five wrappers on the strength of a comment and a filename would repeat today's
mistake in the opposite direction.

Each needs one experiment before it moves: run the program on an mmCIF model and
see whether it reads it. `chainsaw` cannot be tested trivially — it exits on a
missing alignment file before it looks at coordinates — so the experiment is an
i2run run with a `.cif` model rather than a bare invocation. `molrep` already
guards itself (`elif getExt() != '.pdb': convertFormat(...)`), and the
`phasertng` pair guard with an extension test that becomes redundant once
detection is content-based.

Until then the position is unchanged from before this work, and no worse: an
mmCIF with no selection reaches those five programs as mmCIF, exactly as it did.

## When to apply this

Not yet, and not alongside C1.

The pre-C1 i2run baseline is the comparator for C1's impact. These changes alter
what several wrappers hand their programs, so they will move i2run results on
their own account. Landing them into the same interval would leave neither
change measurable.

The order is: finish the C1 triage against the pre-C1 baseline; merge C1; take
the post-C1 run as the new baseline; then apply this work and compare against
*that*. `test_substitute_ligand` ×2 are expected to go green, and any other
movement is this work's to explain.

## What it says about the wider picture

Three of the failure modes here are the same shape as the ones in
[error-handling-remediation.md](error-handling-remediation.md#a-second-axis--nothing-to-discard-in-the-first-place):
a value inferred from something incidental (a filename) rather than from the
thing itself (the bytes), a fallback path that silently does something different
from the main path (copy versus write), and a return code nobody reads. The
audit found four defects — the SubstituteLigand format lie, zanuda's `.tmp`,
five unguarded mmCIF-to-PDB-only-program paths, and a hand-sniffed `isMMCIF`
that misreads a commented mmCIF — of which exactly one had ever been observed,
because the only path anyone tests is the one where the question does not arise.
