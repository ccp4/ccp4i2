# Traps, by stage

Each was met for real; most cost an hour the first time.

## Safety

- **Never point anything at a live home.** Scenarios and `devserver.sh`
  refuse `~/.ccp4i2`, `~/.ccp4i2-django`, `~/.ccp4i2x`. Other sessions hold
  ports 3000/3001 and sometimes 3200/3201 on the live database: the help's
  app runs on 3420/3421.
- **Never run an interactive task in a scenario.** Some tasks' "program" is a
  desktop window: qtpisa opens QtPISA, the coot and ccp4mg tasks open theirs,
  PdbView (`pdbview_edit`) its editor, Lidia's "sketch" input a sketcher. Run
  from a scenario they pop up on the developer's screen (qtpisa did,
  2026-09-30).
- **Never send a scenario's data to an outside service.** The deposition
  task's "Use validation server" is on by default and uploads the structure
  to wwPDB: scenarios pass `--SENDTOVALIDATIONSERVER False`. PDB-REDO and
  BUSTER need accounts or licences: document them without running them.

## Scenarios

- **Make the scenario check its own story.** A run can succeed and show
  nothing: SubstituteLigand "fitted" no ligand because `not (NUT or HOH)`
  matched nothing (the language wants `not (NUT) and not (HOH)`). Assert the
  outcome, and measure against the truth where it is known.
- **`task[1]` counts unrun clones; `task[-1]` skips jobs without the file.**
  Name one particular earlier job by its number: `fileOut=[1].HKLOUT[0]`.
- **A file inside a list item** takes `.../fileOut=` (`pdbItemList/structure/
  fileOut=chainsaw[-1].XYZOUT`, since #712). A path re-imports it as a new
  file of no known origin; `dbFileId=` (`scenario_common.output_file_id`)
  works but needs the database.
- **Pass the dictionary the model was refined with.** A dictionary made from
  SMILES names the atoms differently from a model refined with the monomer
  library's: Refmac finds no restraints and PAIREF failed.
- **A task's own lineage logic**, such as the deposition task tracing a model
  back to its refinement and scaling jobs, is called from the scenario
  (`scenario_deposit.py`), not re-implemented as references.
- **Clone the top-level job.** A pipeline can run a sub-job of another task
  with a page; `clone_last` takes top-level jobs only.
- **Deleting a job deletes the jobs that used its outputs.** Tidy a project by
  job number, children last, and check the ids before deleting: a range one
  too long took three new jobs with it.
- **Deleting a project leaves its directory**: files imported again get
  `_1` names. Remove the directory too before rebuilding a project.
- **A script outside `server/` imports the INSTALLED ccp4i2** from
  ccp4-python's site-packages: run checks from `server/`, or with
  `PYTHONPATH`, or a fix will seem not to work.

## The app

- **Django runs `--noreload`**: `devserver.sh restart django` after any server
  change, and `clear-reports` for reports cached (`report_xml.xml`) before it.
- **The Next dev server dies in long sessions**: "Nothing to click labelled
  View" everywhere means it is gone; `devserver.sh start`.
- **Run `capture.mjs` from a shell that has not sourced CCP4**: its setup puts
  Node 20 first on the PATH ("WebSocket is not defined").
- **Developer mode starts on** in every page load; the capture turns it off.
  A picture showing *Def XML* or *Job container* tabs skipped that.
- **"Application error: a client-side exception"** on a job page was the
  resizable panels' generated ids differing between server and client; the
  layouts now give them fixed ids. If it recurs, `<out>.failed.txt` has it.
- **Viewing a pending job can change it.** A component that writes on mount
  rewrites the job (the selection builder erased clones' selections until
  #702): diff the clone's input_params.xml when a value looks missing.
- **An unrun clone shows what autofill does, and fails to do.** Values filled
  only when a file is picked (cell, wavelength) stay empty on i2run jobs and
  clones: a bug in the page, not the scenario.

## Capture

- **Outline first, picture last.** A full-page screenshot read to learn a
  page's labels costs many times an outline.
- **A shot's section must be a heading the page shows**; a select's value or
  a canvas title is not text (`"Conformational landscape"` in a dropdown).
  Callouts on a plot target its "Plot" selector label.
- **`text` matches the first visible text that starts with it**: "Residue"
  hit "Residues whose…" before the table header. Use a longer prefix or a
  column header with `"closest": "table"`.
- **Collapsed folds** photograph as headings: list them under `"expand"`.
  Unlabelled controls are pressed by CSS selector under `"click"`
  (`button[value="10"]` for the PHIL expert level *All*).
- **Badges sit in the left margin**: two on one row overlap; badge the left
  column and say "beside it".
- **A section that appears only for some choices** cannot be a callout on a
  default job: describe it in the text.
- **Heavy reports need longer to settle** (`"settle": 30000`).
- **A label from the program's own definitions** (a PHIL scope's caption) is
  in no source file: mark the callout `"dynamic": true` and the shot
  `"dynamic_section": true`, or the stamp check calls the page stale.

## Reports

- **A report's text is parsed as XML.** HTML entities (`&Aring;`) are
  undefined there and fail the whole report: write the characters. Escape
  anything taken from a log (`27 < 50`). Embed SVG through
  `ccp4i2.report.svg.inline_svg` (an XML declaration mid-document fails it)
  and `fit_svg` (a fixed-size drawing overflows its box).
- **Build every new or changed report in a unit test**; it catches the above
  before a capture does (`test_pairef_results.py`, `test_areaimol_areas.py`).
- **Graph titles cannot contain spaces** (the loggraph header is split on
  whitespace): use underscores.
- **A score a program did not calculate may be stored as 0.** Phaser keeps
  every score on a solution; the report shows "–" unless the placements say
  it was calculated. Check what a zero means before quoting it.

## Stamps and the build

- **A stamp sees the task's own files**, not shared widgets or a pipeline's
  sub-wrappers (on purpose: a widget change would flag every page). Moving
  code out of a stamped file marks the page stale: read `stamp.py diff`,
  then restamp.
- **`build.sh` always builds from scratch** (`-E`): an incremental build
  undercounts warnings and lowers the ratchet wrongly.
- **JPEG is wrong for screenshots**: 256-colour PNG (`compress_images.py`).
