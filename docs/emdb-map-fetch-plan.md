# Fetching maps from EMDB, aware of the map subtype

A plan for later attention. Written 2026-09-17 against django `f94199b9e`.
**Status:** implemented 2026-09-17 on branch `emdb-map-fetch`: parts 1 to 3
below (the `emdb` download mode, the server-side fetch, half-map pairs);
part 4 (fitted-model fetch, metadata field) is not done.

## What exists

The repository download modes (`ebiPdb`, `rcsbPdb`, `uniprotAFPdb`, `ebiSFs`,
`uniprotFasta`, PDB-REDO) are client-driven. A `downloadModes` qualifier on
the CData class (`core/CCP4ModelData.py`, `core/CCP4XtalData.py`) lists the
repositories a parameter accepts; the file selector's fetch dialog
(`client/renderer/components/task/task-elements/fetch-file-for-param.tsx`)
resolves the identifier through the repository's API, downloads the file in
the browser, and uploads it through `upload_file_param`. The "back end" in the
path is the Next proxy route (`app/api/proxy/pdbe/[...path]/route.ts`,
`.../uniprot/...`), not Django: it exists because the Moorhen pages carry
COEP and PDBe sends no `Cross-Origin-Resource-Policy` header, so a direct
browser fetch fails. Every byte therefore makes a round trip through the
browser, which is fine for a 2 MB mmCIF.

`CMapDataFile` subtypes (`core/CCP4XtalData.py`): 1 normal, 2 difference,
3 anomalous difference, 4 mask, 5 half map. `ImportMap` has a `MAP_SUBTYPE`
selector over the same list; `servalcat` takes `MAPIN1`/`MAPIN2` with
`requiredSubType` 5 and `MAPMASK` with 4.

## What EMDB exposes

The entry API, `https://www.ebi.ac.uk/emdb/api/entry/EMD-NNNNN`, names every
file, and the archive at
`https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-NNNNN/` lays them
out by kind. The kinds map directly onto our subtypes:

| Our subtype | EMDB JSON | Archive path |
|---|---|---|
| 1 normal | `map.file` | `map/emd_NNNNN.map.gz` |
| 5 half map | `interpretation.half_map_list.half_map[i].file` | `other/emd_NNNNN_half_map_{1,2}.map.gz` |
| 4 mask | `interpretation.segmentation_list.segmentation[i].file` | `masks/emd_NNNNN_msk_1.map` |

Availability varies per entry (checked 2026-09-17): EMD-11638 has all three
kinds, EMD-8117 has half maps but no mask, EMD-30210 has only the main map.
The same JSON gives pixel spacing, dimensions, origin, data type, the
author's contour level, the reported resolution, and cross-referenced PDB
ids (`crossreferences.pdb_list`).

Sizes: a 256³ float half map is 67 MB gzipped and 268 MB raw; main maps of
several hundred MB are routine.

## Why the bytes must not go through the browser

- The upload leg goes through the Next proxy, whose middleware body cap on
  a served deployment is 100 MB; staged import (#538, #540) exists to work
  around exactly that, and would make every map fetch a two-hop copy.
- The proxy routes buffer the whole upstream body in memory
  (`arrayBuffer()`), so a 300 MB map is held in the Next process and again
  in the browser.
- The upload path does no gunzip at all (`upload_param.py` has none), and
  EMDB serves maps gzipped; downstream programs want them uncompressed.

So for maps the fetch belongs on the server: download from the archive
straight into the project, gunzip on the way, set the subtype, set the
parameter. That is what "via the back end" should mean here.

## The plan

### 1. Client: an `emdb` download mode

- Add `"downloadModes": ["emdb"]` to `CMapDataFile.Meta.qualifiers` so every
  map parameter's selector offers it.
- Add `/api/proxy/emdb/[...path]/route.ts` mirroring the PDBe route, for the
  entry JSON only (small; the bytes do not come this way).
- In the fetch dialog, an EMDB branch that takes `EMD-NNNNN` or `NNNNN`,
  queries the entry, and lists what the entry actually has, grouped as main
  map, half maps, masks, with sizes and the resolution. Preselect by the
  parameter's `requiredSubType`, or by `MAP_SUBTYPE` on ImportMap.
- The dialog is one long component with a branch per mode; add the EMDB
  branch as its own module (`fetch-modes/emdb.tsx`) and let the others
  follow when next touched.

### 2. Server: fetch into the project

A new job action, `POST jobs/{id}/fetch_repository_file/`, body
`{object_path, repository: "emdb", entry: "EMD-11638", file: "emd_11638_half_map_1.map.gz", sub_type: 5}`:

- Builds the archive URL from the kind (`map/`, `other/`, `masks/`), never
  from a client-supplied URL.
- Streams the download to `CCP4_IMPORTED_FILES` (gunzipping when the name
  ends in `.gz`), records the source checksum and a `FileImport` row with the
  EMDB URL as provenance, exactly as `upload_file_param` does for a body.
- Sets `sub_type` from the request (validated against the parameter's
  `requiredSubType`) and an annotation from the JSON:
  `EMD-11638 half map 1, 0.53 Å/px, 256³, 3.2 Å`.
- Sets the parameter through the same code path as the import-by-path
  branch of `upload_file_param`, so the container, `input_params.xml` and the
  file rows are handled identically.
- Returns the same shape as `upload_file_param`.

The endpoint is repository-agnostic by design: a later `pdbe` repository
could carry structure-factor files the same way and retire the browser round
trip for the large cases.

On Azure the server container needs outbound HTTPS to `ftp.ebi.ac.uk` and
`www.ebi.ac.uk`; confirm before relying on it there.

### 3. Half-map pairs

Half maps only make sense in twos, and our tasks take them as two
parameters. When the dialog is opened on a half-map parameter whose sibling
exists (`MAPIN1` beside `MAPIN2`), offer "half map 1 here, half map 2 into
the sibling" and issue two fetches. Otherwise let the user pick which half.
`ImportMap` imports one file per job; fetching a pair there is two jobs, and
the dialog should say so rather than pretend.

### 4. Extras worth taking while there

- The cross-referenced PDB ids in the entry JSON: offer "also fetch the
  fitted model" through the existing `ebiPdb` path into the sibling
  coordinate parameter when there is one.
- The author contour level and pixel spacing are useful to the Moorhen
  viewer and the molrep_map pipeline; keep them in the annotation now and
  consider a file-metadata field later.

## Tests

- Unit: archive URL construction per kind; gunzip-on-stream; annotation
  from a captured entry JSON (fixture files for EMD-11638, EMD-8117,
  EMD-30210 so the three availability cases are covered without the network).
- API unit: the endpoint refuses an unknown repository, a `sub_type` that
  contradicts `requiredSubType`, and a file name the entry does not list;
  sets the parameter and creates the `FileImport` row (mock the download).
- Client: the dialog lists the right groups for each fixture and preselects
  by `requiredSubType`.

## Order

Server endpoint first (testable without the client, useful to i2run at
once), then the dialog, then the pair handling.
