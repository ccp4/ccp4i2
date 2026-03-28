# Licensing Checklist

CCP4i2 is licensed under the **GNU Lesser General Public License version 3**, modified in accordance with the provisions of the license to address the requirements of UK law. See https://www.ccp4.ac.uk/ccp4license.php for the full text.

## Copyright Holders

Multiple institutions hold copyright over different parts of the codebase. Under UK law (CDPA 1988, s.11(2)), copyright for code written in the course of employment belongs to the employer.

| Area | Primary Copyright Holder(s) |
|------|-----------------------------|
| `client/` (frontend) | Newcastle University |
| `apps/` (compounds) | Newcastle University |
| `server/ccp4i2/` | Newcastle University, University of York |
| `pimple/`, `smartie/`, `report/` | University of York, Newcastle University |
| `wrappers/`, `wrappers2/`, `pipelines/` | Multiple (see below) |

The wrappers/pipelines directories have contributions from Newcastle University, University of York/STFC, Global Phasing Ltd, MRC Laboratory of Molecular Biology, and others. These files do not yet carry copyright headers and will need separate review.

## Checklist for New Code

- [ ] **File headers**: Every new `.py`, `.ts`, and `.tsx` file must include a copyright and license header at the top. See [Header Templates](#header-templates) below.
- [ ] **Copyright holder**: Use your employer's institution name, not your personal name. If unsure, check with your institution's IP/contracts office.
- [ ] **Year**: Use the year(s) of creation/modification. A range (e.g., 2025-2026) is appropriate for files modified over multiple years.

## Checklist for New Dependencies

Before adding any dependency:

- [ ] **Check the license** of the package (look at its `package.json`, `setup.cfg`, or `LICENSE` file).
- [ ] **Verify compatibility** with LGPL v3. The following are compatible:
  - MIT, BSD-2-Clause, BSD-3-Clause, ISC, Apache-2.0, Zlib, Unlicense, CC0
- [ ] **Flag for review** if the dependency uses:
  - **GPL** (any version) -- imposes copyleft on the combined work
  - **AGPL** -- imposes network copyleft (triggered by server-side use)
  - **Non-OSI-approved** or proprietary licenses
  - **No license declared** -- legally means "all rights reserved"
- [ ] **Document the dependency** with its license in the relevant `package.json` or `pyproject.toml`.
- [ ] **Avoid packages with known licensing controversy** (e.g., SheetJS/xlsx was replaced with ExcelJS for this reason).

### Quick Reference: Current Dependencies

All current Python dependencies are MIT or BSD-3-Clause. All current JavaScript dependencies are MIT, BSD-3-Clause, or Apache-2.0.

Notable packages and their licenses:

| Package | License | Notes |
|---------|---------|-------|
| Django | BSD-3-Clause | |
| React / Next.js | MIT | |
| ExcelJS | MIT | Replaced xlsx (SheetJS) due to licensing concerns |
| MUI (@mui/*) | MIT | |
| Electron | MIT | |
| Moorhen | BSD-3-Clause | Bundles Coot (GPL-2.0+) as WASM -- see note below |
| @rdkit/rdkit | BSD-3-Clause | |
| Biopython | Biopython License (BSD-like, OSI-approved) | |

**Moorhen note**: Moorhen bundles compiled WASM from Coot (GPL-2.0+) and Clipper (LGPL-2.1+). As a CCP4-funded project this is not a practical concern, but the interaction between Coot's GPL-2.0+ and our LGPL-3.0 should be formally documented.

## Checklist for Releases

- [ ] **Root LICENSE file** is present and up to date.
- [ ] **`client/package.json`** license field matches root LICENSE (`LGPL-3.0-or-later`).
- [ ] **No files without headers** in the main source directories (run the header check below).
- [ ] **No new dependencies** added without license review.
- [ ] **wrappers/pipelines provenance** has been reviewed (when headers are added to those directories).
- [ ] **Modified LGPL v3 text** is accessible at https://www.ccp4.ac.uk/ccp4license.php.

## Header Templates

### Python (.py)

```python
# Copyright (C) YYYY Institution Name
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
```

Place after any shebang (`#!/...`) or encoding declaration, before imports.

### TypeScript / TSX (.ts, .tsx)

```typescript
/*
 * Copyright (C) YYYY Institution Name
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
```

Place at the very top of the file, before `'use client'` / `'use server'` directives.

### Multiple Copyright Holders

If a file has significant contributions from multiple institutions:

```python
# Copyright (C) 2025-2026 Newcastle University
# Copyright (C) 2025 University of York
```

List the institution with the most contributions first.

## Checking for Missing Headers

To find source files without copyright headers:

```bash
# Python files
grep -rL "Copyright" --include="*.py" server/ccp4i2/ core/ cli/ apps/ pimple/ smartie/ report/

# TypeScript/TSX files
grep -rL "Copyright" --include="*.ts" --include="*.tsx" client/renderer/ client/main/ apps/
```

## Institutional Sign-Off

Before distributing CCP4i2, each contributing institution should formally agree to the LGPL v3 (modified) license for their contributions. Required sign-offs:

- [ ] **Newcastle University** -- IP/contracts office (covers Martin Noble's contributions)
- [ ] **University of York / STFC** -- (covers Paul Bond, Stuart McNicholas, and other York contributors)
- [ ] **Other contributors** -- as identified in wrappers/pipelines review

A Contributor License Agreement (CLA) or Developer Certificate of Origin (DCO) should be considered for future contributors.
