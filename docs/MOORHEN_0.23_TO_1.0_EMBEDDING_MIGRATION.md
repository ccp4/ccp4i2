# Migrating an embedding app from Moorhen 0.23 → 1.0

Notes for a program that **embeds** Moorhen as a React component library (e.g.
CCP4i2, pandda-inspect) — i.e. you import `MoorhenContainer` and friends into
your own app and manage your own Redux store, rather than running Moorhen's
standalone app.

This is written from the worked CCP4i2 upgrade (0.23.1-alpha.0 → 1.0.0-alpha.3).
The version numbers are illustrative; the structural changes are what matter.

---

## TL;DR — what actually breaks

1. **The package has no root entry any more.** `from "moorhen"` no longer
   resolves. Use `from "moorhen/react-lib"`.
2. **`<MoorhenContainer>` must be wrapped in `<MoorhenInstanceProvider>`, and the
   provider now *requires* a `menuSystem` prop** (`new MoorhenMenuSystem()`).
3. **`MoorhenMenuSystem` must be exported from the `react-lib` entry.** Early
   1.0.0-alpha builds did not export it — see the upstream note at the bottom.

Everything else (asset layout, the named action-creator / component exports, the
`moorhen/types/*` type imports, React 19, MUI 7) is unchanged or backward
compatible.

---

## 1. Package entry point: `main` → `exports` map

**0.23** `package.json`:

```jsonc
{
  "main":  "moorhen.js",
  "types": "moorhen.d.ts"
}
```

So `import { MoorhenContainer } from "moorhen"` resolved to the root bundle, and
arbitrary subpaths like `moorhen/types/moorhen` "just worked" (no `exports`
field means no subpath restriction).

**1.0** `package.json` — no `main`/`module`/`types`, an `exports` map instead:

```jsonc
{
  "exports": {
    "./react-lib":          { "types": "./types/src/moorhen.d.ts",            "default": "./moorhen-react-lib.js" },
    "./web-component":       { "types": "./types/src/WebComponent/entry.d.ts", "default": "./moorhen.js" },
    "./web-component/utils": { "types": "./types/src/WebComponent/utils/entry.d.ts", "default": "./MoorhenWebComponentUtils.js" },
    "./vite-plugin":         { "types": "./scripts/VitePlugin.d.ts",           "default": "./scripts/VitePlugin.mjs" }
  }
}
```

Consequences for an embedder:

- **`from "moorhen"` is now invalid** — there is no `.` key in the map. A modern
  bundler (webpack 5 / Next / Vite / esbuild) will fail to resolve it at build
  time. Switch to the React library entry:

  ```diff
  - import { MoorhenContainer, MoorhenMolecule, addMolecule } from "moorhen";
  + import { MoorhenContainer, MoorhenMolecule, addMolecule } from "moorhen/react-lib";
  ```

  As an embedding app you want **`moorhen/react-lib`** (the React component
  library), *not* `moorhen/web-component` (the `<moorhen-app>` custom element).

- **`moorhen/types/*` imports still work** because those `.d.ts` files are still
  physically present in the tarball (`types/moorhen.d.ts`, `types/mgWebGL.d.ts`,
  …). They are type-only and erased before the bundler sees them, so the
  `exports` map never gates them. Leave them as-is:

  ```ts
  import { moorhen } from "moorhen/types/moorhen";
  import { webGL }   from "moorhen/types/mgWebGL";
  ```

### Module resolution gotcha (the one that wasted the most time)

How `moorhen/react-lib` resolves depends on your `tsconfig` `moduleResolution`:

| `moduleResolution` | `moorhen/react-lib` | `moorhen/types/*` |
|---|---|---|
| `bundler` / `node16` / `nodenext` | ✅ honours `exports` map | ❌ **blocked** (not in `exports`) |
| `node` (classic) | ❌ ignores `exports`; looks for `react-lib.js`/`react-lib.d.ts` that don't exist | ✅ resolves the physical `.d.ts` |

So there is **no single `moduleResolution` that makes both work out of the box.**
Two pragmatic fixes:

- **Stay on classic `node`** (what CCP4i2 / Next 15 do) and add a `paths` alias so
  TypeScript can find the react-lib types. The *bundler* (webpack/Vite) resolves
  the actual JS via the `exports` map at build time; TypeScript only needs the
  declarations:

  ```jsonc
  // tsconfig.json — note path is relative to baseUrl
  {
    "compilerOptions": {
      "moduleResolution": "node",
      "baseUrl": ".",
      "paths": {
        "moorhen/react-lib": ["./node_modules/moorhen/types/src/moorhen.d.ts"]
        // adjust the leading ../ if node_modules is in a parent dir
      }
    }
  }
  ```

  (`moorhen/types/*` keeps resolving against the physical files; no alias needed.)

- **Or switch to `bundler` resolution** and rewrite your `moorhen/types/*`
  imports to pull the same symbols from `moorhen/react-lib` (which re-exports the
  `moorhen`/`webGL` namespaces' types). More churn; only worth it if you're
  modernising resolution anyway.

> Tip: the bundler and TypeScript resolve **independently**. A green
> `tsc --noEmit` does not prove the bundler is happy, and vice-versa. Validate
> both: `tsc --noEmit` *and* a real build (`next build` / `vite build`).

---

## 2. `MoorhenInstanceProvider` now requires a `menuSystem`

1.0 moved Moorhen's per-instance state out of a singleton into a React context
(`MoorhenInstanceContext`). `MoorhenContainer` calls `useMoorhenInstance()`
internally, which **throws** if it isn't rendered inside a
`<MoorhenInstanceProvider>`:

```
useMoorhenInstance must be used within a MoorhenInstanceProvider.
```

And the provider's props changed — `menuSystem` is now **required**:

```ts
interface MoorhenInstanceProviderProps {
  children: ReactNode;
  menuSystem: MoorhenMenuSystem;   // <-- new, required
}
```

The provider builds a `new MoorhenInstance(containerRef, menuSystem)` from it.
`MoorhenMenuSystem` has a no-arg constructor that builds the default menu tree,
so an embedder just instantiates one (once, memoised) and passes it:

```diff
+ import { MoorhenInstanceProvider, MoorhenMenuSystem, MoorhenContainer } from "moorhen/react-lib";
  import { useMemo } from "react";

  function MyMoorhenWrapper(props) {
+   // One menu system per provider; memoise so it isn't rebuilt each render.
+   const menuSystem = useMemo(() => new MoorhenMenuSystem(), []);

    return (
-     <MoorhenInstanceProvider>
+     <MoorhenInstanceProvider menuSystem={menuSystem}>
        <MoorhenContainer {...collectedProps} />
      </MoorhenInstanceProvider>
    );
  }
```

### Why not just use `MoorhenProvider`?

1.0 ships a convenience `MoorhenProvider` that wires up the store *and* a default
`menuSystem` for you. **Avoid it if you manage your own Redux store**, for two
reasons:

1. It calls `createMoorhenStore()` internally, so you'd end up with a second
   store and a double `<Provider>`.
2. (In the alpha at least) it calls `createMoorhenStore()` **un-memoised in the
   component body**, recreating the store on every render.

Embedding apps typically build their own store from Moorhen's reducers:

```ts
import { configureStore } from "@reduxjs/toolkit";
import { MoorhenStoreReducers } from "moorhen/react-lib";

export const store = configureStore({
  reducer: MoorhenStoreReducers,
  middleware: (g) => g({ serializableCheck: false }),
});
```

Keep doing that, wrap your tree in your own `<Provider store={store}>`, and feed
`MoorhenInstanceProvider` a `menuSystem` directly (as above). Customising menus?
Use the `MoorhenMenuSystem` API (`addMainMenu`, `addSubmenu`,
`addToExistingSubmenu`, `deleteItemById`) on the instance before/while passing it
in.

---

## 3. Verify your named exports survived

The 1.0 `react-lib` entry exports essentially the same surface an embedder uses
(action creators like `addMolecule` / `showMap` / `setContourLevel`, components
`MoorhenContainer` / `MoorhenMolecule` / `MoorhenMap`, the `MoorhenPanel` type for
`extraSidePanels`, `MoorhenStoreReducers`, …). But the entry was reorganised, so
**diff your actual imports against the new entry** rather than assume:

```bash
# list everything your app imports from "moorhen"
grep -rhoE 'from "moorhen(/react-lib)?"' -A0 src | ...

# confirm each symbol exists in the installed react-lib declarations
grep -nE "MoorhenContainer|addMolecule|setContourLevel|..." \
  node_modules/moorhen/types/src/moorhen.d.ts
```

If a symbol is missing, it has usually just moved sub-module and needs a
different specifier (or, like `MoorhenMenuSystem`, needs to be added to the entry
upstream — see below).

`extraSidePanels` (the typed "drawer" of custom side panels, `Record<string,
MoorhenPanel>`) is unchanged in 1.0 — `MainContainer` still accepts it and
`MoorhenPanel` is still exported. No change needed there.

---

## 4. Assets (no layout change 0.23 → 1.0, but mind stale copies)

The `public/` tree is the **same shape** in 0.23 and 1.0:

```
public/
  CootWorker.js
  coot_env_web.js
  baby-gru/
  manifest.json, robots.txt, logo*.png, splash*.png
  MoorhenAssets/
    moorhen.css
    monomers/  pixmaps/  tutorials/  mathjax/
    wasm/        <-- the .wasm + emscripten glue live here
```

If you copy these into your app's public dir at build time (CCP4i2 uses
`copyfiles`), there's nothing to change. **But** copy tools generally *overwrite,
never prune*, so anything from an older Moorhen lingers. When jumping versions:

- Delete stale top-level glue files from pre-0.23 Moorhen if present in your dest
  (`moorhen.js`, `moorhen64.js`, `moorhen_st.js` — these moved into
  `MoorhenAssets/wasm/` long ago).
- Clear your bundler cache (`.next/`, Vite cache) so no stale module graph
  survives the jump.
- Re-run your asset copy after `npm install` of the new package.

Moorhen is cross-origin-isolation / SharedArrayBuffer dependent (WASM threads).
If you already serve `Cross-Origin-Opener-Policy: same-origin` +
`Cross-Origin-Embedder-Policy: require-corp` (or `credentialless`) for 0.23, that
requirement is unchanged.

---

## 5. Suggested upgrade sequence

```bash
# 1. install the new package
npm install moorhen@<new>      # or: npm install ./moorhen-<new>.tgz for a local build

# 2. rewrite imports
#    from "moorhen"            -> from "moorhen/react-lib"   (value + `import type`)
#    from "moorhen/types/..."  -> unchanged

# 3. tsconfig: add the paths alias (if on classic `node` resolution)

# 4. add `menuSystem={useMemo(() => new MoorhenMenuSystem(), [])}` to every
#    <MoorhenInstanceProvider>

# 5. clean stale assets + bundler cache, re-copy assets

# 6. validate BOTH:
npx tsc --noEmit                # type resolution
npm run build                   # bundler resolution (the exports map)
#    then run the app and confirm the viewer actually mounts (the InstanceProvider
#    throw is a *runtime* error a build won't catch)
```

---

## Upstream note: export `MoorhenMenuSystem` from `react-lib`

Some 1.0.0-alpha builds export `MoorhenInstanceProvider` (which *requires* a
`menuSystem`) but **not** `MoorhenMenuSystem` itself — leaving an embedder with
no way to construct the required argument (unless they use `MoorhenProvider`,
which isn't viable when you manage your own store, see §2). If you build Moorhen
yourself, ensure the react-lib entry (`src/moorhen.ts`) re-exports it:

```ts
export { MoorhenMenuSystem } from "./components/menu-system/MenuSystem";
```

Then `grep MoorhenMenuSystem dist/types/src/moorhen.d.ts` and
`grep -c MoorhenMenuSystem dist/moorhen-react-lib.js` should both be non-empty in
the packed tarball.
