#!/usr/bin/env node
/**
 * Precompress the large WASM/JS assets the Moorhen API route serves.
 *
 * Next compresses what it serves from public/, but not the responses of
 * app/api/moorhen/[...path]/route.ts, so the single largest file on the page
 * went out raw: moorhen64.wasm is 17.8 MB, and brotli takes it to 4.4 MB. On a
 * poor link that is the difference between a viewer that starts and one that
 * does not -- a cold load over mobile measured 50 s for a 1.3 MB archive, at
 * which rate the uncompressed WASM needs five to six minutes.
 *
 * Why build time, not request time: brotli quality 10 takes ~28 s for
 * moorhen.wasm (and quality 11, ~41 s, for about 1% fewer bytes -- not worth
 * it). That cost belongs in the build, paid once, not in a request handler.
 *
 * Why this is cheap in practice: the output is keyed by a hash of the input, so
 * a rebuild with an unchanged Moorhen reuses what is already there and the
 * whole step takes milliseconds. It is paid once per Moorhen bump, which is
 * infrequent, not once per CCP4i2 build.
 *
 * Only the web build needs this (Electron serves these files from disk with no
 * HTTP layer in between), so copy-assets runs it for BUILD_TARGET=web.
 */
const { createHash } = require("crypto");
const {
  existsSync,
  readFileSync,
  writeFileSync,
  mkdirSync,
  statSync,
  readdirSync,
  unlinkSync,
} = require("fs");
const { join, dirname } = require("path");
const zlib = require("zlib");

const PUBLIC_DIR = join(__dirname, "..", "renderer", "public");

// Brotli quality. 10 rather than 11 deliberately: on moorhen.wasm, q11 spends
// 41 s to q10's 22 s and saves about 1% more of the original.
const BROTLI_QUALITY = 10;
const GZIP_LEVEL = 9;

// Don't bother below this: small files are a rounding error on the wire and
// the build time stops being worth it.
const MIN_BYTES = 64 * 1024;

/**
 * Files to compress, relative to public/.
 *
 * Only what the Moorhen route actually serves. The root-level moorhen*.wasm
 * copies are NOT here: the worker loads `${urlPrefix}/wasm/CootWorker.js` and
 * takes everything from MoorhenAssets/wasm/, while the root copies are stale
 * leftovers from an older layout (copy-moorhen-root explicitly excludes
 * *.wasm). Compressing them would add ~50 s for files nothing fetches.
 */
const TARGET_DIRS = ["MoorhenAssets/wasm"];
const TARGET_FILES = ["RDKit_minimal.wasm", "RDKit_minimal.js"];

const COMPRESSIBLE = /\.(wasm|js)$/;

function collect() {
  const out = [];
  for (const dir of TARGET_DIRS) {
    const full = join(PUBLIC_DIR, dir);
    if (!existsSync(full)) continue;
    for (const entry of readdirSync(full, { withFileTypes: true })) {
      if (entry.isFile() && COMPRESSIBLE.test(entry.name)) {
        out.push(join(dir, entry.name));
      }
    }
  }
  for (const file of TARGET_FILES) {
    if (existsSync(join(PUBLIC_DIR, file))) out.push(file);
  }
  return out;
}

function hashOf(buffer) {
  return createHash("sha256").update(buffer).digest("hex").slice(0, 16);
}

/**
 * Write `<file>.br` / `<file>.gz` beside the file, plus a `.hash` stamp.
 *
 * The stamp records which input the siblings were made from, so an unchanged
 * asset is skipped and a changed one is recompressed. Without it a rebuild
 * would either redo a minute of work every time or, worse, leave the previous
 * version's compressed bytes next to this version's file.
 */
function compress(relPath) {
  const source = join(PUBLIC_DIR, relPath);
  const raw = readFileSync(source);
  if (raw.length < MIN_BYTES) return { relPath, skipped: "too small" };

  const hash = hashOf(raw);
  const stampPath = source + ".hash";
  const brPath = source + ".br";
  const gzPath = source + ".gz";

  if (
    existsSync(stampPath) &&
    existsSync(brPath) &&
    existsSync(gzPath) &&
    readFileSync(stampPath, "utf8").trim() === hash
  ) {
    return { relPath, skipped: "unchanged", br: statSync(brPath).size };
  }

  const started = Date.now();
  const br = zlib.brotliCompressSync(raw, {
    params: {
      [zlib.constants.BROTLI_PARAM_QUALITY]: BROTLI_QUALITY,
      [zlib.constants.BROTLI_PARAM_LGWIN]: 24,
      [zlib.constants.BROTLI_PARAM_SIZE_HINT]: raw.length,
    },
  });
  const gz = zlib.gzipSync(raw, { level: GZIP_LEVEL });

  mkdirSync(dirname(brPath), { recursive: true });
  writeFileSync(brPath, br);
  writeFileSync(gzPath, gz);
  writeFileSync(stampPath, hash + "\n");

  return {
    relPath,
    raw: raw.length,
    br: br.length,
    gz: gz.length,
    ms: Date.now() - started,
  };
}

function main() {
  if (process.argv.includes("--force")) {
    // Drop the stamps so everything is redone.
    for (const rel of collect()) {
      const stamp = join(PUBLIC_DIR, rel + ".hash");
      if (existsSync(stamp)) unlinkSync(stamp);
    }
  }

  const files = collect();
  if (files.length === 0) {
    console.log("[precompress] nothing to do (assets not copied yet?)");
    return;
  }

  let rawTotal = 0;
  let brTotal = 0;
  let built = 0;
  const started = Date.now();

  for (const rel of files) {
    const result = compress(rel);
    if (result.skipped === "too small") continue;
    if (result.skipped === "unchanged") {
      brTotal += result.br;
      continue;
    }
    built += 1;
    rawTotal += result.raw;
    brTotal += result.br;
    console.log(
      `[precompress] ${rel}  ${(result.raw / 1048576).toFixed(2)}MB -> ` +
        `${(result.br / 1048576).toFixed(2)}MB br  (${(result.ms / 1000).toFixed(1)}s)`
    );
  }

  if (built === 0) {
    console.log(`[precompress] up to date (${files.length} files)`);
  } else {
    console.log(
      `[precompress] compressed ${built} file(s), ` +
        `${(rawTotal / 1048576).toFixed(1)}MB -> ${(brTotal / 1048576).toFixed(1)}MB, ` +
        `${((Date.now() - started) / 1000).toFixed(1)}s`
    );
  }
}

main();
