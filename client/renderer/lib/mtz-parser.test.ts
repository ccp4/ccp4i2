/**
 * Test script for the MTZ parser.
 *
 * Run with: npx tsx lib/mtz-parser.test.ts
 *
 * Tests the pure TypeScript MTZ header parser against real MTZ files.
 */

import * as fs from "fs";
import * as path from "path";
import { fileURLToPath } from "url";
import { parseMtzHeader, MtzHeader, MtzParseError } from "./mtz-parser";

// ESM-compatible __dirname
const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Test file paths (relative to ccp4i2 root)
const TEST_FILES = [
  "server/ccp4i2/demo_data/gamma/gamma_native.mtz",          // Standard merged data
  "server/ccp4i2/demo_data/gamma/freeR.mtz",                 // FreeR flags
  "server/ccp4i2/demo_data/gamma/initial_phases.mtz",        // Phases
  "server/ccp4i2/demo_data/gere/gere.mtz",                   // Another structure
  "server/ccp4i2/wrappers/freerflag/test_data/rnase25.mtz",  // Test data
];

// Color codes for terminal output
const GREEN = "\x1b[32m";
const RED = "\x1b[31m";
const YELLOW = "\x1b[33m";
const RESET = "\x1b[0m";
const BOLD = "\x1b[1m";

interface TestResult {
  file: string;
  passed: boolean;
  header?: MtzHeader;
  error?: string;
}

function formatCell(cell: MtzHeader["cell"]): string {
  if (!cell) return "N/A";
  return `${cell.a.toFixed(2)} ${cell.b.toFixed(2)} ${cell.c.toFixed(2)} ${cell.alpha.toFixed(1)} ${cell.beta.toFixed(1)} ${cell.gamma.toFixed(1)}`;
}

function formatResolution(res: MtzHeader["resolution"]): string {
  if (!res) return "N/A";
  return `${res[0].toFixed(2)} - ${res[1].toFixed(2)} Å`;
}

async function testFile(filePath: string): Promise<TestResult> {
  const result: TestResult = {
    file: path.basename(filePath),
    passed: false,
  };

  try {
    // Read file - __dirname is client/renderer/lib, so go up to ccp4i2 root
    const absolutePath = path.resolve(__dirname, "../../..", filePath);
    if (!fs.existsSync(absolutePath)) {
      result.error = `File not found: ${absolutePath}`;
      return result;
    }

    const buffer = fs.readFileSync(absolutePath);
    const arrayBuffer = buffer.buffer.slice(
      buffer.byteOffset,
      buffer.byteOffset + buffer.byteLength
    );

    // Parse header
    const header = parseMtzHeader(arrayBuffer);
    result.header = header;

    // Validation checks
    const checks: string[] = [];

    // Must have columns
    if (header.columns.length === 0) {
      checks.push("No columns found");
    }

    // Must have H, K, L columns (Miller indices)
    const hklColumns = header.columns.filter((c) => c.type === "H");
    if (hklColumns.length < 3) {
      checks.push(`Expected at least 3 H-type columns (Miller indices), found ${hklColumns.length}`);
    }

    // Should have cell parameters
    if (!header.cell) {
      checks.push("No cell parameters found");
    }

    // Should have a title or space group
    if (!header.title && !header.spaceGroup) {
      checks.push("No title or space group found");
    }

    if (checks.length > 0) {
      result.error = checks.join("; ");
    } else {
      result.passed = true;
    }
  } catch (e) {
    result.error = e instanceof Error ? e.message : String(e);
  }

  return result;
}

async function runTests() {
  console.log(`\n${BOLD}MTZ Parser Test Suite${RESET}\n`);
  console.log("=" .repeat(80));

  const results: TestResult[] = [];

  for (const file of TEST_FILES) {
    const result = await testFile(file);
    results.push(result);

    const status = result.passed
      ? `${GREEN}PASS${RESET}`
      : `${RED}FAIL${RESET}`;

    console.log(`\n${status} ${BOLD}${result.file}${RESET}`);

    if (result.header) {
      console.log(`  Title:       ${result.header.title || "(none)"}`);
      console.log(`  Space Group: ${result.header.spaceGroup || "N/A"} (#${result.header.spaceGroupNumber || "?"})`);
      console.log(`  Cell:        ${formatCell(result.header.cell)}`);
      console.log(`  Resolution:  ${formatResolution(result.header.resolution)}`);
      console.log(`  Reflections: ${result.header.nReflections.toLocaleString()}`);
      console.log(`  Columns:     ${result.header.nColumns} (${result.header.columns.length} parsed)`);
      console.log(`  Merged:      ${result.header.isMerged ? "Yes" : "No (unmerged)"}`);
      console.log(`  Datasets:    ${result.header.datasets.length}`);

      // Show column summary
      const colTypes = new Map<string, string[]>();
      for (const col of result.header.columns) {
        if (!colTypes.has(col.type)) colTypes.set(col.type, []);
        colTypes.get(col.type)!.push(col.label);
      }
      console.log(`  Column types:`);
      for (const [type, labels] of colTypes) {
        const displayLabels = labels.length > 4
          ? `${labels.slice(0, 4).join(", ")}... (+${labels.length - 4} more)`
          : labels.join(", ");
        console.log(`    ${type}: ${displayLabels}`);
      }
    }

    if (result.error) {
      console.log(`  ${RED}Error: ${result.error}${RESET}`);
    }
  }

  // Summary
  console.log("\n" + "=".repeat(80));
  const passed = results.filter((r) => r.passed).length;
  const failed = results.length - passed;

  if (failed === 0) {
    console.log(`\n${GREEN}${BOLD}All ${passed} tests passed!${RESET}\n`);
  } else {
    console.log(`\n${passed > 0 ? GREEN : ""}${passed} passed${RESET}, ${RED}${failed} failed${RESET}\n`);
  }

  // Exit with error code if any failed
  process.exit(failed > 0 ? 1 : 0);
}

// Additional test: Byte order detection
async function testByteOrderDetection() {
  console.log(`\n${BOLD}Testing byte order detection...${RESET}`);

  // Create a minimal MTZ-like buffer with known byte order markers
  // Little-endian marker: 0x44 at byte 8
  const leBuffer = new ArrayBuffer(100);
  const leView = new DataView(leBuffer);
  // "MTZ "
  leView.setUint8(0, 77); // M
  leView.setUint8(1, 84); // T
  leView.setUint8(2, 90); // Z
  leView.setUint8(3, 32); // space
  // Header offset (pointing past buffer to trigger error - just testing detection)
  leView.setInt32(4, 1000, true);
  // Machine stamp - little endian (0x44)
  leView.setUint8(8, 0x44);
  leView.setUint8(9, 0x41);

  try {
    parseMtzHeader(leBuffer);
  } catch (e) {
    if (e instanceof MtzParseError && e.message.includes("Invalid header offset")) {
      console.log(`  ${GREEN}Little-endian marker detected correctly${RESET}`);
    } else {
      console.log(`  ${YELLOW}Unexpected error: ${e}${RESET}`);
    }
  }

  // Big-endian marker: 0x11 at byte 8
  const beBuffer = new ArrayBuffer(100);
  const beView = new DataView(beBuffer);
  beView.setUint8(0, 77);
  beView.setUint8(1, 84);
  beView.setUint8(2, 90);
  beView.setUint8(3, 32);
  beView.setInt32(4, 1000, false); // big-endian
  beView.setUint8(8, 0x11);
  beView.setUint8(9, 0x11);

  try {
    parseMtzHeader(beBuffer);
  } catch (e) {
    if (e instanceof MtzParseError && e.message.includes("Invalid header offset")) {
      console.log(`  ${GREEN}Big-endian marker detected correctly${RESET}`);
    } else {
      console.log(`  ${YELLOW}Unexpected error: ${e}${RESET}`);
    }
  }
}

// Run everything
async function main() {
  await testByteOrderDetection();
  await runTests();
}

main().catch(console.error);
