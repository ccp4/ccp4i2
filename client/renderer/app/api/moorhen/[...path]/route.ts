import { NextRequest, NextResponse } from "next/server";
import { readFile } from "fs/promises";
import { join } from "path";

/**
 * API route to serve Moorhen static files with proper CORP headers.
 *
 * Next.js doesn't apply middleware or headers() config to files in public/.
 * This route serves those files through the API, ensuring the
 * Cross-Origin-Resource-Policy header is set for COEP compatibility.
 *
 * Usage: /api/moorhen/moorhen.js instead of /moorhen.js
 */

// Content type mapping
const contentTypes: Record<string, string> = {
  ".js": "text/javascript",
  ".wasm": "application/wasm",
  ".data": "application/octet-stream",
  ".gz": "application/gzip",
  ".html": "text/html",
  ".css": "text/css",
  ".json": "application/json",
  ".ts": "text/typescript",
  ".png": "image/png",
  ".ico": "image/x-icon",
  ".svg": "image/svg+xml",
  ".cif": "chemical/x-cif",
  ".pdb": "chemical/x-pdb",
};

export async function GET(
  _request: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  const { path } = await params;
  const filePath = path.join("/");

  // Security: only allow specific file patterns for moorhen resources
  const allowedPatterns = [
    // Root-level moorhen files
    /^moorhen(64)?\.js$/,
    /^moorhen(64)?\.wasm$/,
    /^CootWorker\.js$/,
    /^CootWorker\.ts$/,
    /^coot_env_web\.js$/,
    /^RDKit_minimal\.js$/,
    /^RDKit_minimal\.wasm$/,
    /^index\.d\.ts$/,
    // baby-gru folder - CSS, data files, images, etc.
    /^baby-gru\/.+\.(js|wasm|data|gz|css|json|html|png|ico|svg|cif|pdb)$/,
    // baby-gru subfolders (monomers, pixmaps, tutorials)
    /^baby-gru\/monomers\/.+$/,
    /^baby-gru\/pixmaps\/.+$/,
    /^baby-gru\/tutorials\/.+$/,
  ];

  const isAllowed = allowedPatterns.some((pattern) => pattern.test(filePath));
  if (!isAllowed) {
    return NextResponse.json({ error: "File not allowed" }, { status: 403 });
  }

  // Determine the public directory path
  // In production (standalone build), files are in the public directory relative to cwd
  // In development, cwd is 'client/' but public files are in 'client/renderer/public/'
  const cwd = process.cwd();
  const publicDir =
    process.env.NODE_ENV === "production"
      ? join(cwd, "public")
      : cwd.endsWith("client")
        ? join(cwd, "renderer", "public")
        : join(cwd, "public");

  const fullPath = join(publicDir, filePath);

  try {
    const fileContent = await readFile(fullPath);

    // Determine content type
    const ext = "." + filePath.split(".").pop();
    const contentType = contentTypes[ext] || "application/octet-stream";

    // Build headers - all files need CORP for COEP compatibility
    const headers: Record<string, string> = {
      "Content-Type": contentType,
      "Cross-Origin-Resource-Policy": "cross-origin",
      "Access-Control-Allow-Origin": "*",
      "Cache-Control": "public, max-age=31536000, immutable",
    };

    // Worker scripts (moorhen.js, CootWorker.js) need COEP/COOP headers
    // because they run in a worker context that needs SharedArrayBuffer
    if (filePath.endsWith(".js")) {
      headers["Cross-Origin-Embedder-Policy"] = "require-corp";
      headers["Cross-Origin-Opener-Policy"] = "same-origin";
    }

    return new NextResponse(fileContent, {
      status: 200,
      headers,
    });
  } catch (error) {
    console.error(`Failed to read moorhen file: ${filePath}`, error);
    return NextResponse.json({ error: "File not found" }, { status: 404 });
  }
}
