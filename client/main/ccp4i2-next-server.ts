import path from "path";
import { fileURLToPath } from "node:url";
import { createRequire } from "node:module";
import { Server } from "node:http";
import fs from "node:fs";

// Use createRequire for CJS modules to avoid ESM/CJS interop issues in Electron
const require = createRequire(import.meta.url);
const next = require("next");
const express = require("express");
//import { createProxyMiddleware } from "http-proxy-middleware";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

/**
 * Serves static files from the public directory, handling ASAR archives.
 * For files in report_files/, serves directly from the ASAR-aware fs.
 */
const serveStaticWithAsar = (publicPath: string) => {
  // Log once at setup time
  console.log(`[serveStaticWithAsar] Initialized with publicPath: ${publicPath}`);
  console.log(`[serveStaticWithAsar] publicPath exists: ${fs.existsSync(publicPath)}`);

  // Check if report_files exists
  const reportFilesPath = path.join(publicPath, "report_files", "0.1.0");
  console.log(`[serveStaticWithAsar] report_files path: ${reportFilesPath}`);
  console.log(`[serveStaticWithAsar] report_files exists: ${fs.existsSync(reportFilesPath)}`);

  if (fs.existsSync(reportFilesPath)) {
    try {
      const files = fs.readdirSync(reportFilesPath);
      console.log(`[serveStaticWithAsar] report_files contents: ${files.join(", ")}`);
    } catch (e) {
      console.log(`[serveStaticWithAsar] Error reading report_files: ${e}`);
    }
  }

  return (req: any, res: any, next: any) => {
    // Only handle specific static paths that might have issues
    const staticPaths = ["/report_files/", "/qticons/", "/svgicons/"];
    const isStaticPath = staticPaths.some((p) => req.url.startsWith(p));

    if (!isStaticPath) {
      return next();
    }

    const filePath = path.join(publicPath, req.url.split("?")[0]);
    console.log(`[serveStaticWithAsar] Requested: ${req.url} -> ${filePath}`);

    // Check if file exists (fs is ASAR-aware in Electron)
    const fileExists = fs.existsSync(filePath);
    console.log(`[serveStaticWithAsar] File exists: ${fileExists}`);

    if (fileExists) {
      const stat = fs.statSync(filePath);
      if (stat.isFile()) {
        // Determine content type
        const ext = path.extname(filePath).toLowerCase();
        const contentTypes: Record<string, string> = {
          ".png": "image/png",
          ".jpg": "image/jpeg",
          ".jpeg": "image/jpeg",
          ".gif": "image/gif",
          ".svg": "image/svg+xml",
          ".css": "text/css",
          ".js": "application/javascript",
          ".html": "text/html",
        };
        const contentType = contentTypes[ext] || "application/octet-stream";

        res.setHeader("Content-Type", contentType);
        res.setHeader("Content-Length", stat.size);

        // Stream the file (fs.createReadStream is ASAR-aware)
        const stream = fs.createReadStream(filePath);
        stream.pipe(res);
        return;
      }
    }

    next();
  };
};

/**
 * Starts a Next.js server with custom configurations, including a Content Security Policy (CSP)
 * and optional proxying of API requests to a Django server.
 *
 * @param isDev - A boolean indicating whether the server is running in development mode.
 * @param nextServerPort - The port number on which the Next.js server will listen.
 * @param djangoServerPort - The port number of the Django server to proxy API requests to.
 * @returns A promise that resolves to the created server instance.
 *
 * @remarks
 * - The server sets a Content Security Policy (CSP) header to enhance security.
 * - The `dir` option for the Next.js app is set to the `../renderer` directory relative to the current file.
 * - The server uses Express.js to handle requests and apply middleware.
 * - Proxying to the Django server is currently commented out but can be enabled by uncommenting the relevant code.
 *
 * @example
 * ```typescript
 * const server = await startNextServer(true, 3000, 8000);
 * console.log("Server started successfully");
 * ```
 */
export const startNextServer = async (
  isDev: boolean,
  nextServerPort: number,
  djangoServerPort: number
): Promise<Server> => {
  const nextApp = next({
    dev: isDev,
    dir: path.join(__dirname, "../renderer"), // this is your Next app root
  });

  const handle = nextApp.getRequestHandler();

  await nextApp.prepare();

  // Set the Content Security Policy header
  const csp = {
    defaultSrc: "'self'",
    imgSrc: "'self' data: blob:",
    connectSrc:
      "'self' https://www.ebi.ac.uk https://www.uniprot.org https://pubmed.ncbi.nlm.nih.gov https://raw.githubusercontent.com/MonomerLibrary/monomers/master/ " +
      "https://login.microsoftonline.com https://graph.microsoft.com https://*.microsoftonline.com https://*.microsoft.com " +
      "https://graph.windows.net https://management.azure.com " +
      (process.env.API_BASE_URL || process.env.NEXT_PUBLIC_API_BASE_URL || ""),
    styleSrc:
      "'self' https://cdn.jsdelivr.net 'unsafe-inline' https://fonts.googleapis.com/css2",
    fontSrc:
      "'self' https://cdn.jsdelivr.net 'unsafe-inline' https://fonts.gstatic.com",
    scriptSrc: "'self' https://cdn.jsdelivr.net 'unsafe-inline' 'unsafe-eval'",
    workerSrc: "'self' blob:",
    frameSrc:
      "'self' https://login.microsoftonline.com https://*.microsoftonline.com",
    frameAncestors: "'self'",
  };

  const cspString = [
    `default-src ${csp.defaultSrc}`,
    `img-src ${csp.imgSrc}`,
    `connect-src ${csp.connectSrc}`,
    `style-src ${csp.styleSrc}`,
    `font-src ${csp.fontSrc}`,
    `script-src ${csp.scriptSrc}`,
    `worker-src ${csp.workerSrc}`,
    `frame-src ${csp.frameSrc}`,
    `frame-ancestors ${csp.frameAncestors}`,
  ]
    .join("; ")
    .trim(); // Clean up whitespace

  //const server = createServer((req, res) => handle(req, res));
  const server = express();

  // Proxy Django API requests
  /*
  server.use(
    "/api",
    createProxyMiddleware({
      target: `http://localhost:${djangoServerPort}`, // Your Django server URL
      changeOrigin: true,
      pathRewrite: {
        "^/api": "", // Remove `/api` prefix when forwarding to Django
      },
    })
  );
*/
  server.use((req, res, next_operator) => {
    res.setHeader("Content-Security-Policy", cspString);
    res.setHeader("Cross-Origin-Opener-Policy", "same-origin");
    res.setHeader("Cross-Origin-Embedder-Policy", "require-corp");
    next_operator();
  });

  // Serve static files from public directory
  const publicPath = path.join(__dirname, "../renderer/public");
  console.log(`[Next Server] Static files path: ${publicPath}`);

  // Use custom ASAR-aware middleware for specific paths
  server.use(serveStaticWithAsar(publicPath));

  // Also use express.static as fallback for other static files
  server.use(express.static(publicPath));

  server.use((req, res) => handle(req, res));

  const serverInstance: Server = server.listen(nextServerPort, () => {
    console.log(`Next.js ready on http://localhost:${nextServerPort}`);
  });

  return serverInstance;
};
