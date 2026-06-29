/**
 * Moorhen Scene — ccp4i2 dialect.
 *
 * Extends the portable core file ref (./core.ts) with deployment-coupled
 * references that are meaningful only inside a ccp4i2 deployment:
 *
 *   - relativeUrl  — origin-relative loader URL (e.g. /api/proxy/pdbe/…,
 *                    /api/ccp4i2/download?…). Renamed from the old `path`,
 *                    which read as a filesystem path but actually held these
 *                    origin-relative URLs. Portable WITHIN a deployment only.
 *   - fileId (+projectId)        — symbolic project file.
 *   - job + param (+projectId)   — symbolic job output.
 *
 * These are "recipes for producing a core ref": on a strict-portable export
 * they lower to bundle assets or absolute URLs. See design doc §5.
 *
 * This module is intentionally NOT part of the upstreamable core schema.
 */
import { z } from "zod";

import { CoreFileRefShape, buildScene } from "./core";

/** ccp4i2 file ref: core sources plus the deployment-coupled ones. */
export const Ccp4i2FileRef = CoreFileRefShape.extend({
  relativeUrl: z
    .string()
    .optional()
    .describe("origin-relative URL (/api/…); not portable across deployments"),
  projectId: z.string().optional().describe("UUID; required with fileId or job+param"),
  projectName: z.string().optional().describe("advisory"),
  fileId: z.number().optional(),
  job: z.number().optional().describe("pair with param"),
  param: z.string().optional().describe('job parameter, e.g. "XYZOUT"'),
})
  .strict()
  .superRefine((ref, ctx) => {
    const hasFileId = ref.fileId != null;
    const hasJobParam = ref.job != null && ref.param != null;
    const sources = [
      ref.pdb,
      ref.url,
      ref.bundle,
      ref.cifText,
      ref.relativeUrl,
      hasFileId ? "fileId" : undefined,
      hasJobParam ? "job+param" : undefined,
    ].filter((v) => v != null && v !== "");
    if (sources.length !== 1) {
      ctx.addIssue({
        code: z.ZodIssueCode.custom,
        message:
          "set exactly one of: pdb, url, bundle, cifText, relativeUrl, fileId(+projectId), job+param(+projectId)",
      });
    }
    if ((hasFileId || hasJobParam) && ref.projectId == null) {
      ctx.addIssue({
        code: z.ZodIssueCode.custom,
        path: ["projectId"],
        message: "projectId required with fileId or job+param",
      });
    }
    if (ref.job != null && ref.param == null) {
      ctx.addIssue({
        code: z.ZodIssueCode.custom,
        path: ["param"],
        message: "job requires param",
      });
    }
    if (ref.cifText && ref.kind !== "dictionary") {
      ctx.addIssue({
        code: z.ZodIssueCode.custom,
        path: ["cifText"],
        message: "cifText is only valid on kind: dictionary",
      });
    }
    if (ref.pdb && ref.kind === "dictionary") {
      ctx.addIssue({
        code: z.ZodIssueCode.custom,
        path: ["pdb"],
        message: "pdb is for coordinates, not dictionaries",
      });
    }
  });

/** The ccp4i2-flavoured scene schema (core structure + dialect file refs). */
export const Ccp4i2Scene = buildScene(Ccp4i2FileRef);
export type Ccp4i2Scene = z.infer<typeof Ccp4i2Scene>;
