// Capture the screenshots a help page is illustrated with, from the running app.
//
//   node capture.mjs <shots.json> [--base http://localhost:3420] [--api http://127.0.0.1:3421]
//
// The app must be running against a scratch home holding the page's scenario
// project (see scenario_*.py). Each shot opens a job page in a fresh tab,
// turns developer mode off, clicks through tabs, crops to one section of the
// interface and numbers the fields the page's text refers to, (1), (2)...
//
// Fields are found by their label as the user sees it. If a label changes the
// capture fails, naming it: the page's text then needs looking at too.
//
// Needs Node 22+ (global fetch and WebSocket) and Google Chrome (CHROME=path).
import { spawn } from "node:child_process";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";

const args = process.argv.slice(2);
const opt = (name, fallback) => {
  const i = args.indexOf(name);
  return i >= 0 ? args[i + 1] : fallback;
};
const specPath = args.find((a) => a.endsWith(".json"));
const BASE = opt("--base", "http://localhost:3420");
const API = opt("--api", "http://127.0.0.1:3421");
const CHROME =
  process.env.CHROME ||
  "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";
const spec = JSON.parse(fs.readFileSync(specPath, "utf-8"));
const outDir = path.resolve(path.dirname(specPath), spec.outDir || ".");
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

// ---- resolve the scenario's project and jobs --------------------------------
async function api(p) {
  const r = await fetch(`${API}/api/ccp4i2/${p}`);
  if (!r.ok) throw new Error(`${p}: HTTP ${r.status}`);
  return r.json();
}
const projects = await api("projects/");
const project = projects.find((p) => p.name === spec.project);
if (!project) throw new Error(`No project ${spec.project}: run its scenario first`);
const jobs = (await api("jobs/")).filter((j) => j.project === project.id);
const jobId = (number) => {
  const job = jobs.find((j) => j.number === String(number));
  if (!job) throw new Error(`No job ${number} in ${spec.project}`);
  return job.id;
};

// ---- Chrome over the DevTools protocol ---------------------------------------
const profile = fs.mkdtempSync(path.join(os.tmpdir(), "capture-chrome-"));
const chrome = spawn(CHROME, [
  "--headless=new", "--remote-debugging-port=0", `--user-data-dir=${profile}`,
  "--hide-scrollbars", "--use-angle=swiftshader", "--enable-unsafe-swiftshader",
  "--ignore-gpu-blocklist", "about:blank",
], { stdio: "ignore" });
let port;
for (let i = 0; i < 100 && !port; i++) {
  await sleep(100);
  try {
    port = fs.readFileSync(path.join(profile, "DevToolsActivePort"), "utf-8").split("\n")[0];
  } catch {}
}
if (!port) throw new Error("Chrome did not start");

async function connect(wsUrl) {
  const ws = new WebSocket(wsUrl);
  await new Promise((r, e) => { ws.onopen = r; ws.onerror = e; });
  let id = 0;
  const pending = new Map();
  // The page's own errors, kept so that a crash can be reported, not guessed.
  const pageErrors = [];
  ws.onmessage = (ev) => {
    const m = JSON.parse(ev.data);
    if (m.id && pending.has(m.id)) { pending.get(m.id)(m); pending.delete(m.id); }
    else if (m.method === "Runtime.exceptionThrown") {
      const d = m.params.exceptionDetails;
      pageErrors.push(`exception: ${d.exception?.description || d.text}`);
    } else if (m.method === "Runtime.consoleAPICalled" && m.params.type === "error") {
      pageErrors.push(`console.error: ${m.params.args.map((a) => a.value ?? a.description).join(" ")}`);
    }
  };
  const send = (method, params = {}) => new Promise((r) => {
    const i = ++id;
    pending.set(i, r);
    ws.send(JSON.stringify({ id: i, method, params }));
  });
  const evaluate = async (fn, ...fnArgs) => {
    const m = await Promise.race([
      send("Runtime.evaluate", {
        expression: `(${fn})(...${JSON.stringify(fnArgs)})`,
        awaitPromise: true, returnByValue: true,
      }),
      sleep(60000).then(() => { throw new Error(`Page step timed out: ${fn}`.slice(0, 200)); }),
    ]);
    if (m.result?.exceptionDetails) {
      throw new Error(m.result.exceptionDetails.exception?.description || "page error");
    }
    return m.result?.result?.value;
  };
  return { ws, send, evaluate, pageErrors };
}

// ---- page-side helpers (serialised into the page) ------------------------------
const PAGE_HELPERS = () => {
  const norm = (s) => (s || "").replace(/\s+/g, " ").trim().toLowerCase();
  const visible = (el) => { const r = el.getBoundingClientRect(); return r.width > 0 && r.height > 0; };
  window.__cap = {
    // Buttons, tabs and menu items, by their text.
    click(text) {
      const el = [...document.querySelectorAll("button,[role=tab],[role=menuitem],li")]
        .filter(visible).find((e) => norm(e.innerText) === norm(text));
      if (!el) throw new Error(`Nothing to click labelled "${text}"`);
      el.click();
    },
    // The smallest element whose own text is the label.
    label(text) {
      const hits = [...document.querySelectorAll("label,legend,span,p,div,h6,h5,b")]
        .filter(visible)
        .filter((e) => [...e.childNodes].some((n) => n.nodeType === 3 && norm(n.textContent) === norm(text)));
      if (!hits.length) throw new Error(`No field labelled "${text}"`);
      // A field's own <label>/<legend> before a folder heading of the same
      // text ("Reflections" is both, in several interfaces).
      return hits.find((e) => e.tagName === "LABEL" || e.tagName === "LEGEND") || hits[0];
    },
    // A field: its label's form control, or the row that holds a radio group.
    field(text) {
      const l = this.label(text);
      return l.closest(".MuiFormControl-root") && !l.closest("[role=radiogroup]") && l.tagName === "LABEL"
        ? l.closest(".MuiFormControl-root")
        : (l.closest(".MuiFormControl-root") || l.parentElement);
    },
    // A folder of the interface ("Input data", "Controls"): its whole panel.
    // Or {text, scroll: true}: the scrolling panel that shows that text (a report).
    // No section: the open tab's whole panel, for a first look at a new page.
    section(spec) {
      if (spec == null) {
        return [...document.querySelectorAll("[role=tabpanel]")].find(visible) || document.body;
      }
      if (typeof spec === "object") {
        let el = this.containing(spec.text);
        while (el && !/auto|scroll/.test(getComputedStyle(el).overflowY)) el = el.parentElement;
        if (!el) throw new Error(`No scrolling panel shows "${spec.text}"`);
        return el;
      }
      const l = this.label(spec);
      return l.closest(".MuiAccordion-root,.MuiPaper-root") || l.parentElement.parentElement;
    },
    // The element holding the first visible text that starts with the given
    // text. (Walks text nodes: innerText on every element of a report full of
    // plots forces a layout per element and takes minutes.)
    containing(text) {
      const walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
      for (let n = walker.nextNode(); n; n = walker.nextNode()) {
        if (norm(n.textContent).startsWith(norm(text)) && visible(n.parentElement)) {
          return n.parentElement;
        }
      }
      throw new Error(`Nothing shows "${text}"`);
    },
    // A callout target: a field label, or {text, closest} for report parts
    // (the table holding a cell, the plot holding a title).
    target(spec) {
      if (typeof spec === "string") return this.field(spec);
      if (spec.field) return this.field(spec.field);
      const el = this.containing(spec.text);
      const t = spec.closest ? el.closest(spec.closest) : el;
      if (!t) throw new Error(`"${spec.text}" is in no ${spec.closest}`);
      return t;
    },
    rect(el) { const r = el.getBoundingClientRect(); return { x: r.x, y: r.y, width: r.width, height: r.height }; },
  };
  const style = document.createElement("style");
  style.textContent = "nextjs-portal{display:none!important}";
  document.head.appendChild(style);
};

async function shoot(shot) {
  const { targetId } = (await browser.send("Target.createTarget", { url: "about:blank" })).result;
  const target = (await (await fetch(`http://127.0.0.1:${port}/json/list`)).json())
    .find((t) => t.id === targetId);
  const page = await connect(target.webSocketDebuggerUrl);
  const [w, h] = shot.viewport || spec.viewport || [1400, 1000];
  await page.send("Emulation.setDeviceMetricsOverride", { width: w, height: h, deviceScaleFactor: 2, mobile: false });
  await page.send("Page.enable");
  await page.send("Runtime.enable");
  await page.send("Page.navigate", { url: `${BASE}/ccp4i2/project/${project.id}/job/${jobId(shot.job)}` });
  await sleep(shot.settle || spec.settle || 15000);
  await page.evaluate(PAGE_HELPERS);
  // The page can still be compiling or fetching after the settle time (a
  // restarted server, a recompiled interface): wait for its menu bar rather
  // than fail on "Nothing to click labelled View".
  for (let tries = 0; ; tries++) {
    try {
      // Developer mode starts on in every session: turn it off as a user would.
      await page.evaluate(() => __cap.click("View"));
      break;
    } catch (e) {
      if (tries >= 20) {
        // Keep what the page showed instead: a crashed page has no menu.
        const shotOnFail = await page.send("Page.captureScreenshot", { format: "png" });
        fs.writeFileSync(path.join(outDir, `${shot.out}.failed.png`),
          Buffer.from(shotOnFail.result.data, "base64"));
        fs.writeFileSync(path.join(outDir, `${shot.out}.failed.txt`),
          page.pageErrors.join("\n") + "\n");
        throw e;
      }
      await sleep(3000);
      await page.evaluate(PAGE_HELPERS);
    }
  }
  await sleep(400);
  await page.evaluate(() => __cap.click("Turn Dev Mode Off"));
  await sleep(1200);
  for (const tab of shot.tabs || []) {
    await page.evaluate((t) => __cap.click(t), tab);
    await sleep(2500);
  }
  // "expand": open folds by their heading (a report's collapsed sections).
  for (const heading of shot.expand || []) {
    await page.evaluate((t) => {
      const el = __cap.containing(t);
      const toggle = el.closest("[aria-expanded]") || el.closest("[role=button]") || el;
      if (toggle.getAttribute("aria-expanded") !== "true") toggle.click();
    }, heading);
    await sleep(1500);
  }
  // "click": press controls that have no visible label, by CSS selector
  // (a file row's "Expand options" chevron). Every match is pressed; a
  // selector matching nothing fails the shot, as a missing label does.
  for (const selector of shot.click || []) {
    await page.evaluate((sel) => {
      const els = document.querySelectorAll(sel);
      if (!els.length) throw new Error(`Nothing matches ${sel}`);
      els.forEach((el) => el.click());
    }, selector);
    await sleep(1500);
  }
  // Bring the section into view, then measure it and its fields.
  // (Two steps: an error thrown in a timer callback would never reach us.)
  await page.evaluate((sec) => __cap.section(sec).scrollIntoView({ block: "start" }), shot.section);
  await sleep(800);
  const measured = await page.evaluate((sec, labels, until, from, through) => {
    const el = __cap.section(sec);
    const section = __cap.rect(el);
    // A scrolling panel is often taller than what is in it: end the crop at
    // the bottom of its lowest visible content.
    if (typeof sec === "object") {
      let bottom = section.y;
      for (const d of el.querySelectorAll("*")) {
        const r = d.getBoundingClientRect();
        if (r.width > 0 && r.height > 0 && r.bottom > bottom) bottom = r.bottom;
      }
      section.height = Math.min(section.height, bottom + 12 - section.y);
    }
    // "until": end the crop where that text begins (a report's file lists).
    if (until) section.height = __cap.rect(__cap.containing(until)).y - 12 - section.y;
    // "through": extend the crop to the end of a later folder.
    if (through) {
      const t = __cap.rect(__cap.section(through));
      section.height = t.y + t.height - section.y;
    }
    // "from": start the crop where that text begins (below a report's tab bar).
    if (from) {
      const top = __cap.rect(__cap.containing(from)).y - 4;
      section.height -= top - section.y;
      section.y = top;
    }
    return { section, fields: labels.map((l) => __cap.rect(__cap.target(l))) };
  }, shot.section, shot.callouts || [], shot.until || null, shot.from || null, shot.through || null);
  // Number the fields: a badge in a margin left of the section, level with each.
  const margin = (shot.callouts || []).length ? 48 : 0;
  // A callout may carry its own badge text ("1.1"), as the older pages number.
  const badges = (shot.callouts || []).map((c) => (c && c.badge) || null);
  await page.evaluate((sec, fields, badges) => {
    if (fields.length) {
      // Blank the margin: the job list and its divider sit there.
      const strip = document.createElement("div");
      Object.assign(strip.style, {
        position: "fixed", zIndex: 99998, background: "white",
        left: `${sec.x - 64}px`, width: "64px", top: `${sec.y - 10}px`,
        height: `${sec.height + 20}px`,
      });
      document.body.appendChild(strip);
    }
    fields.forEach((f, i) => {
      const b = document.createElement("div");
      b.textContent = badges[i] || String(i + 1);
      Object.assign(b.style, {
        position: "fixed", zIndex: 99999, left: `${sec.x - 44}px`,
        top: `${f.y + f.height / 2 - 13}px`, minWidth: "26px", height: "26px",
        padding: "0 5px", boxSizing: "border-box",
        borderRadius: "13px", background: "#c62828", color: "white",
        font: "bold 15px/26px Roboto, Arial, sans-serif", textAlign: "center",
        boxShadow: "0 1px 3px rgba(0,0,0,.4)",
      });
      document.body.appendChild(b);
    });
  }, measured.section, measured.fields, badges);
  const pad = 8;
  const left = Math.max(0, measured.section.x - pad - margin);
  const clip = {
    x: left, y: Math.max(0, measured.section.y),
    width: measured.section.x + measured.section.width + pad - left,
    height: Math.min(measured.section.height + pad, h - Math.max(0, measured.section.y)),
    scale: 1,
  };
  const png = (await page.send("Page.captureScreenshot", { format: "png", clip })).result.data;
  const out = path.join(outDir, shot.out);
  fs.writeFileSync(out, Buffer.from(png, "base64"));
  console.log(`${shot.out}: ${JSON.stringify(shot.section)}, ${(shot.callouts || []).length} callouts`);
  page.ws.close();
  await browser.send("Target.closeTarget", { targetId });
}

const version = await (await fetch(`http://127.0.0.1:${port}/json/version`)).json();
const browser = await connect(version.webSocketDebuggerUrl);
try {
  for (const shot of spec.shots) {
    if (!args.includes("--only") || opt("--only") === shot.out) await shoot(shot);
  }
} finally {
  browser.ws.close();
  // Chrome writes to its profile as it exits: wait before removing it.
  const exited = new Promise((r) => chrome.once("exit", r));
  chrome.kill();
  await Promise.race([exited, sleep(5000)]);
  fs.rmSync(profile, { recursive: true, force: true, maxRetries: 5 });
}
