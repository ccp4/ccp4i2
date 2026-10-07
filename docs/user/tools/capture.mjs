// Capture the screenshots a help page is illustrated with, from the running app.
//
//   node capture.mjs <shots.json> [--base http://localhost:3420] [--api http://127.0.0.1:3421]
//
// The app must be running against a scratch home holding the page's scenario
// project (see scenario_*.py). Each shot opens a job page in a fresh tab,
// turns developer mode off (if on), clicks through tabs, crops to one section of the
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

// ---- One capture at a time ---------------------------------------------------
// Several captures at once (subagents drafting pages in parallel) slowed every
// shot five-fold, against one dev server. Queue on a lock directory instead;
// one left by a capture that died is taken over after 20 minutes.
const LOCK = path.join(os.tmpdir(), "ccp4i2-capture.lock");
for (let waited = 0; ; waited += 2000) {
  try {
    fs.mkdirSync(LOCK);
    break;
  } catch {
    let age = 0;
    try { age = Date.now() - fs.statSync(LOCK).mtimeMs; } catch { continue; }
    if (age > 20 * 60 * 1000) { fs.rmSync(LOCK, { recursive: true, force: true }); continue; }
    if (waited % 60000 === 0) console.log("Waiting for another capture to finish...");
    await sleep(2000);
  }
}
const releaseLock = () => fs.rmSync(LOCK, { recursive: true, force: true });
process.on("exit", releaseLock);
for (const signal of ["SIGINT", "SIGTERM"]) process.on(signal, () => process.exit(1));

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
      if (typeof spec === "object" && spec.graph) {
        // {graph: "..."}: one graph and its menu, the menu named by its label
        // or by the graph it shows (as in "plots"; give "plots" first to
        // choose the graph, then name the choice here). The Qt pictures
        // were one graph each; a fold's crop showed a dozen.
        const sq = (t) => (t || "").replace(/\s+/g, " ").trim();
        const roots = [...document.querySelectorAll(".MuiAutocomplete-root")].filter((r) => {
          const label = r.querySelector("label"), input = r.querySelector("input");
          return (label && sq(label.textContent).startsWith(sq(spec.graph))) ||
            (input && sq(input.value).startsWith(sq(spec.graph)));
        });
        const root = roots[spec.nth || 0];
        if (!root) throw new Error(`No graph shows "${spec.graph}"`);
        let el = root.parentElement;
        while (el && !el.querySelector("canvas,svg")) el = el.parentElement;
        if (!el) throw new Error(`No drawing beside the graph menu "${spec.graph}"`);
        return el.closest(".MuiPaper-root") || el;
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
    // Open every collapsed fold around el, outermost first: a graph in a
    // closed fold is in the page, measured at no height, and a crop of it
    // photographs the fold headers in front of it. Returns how many it opened.
    // A graph menu by its label or the graph it shows now ("plots", {graph}).
    graphMenu(name, nth) {
      const sq = (t) => (t || "").replace(/\s+/g, " ").trim();
      const roots = [...document.querySelectorAll(".MuiAutocomplete-root")].filter((r) => {
        const label = r.querySelector("label"), input = r.querySelector("input");
        return (label && sq(label.textContent).startsWith(sq(name))) ||
          (input && sq(input.value).startsWith(sq(name)));
      });
      return roots[nth || 0] || null;
    },
    reveal(el) {
      const folds = [];
      for (let a = el.closest(".MuiAccordion-root"); a; a = a.parentElement && a.parentElement.closest(".MuiAccordion-root")) {
        folds.unshift(a);
      }
      let opened = 0;
      for (const fold of folds) {
        const summary = fold.querySelector(".MuiAccordionSummary-root");
        if (summary && summary.getAttribute("aria-expanded") !== "true") { summary.click(); opened++; }
      }
      return opened;
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
  // Developer mode is off by default since #753; turn it off only if a
  // stored preference left it on (the menu then offers "Turn Dev Mode Off").
  const wasOn = await page.evaluate(() => {
    try { __cap.click("Turn Dev Mode Off"); return true; } catch (e) { return false; }
  });
  if (!wasOn) {
    // Already off: close the View menu, or it covers the top of every shot.
    for (const type of ["keyDown", "keyUp"]) {
      await page.send("Input.dispatchKeyEvent",
        { type, key: "Escape", code: "Escape", windowsVirtualKeyCode: 27 });
    }
  }
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
  // "plots": choose the graph a report's graph group shows, as a user does:
  // open its menu (an Autocomplete labelled with the group's title) and pick
  // the option whose text starts with "choose". {"menu": its label or the
  // graph it shows now, "choose": title, "nth": which of several matches
  // (default 0)}. A group menu picks a table, then its "Plot" menu a graph. The
  // Qt pictures each showed one chosen graph, which a shot could not reach.
  for (const plot of shot.plots || []) {
    const opened = await page.evaluate((pl) => {
      const root = __cap.graphMenu(pl.menu, pl.nth);
      return root ? __cap.reveal(root) : 0;
    }, plot);
    if (opened) await sleep(1500);
    await page.evaluate((pl) => {
      // A menu is named by its label or, better, by what it shows now: the
      // labels are generic ("Plot", "Group of 13 graphs") and there are many.
      const roots = [...document.querySelectorAll(".MuiAutocomplete-root")].filter((r) => {
        const label = r.querySelector("label");
        const input = r.querySelector("input");
        const sq = (t) => t.replace(/\s+/g, " ").trim();
        return (label && sq(label.textContent).startsWith(sq(pl.menu))) ||
          (input && sq(input.value).startsWith(sq(pl.menu)));
      });
      const root = roots[pl.nth || 0];
      if (!root) throw new Error(`No graph menu labelled or showing "${pl.menu}"`);
      if (!root) throw new Error(`"${pl.menu}" is not a graph menu`);
      root.scrollIntoView({ block: "center" });
      const opener = root.querySelector(".MuiAutocomplete-popupIndicator");
      if (opener) opener.click(); else root.querySelector("input").focus();
    }, plot);
    await sleep(800);
    await page.evaluate((pl) => {
      // Graph titles carry runs of spaces ("within   5 sd"): compare squeezed.
      const sq = (t) => t.replace(/\s+/g, " ").trim();
      const options = [...document.querySelectorAll('[role="option"]')];
      const option = options.find((o) => sq(o.textContent).startsWith(sq(pl.choose)));
      if (!option) {
        throw new Error(`Graph menu "${pl.menu}" has no "${pl.choose}"; it offers: ` +
          options.map((o) => o.textContent.trim()).join(" | "));
      }
      option.click();
    }, plot);
    await sleep(2500);
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
  // "outline": the page as text instead of a picture, for finding out what a
  // new page's sections and labels are called (a full-page screenshot read
  // for that costs many times more). Written to <out>.outline.txt and printed.
  if (shot.outline) {
    // Graph menus show only the graph chosen; open each in turn and note what
    // else it offers (for "plots"), then close it.
    const nMenus = await page.evaluate(() =>
      document.querySelectorAll(".MuiAutocomplete-root").length);
    for (let i = 0; i < nMenus; i++) {
      const opened = await page.evaluate((k) => {
        const r = document.querySelectorAll(".MuiAutocomplete-root")[k];
        const opener = r && r.querySelector(".MuiAutocomplete-popupIndicator");
        if (!opener || !r.getClientRects().length) return false;
        r.scrollIntoView({ block: "center" });
        opener.click();
        return true;
      }, i);
      if (!opened) continue;
      await sleep(400);
      await page.evaluate((k) => {
        const r = document.querySelectorAll(".MuiAutocomplete-root")[k];
        const opts = [...document.querySelectorAll('[role="option"]')]
          .map((o) => o.textContent.replace(/\s+/g, " ").trim());
        r.setAttribute("data-cap-options", opts.join(" | "));
        r.querySelector(".MuiAutocomplete-popupIndicator").click();
      }, i);
      await sleep(200);
    }
    const lines = await page.evaluate((sec) => {
      // No section: the job's panel (its fixed id, project/[id]/layout.tsx),
      // not the whole page with the job list beside it.
      const root = sec == null
        ? (document.querySelector('[data-panel-id="project-content"]') || __cap.section(null))
        : __cap.section(sec);
      const text = (el) => (el ? el.innerText || el.textContent || "" : "").replace(/\s+/g, " ").trim();
      const out = [];
      root.querySelectorAll('[role="tab"]').forEach((t) => {
        if (t.offsetParent) out.push(`tab${t.getAttribute("aria-selected") === "true" ? "*" : " "} ${text(t)}`);
      });
      const seen = new Set();
      const walker = document.createTreeWalker(root, NodeFilter.SHOW_ELEMENT);
      for (let el = walker.currentNode; el; el = walker.nextNode()) {
        if ([...seen].some((s) => s.contains(el)) || !el.getClientRects().length) continue;
        const cls = el.className && el.className.baseVal === undefined ? String(el.className) : "";
        if (cls.includes("MuiAccordionSummary-root")) {
          out.push(`section  ${text(el)}`); seen.add(el);
        } else if (el.tagName === "TABLE") {
          const head = [...el.querySelectorAll("thead th, tr:first-child th")].map(text);
          const rows = [...el.querySelectorAll("tbody tr")];
          out.push(`table    ${head.join(" | ")}  (${rows.length} rows)`);
          rows.slice(0, 4).forEach((r) => out.push(`  row    ${[...r.children].map(text).join(" | ")}`));
          seen.add(el);
        } else if (cls.includes("MuiFormControlLabel-root")) {
          // A label with no input is a group's title (a radio group's caption),
          // not an option: printed as a checkbox it read as a blank choice.
          const box = el.querySelector("input");
          if (!box) {
            if (text(el).trim()) out.push(`text     ${text(el)}`);
          } else {
            const kind = box.type === "radio" ? "radio" : "check";
            out.push(`${kind.padEnd(8)} [${box.checked ? "x" : " "}] ${text(el)}`);
          }
          seen.add(el);
        } else if (cls.includes("MuiFormControl-root") || cls.includes("MuiTextField-root")) {
          const label = text(el.querySelector("label"));
          const input = el.querySelector("input, textarea");
          const shown = text(el.querySelector('[role="combobox"], .MuiSelect-select'));
          const value = shown || (input ? input.value : "");
          out.push(`field    ${label || "(no label)"} = ${value.slice(0, 90)}`); seen.add(el);
          const menu = el.closest(".MuiAutocomplete-root");
          const options = menu && menu.getAttribute("data-cap-options");
          if (options) out.push(`  menu   ${options}`);
        } else if (el.tagName === "BUTTON" && text(el) && el.getAttribute("role") !== "tab") {
          out.push(`button   ${text(el)}`); seen.add(el);
        } else if ((el.tagName === "LI" || cls.includes("MuiListItemText-root")) && text(el)) {
          out.push(`item     ${text(el).slice(0, 140)}`); seen.add(el);
        } else if (/^(P|PRE|H[1-6])$/.test(el.tagName) && text(el)) {
          out.push(`text     ${text(el).slice(0, 140)}`); seen.add(el);
        } else if (el.tagName === "SPAN" && !el.children.length && text(el) &&
                   !el.closest("label, button, [role=tab], li, p, table, .MuiFormControl-root")) {
          // A report's own text is a bare span (CCP4i2ReportText): it was
          // missed, so "Number of waters found: 48" read as absent.
          out.push(`text     ${text(el).slice(0, 140)}`); seen.add(el);
        }
      }
      return out;
    }, shot.section ?? null);
    const out = path.join(outDir, `${shot.out}.outline.txt`);
    fs.writeFileSync(out, lines.join("\n") + "\n");
    console.log(`--- ${shot.out} (job ${shot.job}, ${(shot.tabs || []).join(" > ")})`);
    console.log(lines.join("\n"));
    page.ws.close();
    await browser.send("Target.closeTarget", { targetId });
    return;
  }
  // Bring the section into view, then measure it and its fields.
  // (Two steps: an error thrown in a timer callback would never reach us.)
  if (shot.section && shot.section.graph) {
    const opened = await page.evaluate((sec) => {
      const root = __cap.graphMenu(sec.graph, sec.nth);
      return root ? __cap.reveal(root) : 0;
    }, shot.section);
    if (opened) await sleep(1500);
  }
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
    // The clip is in document coordinates, the rectangles in the viewport's:
    // when scrollIntoView scrolled the document itself (a graph card in a
    // report's grid), the picture came from where the card had been.
    return { section, fields: labels.map((l) => __cap.rect(__cap.target(l))),
             scroll: { x: window.scrollX, y: window.scrollY } };
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
  // A graph card sits beside others in a report's grid: padding reaches into
  // its neighbour's text.
  const pad = shot.section && shot.section.graph ? 1 : 8;
  const left = Math.max(0, measured.section.x - pad - margin);
  const clip = {
    x: left + measured.scroll.x, y: Math.max(0, measured.section.y) + measured.scroll.y,
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
