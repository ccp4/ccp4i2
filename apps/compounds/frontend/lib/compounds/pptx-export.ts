/**
 * Browser-side PowerPoint export for the Cards aggregation view.
 *
 * Each card is captured as a PNG via html2canvas and laid out into a
 * .pptx slide deck via pptxgenjs. All client-side — no backend.
 */

import html2canvas from 'html2canvas';
import pptxgen from 'pptxgenjs';
import { pinCaptureFonts } from './html2canvas-fonts';

export interface ExportOptions {
  /** Number of cards per slide. 6 is the sweet spot for 16:9 readability. */
  cardsPerSlide?: number;
  /** Optional first slide containing a single PNG of the scorecard key. */
  keyNode?: HTMLElement | null;
  /** PowerPoint deck title (defaults to "Compound cards"). */
  deckTitle?: string;
  /** Filename — should include the .pptx extension. */
  fileName?: string;
  /** Called after each card is captured. */
  onProgress?: (done: number, total: number) => void;
}

interface CardSnap {
  dataUrl: string;
  /** Source pixel dimensions — used to preserve aspect ratio in the slide. */
  widthPx: number;
  heightPx: number;
}

/**
 * Wait for the moment when card captures will be highest fidelity:
 * web fonts loaded, RDKit-generated structure SVGs present in every
 * card. Polls up to ~5 seconds, then proceeds either way.
 */
async function waitForCardsReady(nodes: HTMLElement[]): Promise<void> {
  if (typeof document !== 'undefined' && document.fonts?.ready) {
    await document.fonts.ready;
  }
  const deadline = Date.now() + 5000;
  while (Date.now() < deadline) {
    const allHaveSvg = nodes.every((n) => n.querySelector('svg') !== null);
    if (allHaveSvg) return;
    await new Promise((r) => setTimeout(r, 100));
  }
}

async function snapNode(node: HTMLElement): Promise<CardSnap | null> {
  try {
    const canvas = await html2canvas(node, {
      backgroundColor: '#ffffff',
      scale: 2,
      logging: false,
      onclone: pinCaptureFonts,
    });
    return {
      dataUrl: canvas.toDataURL('image/png'),
      widthPx: canvas.width,
      heightPx: canvas.height,
    };
  } catch (err) {
    console.error('html2canvas failed for card', err);
    return null;
  }
}

/**
 * Compute the largest (width, height) that fits inside the cell while
 * preserving the source aspect ratio, then centre it within the cell.
 */
function fitInCell(
  srcW: number,
  srcH: number,
  cellW: number,
  cellH: number,
): { w: number; h: number; xOffset: number; yOffset: number } {
  const srcRatio = srcW / srcH;
  const cellRatio = cellW / cellH;
  let w: number;
  let h: number;
  if (srcRatio > cellRatio) {
    // source is wider than cell → constrain by width
    w = cellW;
    h = cellW / srcRatio;
  } else {
    h = cellH;
    w = cellH * srcRatio;
  }
  return {
    w,
    h,
    xOffset: (cellW - w) / 2,
    yOffset: (cellH - h) / 2,
  };
}

/**
 * Snap an array of card-rooted DOM nodes to PNGs and assemble them
 * into a downloadable .pptx file. Triggers the browser download
 * directly via pptxgen's writeFile.
 */
export async function exportCardsToPptx(
  cardNodes: HTMLElement[],
  options: ExportOptions = {},
): Promise<void> {
  const {
    cardsPerSlide = 6,
    keyNode = null,
    deckTitle = 'Compound cards',
    fileName,
    onProgress,
  } = options;

  if (cardNodes.length === 0) {
    throw new Error('No cards to export.');
  }

  await waitForCardsReady(cardNodes);

  // Capture sequentially — html2canvas is heavy enough that running
  // many in parallel can OOM the browser on large card sets.
  const snaps: CardSnap[] = [];
  for (let i = 0; i < cardNodes.length; i++) {
    const snap = await snapNode(cardNodes[i]);
    if (snap) snaps.push(snap);
    // Throttle progress to every 5th card (and the last one) so the
    // toolbar doesn't re-render on every iteration of a 100-card deck.
    if (i % 5 === 4 || i === cardNodes.length - 1) {
      onProgress?.(i + 1, cardNodes.length);
    }
  }

  // Optional key snap (for the lead slide)
  let keySnap: CardSnap | null = null;
  if (keyNode) {
    keySnap = await snapNode(keyNode);
  }

  // Assemble the deck — 16:9 widescreen.
  const pres = new pptxgen();
  pres.title = deckTitle;
  pres.layout = 'LAYOUT_WIDE';
  // pptxgen units are inches; LAYOUT_WIDE is 13.333" × 7.5"
  const slideW = 13.333;
  const slideH = 7.5;
  const margin = 0.3;

  // Lead slide with the scorecard key, if provided.
  if (keySnap) {
    const slide = pres.addSlide();
    const usable = { w: slideW - 2 * margin, h: slideH - 2 * margin };
    const fit = fitInCell(keySnap.widthPx, keySnap.heightPx, usable.w, usable.h);
    slide.addImage({
      data: keySnap.dataUrl,
      x: margin + fit.xOffset,
      y: margin + fit.yOffset,
      w: fit.w,
      h: fit.h,
    });
  }

  // Card slides. Auto-pick a sensible grid for the requested cards/slide:
  //   4 →  2×2,  6 → 3×2,  9 → 3×3, 12 → 4×3
  const grid = gridFor(cardsPerSlide);
  const cellW = (slideW - 2 * margin) / grid.cols;
  const cellH = (slideH - 2 * margin) / grid.rows;

  for (let i = 0; i < snaps.length; i += cardsPerSlide) {
    const slide = pres.addSlide();
    const slot = snaps.slice(i, i + cardsPerSlide);
    for (let j = 0; j < slot.length; j++) {
      const col = j % grid.cols;
      const row = Math.floor(j / grid.cols);
      const cellX = margin + col * cellW;
      const cellY = margin + row * cellH;
      const fit = fitInCell(slot[j].widthPx, slot[j].heightPx, cellW * 0.95, cellH * 0.95);
      slide.addImage({
        data: slot[j].dataUrl,
        x: cellX + (cellW - fit.w) / 2,
        y: cellY + (cellH - fit.h) / 2,
        w: fit.w,
        h: fit.h,
      });
    }
  }

  const finalName = fileName || `compound-cards-${todayIso()}.pptx`;
  await pres.writeFile({ fileName: finalName });
}

function gridFor(cardsPerSlide: number): { cols: number; rows: number } {
  switch (cardsPerSlide) {
    case 4: return { cols: 2, rows: 2 };
    case 6: return { cols: 3, rows: 2 };
    case 9: return { cols: 3, rows: 3 };
    case 12: return { cols: 4, rows: 3 };
    default: {
      const cols = Math.ceil(Math.sqrt(cardsPerSlide));
      const rows = Math.ceil(cardsPerSlide / cols);
      return { cols, rows };
    }
  }
}

function todayIso(): string {
  return new Date().toISOString().slice(0, 10);
}

// ---------------------------------------------------------------------------
// Single-node snapshot export — used by BulletsView for a "screenshot
// of the whole table on one slide" PPTX. Capture the root once, fit it
// into a single LAYOUT_WIDE slide preserving aspect ratio.
//
// For very tall tables (many compounds) the image scales down to fit
// the slide height, narrowing it horizontally. Acceptable for an MVP
// — row-group splitting is a follow-up if chemists routinely export
// >30-compound bullet decks and want the per-row legibility back.
// ---------------------------------------------------------------------------

export interface SnapshotExportOptions {
  /** PowerPoint deck title. */
  deckTitle?: string;
  /** Filename — should include the .pptx extension. */
  fileName?: string;
}

export async function exportSnapshotToPptx(
  node: HTMLElement,
  options: SnapshotExportOptions = {},
): Promise<void> {
  const { deckTitle = 'Bullets table', fileName } = options;

  // Same font-loading guard the cards path uses, plus a SVG presence
  // check across all RDKit thumbnails inside the table.
  await waitForCardsReady([node]);

  const snap = await snapNode(node);
  if (!snap) {
    throw new Error('Failed to capture the bullets table.');
  }

  const pres = new pptxgen();
  pres.title = deckTitle;
  pres.layout = 'LAYOUT_WIDE';
  const slideW = 13.333;
  const slideH = 7.5;
  const margin = 0.3;
  const usable = { w: slideW - 2 * margin, h: slideH - 2 * margin };

  const slide = pres.addSlide();
  const fit = fitInCell(snap.widthPx, snap.heightPx, usable.w, usable.h);
  slide.addImage({
    data: snap.dataUrl,
    x: margin + fit.xOffset,
    y: margin + fit.yOffset,
    w: fit.w,
    h: fit.h,
  });

  const finalName = fileName || `bullets-${todayIso()}.pptx`;
  await pres.writeFile({ fileName: finalName });
}
