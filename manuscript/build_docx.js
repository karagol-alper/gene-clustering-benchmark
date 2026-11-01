// Render a manuscript JSON (produced by build_manuscript.py) to .docx.
// usage: node build_docx.js input.json output.docx [--supplement]
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, ImageRun, Table, TableRow, TableCell,
  WidthType, AlignmentType, BorderStyle, ShadingType, PageBreak, Footer, Header,
  PageNumber, LineNumberRestartFormat, HeadingLevel, PageOrientation, VerticalAlign,
  TableLayoutType,
} = require("docx");

const [, , inFile, outFile, ...flags] = process.argv;
const SUPP = flags.includes("--supplement");
const doc = JSON.parse(fs.readFileSync(inFile, "utf8"));
const FONT = "Times New Roman";
const BODY_PT = SUPP ? 11 : 12;
const LINE = SUPP ? 300 : 480;          // 1.25 / double spacing (240 = single)

function runs(rs, base = {}) {
  return (rs || []).map((r) => new TextRun({
    text: r.t, bold: r.b || base.bold, italics: r.i || base.italics,
    superScript: r.sup || false, subScript: r.sub || false,
    font: r.mono ? "Courier New" : FONT,
    size: Math.round(2 * (r.size || base.size || BODY_PT)),
    color: base.color,
  }));
}

function para(b) {
  const st = b.style || "body";
  const common = { spacing: { line: LINE, after: 0, before: 0 } };
  switch (st) {
    case "title":
      return new Paragraph({ alignment: AlignmentType.LEFT, spacing: { line: 360, after: 240 },
        children: runs(b.runs, { bold: true, size: SUPP ? 15 : 16 }) });
    case "authors":
      return new Paragraph({ spacing: { line: 360, after: 120 }, children: runs(b.runs) });
    case "affil":
      return new Paragraph({ spacing: { line: 300, after: 40 }, children: runs(b.runs, { size: BODY_PT - 1 }) });
    case "abstract":
      return new Paragraph({ ...common, spacing: { line: LINE, after: 120 }, alignment: AlignmentType.JUSTIFIED,
        children: runs(b.runs) });
    case "first":
      return new Paragraph({ ...common, alignment: AlignmentType.JUSTIFIED, children: runs(b.runs) });
    case "body":
      return new Paragraph({ ...common, alignment: AlignmentType.JUSTIFIED,
        indent: { firstLine: 360 }, children: runs(b.runs) });
    case "bodyflush":
      return new Paragraph({ ...common, alignment: AlignmentType.JUSTIFIED, spacing: { line: LINE, before: 120 },
        children: runs(b.runs) });
    case "caption":
      return new Paragraph({ spacing: { line: SUPP ? 260 : 360, before: 120, after: 160 },
        alignment: AlignmentType.JUSTIFIED, children: runs(b.runs, { size: BODY_PT - 1 }) });
    case "tablecaption":
      return new Paragraph({ spacing: { line: 276, before: 240, after: 100 }, keepNext: true,
        children: runs(b.runs, { size: BODY_PT - 1 }) });
    case "tablenote":
      return new Paragraph({ spacing: { line: 252, before: 60, after: 200 }, alignment: AlignmentType.JUSTIFIED,
        children: runs(b.runs, { size: BODY_PT - 2.5 }) });
    case "ref":
      return new Paragraph({ spacing: { line: SUPP ? 260 : 360, after: 60 }, indent: { left: 480, hanging: 480 },
        children: runs(b.runs, { size: BODY_PT - 1 }) });
    case "eq":
      return new Paragraph({ spacing: { line: 300, before: 80, after: 80 }, alignment: AlignmentType.CENTER,
        children: runs(b.runs) });
    default:
      return new Paragraph({ ...common, children: runs(b.runs) });
  }
}

function heading(b) {
  const lvl = { h1: HeadingLevel.HEADING_1, h2: HeadingLevel.HEADING_2, h3: HeadingLevel.HEADING_3 }[b.type];
  const size = { h1: 14, h2: 12.5, h3: 12 }[b.type];
  return new Paragraph({
    heading: lvl, keepNext: true,
    spacing: { before: b.type === "h1" ? 360 : 240, after: 120, line: 276 },
    children: [new TextRun({ text: b.text, bold: b.type !== "h3", italics: b.type === "h3",
      font: FONT, size: 2 * size, color: "000000" })],
  });
}

const thin = { style: BorderStyle.SINGLE, size: 4, color: "000000" };
const none = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };

function table(b, pageWidthDxa) {
  const ncol = b.header[0].length;
  const total = pageWidthDxa;
  let widths = b.widths ? b.widths.slice() : Array(ncol).fill(1);
  const sum = widths.reduce((a, c) => a + c, 0);
  widths = widths.map((w) => Math.floor((w / sum) * total));
  widths[widths.length - 1] += total - widths.reduce((a, c) => a + c, 0);
  const fs = b.font_size || 8.5;
  const mkCell = (content, i, isHead, isLastHead, isLastRow, shade) => new TableCell({
    width: { size: widths[i], type: WidthType.DXA },
    margins: { top: 30, bottom: 30, left: 60, right: 60 },
    verticalAlign: VerticalAlign.CENTER,
    shading: shade ? { fill: "F2F2F2", type: ShadingType.CLEAR, color: "auto" } : undefined,
    borders: {
      top: isHead && content.__firstHead ? thin : none,
      bottom: (isLastHead || isLastRow) ? thin : none,
      left: none, right: none,
    },
    children: [new Paragraph({
      alignment: (b.align && b.align[i] === "l") || (!b.align && i === 0) ? AlignmentType.LEFT : AlignmentType.CENTER,
      spacing: { line: 240, before: 0, after: 0 },
      children: runs(content, { size: fs, bold: isHead }),
    })],
  });
  const rows = [];
  b.header.forEach((hr, hi) => {
    rows.push(new TableRow({ tableHeader: true, cantSplit: true, children: hr.map((c, i) => {
      const cell = mkCell(c, i, true, hi === b.header.length - 1, false, false);
      return cell;
    }) }));
  });
  b.rows.forEach((r, ri) => {
    const shade = (b.shade_rows || []).includes(ri);
    rows.push(new TableRow({ cantSplit: true, children: r.map((c, i) =>
      mkCell(c, i, false, false, ri === b.rows.length - 1, shade)) }));
  });
  // top border on first header row
  const t = new Table({
    columnWidths: widths, width: { size: total, type: WidthType.DXA },
    layout: TableLayoutType.FIXED,
    borders: { top: thin, bottom: thin, left: none, right: none, insideHorizontal: none, insideVertical: none },
    rows,
  });
  return t;
}

function figure(b, maxWidthIn) {
  const buf = fs.readFileSync(b.path);
  // read PNG size
  const w = buf.readUInt32BE(16), h = buf.readUInt32BE(20);
  const widthIn = Math.min(b.width_in || maxWidthIn, maxWidthIn);
  const wpx = Math.round(widthIn * 96), hpx = Math.round(wpx * h / w);
  return new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 120, after: 60 }, keepNext: true,
    children: [new ImageRun({ type: "png", data: buf, transformation: { width: wpx, height: hpx },
      altText: { title: b.alt || "figure", description: b.alt || "figure", name: path.basename(b.path) } })] });
}

// split blocks into sections (portrait/landscape)
const sections = [];
let cur = { landscape: false, children: [] };
const PORTRAIT_W = 12240 - 2 * 1440, LANDSCAPE_W = 15840 - 2 * 1440;
for (const b of doc.blocks) {
  if (b.type === "section") {
    if (cur.children.length) sections.push(cur);
    cur = { landscape: !!b.landscape, children: [] };
    continue;
  }
  const pw = cur.landscape ? LANDSCAPE_W : PORTRAIT_W;
  if (b.type === "h1" || b.type === "h2" || b.type === "h3") cur.children.push(heading(b));
  else if (b.type === "para") cur.children.push(para(b));
  else if (b.type === "table") cur.children.push(table(b, pw));
  else if (b.type === "figure") cur.children.push(figure(b, (pw / 1440)));
  else if (b.type === "pagebreak") cur.children.push(new Paragraph({ children: [new PageBreak()] }));
  else if (b.type === "spacer") cur.children.push(new Paragraph({ spacing: { after: b.after || 120 }, children: [] }));
}
sections.push(cur);

const footer = new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER,
  children: [new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 20 })] })] });
const header = new Header({ children: [new Paragraph({ alignment: AlignmentType.RIGHT,
  children: [new TextRun({ text: doc.meta.running_head || "", font: FONT, size: 18, color: "595959" })] })] });

const document = new Document({
  creator: doc.meta.creator || "",
  title: doc.meta.title || "",
  styles: {
    default: { document: { run: { font: FONT, size: 2 * BODY_PT } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 28, bold: true, font: FONT, color: "000000" }, paragraph: { spacing: { before: 360, after: 120 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 25, bold: true, font: FONT, color: "000000" }, paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 1 } },
      { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, italics: true, font: FONT, color: "000000" }, paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 2 } },
    ],
  },
  sections: sections.map((s, i) => ({
    properties: {
      page: {
        size: s.landscape
          ? { width: 12240, height: 15840, orientation: PageOrientation.LANDSCAPE }
          : { width: 12240, height: 15840 },
        margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 },
      },
      lineNumbers: (!SUPP && !s.landscape && doc.meta.line_numbers) ? { countBy: 1, restart: LineNumberRestartFormat.CONTINUOUS, distance: 360 } : undefined,
    },
    headers: { default: header },
    footers: { default: footer },
    children: s.children,
  })),
});

Packer.toBuffer(document).then((buf) => {
  fs.writeFileSync(outFile, buf);
  console.log(`wrote ${outFile}`);
});
