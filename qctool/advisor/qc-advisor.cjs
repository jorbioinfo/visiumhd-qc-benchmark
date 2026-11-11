#!/usr/bin/env node
const fs = require("node:fs");
const path = require("node:path");
const {parseArgs} = require("node:util");
const {createGunzip} = require("node:zlib");
const benchmark = require("./benchmark.js");
const {countBarcodes, readUmis} = require("./matrix.js");
const {sectionSummary, showRecommendation} = require("./advisor.js");

function mexFile(directory, name) {
  for (const suffix of [".gz", ""]) {
    const file = path.join(directory, name + suffix);
    if (fs.statSync(file, {throwIfNoEntry: false})?.isFile()) return file;
  }
  throw new Error(`Missing ${name}.gz or ${name} in ${directory}`);
}

async function* chunks(file) {
  const source = fs.createReadStream(file);
  const stream = file.endsWith(".gz") ? createGunzip() : source;
  if (stream !== source) {
    source.on("error", error => stream.destroy(error));
    source.pipe(stream);
  }
  try {
    yield* stream;
  } finally {
    source.destroy();
  }
}

async function main() {
  const {values: args} = parseArgs({options: {
    input: {type: "string"},
    "damage-um": {type: "string"},
    edge: {type: "boolean"},
    out: {type: "string", default: "qc-advisor-report.html"},
    help: {type: "boolean", short: "h"},
  }});
  if (args.help) {
    console.log(`Visium HD QC Advisor (Node.js 16.20+)

Usage
  node qc-advisor.cjs --input DIRECTORY [--damage-um MICROMETRES] [--edge] [--out FILE]

DIRECTORY is the gene-expression MEX directory from Space Ranger:
  outs/binned_outputs/square_008um/filtered_feature_bc_matrix/

Use matrix.mtx.gz and barcodes.tsv.gz, or their decompressed counterparts.
Leave --damage-um out when damage size is unknown. The HTML report contains
the recommended QC pipeline and benchmark trade-offs. The default output
is qc-advisor-report.html. Choose a new path if that file already exists.`);
    return;
  }
  if (!args.input) throw new Error("Supply --input with the 8 µm filtered MEX directory. Run --help for usage.");
  const diameter = args["damage-um"] === undefined ? null : Number(args["damage-um"]);
  if (diameter !== null && (!Number.isFinite(diameter) || diameter <= 0)) {
    throw new Error("--damage-um must be a positive number.");
  }

  const matrix = mexFile(args.input, "matrix.mtx");
  const bins = await countBarcodes(chunks(mexFile(args.input, "barcodes.tsv")));
  const summary = await readUmis(chunks(matrix), bins);
  const body = showRecommendation(benchmark, summary.mean, diameter, args.edge);
  const report = `<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Visium HD QC Advisor report</title>
<style>
body {font: 1rem/1.55 system-ui, sans-serif; color: #172f31; background: #f3f6f2; max-width: 65rem; padding: 1.2rem; margin: auto}
.card {background: white; border: 1px solid #cddbd6; border-radius: .6rem; padding: 1rem 1.3rem; margin: 1rem 0}
li {margin: .6rem 0} li span {display: block} .small {font-size: .9rem}
.eyebrow {color: #135b51; font-weight: bold} .note {border-left: .25rem solid #135b51; padding-left: .8rem}
.scroll {overflow-x: auto} table {border-collapse: collapse; width: 100%} th,td {padding: .45rem; text-align: right; border-bottom: 1px solid #cddbd6}
th:first-child {text-align: left} [data-selected] {background: #eaf3ef; font-weight: bold}
</style></head><body>
<h1>Visium HD QC Advisor report</h1>
${sectionSummary(summary)}
<p>Damage diameter ${diameter === null ? "unknown" : `${diameter} µm`}${args.edge ? " · Suspected tissue-edge damage" : ""}</p>
${body}
<p class="small">The pipeline follows the Visium HD benchmark guideline in Figure 6.</p>
</body></html>\n`;
  fs.writeFileSync(args.out, report, {flag: "wx"});
  console.log(`Wrote ${args.out} (${summary.mean.toPrecision(8)} mean UMIs per 8 µm bin).`);
}

main().catch(error => {
  console.error(error.code === "EEXIST" ? "The output report exists. Choose a new --out path." : error.message);
  process.exitCode = 1;
});
