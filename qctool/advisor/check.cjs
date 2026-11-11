// Run with node check.cjs from any directory.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const os = require("node:os");
const {gzipSync} = require("node:zlib");
const {spawnSync} = require("node:child_process");
const data = require("./benchmark.js");
const {choosePipeline, showRecommendation} = require("./advisor.js");
const {countBarcodes, readUmis} = require("./matrix.js");

const choose = (depth, diameter) => choosePipeline(data.settings, depth, diameter);
for (const [depth, bin] of [[100, 8], [99.9, 16], [25, 16], [24.9, 32], [12, 32], [2, 64], [1, 64]]) {
  assert.equal(choose(depth, null).mad, bin);
}
assert.equal(choose(1, null).molecules, 64);

for (const [depth, diameter, bin, neighbours, region] of [
  [100, 44.9, 8, 192, null], [99.9, 44.9, 16, 48, null],
  [100, 45, 8, 768, 8], [99.9, 45, 32, 48, 16],
  [100, 130, 64, 48, 8], [12, 130, 64, 48, 16],
  [100, 360, 128, 48, 8], [12, 512, 128, 48, 16],
]) {
  const plan = choose(depth, diameter);
  assert.deepEqual([plan.spot.bin, plan.spot.k, plan.region], [bin, neighbours, region]);
}
for (const [depth, diameter] of [[0, null], [Infinity, 64], [NaN, 64], [100, 0], [100, -1], [100, NaN]]) {
  assert.throws(() => choose(depth, diameter), /Enter a positive/);
}

const render = (depth, diameter, edge = false) => showRecommendation(data, depth, diameter, edge);
assert.match(render(100, 64), /Either test.*99\.5%.*81\.9%.*3\.0%/s);
assert.match(render(12, 64, true), /tissue boundary.*check whether SpatialArtifacts/s);
assert.match(render(100, 360), /Damage class ≥ 360 µm/);
assert.match(render(12, 512), /512 µm.*reports SpotSweeper alone/s);
assert.match(render(1, null), /below the 100-molecule target.*below the benchmark minimum/s);
assert.doesNotMatch(render(11.999874, 64), /below the benchmark minimum/);

async function checkFiles() {
  const barcodeText = "s_008um_00001_00001-1\ns_008um_00001_00002-1";
  const mtx = "%%MatrixMarket matrix coordinate integer general\r\n% test\r\n2 2 3\r\n1 1 80\r\n2 1 50\r\n2 2 70";
  const chunks = (...parts) => parts.map(part => Buffer.from(part));
  const bins = await countBarcodes(chunks(barcodeText.slice(0, 12), barcodeText.slice(12)));
  assert.deepEqual(await readUmis(chunks(mtx.slice(0, 37), mtx.slice(37, 69), mtx.slice(69)), bins),
    {bins: 2, umis: 200, mean: 100});
  await assert.rejects(countBarcodes(chunks(barcodeText.replace("008um", "016um"))), /Use the 8 µm/);
  await assert.rejects(readUmis(chunks(mtx), 3), /different numbers of bins/);
  await assert.rejects(readUmis(chunks(mtx.replace(/2 2 70$/, "")), bins), /incomplete/);
  await assert.rejects(readUmis(chunks(mtx.replace("2 2 70", "2 2 -70")), bins), /invalid count/);

  const directory = fs.mkdtempSync(path.join(os.tmpdir(), "qc-advisor-check-"));
  try {
    fs.writeFileSync(path.join(directory, "matrix.mtx.gz"), gzipSync(mtx));
    fs.writeFileSync(path.join(directory, "barcodes.tsv.gz"), gzipSync(barcodeText));
    const output = path.join(directory, "report.html");
    const args = [path.join(__dirname, "qc-advisor.cjs"), "--input", directory,
      "--damage-um", "64", "--edge", "--out", output];
    const run = spawnSync(process.execPath, args, {encoding: "utf8"});
    assert.equal(run.status, 0, run.stderr);
    const report = fs.readFileSync(output, "utf8");
    assert.match(report, /100 mean UMIs per 8 µm bin.*2 bins.*200 total UMIs/s);
    assert.match(report, /SpotSweeper \+ SpatialArtifacts.*768 neighbours.*Either test/s);
    assert.equal(spawnSync(process.execPath, args).status, 1);
    assert.equal(fs.readFileSync(output, "utf8"), report);
  } finally {
    fs.rmSync(directory, {recursive: true, force: true});
  }
}

checkFiles().then(() => console.log("QC advisor checks passed."));
