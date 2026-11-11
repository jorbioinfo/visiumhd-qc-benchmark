async function readLines(chunks, visit) {
  const decoder = new TextDecoder("utf-8", {fatal: true});
  let tail = "";
  // ponytail: keep only one partial line and one decoded chunk. A larger metadata
  // line would need a streaming comment skipper instead of raising this limit.
  const maxLine = 1024 * 1024;
  for await (const chunk of chunks) {
    const text = tail + decoder.decode(chunk, {stream: true});
    let start = 0;
    for (let end = text.indexOf("\n"); end !== -1; end = text.indexOf("\n", start)) {
      if (end - start > maxLine) throw new Error("A line in the input file is too long.");
      visit(text.slice(start, end).trim());
      start = end + 1;
    }
    tail = text.slice(start);
    if (tail.length > maxLine) throw new Error("A line in the input file is too long.");
  }
  tail += decoder.decode();
  if (tail.trim()) visit(tail.trim());
}

async function countBarcodes(chunks) {
  let bins = 0;
  await readLines(chunks, line => {
    if (!/^s_008um_\d+_\d+-\d+$/.test(line)) {
      throw new Error("Use the 8 µm barcodes.tsv file from square_008um/filtered_feature_bc_matrix.");
    }
    bins += 1;
  });
  if (!bins) throw new Error("The barcodes file is empty.");
  return bins;
}

async function readUmis(chunks, bins) {
  let hasHeader = false;
  let features, expected, entries = 0, umis = 0;
  await readLines(chunks, line => {
    if (!hasHeader) {
      if (!/^%%MatrixMarket\s+matrix\s+coordinate\s+integer\s+general$/i.test(line)) {
        throw new Error("Use a Space Ranger Matrix Market file (matrix.mtx or matrix.mtx.gz).");
      }
      hasHeader = true;
      return;
    }
    if (!line || line.startsWith("%")) return;
    const fields = /^(\d+)[ \t]+(\d+)[ \t]+(\d+)$/.exec(line);
    if (!fields) throw new Error("The matrix contains an invalid count or dimension.");
    const row = Number(fields[1]), column = Number(fields[2]), count = Number(fields[3]);
    if (features === undefined) {
      if (column !== bins) throw new Error("The matrix and barcodes file contain different numbers of bins.");
      if (!Number.isSafeInteger(row) || row < 1 || !Number.isSafeInteger(count)) {
        throw new Error("The matrix dimensions are invalid.");
      }
      features = row;
      expected = count;
      return;
    }
    if (row < 1 || row > features || column < 1 || column > bins || ++entries > expected) {
      throw new Error("The matrix entries do not match its declared dimensions.");
    }
    umis += count;
  });
  if (expected === undefined || entries !== expected) {
    throw new Error("The matrix is incomplete or its entry count is wrong.");
  }
  if (!Number.isSafeInteger(umis) || umis <= 0) throw new Error("The matrix must contain a positive, valid UMI total.");
  return {bins, umis, mean: umis / bins};
}

module.exports = {countBarcodes, readUmis};
