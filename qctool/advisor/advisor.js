function choosePipeline(settings, depth, diameter) {
  if (!Number.isFinite(depth) || depth <= 0 ||
      (diameter !== null && (!Number.isFinite(diameter) || diameter <= 0))) {
    throw new RangeError("Enter a positive mean UMI count and, if known, a positive damage diameter.");
  }

  if (diameter === null) {
    const {bins, min_molecules_per_bin: target} = settings.mad;
    const bin = bins.find(size => depth * (size / 8) ** 2 >= target) ?? bins[bins.length - 1];
    return {mad: bin, molecules: depth * (bin / 8) ** 2};
  }

  const band = settings.classes.find(c =>
    diameter >= c.min_um && (c.max_um === null || diameter < c.max_um));
  const high = depth >= settings.depth_cut_umi_per_bin;
  const tier = high ? "hi" : "lo";
  return {
    spot: band[tier],
    region: diameter >= settings.spatialartifacts.min_um ? settings.spatialartifacts[tier] : null,
    referenceDiameter: [32, 64, 256, 512][settings.classes.indexOf(band)],
    referenceDepth: high ? 100 : 12,
    band,
  };
}

function sectionSummary({bins, umis, mean}) {
  const count = new Intl.NumberFormat("en", {maximumSignificantDigits: 12});
  return `<p class="small"><strong>${count.format(mean)} mean UMIs per 8 µm bin</strong> · ${count.format(bins)} bins · ${count.format(umis)} total UMIs</p>`;
}

function showRecommendation(data, depth, diameter, edge) {
  const {settings, mad, lqr} = data;
  const plan = choosePipeline(settings, depth, diameter);
  const binText = bin => bin === 8 ? "8 µm native bins" : `${bin} µm bins with A&P`;
  const percent = value => `${value.toFixed(1)}%`;
  const target = settings.mad.min_molecules_per_bin;
  let steps, tradeoff, heading, evidence, classText;

  if (plan.mad) {
    heading = "Section-wide MAD threshold";
    const count = new Intl.NumberFormat("en", {maximumFractionDigits: 1}).format(plan.molecules);
    const why = plan.molecules >= target
      ? `This is the smallest supported bin expected to hold about ${target} molecules or more at this depth (roughly ${count} UMIs per bin).`
      : `This is the largest tested bin. At this depth it holds roughly ${count} UMIs per bin, below the ${target}-molecule target.`;
    steps = `<li><strong>MAD threshold · ${binText(plan.mad)}</strong><span>${why}</span></li>`;
    tradeoff = "At low depth, pooling counts stabilises the false-positive rate. Coarser bins can dilute small areas of signal loss and remove more intact tissue along the edge. The guideline also offers 64 µm bins when a stable false-positive rate across sections is the priority.";
    const rows = [100, 50, 25, 12].map(d => {
      const [n, rates] = mad[d];
      return `<tr><th scope="row">${d}</th>${rates.map((v, i) => `<td${settings.mad.bins[i] === plan.mad ? ' data-selected' : ''}>${percent(v)}</td>`).join("")}<td>${n}</td></tr>`;
    }).join("");
    evidence = `<p>Median share of tissue removed on the uninjected condition in probe-based sections (Figure 3 and Supplementary Table S1). The selected bin is highlighted.</p>
      <div class="scroll" role="region" aria-label="MAD benchmark rates" tabindex="0"><table>
      <thead><tr><th scope="col">UMIs per 8 µm bin</th>${settings.mad.bins.map(b => `<th scope="col">${b} µm</th>`).join("")}<th scope="col">Sections</th></tr></thead>
      <tbody>${rows}</tbody></table></div>`;
  } else {
    const {spot, region, referenceDiameter, referenceDepth, band} = plan;
    const lower = band.min_um > 0 ? `≥ ${band.min_um} µm` : "";
    const upper = band.max_um === null ? "" : `&lt; ${band.max_um} µm`;
    const depthTier = depth >= settings.depth_cut_umi_per_bin ? "≥" : "&lt;";
    classText = `Damage class ${[lower, upper].filter(Boolean).join(" and ")} · ${depthTier} ${settings.depth_cut_umi_per_bin} UMIs per 8 µm bin`;
    heading = region ? "SpotSweeper + SpatialArtifacts" : "SpotSweeper";
    const why = spot.k > 48
      ? "A wider neighbourhood covers enough tissue while keeping native bins, so small damage is less diluted."
      : "Aggregation gives each bin more molecules and lets its neighbourhood span more tissue than the expected damage.";
    steps = `<li><strong>SpotSweeper · ${binText(spot.bin)} · ${spot.k} neighbours</strong><span>${why}</span></li>`;
    if (region) {
      steps += `<li><strong>SpatialArtifacts · ${binText(region)}</strong><span>At this diameter, the region test can recover broad signal loss that SpotSweeper may miss. Remove the bins flagged by either test.</span></li>`;
    }
    tradeoff = region
      ? "Combining the tests can recover more damage and can also remove more uninjected tissue than either test alone. SpatialArtifacts looks for low counts. SpotSweeper provides the mitochondrial check."
      : "Keeping the test local limits how much uninjected tissue is removed. At low depth, recall of small damage can still be limited.";
    const observations = lqr[`${referenceDiameter}/${referenceDepth}`];
    const labels = {spot: "SpotSweeper", region: "SpatialArtifacts", union: "Either test"};
    const rows = Object.entries(observations).map(([arm, [signal, mito, fpr, n]]) =>
      `<tr${arm === "union" || (!region && arm === "spot") ? ' data-selected' : ''}><th scope="row">${labels[arm]}</th>
      <td>${percent(signal)}</td><td>${percent(mito)}</td><td>${percent(fpr)}</td><td>${n}</td></tr>`
    ).join("");
    const scope = referenceDiameter === 512
      ? "The 512 µm benchmark placed LQRs inside the tissue. Supplementary Table S9 reports SpotSweeper alone at this diameter."
      : "These rows use LQRs injected inside the tissue.";
    evidence = `<p>Probe-based section medians for ${referenceDiameter} µm synthetic low-quality regions (LQRs) at ${referenceDepth} UMIs per 8 µm bin (Supplementary Table S9). ${scope}</p>
      <div class="scroll" role="region" aria-label="Recall and false-positive rates" tabindex="0"><table>
      <thead><tr><th scope="col">Test</th><th scope="col">Signal-loss recall</th><th scope="col">Mito. elevation recall</th><th scope="col">Uninjected tissue removed</th><th scope="col">Sections</th></tr></thead>
      <tbody>${rows}</tbody></table></div>`;
  }

  const usesAP = (plan.mad ?? plan.spot.bin) > 8 || plan.region > 8;
  const nominalDepth = Math.round(depth);
  const depthNote = nominalDepth < 12
    ? "This depth is below the benchmark minimum of 12 UMIs per 8 µm bin."
    : nominalDepth > 100 ? "The recall benchmark used depths up to 100 UMIs per 8 µm bin." : "";
  const review = edge
    ? "Inspect removals along the tissue boundary. Neighbourhood recall can drop at the edge, and coarser bins can flag intact boundary tissue. Compare marker genes in removed and matched kept bins."
    : "Inspect where the removed bins lie and compare their marker genes with matched kept bins.";
  const lowRegion = plan.region && depth < settings.depth_cut_umi_per_bin
    ? " At low depth, check whether SpatialArtifacts flagged any bins." : "";

  return `<div class="card">
      <p class="eyebrow">Recommended pipeline</p>
      <h2 id="result-title" tabindex="-1">${heading}</h2>
      ${classText ? `<p class="small">${classText}</p>` : ""}
      <ol>${steps}</ol>
      ${usesAP ? "<p class='small'>With aggregation and projection (A&P), sum native bins into the larger bins, run QC there, then project each flag back to the 8 µm bins inside.</p>" : ""}
      <p class="note"><strong>After QC</strong><br>${review}${lowRegion}</p>
    </div>
    <div class="card">
      <h3>What you gain and what you trade off</h3>
      <p>${tradeoff}</p>
      ${evidence}
      <p class="small">Uninjected tissue removed is the benchmark false-positive rate. Use the medians as reference values for the displayed benchmark conditions. ${depthNote}</p>
    </div>`;
}

module.exports = {choosePipeline, sectionSummary, showRecommendation};
