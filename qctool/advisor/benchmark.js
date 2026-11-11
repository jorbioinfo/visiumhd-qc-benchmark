// Settings from Figure 6. Probe-based medians from Supplementary Table S9 and guideline_numbers.json.
const benchmark = {
  settings: {
    "generated": "2026-11-09 10:46",
    "depth_cut_umi_per_bin": 100,
    "classes": [{"label": "< 45 µm", "min_um": 0, "max_um": 45, "hi": {"bin": 8, "k": 192, "prov": false}, "lo": {"bin": 16, "k": 48, "prov": false}}, {"label": "45–130 µm", "min_um": 45, "max_um": 130, "hi": {"bin": 8, "k": 768, "prov": false}, "lo": {"bin": 32, "k": 48, "prov": false}}, {"label": "130–360 µm", "min_um": 130, "max_um": 360, "hi": {"bin": 64, "k": 48, "prov": false}, "lo": {"bin": 64, "k": 48, "prov": false}}, {"label": "> 360 µm", "min_um": 360, "max_um": null, "hi": {"bin": 128, "k": 48, "prov": false}, "lo": {"bin": 128, "k": 48, "prov": false}}],
    "mad": {"hi": 8, "lo": 32, "min_molecules_per_bin": 100, "bins": [8, 16, 32, 64]},
    "spatialartifacts": {"hi": 8, "lo": 16, "min_um": 45},
  },
  mad: {
    "100": [28, [7.3, 7.5, 8.1, 9.0]],
    "50": [33, [9.4, 6.8, 7.5, 7.9]],
    "25": [33, [21.6, 7.0, 7.4, 7.5]],
    "12": [33, [33.2, 9.5, 6.5, 7.0]],
  },
  lqr: {
    "32/100": {"spot": [84.5, 85.1, 0.5, 28]},
    "32/12": {"spot": [33.5, 46.8, 0.4, 33]},
    "64/100": {"spot": [81.4, 81.1, 0.7, 28], "region": [98.5, 0.4, 1.7, 28], "union": [99.5, 81.9, 3.0, 28]},
    "64/12": {"spot": [35.9, 78.0, 0.8, 33], "region": [67.2, 0.6, 1.6, 33], "union": [70.8, 81.5, 2.8, 33]},
    "256/100": {"spot": [53.9, 83.2, 1.4, 28], "region": [99.7, 0.5, 1.7, 28], "union": [99.7, 83.4, 3.2, 28]},
    "256/12": {"spot": [45.8, 72.4, 1.0, 33], "region": [91.0, 0.8, 1.6, 33], "union": [92.1, 74.5, 2.9, 33]},
    "512/100": {"spot": [59.7, 80.0, 2.1, 27]},
    "512/12": {"spot": [51.7, 84.3, 1.5, 32]},
  },
};

module.exports = benchmark;
