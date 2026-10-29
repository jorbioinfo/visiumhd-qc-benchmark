import math
import os

OUT = os.path.join(os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "methods_diagrams")
W = 800
SANS = "Arial, Helvetica, 'Liberation Sans', sans-serif"
MONO = "Menlo, Consolas, 'Liberation Mono', monospace"
C = dict(paper="#ffffff", paper2="#f4f4f2", ink="#0b0b0b", muted="#52514e", soft="#6b6b6b",
         rule="rgba(11,11,11,0.12)", rule_solid="#d8d7d2", accent="#c8501a",
         accent_tint="rgba(235,104,52,0.10)", link="#2368c0", link_tint="rgba(42,120,214,0.08)",
         grey_tint="rgba(11,11,11,0.04)")
KIND = {"backend": ("#ffffff", C["ink"], 1),
        "focal": (C["accent_tint"], C["accent"], 1.2),
        "native": (C["link_tint"], C["link"], 1.2),
        "store": ("rgba(11,11,11,0.05)", C["muted"], 1),
        "input": ("rgba(82,81,78,0.10)", C["soft"], 1)}
NAME, SUB, LINE = 16, 12, 16
DMG, DMG_TINT = C["ink"], "rgba(11,11,11,0.20)"
TISSUE = "rgba(82,81,78,0.22)"
SERIES = ("#7c8f6f", "#5e7a9b", "#b8915a")

def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def box(x, y, w, h, name, subs=(), kind="backend", node=""):
    fill, stroke, sw = KIND[kind]
    n = len(subs)
    block = NAME + (8 + LINE * n if n else 0)
    top = y + (h - block) / 2
    cx = x + w / 2
    out = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" fill="{C["paper"]}"/>',
           f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" fill="{fill}" '
           f'stroke="{stroke}" stroke-width="{sw}" data-node="{node}"/>',
           f'<text x="{cx}" y="{top + 13:.0f}" fill="{C["ink"]}" font-size="{NAME}" '
           f'font-weight="700" font-family="{SANS}" text-anchor="middle" '
           f'data-in="{node}">{esc(name)}</text>']
    for i, s in enumerate(subs):
        fam, txt = (MONO, s[5:]) if s.startswith("mono:") else (SANS, s)
        out.append(f'<text x="{cx}" y="{top + NAME + 8 + LINE * i + 10:.0f}" fill="{C["muted"]}" '
                   f'font-size="{SUB}" font-family="{fam}" text-anchor="middle" '
                   f'data-in="{node}">{esc(txt)}</text>')
    return "\n".join(out)

def path(points, color=None, marker="arrow", dash=None, r=8):
    color = color or C["muted"]
    d = [f"M {points[0][0]} {points[0][1]}"]
    for i in range(1, len(points) - 1):
        (x0, y0), (x1, y1), (x2, y2) = points[i - 1], points[i], points[i + 1]
        ux, uy = (x1 - x0) / (abs(x1 - x0) + abs(y1 - y0)), (y1 - y0) / (abs(x1 - x0) + abs(y1 - y0))
        vx, vy = (x2 - x1) / (abs(x2 - x1) + abs(y2 - y1)), (y2 - y1) / (abs(x2 - x1) + abs(y2 - y1))
        sweep = 1 if ux * vy - uy * vx > 0 else 0
        d.append(f"L {x1 - ux * r:g} {y1 - uy * r:g}")
        d.append(f"A {r} {r} 0 0 {sweep} {x1 + vx * r:g} {y1 + vy * r:g}")
    d.append(f"L {points[-1][0]} {points[-1][1]}")
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    mk = f' marker-end="url(#{marker})"' if marker else ""
    return (f'<path d="{" ".join(d)}" fill="none" stroke="{color}" stroke-width="1.2"'
            f'{extra}{mk} data-connector="1"/>')

def text(x, y, s, size=SUB, color=None, weight=400, anchor="start", italic=False, fam=SANS,
         node=""):
    st = ' font-style="italic"' if italic else ""
    dn = f' data-in="{node}"' if node else ""
    return (f'<text x="{x:g}" y="{y:g}" fill="{color or C["ink"]}" font-size="{size}" '
            f'font-weight="{weight}" font-family="{fam}" text-anchor="{anchor}"{st}{dn}>{esc(s)}</text>')

LETTER = 9 / 72 * 2.54 / 16 * W

def letter(x, y, s):
    return (f'<text x="{x:g}" y="{y:g}" fill="{C["ink"]}" font-size="{LETTER:.1f}" '
            f'font-weight="700" font-stretch="normal" font-family="DejaVu Sans">{esc(s)}</text>')

def rect(x, y, w, h, fill="none", stroke="none", sw=1, rx=0, extra=""):
    return (f'<rect x="{x:g}" y="{y:g}" width="{w:g}" height="{h:g}" rx="{rx}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{sw}"{extra}/>')

def grid(x, y, nr, nc, s, fill=lambda i, j: C["paper"], stroke=None, sw=0.8):
    return [rect(x + j * s, y + i * s, s, s, fill(i, j), stroke or C["rule_solid"], sw)
            for i in range(nr) for j in range(nc)]

def cross(cx, cy, h, color=None, sw=1.2):
    color = color or C["ink"]
    return (f'<path d="M {cx - h:g} {cy - h:g} L {cx + h:g} {cy + h:g} M {cx - h:g} {cy + h:g} '
            f'L {cx + h:g} {cy - h:g}" stroke="{color}" stroke-width="{sw}" stroke-linecap="round"/>')

def mini_arrow(x1, x2, y, color=None, marker="arrow"):
    return (f'<path d="M {x1:g} {y:g} L {x2:g} {y:g}" fill="none" stroke="{color or C["muted"]}" '
            f'stroke-width="1.2" marker-end="url(#{marker})"/>')

def shell(x, y, w, h, name, kind="backend", node="", name_x=None):
    fill, stroke, sw = KIND[kind]
    anchor, nx = ("middle", x + w / 2) if name_x is None else ("start", name_x)
    return [rect(x, y, w, h, C["paper"], rx=6),
            f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{sw}" data-node="{node}"/>',
            text(nx, y + 26, name, NAME, weight=700, anchor=anchor, node=node)][:3 if name else 2]

def markers():
    return "\n".join(
        f'<marker id="{i}" markerWidth="8" markerHeight="6" refX="7" refY="3" orient="auto">'
        f'<polygon points="0 0, 8 3, 0 6" fill="{c}"/></marker>'
        for i, c in (("arrow", C["muted"]), ("arrow-accent", C["accent"]), ("arrow-link", C["link"])))

def page(slug, eyebrow, title, desc, h, body, caption, defs=""):
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{esc(title)}</title>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    :root {{ --paper: {C["paper"]}; --ink: {C["ink"]}; --muted: {C["muted"]}; --sans: {SANS}; }}
    body {{ font-family: var(--sans); background: var(--paper); color: var(--ink);
           min-height: 100vh; display: flex; align-items: center; justify-content: center;
           padding: 3rem 2rem; }}
    .frame {{ max-width: 860px; width: 100%; }}
    .diagram-container {{ width: 100%; overflow-x: auto; }}
    .eyebrow {{ font-size: 0.72rem; font-weight: 700; letter-spacing: 0.08em;
               text-transform: uppercase; color: var(--muted); margin-bottom: 0.5rem; }}
    h1 {{ font-size: 1.6rem; font-weight: 400; line-height: 1.2; margin-bottom: 1.5rem;
          text-wrap: balance; }}
    .caption {{ font-size: 0.9rem; line-height: 1.5; color: var(--muted); max-width: 65ch;
               margin-top: 1.25rem; }}
    .caption b {{ color: var(--ink); }}
    svg {{ width: 100%; min-width: {W}px; display: block; }}
    @media print {{
      .diagram-container {{ overflow-x: visible; }}
      svg {{ min-width: 0; }}
    }}
  </style>
</head>
<body>
  <div class="frame">
    <p class="eyebrow">{esc(eyebrow)}</p>
    <h1>{esc(title)}</h1>
    <div class="diagram-container">
      <svg viewBox="0 0 {W} {h}" xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="{slug}-title {slug}-desc">
        <title id="{slug}-title">{esc(title)}</title>
        <desc id="{slug}-desc">{esc(desc)}</desc>
        <defs>
{markers()}
{defs}
        </defs>
        <rect width="100%" height="100%" fill="{C["paper"]}"/>
{body}
      </svg>
    </div>
    <p class="caption">{caption}</p>
  </div>
</body>
</html>
"""

def overview():
    b = []
    b += [path([(256, 114), (292, 114)]), path([(508, 114), (544, 114)]),
          path([(652, 188), (652, 220)]), path([(400, 368), (400, 400)]),
          path([(148, 456), (148, 488)]), path([(400, 456), (400, 488)]),
          path([(652, 456), (652, 488)]),
          path([(148, 640), (148, 656), (340, 656), (340, 672)]),
          path([(400, 640), (400, 672)]),
          path([(652, 640), (652, 656), (460, 656), (460, 672)])]

    b += shell(40, 40, 216, 148, "Visium HD sections", "input", "sections")
    blobs = ("M {x} {y} m 12 16 c 6 -10 24 -8 28 2 c 6 12 0 26 -12 24 c -12 4 -22 -6 -20 -16 "
             "c 1 -5 1 -7 4 -10 z",
             "M {x} {y} m 10 24 c 0 -14 20 -18 30 -10 c 8 8 4 26 -10 30 c -12 2 -20 -6 -20 -20 z")
    for k, (sx, lab) in enumerate(((74, "probe-based"), (170, "3′ polyA"))):
        b.append(rect(sx, 80, 52, 52, C["paper"], C["soft"], 1, rx=2))
        b.append(f'<path d="{blobs[k].format(x=sx, y=80)}" fill="{TISSUE}" stroke="{C["soft"]}" '
                 f'stroke-width="0.8"/>')
        b.append(text(sx + 26, 148, lab, color=C["muted"], anchor="middle", node="sections"))
    b.append(text(148, 172, "public, human and mouse", color=C["muted"], anchor="middle",
                  node="sections"))

    b += shell(292, 40, 216, 148, "Depth control", node="depth")
    for k, dep in enumerate((100, 50, 25, 12)):
        h = round(52 * dep / 100)
        b.append(rect(328 + 40 * k, 132 - h, 24, h, C["muted"]))
        b.append(text(340 + 40 * k, 148, str(dep), color=C["muted"], anchor="middle", node="depth"))
    b.append(text(400, 172, "UMIs per bin", color=C["muted"], anchor="middle", node="depth"))

    b += shell(544, 40, 216, 148, "Synthetic LQRs", node="lqr")
    for cx, r, lab in ((588, 3.25, "32"), (616, 6.5, "64"), (668, 26, "256 µm")):
        b.append(f'<circle cx="{cx}" cy="106" r="{r}" fill="{DMG_TINT}" stroke="{DMG}" '
                 f'stroke-width="1"/>')
        b.append(text(cx, 148, lab, color=C["muted"], anchor="middle", node="lqr"))
    b.append(rect(716, 96, 20, 20, TISSUE, rx=3))
    b.append(text(726, 148, "none", color=C["muted"], anchor="middle", node="lqr"))
    b.append(text(652, 172, "signal loss or mito. elevation", color=C["muted"], anchor="middle",
                  node="lqr"))

    b += shell(40, 220, 720, 148, "Each QC method, run two ways", node="qc")
    b.append(text(400, 266, "MAD, SpotSweeper, miQC, SpatialArtifacts", color=C["muted"],
                  anchor="middle", node="qc"))
    b.append(rect(56, 280, 328, 76, C["link_tint"], rx=6))
    b.append(rect(416, 280, 328, 76, C["accent_tint"], rx=6))
    b.append(text(68, 298, "Native", weight=700, node="qc"))
    b.append(text(428, 298, "A&D", weight=700, node="qc"))
    nat = {(1, 1), (1, 2), (3, 3)}
    b += grid(156, 290, 5, 5, 8)
    b.append(mini_arrow(200, 236, 310, C["link"], "arrow-link"))
    b.append(text(218, 302, "QC", color=C["muted"], anchor="middle", node="qc"))
    b += grid(240, 290, 5, 5, 8, lambda i, j: C["link"] if (i, j) in nat else C["paper"])
    b.append(text(176, 346, "8 µm bins", color=C["muted"], anchor="middle", node="qc"))
    b.append(text(260, 346, "8 µm verdict", color=C["muted"], anchor="middle", node="qc"))
    b += grid(472, 290, 4, 4, 10)
    b.append(mini_arrow(516, 556, 310, C["accent"], "arrow-accent"))
    b.append(text(536, 302, "sum", color=C["muted"], anchor="middle", node="qc"))
    b += grid(560, 290, 2, 2, 20, lambda i, j: C["accent"] if (i, j) == (0, 1) else C["paper"])
    b.append(mini_arrow(604, 644, 310, C["accent"], "arrow-accent"))
    b.append(text(624, 302, "project", color=C["muted"], anchor="middle", node="qc"))
    b += grid(648, 290, 4, 4, 10,
              lambda i, j: C["accent"] if i < 2 and j >= 2 else C["paper"])
    b.append(text(492, 346, "8 µm bins", color=C["muted"], anchor="middle", node="qc"))
    b.append(text(580, 346, "coarse QC", color=C["muted"], anchor="middle", node="qc"))
    b.append(text(668, 346, "8 µm verdict", color=C["muted"], anchor="middle", node="qc"))

    b += shell(40, 400, 720, 56, "", "store", "verdict")
    b.append(text(60, 434, "Verdict on the same 8 µm bins", NAME, weight=700, node="verdict"))
    for k in range(14):
        x = 388 + 18 * k
        b.append(rect(x, 421, 14, 14, C["paper"], C["rule_solid"], 0.8))
        if k in (2, 5, 6, 10):
            b.append(cross(x + 7, 428, 4))
    b.append(rect(660, 412, 12, 12, C["paper"], C["rule_solid"], 0.8))
    b.append(cross(666, 418, 3.5))
    b.append(text(680, 422, "removed", color=C["muted"], node="verdict"))
    b.append(rect(660, 432, 12, 12, C["paper"], C["rule_solid"], 0.8))
    b.append(text(680, 442, "kept", color=C["muted"], node="verdict"))

    b += shell(40, 488, 216, 152, "Scoring", node="score")
    inj = {(i, j) for i in range(7) for j in range(7) if (i - 3) ** 2 + (j - 3) ** 2 <= 4}
    missed, fp = {(1, 3), (3, 5), (5, 3)}, {(0, 6), (6, 0)}
    b += grid(113, 524, 7, 7, 10, lambda i, j: DMG_TINT if (i, j) in inj else C["paper"])
    b.append(f'<circle cx="148" cy="559" r="20" fill="none" stroke="{DMG}" stroke-width="1" '
             f'stroke-dasharray="3,2"/>')
    for i, j in (inj - missed) | fp:
        b.append(cross(113 + 10 * j + 5, 524 + 10 * i + 5, 3))
    b.append(text(148, 622, "recall and cost", color=C["muted"], anchor="middle", node="score"))

    b += shell(292, 488, 216, 152, "Marker check", node="marker")
    b.append(rect(356, 548, 24, 30, C["paper"], C["muted"], 1))
    b.append(rect(420, 528, 24, 50, C["muted"]))
    b.append(text(400, 560, "÷", NAME, weight=700, anchor="middle", node="marker"))
    b.append(text(368, 594, "removed", color=C["muted"], anchor="middle", node="marker"))
    b.append(text(432, 594, "kept", color=C["muted"], anchor="middle", node="marker"))

    b += shell(544, 488, 216, 152, "Clustering agreement", node="down")
    m100 = lambda i, j: SERIES[0] if j < 2 else (SERIES[1] if i < 3 else SERIES[2])
    diff = {(2, 2), (3, 1), (4, 3)}
    m12 = lambda i, j: SERIES[2] if (i, j) in diff else m100(i, j)
    b += grid(590, 526, 5, 5, 10, m100, C["paper"], 0.6)
    b += grid(664, 526, 5, 5, 10, m12, C["paper"], 0.6)
    b.append(text(652, 555, "ARI", weight=700, anchor="middle", node="down"))
    b.append(text(615, 592, "100", color=C["muted"], anchor="middle", node="down"))
    b.append(text(689, 592, "12", color=C["muted"], anchor="middle", node="down"))

    b += shell(220, 672, 360, 88, "Section-level inference", node="infer")
    for k in range(9):
        cx = 320 + 20 * k
        hollow = k in (3, 7)
        b.append(f'<circle cx="{cx}" cy="716" r="5" fill="{C["paper"] if hollow else C["ink"]}" '
                 f'stroke="{C["ink"]}" stroke-width="1"/>')

    cap = ("<b>Study design.</b> Public Visium HD sections, probe-based and 3′ polyA, are "
           "thinned to four common depths, and synthetic low-quality regions (LQRs) of three "
           "diameters are added once to the 8 µm counts of each run (one section at one "
           "depth). Every QC method is run natively on the 8 µm grid and through "
           "aggregate-and-deconvolve (A&amp;D), which sums bins, runs the method on the coarse "
           "bins and projects its verdict back to 8 µm, so both routes return a verdict on the "
           "same 8 µm bins. The verdict is judged three ways: against the injected regions "
           "(recall and nominal false-positive rate), by the marker share of removed against "
           "kept bins on undamaged tissue matched for tile and library size, and by the "
           "adjusted Rand index (ARI) between clusterings of one section at 100 and 12 UMIs "
           "per bin. Sections are the unit of inference, and every test is made within one "
           "chemistry.")
    return page("methods-overview", "Materials and Methods · study design",
                "From public sections to a section-level test",
                "Workflow: public Visium HD sections are depth-controlled and given synthetic "
                "damage, every QC method is run natively and through A&D onto the same 8 µm "
                "bins, and the verdict is judged by scoring, a marker check and clustering "
                "agreement before a section-level test.",
                800, "\n".join(b), cap)

def nearest(n, ci, cj, k):
    cells = [(i, j) for i in range(n) for j in range(n) if (i, j) != (ci, cj)]
    cells.sort(key=lambda c: ((c[0] - ci) ** 2 + (c[1] - cj) ** 2, math.atan2(c[0] - ci, c[1] - cj)))
    return set(cells[:k])

def routes():
    b = []
    GREEN, GREEN_TINT = "#1baf7a", "rgba(27,175,122,0.08)"
    FOOT = {"native": "rgba(42,120,214,0.40)", "ad": "rgba(235,104,52,0.45)",
            "wide": "rgba(27,175,122,0.45)"}
    zones = (("Native", 40, C["link_tint"], "native", C["link"]),
             ("A&D", 284, C["accent_tint"], "ad", C["accent"]),
             ("Widened", 528, GREEN_TINT, "wide", GREEN))
    for name, x, tint, _, _ in zones:
        b.append(rect(x, 104, 232, 336, tint, rx=8))
        b.append(text(x + 12, 128, name, 14, weight=700))
    b += [path([(156, 88), (156, 140)]), path([(400, 88), (400, 140)]), path([(644, 88), (644, 140)]),
          path([(156, 420), (156, 456)]), path([(400, 420), (400, 456)]), path([(644, 420), (644, 456)])]
    b += [box(40, 40, 720, 48, "Injected 8 µm counts", kind="store", node="src"),
          box(40, 456, 720, 48, "Verdict on the same 8 µm bins", kind="store", node="verdict")]

    def fine(x0, k, key):
        foot = nearest(18, 8, 8, k)
        return grid(x0, 140, 18, 18, 8, lambda i, j: C["ink"] if (i, j) == (8, 8)
                    else (FOOT[key] if (i, j) in foot else C["paper"]), sw=0.5)
    b += fine(84, 48, "native")
    foot = nearest(9, 4, 4, 48)
    b += grid(328, 140, 9, 9, 16, lambda i, j: C["ink"] if (i, j) == (4, 4)
              else (FOOT["ad"] if (i, j) in foot else C["paper"]), C["soft"], 0.8)
    for k in range(9):
        b.append(f'<line x1="{328 + 16 * k + 8}" y1="140" x2="{328 + 16 * k + 8}" y2="284" '
                 f'stroke="rgba(11,11,11,0.12)" stroke-width="0.4"/>')
        b.append(f'<line x1="328" y1="{140 + 16 * k + 8}" x2="472" y2="{140 + 16 * k + 8}" '
                 f'stroke="rgba(11,11,11,0.12)" stroke-width="0.4"/>')
    b += fine(572, 192, "wide")
    for cx, l1 in ((156, "48 nearest 8 µm bins"), (400, "48 nearest 16 µm bins"),
                   (644, "192 nearest 8 µm bins")):
        b.append(text(cx, 308, l1, anchor="middle"))

    flag = lambda col, cells: (lambda i, j: col if (i, j) in cells else C["paper"])
    b += grid(132, 344, 4, 4, 12, flag(C["link"], {(1, 2)}), sw=0.8)
    b += grid(336, 344, 2, 2, 24, flag(C["accent"], {(0, 1)}), sw=0.8)
    b.append(mini_arrow(388, 410, 368, C["accent"], "arrow-accent"))
    b += grid(416, 344, 4, 4, 12, flag(C["accent"], {(0, 2), (0, 3), (1, 2), (1, 3)}), sw=0.8)
    b += grid(620, 344, 4, 4, 12, flag(GREEN, {(1, 2)}), sw=0.8)
    for cx, lab in ((156, "verdict per 8 µm bin"), (400, "projected to 2 × 2 bins"),
                    (644, "verdict per 8 µm bin")):
        b.append(text(cx, 412, lab, color=C["muted"], anchor="middle"))

    cap = ("<b>Three routes to the same 8 µm verdict, shown for SpotSweeper at the 16 µm "
           "step.</b> Each panel draws one bin (black) and the bins it is compared against, at "
           "the same physical scale. Natively, a bin is compared with its 48 nearest 8 µm "
           "bins. A&amp;D sums the counts of 2 × 2 bins into 16 µm bins, compares each with its "
           "48 nearest 16 µm bins, and projects the verdict back with "
           "<code>barcode_based_deconvolve</code>, so the four 8 µm bins of a coarse bin share "
           "one verdict. The widened arm stays on the 8 µm grid and compares each bin with its "
           "192 nearest bins, the same area as 48 bins of 16 µm. A&amp;D and the widened arm "
           "therefore look across the same area and differ only in whether counts are pooled "
           "before the test. The 32 and 64 µm steps scale the same way (768 and 3072 "
           "neighbours). The widened arms were run on a subset of sections, at the deepest and "
           "shallowest levels.")
    return page("methods-routes", "Materials and Methods · arms",
                "Aggregating counts versus widening the neighbourhood",
                "Three routes from the injected 8 µm counts to a verdict on the same 8 µm bins, "
                "drawn as the area each bin is compared against: 48 native bins, 48 summed "
                "16 µm bins, or 192 native bins covering the same area as A&D.",
                544, "\n".join(b), cap)

def lqr():
    b, d = [], []
    DY = 132
    for x, y, t in ((40, 56, "a"), (416, 56, "b"), (40, 332, "c"), (40, 404 + DY, "d"),
                    (416, 404 + DY, "e")):
        b.append(letter(x, y, t))
    gx, gy, s = 52, 80, 20
    cxb, cyb = gx + 4.5 * s, gy + 4.5 * s
    for i in range(9):
        for j in range(9):
            inside = (i - 4) ** 2 + (j - 4) ** 2 <= 4
            fill = DMG_TINT if inside else C["paper"]
            stroke = DMG if inside else C["rule_solid"]
            b.append(f'<rect x="{gx + j * s}" y="{gy + i * s}" width="{s}" height="{s}" '
                     f'fill="{fill}" stroke="{stroke}" stroke-width="{1 if inside else 0.8}"/>')
    b.append(f'<circle cx="{cxb}" cy="{cyb}" r="{2 * s}" fill="none" stroke="{C["ink"]}" '
             f'stroke-width="1.2" stroke-dasharray="4,3"/>')
    for i in range(9):
        for j in range(9):
            if (i - 4) ** 2 + (j - 4) ** 2 <= 4:
                b.append(f'<circle cx="{gx + j * s + s / 2}" cy="{gy + i * s + s / 2}" r="1.6" '
                         f'fill="{DMG}"/>')
    tx = 248
    b.append(f'<circle cx="256" cy="92" r="7" fill="none" stroke="{C["ink"]}" stroke-width="1.2" '
             f'stroke-dasharray="3,2"/>')
    b.append(text(272, 96, "circle, 32 µm across", color=C["muted"]))
    b.append(f'<circle cx="256" cy="116" r="2" fill="{DMG}"/>')
    b.append(text(272, 120, "bin centre", color=C["muted"]))
    b.append(rect(250, 134, 12, 12, DMG_TINT, DMG, 1))
    b.append(text(272, 144, "bin in the region", color=C["muted"]))
    for k, (dm, nb) in enumerate((("32 µm", "13 bins"), ("64 µm", "49 bins"),
                                  ("256 µm", "797 bins"))):
        b.append(text(tx, 196 + 18 * k, dm, weight=700))
        b.append(text(tx + 56, 196 + 18 * k, nb, color=C["muted"]))
    b.append(text(52, 290, "8 µm bin", color=C["muted"]))
    b.append(f'<line x1="112" y1="286" x2="132" y2="286" stroke="{C["rule_solid"]}" stroke-width="1"/>')
    b.append(f'<rect x="136" y="278" width="16" height="16" fill="{C["paper"]}" '
             f'stroke="{C["rule_solid"]}" stroke-width="0.8"/>')

    GREYDOT = "rgba(11,11,11,0.35)"
    def molecules(x0, y0, n, sp=8):
        out = [rect(x0, y0, 14 + 9 * sp, 14 + 9 * sp, C["paper"], C["rule_solid"], 0.8, rx=2)]
        for k in range(n):
            i, j = divmod(k, 10)
            out.append(f'<circle cx="{x0 + 7 + sp * j}" cy="{y0 + 7 + sp * i}" r="2" '
                       f'fill="{GREYDOT}"/>')
        return out
    b += molecules(424, 116, 100)
    b.append(mini_arrow(520, 556, 159))
    b.append(text(538, 151, "≈ 1%", color=C["muted"], anchor="middle"))
    b += molecules(566, 116, 1)
    b.append(text(467, 222, "before", color=C["muted"], anchor="middle"))
    b.append(text(609, 222, "after", color=C["muted"], anchor="middle"))

    OTHER, U, BH = "rgba(11,11,11,0.22)", 0.6, 10
    def bar(x, y, oth, mt, add=0):
        out = [rect(x, y, oth * U, BH, OTHER)]
        if mt:
            out.append(rect(x + oth * U, y, mt * U, BH, C["muted"]))
        if add:
            out.append(rect(x + (oth + mt) * U, y, add * U, BH, DMG))
        return out
    t0, y0 = 356, 388
    tgt = 10 / 9
    b.append(text(40, t0, "1  Baseline", weight=700))

    for k, mt in enumerate((2, 6, 10, 14, 18)):
        y = y0 + 16 * k
        b += bar(48, y, 100 - mt, mt)
        b.append(text(116, y + 9, f"{mt}%", color=C["muted"]))
    for k in range(3):
        b.append(f'<circle cx="{70 + 8 * k}" cy="{y0 + 86}" r="1.5" fill="{C["muted"]}"/>')
    bx, ya, yb = 148, y0, y0 + 74
    ym = (ya + yb) / 2
    b.append(f'<path d="M {bx} {ya} L {bx + 6} {ya} L {bx + 6} {yb} L {bx} {yb} '
             f'M {bx + 6} {ym} L {bx + 12} {ym}" fill="none" stroke="{C["muted"]}" '
             f'stroke-width="1"/>')
    b.append(text(168, ym - 6, "mean", color=C["muted"]))
    b.append(text(168, ym + 14, "10%", NAME, weight=700))
    b.append(mini_arrow(212, 244, ym))
    b.append(text(252, t0, "2  Target", weight=700))

    b += bar(252, y0, 90, 10)
    b.append(text(320, y0 + 9, "MT : other = 10 : 90", color=C["muted"]))
    xs = 252
    for k, (tok, lab) in enumerate((("10%", "to odds"), ("0.11", "× 10"), ("1.11", "to %"),
                                    ("53%", None))):
        big = lab is None
        w = 34 if big else 24
        b.append(text(xs, ym + (6 if big else 4), tok, NAME if big else SUB, weight=700 if big else 400))
        xs += w + 6
        if lab:
            b.append(mini_arrow(xs, xs + 26, ym))
            b.append(text(xs + 13, ym - 10, lab, color=C["muted"], anchor="middle"))
            xs += 26 + 10
    b += bar(252, y0 + 56, 90, 10, 90)
    b.append(text(374, y0 + 65, "MT : other = 100 : 90", color=C["muted"]))

    b.append(mini_arrow(xs + 2, 510, ym))
    b.append(text(520, t0, "3  Every bin in the region", weight=700))

    for k, mt in enumerate((2, 10, 18)):
        y = y0 + 16 * k
        b += bar(528, y, 100 - mt, mt, tgt * (100 - mt) - mt)
        b.append(text(660, y + 9, f"{mt}% → 53%", color=C["muted"]))
    b += bar(528, y0 + 58, 40, 60)
    b.append(text(660, y0 + 67, "60%, unchanged", color=C["muted"]))
    for x, col, lab in ((40, OTHER, "other genes"), (148, C["muted"], "MT before"),
                        (256, DMG, "MT added")):
        b.append(rect(x, y0 + 108, 10, 10, col))
        b.append(text(x + 14, y0 + 117, lab, color=C["muted"]))

    mx, my = 40, 428 + DY
    cols = (("32 µm", 2), ("64 µm", 4), ("256 µm", 16))
    b.append(text(mx + 152 + 100, my + 20, "region diameter", color=C["muted"], anchor="middle"))
    for k, (lab, _) in enumerate(cols):
        b.append(text(mx + 152 + 72 * k + 28, my + 40, lab, weight=700, anchor="middle"))
    for r, (lab, col) in enumerate((("Signal loss", DMG),
                                    ("Mito elevation", DMG))):
        yy = my + 56 + 44 * r
        b.append(text(mx, yy + 22, lab, weight=700))
        for k, (_, rad) in enumerate(cols):
            ccx = mx + 152 + 72 * k + 28
            b.append(f'<rect x="{ccx - 28}" y="{yy}" width="56" height="36" rx="4" '
                     f'fill="{C["paper"]}" stroke="{C["rule_solid"]}" stroke-width="0.8"/>')
            b.append(f'<circle cx="{ccx}" cy="{yy + 18}" r="{rad}" fill="{DMG_TINT}" '
                     f'stroke="{col}" stroke-width="1"/>')
    yy = my + 56 + 88
    b.append(text(mx, yy + 22, "No damage", weight=700))
    b.append(f'<rect x="{mx + 152}" y="{yy}" width="200" height="36" rx="4" fill="{C["paper"]}" '
             f'stroke="{C["rule_solid"]}" stroke-width="0.8"/>')
    b.append(text(mx + 252, yy + 22, "none", color=C["muted"], anchor="middle"))

    tx0, ty0, tw, th = 440, 440 + DY, 220, 152
    d.append(f'<clipPath id="lqr-tissue"><rect x="{tx0}" y="{ty0}" width="{tw}" height="{th}" '
             f'rx="28"/></clipPath>')
    b.append(f'<rect x="{tx0}" y="{ty0}" width="{tw}" height="{th}" rx="28" '
             f'fill="rgba(82,81,78,0.10)" stroke="{C["soft"]}" stroke-width="1"/>')
    b.append(f'<circle cx="{tx0 + 76}" cy="{ty0 + 80}" r="28" fill="{DMG_TINT}" '
             f'stroke="{DMG}" stroke-width="1.2"/>')
    ecx, ecy = tx0 + tw, ty0 + 76
    b.append(f'<circle cx="{ecx}" cy="{ecy}" r="28" fill="none" stroke="{DMG}" '
             f'stroke-width="1" stroke-dasharray="4,3"/>')
    b.append(f'<circle cx="{ecx}" cy="{ecy}" r="28" fill="{DMG_TINT}" '
             f'stroke="{DMG}" stroke-width="1.2" clip-path="url(#lqr-tissue)"/>')
    b.append(text(tx0 + 76, ty0 + 132, "interior", color=C["muted"], anchor="middle"))
    b.append(text(ecx + 36, ecy - 12, "edge", color=C["muted"]))
    b.append(text(tx0 + 16, ty0 + 24, "tissue", color=C["muted"]))
    cap = ("<b>Synthetic low-quality regions (LQRs).</b> (a) A region is the set of whole 8 µm "
           "bins whose centres fall inside a circle; a 32 µm circle centred on a bin takes 13 "
           "bins. (b) Signal loss keeps each count with probability 0.01. (c) Mitochondrial "
           "(MT) elevation, drawn on bins of 100 UMIs. (1) The section's baseline is the mean "
           "per-bin MT fraction over all its tissue bins at native depth, here 10%. Thinning "
           "does not change this mean, so all depths of a section share one target. (2) The "
           "target is set on the odds scale, the MT molecules counted against all other "
           "molecules. The baseline of 10% is 10 against 90, odds 0.11; ten times that is "
           "1.11, or 100 against 90, which converts back to a fraction of 53%. The target is "
           "never below 25%. (3) "
           "Every bin in the region receives the MT counts that bring it to the target, "
           "spread over MT- genes by their abundance in the section, so a bin that starts low "
           "receives more than one that starts high, and a bin already above the target "
           "receives none. The number added is drawn from a Poisson distribution around that "
           "need, so bins end near the target rather than on it; bars show the expected "
           "values. (d) Every run carries seven core conditions, each with 25 non-overlapping "
           "regions of one diameter. (e) Core regions are centred on tissue; a subset of runs "
           "adds edge regions centred within one radius of the tissue boundary, scored on "
           "their on-tissue part. Regions are added once, to the 8 µm counts, and the 16 to "
           "64 µm views are sums of those bins.")
    return page("methods-lqr", "Materials and Methods · synthetic damage",
                "Synthetic low-quality regions",
                "How synthetic damage is built: whole 8 µm bins inside a circle, signal loss, "
                "mitochondrial elevation set from each section's baseline, seven core "
                "conditions, and interior or edge placement.",
                676 + DY, "\n".join(b), cap, "\n".join(d))

def downstream():
    b = []
    base = lambda i, j: SERIES[0] if j < 2 else (SERIES[1] if i < 3 else SERIES[2])
    def cmap(x, y, diff=()):
        return grid(x, y, 5, 5, 10, lambda i, j: SERIES[2] if (i, j) in diff else base(i, j),
                    C["paper"], 0.6)

    b.append(letter(40, 30, "a"))
    for cx, lab in ((88, "one section"), (216, "thin"), (336, "same QC arm"),
                    (456, "pool to 32 µm"), (576, "Leiden clusters"), (716, "compare")):
        b.append(text(cx, 64, lab, weight=700, anchor="middle"))
    b += [path([(120, 184), (140, 184), (140, 136), (188, 136)]),
          path([(120, 208), (156, 208), (156, 256), (188, 256)])]
    for y in (136, 256):
        b += [path([(244, y), (308, y)]), path([(364, y), (428, y)]), path([(484, y), (548, y)])]
    b += [path([(602, 136), (632, 136), (632, 184), (672, 184)]),
          path([(602, 256), (648, 256), (648, 208), (672, 208)])]

    b.append(rect(56, 164, 64, 64, C["paper"], C["soft"], 1, rx=2))
    b.append(f'<path d="M 64 196 c 0 -18 22 -26 38 -18 c 12 6 12 26 -2 34 c -14 8 -36 4 -36 -16 z" '
             f'fill="{TISSUE}" stroke="{C["soft"]}" stroke-width="0.8"/>')
    b.append(text(88, 246, "undamaged", color=C["muted"], anchor="middle"))

    for y, h in ((112, 48), (232, 6)):
        b.append(rect(204, y, 24, 48, C["paper"], C["rule_solid"], 0.8))
        b.append(rect(204, y + 48 - h, 24, h, C["muted"]))
    b.append(text(216, 176, "100 UMIs/bin", color=C["muted"], anchor="middle"))
    b.append(text(216, 296, "12 UMIs/bin", color=C["muted"], anchor="middle"))

    for y, rem in ((112, {(1, 4), (4, 1)}), (232, {(0, 2), (1, 4), (2, 1), (4, 1), (4, 4)})):
        b += grid(312, y, 6, 6, 8)
        for i, j in rem:
            b.append(cross(312 + 8 * j + 4, y + 8 * i + 4, 2.6))
    b.append(text(336, 176, "× removed", color=C["muted"], anchor="middle"))

    for y in (112, 232):
        b += grid(432, y, 3, 3, 16, lambda i, j: "rgba(11,11,11,0.16)")
        for k in range(3):
            b.append(f'<line x1="{432 + 16 * k + 8}" y1="{y}" x2="{432 + 16 * k + 8}" y2="{y + 48}" '
                     f'stroke="rgba(11,11,11,0.18)" stroke-width="0.4"/>')
            b.append(f'<line x1="432" y1="{y + 16 * k + 8}" x2="480" y2="{y + 16 * k + 8}" '
                     f'stroke="rgba(11,11,11,0.18)" stroke-width="0.4"/>')
    b.append(text(456, 176, "survivors only", color=C["muted"], anchor="middle"))

    b += cmap(552, 111)
    b += cmap(552, 231, {(2, 2), (3, 1), (4, 3)})

    b.append(box(672, 168, 88, 56, "ARI", ["bin by bin"], "store", "ari"))

    b.append(letter(40, 340, "b"))
    b += cmap(56, 356)
    b.append(text(128, 385, "ARI", weight=700, anchor="middle"))
    b += cmap(150, 356, {(1, 2), (4, 0)})
    b.append(text(81, 422, "100", color=C["muted"], anchor="middle"))
    b.append(text(175, 422, "12", color=C["muted"], anchor="middle"))
    b.append(text(216, 389, "no filter", weight=700))
    b += cmap(432, 356)
    b.append(text(504, 385, "ARI", weight=700, anchor="middle"))
    b += cmap(526, 356, {(2, 1)})
    b.append(text(457, 422, "seed 0", color=C["muted"], anchor="middle"))
    b.append(text(551, 422, "seed 1", color=C["muted"], anchor="middle"))
    b.append(text(592, 389, "second seed", weight=700))

    cap = ("<b>Downstream consistency.</b> One section without injected damage is thinned to "
           "100 and 12 UMIs per bin, so tissue, cells and bins are identical and only the "
           "library differs. The same QC arm is applied at each depth, the surviving 8 µm bins "
           "are pooled into 32 µm bins and clustered (scanpy, 2,000 highly variable genes, 30 "
           "principal components, 15 neighbours, Leiden resolution 1.0), and the two "
           "clusterings are compared bin by bin with the adjusted Rand index (ARI) over the "
           "32 µm bins present at both depths. Barcodes are shared, so nothing has to be "
           "registered. The no-filter arm gives the agreement that depth alone allows, which a "
           "QC arm has to beat, and reclustering the same survivors with a new seed gives the "
           "floor below which a difference is the clustering algorithm's.")
    return page("methods-downstream", "Materials and Methods · downstream consistency",
                "Does the QC choice change what the clustering finds?",
                "One section is thinned to 100 and 12 UMIs per bin, filtered by the same QC arm, "
                "pooled to 32 µm and clustered, and the two clusterings are compared by adjusted "
                "Rand index against a no-filter reference and a reseeding floor.",
                464, "\n".join(b), cap)

def marker():
    b = []
    UP1, UP2 = 72, 144
    for x, y, t in ((40, 56, "a"), (416, 56, "b"), (40, 356 - UP1, "c"), (416, 356 - UP1, "d"),
                    (40, 628 - UP2, "e")):
        b.append(letter(x, y, t))

    for k in range(12):
        h = round(96 * 0.84 ** k)
        mt = k in (2, 6)
        b.append(rect(56 + 18 * k, 200 - h, 12, h, C["rule_solid"] if mt else C["muted"]))
        if mt:
            b.append(text(62 + 18 * k, 216, "MT", color=C["muted"], anchor="middle"))
    b.append(text(279, 196, "…", color=C["muted"], anchor="middle"))
    for m, h in enumerate((12, 10, 8)):
        b.append(rect(292 + 18 * m, 200 - h, 12, h, C["paper"], C["muted"], 0.8))
    b.append(text(316, 216, "other genes", color=C["muted"], anchor="middle"))
    b.append(f'<path d="M 56 222 L 56 228 L 266 228 L 266 222" fill="none" stroke="{C["ink"]}" '
             f'stroke-width="1"/>')
    b.append(text(161, 244, "marker panel", color=C["ink"], anchor="middle"))

    LIGHT = "rgba(11,11,11,0.08)"
    for cx, h, share, lab in ((476, 104, 0.6, "intact"), (576, 40, 0.6, "intact, shallow"),
                              (676, 40, 0.2, "damaged")):
        dark = round(h * share)
        b.append(rect(cx - 18, 200 - h, 36, h - dark, LIGHT, C["rule_solid"], 0.8))
        b.append(rect(cx - 18, 200 - dark, 36, dark, C["muted"]))
        b.append(text(cx, 200 - h - 6, f"{round(share * 100)}%", weight=700, anchor="middle"))
        b.append(text(cx, 216, lab, color=C["muted"], anchor="middle"))
    b.append(rect(432, 234, 12, 12, C["muted"]))
    b.append(text(450, 244, "panel genes", color=C["muted"]))
    b.append(rect(552, 234, 12, 12, LIGHT, C["rule_solid"], 0.8))
    b.append(text(570, 244, "other genes", color=C["muted"]))

    b.append(f'<g transform="translate(0,-{UP1})">')
    b.append(rect(56, 376, 120, 120, TISSUE, C["soft"], 0.8, rx=10))
    for k in (1, 2):
        b.append(f'<line x1="{56 + 40 * k}" y1="376" x2="{56 + 40 * k}" y2="496" '
                 f'stroke="{C["soft"]}" stroke-width="0.8"/>')
        b.append(f'<line x1="56" y1="{376 + 40 * k}" x2="176" y2="{376 + 40 * k}" '
                 f'stroke="{C["soft"]}" stroke-width="0.8"/>')
    b.append(rect(96, 416, 40, 40, "none", C["ink"], 1.6))
    b.append(mini_arrow(180, 204, 436))
    shade = ("rgba(11,11,11,0.05)", "rgba(11,11,11,0.16)", "rgba(11,11,11,0.30)")
    b += grid(208, 376, 8, 8, 15, lambda i, j: shade[(i * i + 3 * j + i * j) % 3], C["paper"], 0.8)
    b.append(rect(208, 376, 120, 120, "none", C["ink"], 1.6))
    for i, j in ((1, 2), (2, 6), (4, 1), (5, 5), (6, 3), (3, 4)):
        b.append(cross(208 + 15 * j + 7.5, 376 + 15 * i + 7.5, 4))
    b.append(text(116, 512, "section, 256 µm tiles", color=C["muted"], anchor="middle"))
    b.append(text(268, 512, "one tile", color=C["muted"], anchor="middle"))

    b.append(text(448, 382, "marker share", color=C["muted"]))
    b.append(f'<path d="M 448 392 L 448 480 L 640 480" fill="none" stroke="{C["rule_solid"]}" '
             f'stroke-width="1"/>')
    kept = (404, 412, 400, 420, 408, 416)
    rem = (436, 430, 444, 440, 420, 456)
    mk, mr = sum(kept) / 6, sum(rem) / 6
    b.append(f'<line x1="452" y1="{mk:.1f}" x2="640" y2="{mk:.1f}" stroke="{C["muted"]}" '
             f'stroke-width="0.8" stroke-dasharray="4,3"/>')
    b.append(f'<line x1="452" y1="{mr:.1f}" x2="640" y2="{mr:.1f}" stroke="{C["ink"]}" '
             f'stroke-width="0.8" stroke-dasharray="4,3"/>')
    for k, (yk, yr) in enumerate(zip(kept, rem)):
        cx = 468 + 32 * k
        b.append(f'<line x1="{cx}" y1="{yk}" x2="{cx}" y2="{yr}" stroke="{C["rule_solid"]}" '
                 f'stroke-width="1"/>')
        b.append(f'<circle cx="{cx}" cy="{yk}" r="4" fill="{C["muted"]}"/>')
        b.append(f'<circle cx="{cx}" cy="{yr}" r="4" fill="{C["paper"]}" stroke="{C["ink"]}" '
                 f'stroke-width="1.2"/>')
    b.append(text(544, 496, "units (tile × stratum)", color=C["muted"], anchor="middle"))
    b.append(f'<circle cx="662" cy="408" r="4" fill="{C["muted"]}"/>')
    b.append(text(672, 412, "kept", color=C["muted"]))
    b.append(f'<circle cx="662" cy="428" r="4" fill="{C["paper"]}" stroke="{C["ink"]}" '
             f'stroke-width="1.2"/>')
    b.append(text(672, 432, "removed", color=C["muted"]))
    b.append(text(672, 452, "- - mean", color=C["muted"]))

    b.append("</g>")
    b.append(f'<g transform="translate(0,-{UP2})">')
    for k, (x, w, lab) in enumerate(((40, 200, "256 µm tiles, 20 strata"),
                                     (288, 200, "512 µm tiles, 10 strata"),
                                     (536, 224, "1,024 µm tiles, 5 strata"))):
        b.append(rect(x, 644, w, 32, C["paper"], C["ink"] if k == 0 else C["rule_solid"], 1, rx=4))
        b.append(text(x + w / 2, 664, lab, weight=700 if k == 0 else 400, anchor="middle"))
    b.append(mini_arrow(244, 284, 660))
    b.append(mini_arrow(492, 532, 660))

    b.append("</g>")
    cap = ("<b>Marker check.</b> On the undamaged condition of every run, the bins each QC arm "
           "removes are compared with the bins it keeps, following the marker validation of "
           "SpatialArtifacts with two changes, a panel taken from the data and a comparison "
           "matched for location and library size. (a) The marker panel is the 20 most "
           "abundant non-mitochondrial genes of the section, taken from the data rather than "
           "chosen by hand. (b) The marker share of a bin is its panel UMIs divided by all its "
           "UMIs. Damaged tissue loses its own RNA and keeps ambient RNA, so its share falls, "
           "whereas an intact bin that is merely shallow keeps its share. (c) Removed and kept "
           "bins are compared only within one 256 µm tile and one of 20 library-size strata, "
           "so a filter that removes shallow bins, or one kind of tissue, is not credited with "
           "finding damage. (d) The median shares of the two groups are paired across units, "
           "and the retention ratio is the mean of removed medians over the mean of kept "
           "medians, tested by the Wilcoxon signed-rank test. A ratio below 1.00 means part of "
           "the nominal false-positive rate was real damage. (e) When an arm removes too few "
           "bins for five units to qualify, the match is coarsened along a ladder fixed in "
           "advance.")
    return page("methods-marker", "Supplementary S1 · marker check",
                "Were the removed bins really damaged?",
                "Marker check: a data-driven panel of the 20 most abundant non-mitochondrial "
                "genes, the panel share of each bin, removed and kept bins compared within one "
                "tile and library-size stratum, and the retention ratio of their medians.",
                740 - UP2 - 24, "\n".join(b), cap)

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for name, fn in (("methods-overview", overview), ("methods-routes", routes),
                     ("methods-lqr", lqr), ("methods-downstream", downstream),
                     ("methods-marker", marker)):
        p = os.path.join(OUT, name + ".html")
        html = fn()
        open(p, "w").write(html)
        i, j = html.index("<svg "), html.index("</svg>") + len("</svg>")
        open(os.path.join(OUT, name + ".svg"), "w").write(
            '<?xml version="1.0" encoding="UTF-8"?>\n' + html[i:j].replace(
                f'role="img"', f'width="{W}" role="img"', 1) + "\n")
        print("wrote", p, "and .svg")
