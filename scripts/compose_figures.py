#!/usr/bin/env python3
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager

import figstyle as FS

FIG = "/work/benchmark_grid/figures"
FONT = font_manager.findfont(font_manager.FontProperties(family=FS.FONT, weight="bold"))
SETS = {
    "figS6": [f"recall_fpr_{m}.png" for m in ("mad", "miqc", "spotsweeper", "spatialartifacts")],
    "figS1_recall_mito": [f"recall_mito_{c}.png" for c in ("probe", "3prime")],
    "figS2_recall_signal_loss": [f"recall_signal_loss_{c}.png" for c in ("probe", "3prime")],
    "figS4_false_positive_rate": [f"false_positive_rate_{c}.png" for c in ("probe", "3prime")],
    "figS3_marker_retention": [f"marker_retention_{c}.png" for c in ("probe", "3prime")],
}
PAGES = [f"demo_downstream/Visium_HD_{s}.png" for s in
                         ("Human_Colon_Cancer", "Human_Kidney_FFPE", "Human_Lung_Cancer",
                          "3prime_Human_Ovarian_Cancer")]

def lettered(im, letter, width):
    px = FS.letter_px(width)
    font = ImageFont.truetype(FONT, int(round(px)))
    box = ImageDraw.Draw(Image.new("L", (1, 1))).textbbox((0, 0), letter, font=font)
    lw, lh = box[2] - box[0], box[3] - box[1]
    pad = int(round(0.6 * px))
    x0 = (width - im.width) // 2
    corner = np.asarray(im.convert("L"))[:lh + 2 * pad, :lw + 2 * pad]
    strip = 0 if corner.min() > 245 else lh + 2 * pad
    out = Image.new("RGB", (width, im.height + strip), "white")
    out.paste(im.convert("RGB"), (x0, strip))
    ImageDraw.Draw(out).text((pad - box[0], pad - box[1]), letter, font=font, fill="#0b0b0b")
    return out

def main():
    for name, files in SETS.items():
        ims = [Image.open(os.path.join(FIG, f)) for f in files]
        width = max(im.width for im in ims)
        parts = [lettered(im, "abcdefgh"[k], width) for k, im in enumerate(ims)]
        if name == "figS6":
            for letter, source, part in zip("abcd", files, parts):
                part.save(os.path.join(FIG, f"figS6{letter}_{source.removeprefix('recall_fpr_')}"),
                          dpi=(600, 600))
            continue
        out = Image.new("RGB", (width, sum(p.height for p in parts)), "white")
        y = 0
        for p in parts:
            out.paste(p, (0, y)); y += p.height
        out.save(os.path.join(FIG, name + ".png"), dpi=(600, 600))
        print("wrote", name, out.size, f"{16 * out.height / out.width:.1f} cm high at 16 cm")
    for letter, f in zip("abcd", PAGES):
        source = os.path.join(FIG, f)
        if not os.path.exists(source):
            print("not drawn, skipped:", f)
            continue
        with Image.open(source) as im:
            out = lettered(im, letter, im.width)
        name = f"figS5{letter}_clustering_maps"
        out.save(os.path.join(FIG, name + ".png"), dpi=(600, 600))
        print("wrote", name, out.size)


if __name__ == "__main__":
    main()
