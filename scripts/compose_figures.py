import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager

import figstyle as FS

FIG = "/work/benchmark_grid/figures"
FONT = font_manager.findfont(font_manager.FontProperties(family=FS.FONT, weight="bold"))
SETS = {
    "figS1_cost_recall": [f"fig1_{m}.png" for m in ("mad", "miqc", "spotsweeper")],
    "figS2_recall_mito": [f"fig1_1_recall_mito_{c}.png" for c in ("probe", "3prime")],
    "figS3_recall_signal_loss": [f"fig1_2_recall_signal_loss_{c}.png" for c in ("probe", "3prime")],
    "figS8_cost_grid": [f"fig1_3_false_positive_{c}.png" for c in ("probe", "3prime")],
    "figS9_cost_marker": [f"fig2_discarded_tissue_{c}.png" for c in ("probe", "3prime")],
}
PAGES = {"figS10_maps": [f"demo_downstream/Visium_HD_{s}.png" for s in
                         ("Human_Colon_Cancer", "Human_Kidney_FFPE", "Human_Lung_Cancer",
                          "3prime_Human_Ovarian_Cancer")]}

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
        out = Image.new("RGB", (width, sum(p.height for p in parts)), "white")
        y = 0
        for p in parts:
            out.paste(p, (0, y)); y += p.height
        out.save(os.path.join(FIG, name + ".png"), dpi=(600, 600))
        print("wrote", name, out.size, f"{16 * out.height / out.width:.1f} cm high at 16 cm")
    for name, files in PAGES.items():
        for k, f in enumerate(files):
            if not os.path.exists(os.path.join(FIG, f)):
                print("not drawn, skipped:", f); continue
            im = Image.open(os.path.join(FIG, f))
            out = lettered(im, "abcdefgh"[k], im.width)
            out.save(os.path.join(FIG, f"{name}_{'abcd'[k]}.png"), dpi=(600, 600))
            print("wrote", f"{name}_{'abcd'[k]}", out.size)

if __name__ == "__main__":
    main()
