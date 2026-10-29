PRINT_CM = 16.0
LETTER_PT = 9.0
FONT = "DejaVu Sans"

def letter_size(width_in, print_cm=PRINT_CM):
    return LETTER_PT * width_in * 2.54 / print_cm

def panel_letter(target, x, y, letter, width_in, print_cm=PRINT_CM):
    target.text(x, y, letter, fontsize=letter_size(width_in, print_cm), fontweight="bold",
                family=FONT, ha="left", va="top", color="#0b0b0b")

def letter_px(width_px, print_cm=PRINT_CM):
    return LETTER_PT / 72 * width_px / (print_cm / 2.54)

def letter_axes(fig, ax, letter, width_in=None):
    width_in = width_in or fig.get_figwidth()
    fig.canvas.draw()
    bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    h = letter_size(width_in) / 72 / fig.get_figheight()
    panel_letter(fig, bb.x0, bb.y1 + 1.2 * h, letter, width_in)
