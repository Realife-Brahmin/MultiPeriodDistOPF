"""Compose 5 T=24 rho-sweep PNGs into a 2x3 grid with per-panel titles."""
from PIL import Image, ImageDraw, ImageFont
from pathlib import Path

REVIEW = Path(__file__).parent / "pngs_review" / "ieee2552C_1ph"
RHOS = [4000, 8000, 16000, 24000, 32000]
files = [REVIEW / f"T024_rho{r:06d}_convergence.png" for r in RHOS]

PANEL_W = 720  # target width per panel; height scales proportionally
imgs = []
for f in files:
    im = Image.open(f)
    scale = PANEL_W / im.size[0]
    new_h = int(round(im.size[1] * scale))
    imgs.append(im.resize((PANEL_W, new_h), Image.LANCZOS))
W, H = imgs[0].size

cols, rows = 3, 2
title_h = 40
header_h = 60
pad = 10

cell_w = W
cell_h = H + title_h

canvas_w = cols * cell_w + (cols + 1) * pad
canvas_h = header_h + rows * cell_h + (rows + 1) * pad

canvas = Image.new("RGB", (canvas_w, canvas_h), "white")
draw = ImageDraw.Draw(canvas)

# Pick a font (fall back gracefully).
def get_font(size):
    for name in ("C:/Windows/Fonts/cmunrm.ttf",
                 "C:/Windows/Fonts/times.ttf",
                 "C:/Windows/Fonts/arial.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()

font_header = get_font(28)
font_panel = get_font(22)

header = "ieee2552C  T=24  rho sweep  (eps_pri=1e-4, eps_dual=5e-4)"
hw, hh = draw.textbbox((0, 0), header, font=font_header)[2:]
draw.text(((canvas_w - hw) // 2, (header_h - hh) // 2), header,
          fill="black", font=font_header)

positions = [(r, c) for r in range(rows) for c in range(cols)]
for (r, c), img, rho in zip(positions, imgs, RHOS):
    x = pad + c * (cell_w + pad)
    y = header_h + pad + r * (cell_h + pad)
    title = f"rho_init = {rho}"
    tw, th = draw.textbbox((0, 0), title, font=font_panel)[2:]
    draw.text((x + (cell_w - tw) // 2, y + (title_h - th) // 2),
              title, fill="black", font=font_panel)
    canvas.paste(img, (x, y + title_h))

out = REVIEW / "T024_rho_sweep_composite.png"
canvas.save(out, "PNG", optimize=True)
print(f"Saved: {out}  ({canvas_w}x{canvas_h})")
