"""Generate the package-tutorial badge (assets/logo.svg).

Circular badge (circle = learning material) with an inverted teal & cream
palette, used for Min-Yao Jhu's package-tutorial repositories (book notes
use the dark notebook badge instead). A package hexagon is assembled from
pieces, with the last piece arriving: learning a package step by step.
Change TITLE to reuse it for another package.

Usage: python3 assets/make_logo.py
"""
import math
from pathlib import Path

TITLE = "hdWGCNA"

CREAM, TEAL_DK, TEAL, TEAL_LT, CORAL, GOLD, LILAC = "#F3E9D2", "#0F3B3A", "#2E8B83", "#7ED0C4", "#E9967A", "#F2C14E", "#B9A3D9"
FONT = 'font-family="Helvetica, Arial, sans-serif"'


def logo(title=TITLE):
    b = [f'<circle cx="120" cy="120" r="114" fill="{CREAM}" stroke="{TEAL_DK}" stroke-width="7"/>']
    cx, cy, r = 116, 100, 50
    v = [(cx + r * math.cos(math.radians(60 * i - 90)), cy + r * math.sin(math.radians(60 * i - 90))) for i in range(6)]
    for i, c in enumerate([TEAL_DK, TEAL, TEAL_DK, TEAL, TEAL_DK]):          # five pieces in place
        a, z = v[i], v[(i + 1) % 6]
        b.append(f'<polygon points="{cx},{cy} {a[0]:.1f},{a[1]:.1f} {z[0]:.1f},{z[1]:.1f}" fill="{c}" '
                 f'stroke="{CREAM}" stroke-width="3" stroke-linejoin="round"/>')
    a, z = v[5], v[0]                                                       # last piece arriving
    dx, dy = 22, -10
    b.append(f'<polygon points="{cx + dx},{cy + dy} {a[0] + dx:.1f},{a[1] + dy:.1f} {z[0] + dx:.1f},{z[1] + dy:.1f}" '
             f'fill="{CORAL}" stroke="{CREAM}" stroke-width="3" stroke-linejoin="round"/>')
    b.append(f'<path d="M{cx + 30},{cy - 34} q10,6 4,16" fill="none" stroke="{CORAL}" stroke-width="2.4" stroke-linecap="round"/>')
    for (x, y, c) in [(cx - 18, cy + 8, TEAL_LT), (cx + 4, cy + 26, GOLD), (cx - 24, cy - 14, GOLD),
                      (cx + 10, cy + 6, TEAL_LT), (cx - 6, cy - 30, LILAC)]:   # cells
        b.append(f'<circle cx="{x}" cy="{y}" r="4.5" fill="{c}"/>')
    b.append(f'<text x="120" y="182" text-anchor="middle" {FONT} font-size="21" font-weight="bold" fill="{TEAL_DK}">{title}</text>')
    b.append(f'<text x="120" y="201" text-anchor="middle" {FONT} font-size="11" letter-spacing="2.2" fill="{TEAL}">TUTORIAL NOTES</text>')
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 240 240" width="240" height="240" role="img" '
            f'aria-label="{title} tutorial notes badge">\n<defs></defs>\n' + "\n".join(b) + "\n</svg>\n")


if __name__ == "__main__":
    out = Path(__file__).with_name("logo.svg")
    out.write_text(logo())
    print("wrote", out)
