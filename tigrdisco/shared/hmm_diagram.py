"""Domain diagram for the Excel reports."""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

_CANDIDATES = {
    "regular": ["C:/Windows/Fonts/arial.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"],
    "bold": [
        "C:/Windows/Fonts/arialbd.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    ],
}
_cache = {}


def _font(kind, size):
    key = (kind, size)
    if key not in _cache:
        for path in _CANDIDATES[kind]:
            try:
                _cache[key] = ImageFont.truetype(path, size)
                break
            except OSError:
                continue
        else:
            _cache[key] = ImageFont.load_default()
    return _cache[key]


def render(protein_length, intervals, protein_coverage, label, output):
    if not intervals:
        return False
    width, height = 1500, 285
    left, right, y = 30, 1410, 130
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw.text(
        (left, 18),
        f"Protein length: {protein_length} aa",
        fill="#333333",
        font=_font("regular", 22),
    )
    draw.line((left, y, right, y), fill="#222222", width=5)

    def x(pos):
        return left + (pos - 1) / max(1, protein_length - 1) * (right - left)

    for start, end in intervals:
        x0, x1 = x(start), x(end)
        draw.rounded_rectangle(
            (x0, y - 31, x1, y + 31), radius=8, fill="#3478C0", outline="#174A7E", width=2
        )
        if x1 - x0 >= 150:
            draw.text(((x0 + x1) / 2, y), label, fill="white", font=_font("bold", 22), anchor="mm")
    for pos in sorted(set((1, *range(50, protein_length, 50), protein_length))):
        tx = x(pos)
        draw.line((tx, y + 35, tx, y + 47), fill="black", width=2)
        draw.text((tx, y + 56), str(pos), fill="black", font=_font("regular", 17), anchor="ma")
    text = ", ".join(f"aa {a}\u2013{b}" for a, b in intervals)
    draw.text(
        (left, 238),
        f"Protein match: {text}  |  protein coverage: {protein_coverage:.1%}",
        fill="black",
        font=_font("regular", 21),
    )
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    image.save(output)
    return True
