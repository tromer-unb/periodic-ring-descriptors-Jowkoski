from pathlib import Path
import itertools
import subprocess
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"

def overlaps(pdf):
    with tempfile.NamedTemporaryFile(suffix=".html") as tmp:
        subprocess.run(
            ["pdftotext", "-bbox", str(pdf), tmp.name],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        root = ET.parse(tmp.name).getroot()
    words = []
    for node in root.iter():
        if not node.tag.endswith("word"):
            continue
        a = node.attrib
        words.append((
            float(a["xMin"]), float(a["yMin"]),
            float(a["xMax"]), float(a["yMax"]),
            "".join(node.itertext()),
        ))
    hits = []
    for left, right in itertools.combinations(words, 2):
        # pdftotext reports intentional superscript pieces of scientific-notation
        # tick labels (for example 10^-6) as overlapping word boxes. These are
        # typographic composition, not collisions between independent labels.
        math_tokens = {"10", "×", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"}
        if left[4] in math_tokens and right[4] in math_tokens:
            continue
        dx = min(left[2], right[2]) - max(left[0], right[0])
        dy = min(left[3], right[3]) - max(left[1], right[1])
        if dx > 0.5 and dy > 0.5 and dx * dy > 1.0:
            hits.append((dx * dy, left[4], right[4]))
    return sorted(hits, reverse=True)

failed = False
for pdf in sorted(FIG.glob("fig*.pdf")):
    hits = overlaps(pdf)
    print(f"{pdf.name}: {len(hits)} text-overlap pairs")
    if hits:
        failed = True
        for hit in hits[:10]:
            print("   ", hit)

if failed:
    raise SystemExit("Figure text-overlap audit failed.")
print("Figure text-overlap audit passed.")