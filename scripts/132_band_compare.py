from PIL import Image
import subprocess, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
D = "figures/Fig8_reproducibility.png"
R = "out/_layout_repro/Fig8_reproducibility.png"
for tag, p in (("delivered", D), ("reproduced", R)):
    im = Image.open(p)
    band = im.crop((0, 150, im.width, 380))
    out = "out/_band_%s.png" % tag
    band.save(out)
    print("=" * 60)
    print(tag, im.size, "-> band", band.size)
    r = subprocess.run(["tesseract", out, "stdout", "-l", "eng", "--psm", "6"],
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace")
    print(r.stdout.strip()[:600])
