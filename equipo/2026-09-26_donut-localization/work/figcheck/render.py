"""Render paper/figures/<name>.pdf to PNG next to this file (PyMuPDF)."""
import os, sys, fitz
here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, "..", "..", "..", ".."))
for name in sys.argv[1:]:
    doc = fitz.open(os.path.join(root, "paper", "figures", name + ".pdf"))
    pix = doc[0].get_pixmap(dpi=200)
    out = os.path.join(here, name + ".png"); pix.save(out)
    print(out, pix.width, pix.height, "page size in:", doc[0].rect.width / 72, doc[0].rect.height / 72)
