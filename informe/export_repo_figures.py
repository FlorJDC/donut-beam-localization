# -*- coding: utf-8 -*-
"""Rasterize the manuscript figures reused in the short report (needs PyMuPDF: pip install pymupdf)."""
import os
import pymupdf

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for name in ("fig3_crb_maps", "fig7_zero_depth"):
    doc = pymupdf.open(os.path.join(ROOT, "paper", "figures", name + ".pdf"))
    out = os.path.join(ROOT, "informe", "figures", "repo_" + name + ".png")
    doc[0].get_pixmap(dpi=220).save(out)
    print("wrote", out)
