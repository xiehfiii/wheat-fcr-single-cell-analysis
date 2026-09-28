"""Recover an English-only, page-indexed code archive from the author's PDF.

This preserves code-like text for audit. PDF wrapping, font overlays and image-only
pages can prevent exact reconstruction; the archive is not an executable pipeline.
Run curated scripts in the main analysis folders for reproducible analyses.
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

import pdfplumber


SECTIONS = [
    (1, 4, "00_reference_and_environment", "Reference preparation and environments"),
    (5, 19, "01_quality_control", "Ambient RNA, QC and doublet filtering"),
    (20, 27, "02_atlas_clustering", "Atlas integration, clustering and marker discovery"),
    (28, 52, "03_cell_annotation", "Cell-type annotation and UMAP"),
    (53, 79, "04_differential_expression", "Pseudobulk contrast and enrichment"),
    (80, 86, "05_rna_velocity", "RNA velocity"),
    (87, 90, "06_pseudotime", "Monocle2 pseudotime and BEAM"),
    (91, 121, "07_scenic", "SCENIC gene-regulatory network"),
    (122, 125, "08_virtual_knockout", "scTenifoldKnk virtual knockout"),
    (126, 141, "09_coexpression", "hdWGCNA coexpression"),
    (142, 192, "10_cell_communication", "Plant ligand-receptor and CellChat"),
]

NUMBERED = re.compile(r"^\s{0,20}(\d{1,4})\s+(.*)$")
CJK = re.compile(r"[\u3400-\u9fff\uf900-\ufaff]")


def page_code(page: pdfplumber.page.Page) -> list[str]:
    """Collect numbered code lines, joining visual line wraps when detectable."""
    lines = (page.extract_text(layout=True) or "").splitlines()
    recovered: list[str] = []
    active = False
    for raw in lines:
        numbered = NUMBERED.match(raw)
        if numbered:
            code = numbered.group(2).rstrip()
            if CJK.search(code):
                # Chinese prose in comments/plot messages is not a code token.
                # Keep the source line in the audit count but remove untranslated
                # wording; curated, executable code is maintained separately.
                code = CJK.sub("", code)
            recovered.append(code)
            active = True
            continue
        if not active or not raw.strip() or CJK.search(raw):
            continue
        if raw.strip() in {"R", "Bash", "Python", "1", "2"}:
            continue
        # Continuations are indented to the same code column but unnumbered.
        if len(raw) - len(raw.lstrip()) < 10 or not recovered:
            continue
        fragment = raw.strip()
        if fragment.startswith("#") and recovered[-1].strip().startswith("#"):
            recovered.append(fragment)
            continue
        separator = " " if recovered[-1].endswith(" ") else ""
        recovered[-1] = recovered[-1].rstrip() + separator + fragment
    return recovered


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdf", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    logging.getLogger("pdfminer").setLevel(logging.ERROR)
    with pdfplumber.open(args.pdf) as pdf:
        if len(pdf.pages) != 192:
            raise ValueError(f"Expected 192 pages, found {len(pdf.pages)}")
        for first, last, stem, title in SECTIONS:
            blocks = [
                f"# {title}",
                "",
                "PDF transcription only. The PDF may wrap code or omit glyphs;",
                "use the curated scripts for executable analysis.",
                "",
            ]
            for page_no in range(first, last + 1):
                code = page_code(pdf.pages[page_no - 1])
                if not code:
                    continue
                blocks.extend([f"## Source PDF page {page_no}", "", "```text", *code, "```", ""])
            out = args.output_dir / f"{stem}.md"
            out.write_text("\n".join(blocks), encoding="utf-8")
            print(f"{out.name}: {sum(len(x) for x in blocks)} characters")


if __name__ == "__main__":
    main()
