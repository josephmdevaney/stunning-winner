#!/usr/bin/env python3
"""Build epigenetics/data/imprinted_dmrs.json from the published table of human germline imprinted DMRs.

Source: "Multi-locus imprinting disturbance (MLID): interim joint statement for clinical and molecular
diagnosis", Clinical Epigenetics 2024 (PMC11295890), Table 1 -- GRCh38 coordinates, the number of CpGs,
a clinical-association category, and the imprint origin of each DMR. The article is open access (CC BY);
the table is fetched from Europe PMC's full-text API and re-parsed, so the panel can be rebuilt and
audited rather than hand-copied.

Two things are added from other published sources, and each entry says which:
  * `methylated_allele` (maternal / paternal): definitional for oocyte- and sperm-derived germline DMRs;
    for secondary DMRs it comes from the Ivashchenko et al. 2025 Genome Medicine supplementary list
    (SECONDARY_ALLELE below, checked by interval overlap), and is null where neither source says.
  * `core` windows for the clinical loci: CpG-island-level windows from the HG002 imprinting analysis
    this project was seeded from (CORE_WINDOWS below).

    python backend/scripts/reference_data/build_imprinted_dmrs.py
"""

from __future__ import annotations

import json
import re
import sys
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

PMCID = "PMC11295890"
URL = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/fullTextXML"
OUT = Path(__file__).resolve().parents[2] / "epigenetics" / "data" / "imprinted_dmrs.json"

# Disorder abbreviations used in the source table -> names.
DISORDERS = {
    "TNDM": "Transient neonatal diabetes mellitus type 1",
    "SRS": "Silver-Russell syndrome",
    "BWS": "Beckwith-Wiedemann syndrome",
    "TS14": "Temple syndrome",
    "KOS14": "Kagami-Ogata syndrome",
    "PWS": "Prader-Willi syndrome",
    "AS": "Angelman syndrome",
    "PHP": "Pseudohypoparathyroidism type 1B",
}

# Parent-of-origin of methylation for secondary DMRs (Ivashchenko 2025 supplementary Table S3, via Joshi 2016),
# keyed by a published DMR name in the Table 1 above. Everything else is derived from `imprint origin`.
SECONDARY_ALLELE = {
    "MEG3:TSS-DMR": "paternal",
    "MAGEL2:TSS-DMR": "maternal",
    "NDN:TSS-DMR": "maternal",
    "SNRPN:alt-TSS-DMR": "maternal",
    "ZNF597:TSS-DMR": "paternal",
    "GNAS-NESP:TSS-DMR": "paternal",
    "GNAS-A/B:TSS-DMR": "maternal",
}

# Tighter CpG-island windows (0-based half-open) for the clinical loci, from the HG002 imprinting analysis.
CORE_WINDOWS = {
    "SNURF:TSS-DMR": ("chr15", 24954888, 24955906),
    "H19/IGF2:IG-DMR & TSS-DMR": ("chr11", 1998200, 2003508),
    "KCNQ1OT1:TSS-DMR": ("chr11", 2699180, 2700857),
    "PLAGL1:alt-TSS-DMR": ("chr6", 144007779, 144008709),
    "MEG3:TSS-DMR": ("chr14", 100825706, 100826764),
    "MEST:alt-TSS-DMR": ("chr7", 130490898, 130493270),
    "GRB10:alt-TSS-DMR": ("chr7", 50782055, 50783174),
}


# Where two published sources give different coordinates for the same DMR, both are kept here.
COORDINATE_NOTES = {
    "MEG3/DLK1:IG-DMR": (
        "The 2022 multi-locus testing consensus (Clin Epigenet, PMC9641836, Table 2) gives chr14:100,811,001-100,811,037 "
        "(determined by a different study); this 2024 table gives chr14:100,809,089-100,811,721. The wider 2024 "
        "interval is used."
    ),
}


def slug(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")


def parse_coords(text: str) -> tuple[str, int, int]:
    m = re.match(r"(chr\w+):([\d,]+)[–-]([\d,]+)", text.replace("–", "–").strip())
    if not m:
        raise ValueError(f"unparseable coordinates: {text!r}")
    chrom, a, b = m.group(1), int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    return chrom, a - 1, b  # published 1-based inclusive -> 0-based half-open


def origin_allele(origin: str) -> str | None:
    if origin == "Oocyte_gDMR":
        return "maternal"
    if origin == "Sperm_gDMR":
        return "paternal"
    return None


def main() -> int:
    with urllib.request.urlopen(URL, timeout=60) as resp:
        root = ET.fromstring(resp.read())
    table = next(
        tw for tw in root.iter("table-wrap")
        if " ".join("".join(tw.find("label").itertext()).split()) == "Table 1"
    )
    rows = []
    for tr in table.iter("tr"):
        cells = [" ".join("".join(c.itertext()).split()) for c in tr if c.tag in ("td", "th")]
        if len(cells) == 7 and cells[3].startswith("chr"):
            rows.append(cells)
    if not (45 <= len(rows) <= 60):
        print(f"unexpected number of DMR rows: {len(rows)}", file=sys.stderr)
        return 1

    regions = []
    for name, aliases, disorders, loc, n_cpgs, typ, origin in rows:
        chrom, start, end = parse_coords(loc)
        category = typ.rstrip("*")
        allele = origin_allele(origin)
        allele_source = "imprint origin" if allele else None
        if allele is None and name in SECONDARY_ALLELE:
            allele, allele_source = SECONDARY_ALLELE[name], "Ivashchenko 2025 Table S3 (Joshi 2016)"
        codes = [d.strip() for d in re.split(r"[/,;]", disorders) if d.strip()]
        entry = {
            "id": slug(name),
            "name": name,
            "aliases": [a.strip() for a in aliases.split(",") if a.strip()],
            "disorders": [DISORDERS[c] for c in codes if c in DISORDERS],
            # molecular subtypes / extra associations the table names (e.g. "upd(7)mat", "MLID (ZFP57)")
            "disorder_notes": [c for c in codes if c not in DISORDERS],
            "chrom": chrom, "start": start, "end": end,
            "n_cpgs_published": int(n_cpgs),
            "category": category,              # CA clinically associated | NC non-clinical | X polymorphic/tissue-variable
            "primary_diagnostic": typ.endswith("*"),
            "imprint_origin": origin,
            "methylated_allele": allele,
            "methylated_allele_source": allele_source,
        }
        for key, (cc, cs, ce) in CORE_WINDOWS.items():
            if name == key:
                entry["core"] = {"chrom": cc, "start": cs, "end": ce}
        if name in COORDINATE_NOTES:
            entry["coordinate_notes"] = COORDINATE_NOTES[name]
        regions.append(entry)

    doc = {
        "source": {
            "title": "Multi-locus imprinting disturbance (MLID): interim joint statement for clinical and molecular diagnosis",
            "journal": "Clinical Epigenetics 2024", "pmcid": PMCID, "table": "Table 1", "genome_build": "GRCh38",
            "license": "CC BY 4.0 (open access)",
            "coordinates": "start is 0-based, end exclusive (bed convention); the paper prints 1-based inclusive",
        },
        "regions": regions,
    }
    OUT.write_text(json.dumps(doc, indent=2) + "\n")
    cats = {c: sum(1 for r in regions if r["category"] == c) for c in ("CA", "NC", "X")}
    print(f"wrote {len(regions)} DMRs to {OUT} | categories {cats} | "
          f"allele known for {sum(1 for r in regions if r['methylated_allele'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
