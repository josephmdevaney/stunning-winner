"""Runs-of-Homozygosity (ROH) tool: autozygosity, consanguinity and UPD screening from
real `bcftools roh` output, plus coverage-based sex inference.

Self-contained: nothing here imports from the RepeatScope package, so this folder can be
lifted into another application as-is. Entry points:

    analysis.py        parsing, region classification, F-value, AR candidates, UPD flags
    reference.py       GRCh38 geometry, curated gene model, thresholds
    sex_inference.py   chrX/chrY coverage-ratio sex call from a mosdepth summary
    router.py          FastAPI routes, mounted under /api/roh
"""
