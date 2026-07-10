"""MitoLR: mtDNA variant review for PacBio HiFi long-read sequencing.

Sibling product to RepeatScope. Shares the same FastAPI app and visual system,
but presents a genome-wide mitochondrial variant call table (vs. the rCRS
reference), per-variant heteroplasmy/read-level detail, and a sequencing QC
dashboard rather than nuclear tandem-repeat genotypes.
"""

__version__ = "0.1.0"
