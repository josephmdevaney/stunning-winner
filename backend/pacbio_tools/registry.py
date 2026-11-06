"""The PacBio GitHub tools this suite reads output from, runs, or sits next to.

`used` tools appear in this suite's data flow (their version can be read from the
files the suite loads, or from the script it vendors). The rest are PacBio tools
a lab would likely run alongside them; they're listed so "is everything on the
newest release?" can be answered for the whole toolchain in one place, but the
suite doesn't read their output, so they have no observed version.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Tool:
    key: str
    name: str
    repo: str          # owner/name on GitHub
    role: str          # what it does in a long-read pipeline
    used: bool         # part of this suite's data flow
    suite_use: str     # how this suite relates to it ("" when not used)


OWNER = "PacificBiosciences"

TOOLS: tuple[Tool, ...] = (
    Tool("paraphase", "Paraphase", f"{OWNER}/paraphase",
         "Phases haplotypes and calls copy number for highly similar gene families (SMN1/2, GBA, CYP21A2/RCCX, HBA, ...).",
         True, "Its JSON and per-region VCFs feed the Paraphase suite and the Carrier panel; its annotation script is vendored and run."),
    Tool("trgt", "TRGT", f"{OWNER}/trgt",
         "Genotypes tandem repeats and reports motif counts, interruptions and methylation.",
         True, "Its VCFs are the input to RepeatScope, PureTarget and the repeat loci of Carrier."),
    Tool("sawfish", "sawfish", f"{OWNER}/sawfish",
         "Calls structural variants from HiFi alignments.",
         True, "Large HBA1/2 deletions in the Carrier panel."),
    Tool("mitorsaw", "MitoRSaw", f"{OWNER}/MitoRSaw",
         "Calls mitochondrial variants, including heteroplasmy, from HiFi reads.",
         True, "Its VCFs are the input to MitoLR."),
    Tool("ptcp", "PTCP (PureTarget Carrier Pipeline)", f"{OWNER}/ptcp",
         "Runs alignment, TRGT, Paraphase, sawfish and QC for the PureTarget carrier panel.",
         True, "Produces the files the Carrier panel reads."),
    Tool("pbmm2", "pbmm2", f"{OWNER}/pbmm2",
         "Aligns HiFi reads (a minimap2 wrapper); the alignment step ahead of every tool above.",
         True, "Upstream of all of it; its version lives in the BAM header, which this suite doesn't read."),
    Tool("hiphase", "HiPhase", f"{OWNER}/HiPhase",
         "Jointly phases small variants, structural variants and tandem repeats.", False, ""),
    Tool("pb-cpg-tools", "pb-CpG-tools", f"{OWNER}/pb-CpG-tools",
         "Summarizes CpG methylation from HiFi reads.", False, ""),
    Tool("pbsv", "pbsv", f"{OWNER}/pbsv",
         "Calls structural variants from HiFi alignments (the older SV caller).", False, ""),
    Tool("hificnv", "HiFiCNV", f"{OWNER}/hificnv",
         "Calls copy number variants from HiFi read depth.", False, ""),
    Tool("pb-starphase", "pb-StarPhase", f"{OWNER}/pb-StarPhase",
         "Calls pharmacogene star alleles and diplotypes.", False, ""),
    Tool("pbaa", "pbaa", f"{OWNER}/pbaa",
         "Clusters HiFi amplicon reads into consensus sequences.", False, ""),
    Tool("pbfusion", "pbfusion", f"{OWNER}/pbfusion",
         "Detects gene fusions from Iso-Seq reads.", False, ""),
    Tool("pbtk", "pbtk", f"{OWNER}/pbtk",
         "Utilities for PacBio BAM and index files (pbindex, bam2fasta, ...).", False, ""),
)

BY_KEY = {t.key: t for t in TOOLS}
