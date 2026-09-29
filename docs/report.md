# Final report

## The question

Does spot-level gene expression across tissue sections reveal how a benign
mole (nevus) transitions into melanoma, and how the immune system interacts
with that process spatially? Dataset: **GSE298774** (Kreuger et al., *Cancer
Research*, 2026) - 10x Genomics Visium spatial transcriptomics, 19 tissue
sections from 18 patients, each containing both nevus and melanoma-in-situ
components.

This is an exploratory computational reanalysis of public data. It makes no
clinical or diagnostic claims, and no mutation calls are made anywhere - the
assay measures RNA, not DNA.

## Data and pipeline

The deposit had real structural problems that had to be resolved before any
analysis was trustworthy: two gene probe panels in circulation, several
slides sharing one matrix across sibling tissue sections (risking
double-counted spots), and one library (A2) missing its barcode file
entirely. A custom loader (`src/geo_loading.py`, unit-tested) resolved all
of this - including recovering A2's barcodes from a sibling section and
correctly zero-filling 42 spots with no detectable transcripts rather than
silently dropping them. All 19 libraries were validated end-to-end against
their own declared spot counts before any analysis began (`notebooks/01`).

No spot-level nevus/melanoma labels exist in this deposit. The published
gene signatures from the source paper's supplement were used only as an
external prior for later cross-checks, never to "discover" nevus-vs-melanoma
differences in the same data that produced those signatures - that would be
circular.

## What was found

**Clustering (26 groups -> 23 usable, `notebooks/02`).** Spots were grouped
into 26 clusters by expression similarity. Three were flagged and excluded
as batch/technical artefacts (one library alone accounted for over 80% of
each). The remaining 23 span most of the 18 patients - real, shared
structure, not noise.

**Cell-type character (`notebooks/03`).** Marker-gene scoring identified a
clear melanocytic cluster, a clear smooth-muscle cluster, several
keratinocyte-like clusters, and roughly ten clusters honestly left
"ambiguous" rather than forced into a label. As an independent check, the
genes with the strongest spatial patterning (found with no cluster label as
input) were dominated by exactly the marker genes expected - melanocytic
(DCT, PMEL, TYRP1, MLANA), keratinocyte (KRT1, KRTDAP), immune (IGKC,
S100A8) - good evidence the clustering reflects real biology, not noise
(`notebooks/04`).

**Tumour-immune architecture (`notebooks/05`).** The initial framing - "one
melanocytic cluster, one immune cluster" - was too simple. A cluster first
labelled "the immune cluster" turned out to be elevated on melanocytic
markers *and* all three immune panels: not a separate immune region, but
tumour tissue with immune cells in it. Rebuilt as four compartments - a cold
tumour core, an immune-infiltrated tumour margin, non-tumour immune tissue,
and everything else - the data show infiltrated tissue tends to sit within
a few cell-spots of the tumour core, a pattern checked per-section before
being trusted: an early pooled version of this analysis suggested a third,
more distant immune layer, which turned out to be one section's idiosyncrasy
(99.3% of that signal came from a single library) and was retracted from the
conclusions rather than kept as a headline finding.

**The epidermal/keratinocyte finding (`notebooks/06`, `notebooks/09`) - the
strongest single result in this project.** A direct comparison between the
cold tumour core and the infiltrated margin initially flagged 90% of all
genes as "different" - an implausible number that traced to a ~23x
sequencing-depth gap between the two regions, not biology, and was corrected
for before drawing any conclusion. Once corrected, the genes that actually
distinguished the cold tumour core were not melanocyte genes at all - they
were skin/epidermis differentiation genes (keratins, filaggrin,
desmosome/cornified-envelope machinery). This unplanned result was chased to
a specific, falsifiable hypothesis: the "cold tumour" compartment is not one
uniform population, but spans a real gradient from epidermis-proximal
(keratinocyte-admixed, closer to normal skin) to purely melanocytic (closer
to where immune infiltration begins).

Three independent lines of evidence supported this:
- **Heterogeneity:** melanocytic and cornified-envelope scores anti-correlate
  within the cold tumour core (r=-0.45); sub-clustering splits it into a
  clearly keratinocyte-dominant group (~25% of spots) and a melanocytic-
  dominant group with its own internal spread.
- **Distance:** the cornified-envelope score rises, and the melanocytic
  score falls, with graph-distance from the infiltrated margin (p<1e-44
  pooled), and in the predicted direction in 11 of 15 individually-checked
  libraries - not a single outlier section driving a pooled illusion.
- **Direct visual confirmation:** plotting these spots on the actual
  histology images shows them tracing a thin line along the tissue's outer
  edge - exactly where epidermis sits anatomically - with the
  melanocytic/keratinocyte gradient running visibly along that same edge in
  two independent sections.

Read together, this reframes the "cold tumour" cluster as an
**epidermis-proximal population spanning a junctional-to-invasive
gradient** - which maps directly onto the histological definition of the
nevus-to-melanoma transition this project is built around. This is
co-expression and spatial-distance evidence, not a pathologist's layer
annotation, and is reported at that resolution.

**Pathway signals, separated from the same depth confound.** After
correcting for the sequencing-depth gap: a melanoma-specific MAPK-pathway
signature is genuinely elevated in the cold tumour core (invisible in the
uncorrected comparison, masked by the depth gap); a generic RAS-pathway
signature is elevated in the infiltrated margin, though likely partly
reflecting immune-cell content rather than tumour-intrinsic biology; and an
apparent PI3K/AKT signal mostly evaporated once depth was accounted for - it
had been the confound, not a real difference.

**External validation (`notebooks/07`).** A second, independent dataset
(GSE300445: 4 melanoma sections, different patients, different study, no
nevus component) was checked for the one finding sturdy enough to travel:
the infiltrated-margin-near-core pattern. Only one of the four sections had
enough of the relevant tissue types to test it at all - and in that section,
the same pattern held, even more sharply. Reported as one corroborating data
point from an independent study, not a validated cohort effect; the other
three sections either lacked the infiltrated compartment entirely or lacked
both compartments. Writing the loader for this dataset's different file
layout also surfaced and fixed a real bug in the shared loading code (a
hard-coded assumption that pixel coordinates were always integers, which
broke on this dataset's genuine sub-pixel float values).

**A small ML check (`notebooks/08`).** A simple classifier (logistic
regression, cross-validated by patient to avoid leakage) separated cold
from infiltrated tumour spots with 0.989 AUC - which looked excellent until
a baseline using *only* sequencing depth scored 0.894 almost by itself,
confirming most of that performance was the same depth confound found in
notebook 06. The genuine signal, once depth was removed, was more modest
(0.665 AUC) but built on the same real genes notebook 06 had already found -
consistency across two independent analyses, not a new discovery.

## What this project is careful not to claim

No diagnostic or clinical conclusions. No mutation calls from expression
data. No claim that any spot "is" nevus or melanoma tissue - only that it's
enriched for expression associated with published nevus/melanoma
signatures. No pretending a single-section pattern is a cohort-wide
finding - every pooled result was checked against its individual sections
before being trusted.

## Methodological rigor - the self-audit trail

This project's most distinctive feature is arguably not any single
biological result but the record of catching and correcting real analytical
mistakes before they became false conclusions, documented in place rather
than quietly fixed:

- A loader bug that would have double-counted spots on shared-matrix slides
- Library A2's missing barcode file, recovered and cross-validated rather
  than assumed
- Two independent instances of a marker-score "elevated everywhere, not just
  here" artefact (notebook 03), both caught by inspecting the whole score
  table rather than each cluster's top value in isolation
- A spatial-statistics library (squidpy) whose output was silently a ratio
  rather than the probability its parameter name implied, caught by reading
  the library's source rather than trusting its docs (notebook 04)
- A Windows-specific multiprocessing crash in that same library, fixed by
  identifying the correct backend rather than working around the symptom
- A sequencing-depth confound that surfaced in two independent forms - a
  differential-expression comparison (notebook 06) and a predictive model
  (notebook 08) - both times quantified and corrected rather than reported
  raw
- An over-generalised spatial pattern retracted after a per-library
  breakdown showed it was driven by a single outlier section (notebook 05),
  and the same discipline applied prospectively in every subsequent
  distance-based analysis (notebooks 07, 09)
- A second loader bug (pixel-coordinate typing) surfaced only when the
  pipeline was run against a second, independent dataset (notebook 07)

## Reproducing this report

Every finding above is backed by an executed, version-controlled notebook
(`notebooks/01` through `notebooks/09`), unit-tested loading code
(`tests/test_geo_loading.py`), and the result tables/figures each notebook
produced (`results/`). See `README.md` for setup and `data/README.md` for
how to obtain the source data (not redistributed in this repository).
