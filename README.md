# Melanoma Spatial Genomics

[![tests](https://github.com/aakibinesar/melanoma-spatial-genomics/actions/workflows/tests.yml/badge.svg)](https://github.com/aakibinesar/melanoma-spatial-genomics/actions/workflows/tests.yml)

Computational analysis of tumour-immune architecture in melanoma using
spatial transcriptomics.

**Question:** how does the spatial organisation of transcriptional programmes
and tumour-immune interactions differ across tissue regions in nevus-associated
melanoma, and can computational analysis recover biologically meaningful
spatial patterns linked to progression from nevus to melanoma?

This is an exploratory computational analysis of public data. It makes no
clinical or diagnostic claims.

## Status

Complete. All stages below have been run; see [`docs/report.md`](docs/report.md)
for the full write-up of methods and findings.

## Data

Public 10x Genomics Visium (FFPE, CytAssist) data of nevus-associated melanoma
from GEO accession **GSE298774** (19 tissue sections on 11 slides). Data are
not redistributed here - see [`data/README.md`](data/README.md) for download
instructions and [`docs/dataset_description.md`](docs/dataset_description.md)
for what the deposit does and does not contain.

## Notebooks

1. `01_dataset_exploration` - data inspection and dictionary
2. `02_qc_clustering` - quality control, normalisation, clustering
3. `03_de_annotation_enrichment` - differential expression, marker
   annotation, pathway/gene-set enrichment
4. `04_spatial_statistics` - Squidpy neighbourhood enrichment, co-occurrence,
   spatial autocorrelation
5. `05_tumour_immune_architecture` - compartment definitions, tumour-immune
   spatial architecture
6. `06_pathway_interpretation` - melanoma driver-pathway activity
   (expression-based only; no mutation calling)
7. `07_external_validation` - checking findings against a second, independent
   dataset (GSE300445)
8. `08_ml_addon` - a small supervised-learning check, cross-validated by
   patient
9. `09_epidermal_signal` - following up an unplanned finding from notebook 6
   to a specific, testable spatial hypothesis

Throughout: observed, inferred and hypothesised statements are kept separate,
negative results are reported, and gene-signature scores are described as
"regions enriched for X-associated expression", not as verified cell types.

## Repository structure

```text
melanoma-spatial-genomics/
├── README.md
├── LICENSE
├── docs/          dataset description, background, limitations, report
├── data/          download instructions (data itself is not committed)
├── notebooks/     numbered analysis notebooks
├── src/           reusable loading / QC / spatial / plotting code
├── tests/         unit tests for src/ (run in CI - see .github/workflows/)
├── results/       figures and tables
└── references/    reading list
```

`requirements-test.txt` is a small subset of `requirements.txt` used by CI to
run the unit tests quickly, without installing the full scanpy/squidpy stack.

## Reproducing

```bash
python -m venv .venv
.venv/Scripts/activate          # Windows; use source .venv/bin/activate on Linux/macOS
pip install -r requirements.txt
```

Then download the data as described in `data/README.md` and run the notebooks
in order.

## License

Code: MIT (see `LICENSE`). Data: see the original GEO deposit and publication.
