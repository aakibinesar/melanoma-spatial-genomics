# Melanoma Spatial Genomics

Computational analysis of tumour-immune architecture in melanoma using
spatial transcriptomics.

**Question:** how does the spatial organisation of transcriptional programmes
and tumour-immune interactions differ across tissue regions in nevus-associated
melanoma, and can computational analysis recover biologically meaningful
spatial patterns linked to progression from nevus to melanoma?

This is an exploratory computational analysis of public data. It makes no
clinical or diagnostic claims.

## Status

In progress - currently at the dataset-understanding stage. Nothing below the
"Analysis plan" heading has been run yet; results will be added as they exist.

## Data

Public 10x Genomics Visium (FFPE, CytAssist) data of nevus-associated melanoma
from GEO accession **GSE298774** (19 tissue sections on 11 slides). Data are
not redistributed here - see [`data/README.md`](data/README.md) for download
instructions and [`docs/dataset_description.md`](docs/dataset_description.md)
for what the deposit does and does not contain.

## Analysis plan

1. Data inspection and dictionary (what each file is, how libraries map to slides)
2. Quality control, documented per sample and per threshold
3. Normalisation, dimensionality reduction, clustering (clusters stay
   unlabelled until markers are examined)
4. Spatial visualisation of clusters, genes and gene-set scores
5. Differential expression and marker analysis
6. Biological annotation using published gene signatures
7. Pathway / gene-set enrichment
8. Spatial neighbourhood analysis (Squidpy): neighbourhood enrichment,
   co-occurrence, spatial autocorrelation
9. Tumour-immune architecture
10. Pathway-level interpretation of melanoma driver biology (expression-based
    only; no mutation calling)

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
├── results/       figures and tables
└── references/    reading list
```

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
