# Data

Nothing in `data/` is committed to git. The dataset is public and is
downloaded by following the steps below.

## Primary dataset: GEO GSE298774

"Spatial gene expression and microenvironmental changes in the transition of
nevus to melanoma" (Leiden University Medical Center).

- **GEO series:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298774
- **Assay:** 10x Genomics Visium (FFPE, CytAssist, probe-based), human skin
- **Design (per GEO):** slides containing both nevus and melanoma tissue
- **Samples:** 11 GEO samples (GSM9022963 - GSM9022973); each is one Visium
  slide carrying 1-3 patient tissue sections ("libraries"), 19 libraries in total
- **Processing (per GEO):** Space Ranger v2.1.1, probe set GRCh38-2020-A. GEO
  says the depositors then kept spots with > 500 genes and < 30 %
  mitochondrial counts in their own analysis; the deposited counts themselves
  are **not** filtered (see `docs/dataset_description.md`)
- **Publication:** Kreuger et al., *Cancer Research* 86(17):4400-4413 (2026),
  doi:10.1158/0008-5472.CAN-25-2934 (PubMed 42258703)

### Files

Deposited as a single archive of 150 files (about 1.4 GB):

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298774/suppl/GSE298774_RAW.tar
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298774/suppl/GSE298774_A2_features.tsv.gz
```

Download and unpack into `data/raw/`:

```bash
mkdir -p data/raw && cd data/raw
curl -L -C - -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298774/suppl/GSE298774_RAW.tar
curl -L -C - -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298774/suppl/GSE298774_A2_features.tsv.gz
tar -xf GSE298774_RAW.tar
```

### Published supplementary tables (patient mapping, signature gene lists)

Public, CC BY, about 0.8 MB - Figshare doi:10.1158/0008-5472.33414485:

```bash
mkdir -p data/raw/supplementary && cd data/raw/supplementary
curl -L -o can-25-2934_supplementary_tables_suppst1.xlsx https://ndownloader.figshare.com/files/68129117
# md5 should be 64c1dbadae6ed7509fdc2362bffe60c7
```

## Companion data (not used)

The same study deposited imaging mass cytometry data on Zenodo
(doi:10.5281/zenodo.15582615, about 26 GB, CC BY 4.0). It is not part of this
analysis.

## Citation

If you use the data, cite the study above and GEO accession GSE298774.
