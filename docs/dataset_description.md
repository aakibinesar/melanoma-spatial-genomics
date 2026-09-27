# Dataset description - GEO GSE298774

Status: metadata, file listing and the published supplementary tables were read
first; the numbers in "What the files contain" were then measured on the
downloaded files during exploration. `notebooks/01_dataset_exploration.ipynb`
regenerates them from scratch and is the source of truth if the two differ.

## Provenance

| Field | Value |
|---|---|
| Title | Spatial gene expression and microenvironmental changes in the transition of nevus to melanoma |
| Accession | GSE298774 (BioProject PRJNA1271082) |
| Organism / tissue | Homo sapiens / skin |
| Assay | 10x Genomics Visium spatial gene expression, FFPE, CytAssist, probe-based |
| Platform | GPL24676 (Illumina NovaSeq 6000) |
| Depositing group | Leiden University Medical Center, Medical Oncology |
| Stated design | "Slides with both nevus and melanoma were used for 10X Visium spatial transcriptomics" |
| Publication | Kreuger et al., Cancer Research 86(17):4400-4413 (2026), doi:10.1158/0008-5472.CAN-25-2934 (closed access; abstract and supplementary tables read) |
| Supplementary tables | Public, CC BY: Figshare doi:10.1158/0008-5472.33414485 (Excel, Tables S1-S10) |
| Public since | 2 Jun 2026 (GEO record updated 1 Sep 2026) |

GEO states the depositors processed data with Space Ranger v2.1.1 (probe set
GRCh38-2020-A), then with the R package *semla*, "retaining spots with more than
500 genes and less than 30 % mitochondrial counts". **That describes their
downstream analysis, not the deposited counts** - see "Not filtered" below.

## Patients and lesions (from published Supplementary Table S1)

The 19 libraries come from **18 patients**. Library id = patient number
(A2 = patient 2, B13 = patient 13, ...). Patient 38 has two libraries
(`B38ex` and `B38_2` in the file names; S1 lists only "B38", so that both are
patient 38 is an inference). All 18 lesions are recorded as Low-CSD melanoma
with **both a melanoma in situ component and a dermal nevus component**, so
every lesion is expected to contain both tissue types. PRAME IHC is
negative in nevus and weak/positive in melanoma, and p16 shows a checkerboard
pattern in nevus versus aberrant in melanoma - independent, non-transcriptomic
context for interpreting results. Because samples are patients, statistical
comparisons must treat patient (not spot, not library) as the unit of
replication.

| GEO sample | Libraries (patient) |
|---|---|
| GSM9022963 | A2 (2), A17 (17) |
| GSM9022964 | A7 (7) |
| GSM9022965 | A21 (21) |
| GSM9022966 | A16 (16), A23 (23) |
| GSM9022967 | A22 (22), A24 (24) |
| GSM9022968 | A34 (34) |
| GSM9022969 | B1 (1), B38ex (38) |
| GSM9022970 | B4 (4) |
| GSM9022971 | B8 (8), B13 (13) |
| GSM9022972 | B14 (14), B15 (15), B38_2 (38) |
| GSM9022973 | D3 (3), D20 (20) |

## What the files contain (measured)

- 150 files, 1.4 GB, no corrupt archives. Flat `<GSM>_<library>_<file>.gz`
  naming, not a Space Ranger directory layout.
- Matrices are genes x spots Matrix Market. `tissue_positions.csv` has a
  header row and lists **every spot of the slide grid** (14,336 spots on the
  11 mm CytAssist area, 4,992 on the 6.5 mm area) with an `in_tissue` flag.
- **A library is defined by its `in_tissue` spots, not by its matrix.**
  In-tissue barcode sets never overlap between libraries on a slide, but:
  - For B1/B38ex, B8/B13, B14/B15/B38_2 and D3/D20 the matrix file is
    *byte-identical* across the slide's libraries (one slide-level matrix,
    duplicated).
  - For A2/A17 and A22/A24 the matrices carry identical data apart from
    all-zero columns (checked entry by entry after aligning barcodes).
  - Only A16/A23 (and single-library slides) ship a matrix that already
    contains just that library's spots.
  Loading each library's matrix whole would double-count spots.
- **Library A2 has no `barcodes.tsv.gz`/`features.tsv.gz`.** Its matrix has
  13,916 columns (A17's has 13,887); of those, exactly 13,750 are non-zero in
  both, and are entry-by-entry identical between the two matrices, so their
  barcode identities are recovered from A17. An all-zero column carries no
  fingerprint to match, so this only recovers 3,325 of A2's 3,367 declared
  in-tissue barcodes directly. The remaining 42 are barcodes A2's own
  `tissue_positions.csv` declares in-tissue but which never appear among the
  13,750 recovered identities - by elimination over the shared barcode
  universe, they have no detected transcript in A2, so they are loaded as
  explicit all-zero rows (flagged `expression_zero_filled` in `obs`) instead
  of being silently dropped. All 3,367 declared in-tissue spots are loaded.
  The series-level `GSE298774_A2_features.tsv.gz` is identical to the
  37,082-feature panel. Implemented and tested in `src/geo_loading.py`.
- **Two feature panels:** 37,082 features in 15 libraries and 18,085 in three
  (A21, A16, A23); the smaller is a strict subset of the larger. All are
  "Gene Expression". 10 (large panel) and 3 (small panel) gene *symbols* are
  duplicated, so gene ids, not symbols, are the key.

### Not filtered

The deposited counts are unfiltered. Many in-tissue spots fall below the
500-gene threshold GEO mentions (in-tissue spots with < 500 genes: B8 4,902 of
5,780; A34 2,520 of 4,151; A16 2,424 of 5,195; B13 2,073 of 3,614; A17 1,410 of
2,274). Quality control is therefore meaningful, and the filter has to be
applied (and justified) here.

### Depth varies about 70-fold between libraries

Median genes detected per in-tissue spot: B8 70, B14 147, B15 151, A17 244,
A34 300, B13 402, A16 554, B38_2 533, D3 1,622, A22 1,546, A23 1,099, B4 1,199,
A21 2,129, B38ex 3,469, A24 3,400, B1 5,075, D20 5,637. Median mitochondrial
fraction is low everywhere (2.4-5.9 %) and almost no in-tissue spot exceeds
30 %. Between-library depth differences are a batch/quality confound that any
cross-library comparison must handle.

### Library A2 looks like background

A2's in-tissue spots (all 3,367, including the 42 zero-filled ones above)
have a median of 35 genes, lower than the 68 of spots outside all tissue on
the same slide, and 42 of them have literally zero detected genes. Either A2
failed, or its coordinates do not match its matrix. This is checked with a
spatial overlay in notebook 01 before A2 is used for anything.

## Labels: what is and is not available

- **No spot- or region-level nevus/melanoma label is deposited**, in GEO, in
  the companion Zenodo record (imaging mass cytometry only, 26 GB), or in any
  of the ten supplementary tables (read in full). The paper's own annotations
  are in its figures and methods, which are behind a paywall. No preprint or
  code repository was found.
- What *is* published, and usable: (a) the patient/library mapping above;
  (b) Table S4, about 1,975 genes differing between nevus and "melanoma 1"
  spots, with log2 fold change and fraction of spots expressing; (c) Table S5,
  genes with on/off differences; (d) S8, transcription factors; (e) S6/S7,
  enriched Reactome pathways; (f) S9/S10, spatially co-expressed
  ligand-receptor pairs. These gene lists were derived by the authors from
  these same data.
- **Consequence for the analysis:** the published signatures can be used as an
  *external prior* to score spots, but they cannot be used to "discover" a
  nevus-versus-melanoma difference in the same data (circular). Every result
  that depends on them is labelled as such. The authors could be asked for the
  spot annotations directly (contact on the GEO record).

## Data dictionary

Regenerated by `notebooks/01_dataset_exploration.ipynb` into
`results/tables/library_inventory.csv`: per library, patient, slide, spots in
tissue, genes detected, median genes/UMI/mitochondrial fraction, feature panel,
and whether barcodes were recovered.
