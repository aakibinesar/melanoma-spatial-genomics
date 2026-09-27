"""Load the GEO GSE298774 Visium deposit into AnnData objects, one per library.

The deposit is a flat archive of ``<GSM>_<library>_<file>.gz`` files rather than
a Space Ranger directory layout, and it has structure that a naive loader gets
wrong (verified against the real files, see docs/dataset_description.md):

* A GEO sample (slide) can hold several tissue sections ("libraries"). For most
  multi-library slides every library ships the *same slide-level matrix*, and
  the libraries differ only in which spots are flagged ``in_tissue`` in their
  ``tissue_positions.csv``. Loading each library's matrix whole would
  double-count spots, so a library is defined by its in-tissue spots.
* Some slides ship one matrix per library instead (barcode sets are disjoint
  and the matrix already only contains that library's spots).
* One library (A2) has no barcodes/features file. Its matrix has identical
  non-zero columns to a sibling library's matrix, so barcodes are recovered
  from that sibling (and the equality is checked, not assumed).
* Two gene panels exist (37,082 and 18,085 features).
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

_MATRIX_RE = re.compile(r"^(GSM\d+)_(.+)_matrix\.mtx\.gz$")
_POS_COLS = ["in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]


def _lines(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n") for line in fh]


def list_libraries(raw_dir: str | Path) -> pd.DataFrame:
    """One row per library found in ``raw_dir``: GEO sample (slide) and library id."""
    raw_dir = Path(raw_dir)
    rows = []
    for f in sorted(raw_dir.glob("GSM*_matrix.mtx.gz")):
        m = _MATRIX_RE.match(f.name)
        gsm, lib = m.groups()
        pre = raw_dir / f"{gsm}_{lib}_"
        rows.append({
            "gsm": gsm,
            "library": lib,
            "has_barcodes": Path(str(pre) + "barcodes.tsv.gz").exists(),
            "has_features": Path(str(pre) + "features.tsv.gz").exists(),
        })
    return pd.DataFrame(rows).sort_values(["gsm", "library"]).reset_index(drop=True)


def read_positions(path: str | Path) -> pd.DataFrame:
    """tissue_positions.csv(.gz) -> DataFrame indexed by barcode (with header row)."""
    with gzip.open(path, "rt") as fh:
        rows = list(csv.reader(fh))
    header, body = rows[0], rows[1:]
    if header[0] != "barcode":
        raise ValueError(f"unexpected tissue_positions header: {header}")
    df = pd.DataFrame(body, columns=header).set_index("barcode")
    return df[_POS_COLS].astype(int)


def read_features(path: str | Path) -> pd.DataFrame:
    rows = [line.split("\t") for line in _lines(Path(path))]
    return pd.DataFrame(rows, columns=["gene_id", "gene_symbol", "feature_type"])


def read_matrix(path: str | Path) -> sp.csc_matrix:
    """Matrix Market file -> genes x spots CSC matrix (as stored in the deposit)."""
    return scipy.io.mmread(str(path)).tocsc()


def infer_missing_barcodes(matrix: sp.csc_matrix, sib_matrix: sp.csc_matrix,
                           sib_barcodes: list[str]) -> tuple[list[int], list[str]]:
    """Recover barcodes for a matrix that has no barcodes file.

    Valid only when the matrix carries the same data as a sibling's matrix apart
    from all-zero columns. Returns (column indices, barcodes) for the non-zero
    columns, and raises if the two matrices are not the same data - the
    equivalence is verified entry by entry rather than assumed.
    """
    nz = np.flatnonzero(np.asarray(matrix.sum(axis=0)).ravel())
    sib_nz = np.flatnonzero(np.asarray(sib_matrix.sum(axis=0)).ravel())
    if len(nz) != len(sib_nz) or (matrix[:, nz] != sib_matrix[:, sib_nz]).nnz != 0:
        raise ValueError("matrix is not identical to the sibling's on its non-zero columns; "
                         "cannot infer barcodes")
    return list(nz), [sib_barcodes[i] for i in sib_nz]


def _read_image(path: Path):
    from PIL import Image
    with gzip.open(path, "rb") as fh:
        img = Image.open(io.BytesIO(fh.read()))
        img.load()
    return np.asarray(img.convert("RGB"))


def _slide_matrices(raw_dir: Path, gsm: str, libs: list[str]) -> dict[str, dict]:
    """Load each library's matrix/barcodes/features once, inferring barcodes if absent."""
    out = {}
    for lib in libs:
        pre = raw_dir / f"{gsm}_{lib}_"
        out[lib] = {
            "matrix": read_matrix(str(pre) + "matrix.mtx.gz"),
            "barcodes": _lines(Path(str(pre) + "barcodes.tsv.gz")) if Path(str(pre) + "barcodes.tsv.gz").exists() else None,
            "features": read_features(str(pre) + "features.tsv.gz") if Path(str(pre) + "features.tsv.gz").exists() else None,
            "inferred_barcodes": False,
        }
    complete = [l for l in libs if out[l]["barcodes"] is not None]
    for lib in libs:
        if out[lib]["barcodes"] is None:
            if not complete:
                raise ValueError(f"{gsm}/{lib}: no barcodes and no sibling to infer them from")
            sib = complete[0]
            cols, bcs = infer_missing_barcodes(out[lib]["matrix"], out[sib]["matrix"], out[sib]["barcodes"])
            out[lib]["matrix"], out[lib]["barcodes"], out[lib]["inferred_barcodes"] = out[lib]["matrix"][:, cols], bcs, True
    return out


def load_slide(raw_dir: str | Path, gsm: str, only_in_tissue: bool = True,
               load_images: bool = True, a2_features_path: str | Path | None = None
               ) -> dict[str, ad.AnnData]:
    """Return {library: AnnData} for one GEO sample (slide).

    obs: library, gsm, barcode, in_tissue, array_row, array_col, matrix_shared
    obsm['spatial']: (x, y) = (pxl_col_in_fullres, pxl_row_in_fullres)
    uns['spatial'][library]: scalefactors (+ hires/lowres images if requested)
    """
    raw_dir = Path(raw_dir)
    libs = sorted(list_libraries(raw_dir).query("gsm == @gsm")["library"])
    mats = _slide_matrices(raw_dir, gsm, libs)
    # a library missing features.tsv uses the series-level A2 features file
    for lib in libs:
        if mats[lib]["features"] is None:
            if a2_features_path is None:
                a2_features_path = next(raw_dir.glob("GSE*_features.tsv.gz"), None)
            if a2_features_path is None:
                raise FileNotFoundError(f"{gsm}/{lib}: no features file")
            mats[lib]["features"] = read_features(a2_features_path)

    out = {}
    for lib in libs:
        pre = raw_dir / f"{gsm}_{lib}_"
        pos = read_positions(str(pre) + "tissue_positions.csv.gz")
        m = mats[lib]
        barcodes = m["barcodes"]
        keep = np.array([b in pos.index and (pos.loc[b, "in_tissue"] == 1 or not only_in_tissue)
                         for b in barcodes])
        X = m["matrix"][:, np.flatnonzero(keep)].T.tocsr().astype(np.float32)
        kept_bc = [b for b, k in zip(barcodes, keep) if k]

        # A library with inferred barcodes only gets identities for its non-zero
        # columns (an all-zero column has no fingerprint to match against a
        # sibling). Its own tissue_positions.csv can still declare in-tissue
        # barcodes that never appear among those non-zero identities; by
        # elimination over the shared barcode universe, such a barcode has no
        # detected transcript in this library, so it is added as an explicit
        # all-zero row rather than silently dropped.
        zero_bc = []
        if m["inferred_barcodes"] and only_in_tissue:
            recovered = set(barcodes)
            zero_bc = [b for b in pos.index[pos["in_tissue"] == 1] if b not in recovered]
        if zero_bc:
            X = sp.vstack([X, sp.csr_matrix((len(zero_bc), X.shape[1]), dtype=np.float32)]).tocsr()
            kept_bc = kept_bc + zero_bc

        obs = pos.loc[kept_bc].copy()
        obs.insert(0, "barcode", kept_bc)
        obs.insert(0, "library", lib)
        obs.insert(0, "gsm", gsm)
        obs["in_tissue"] = obs["in_tissue"].astype(bool)
        obs["barcodes_inferred"] = m["inferred_barcodes"]
        obs["expression_zero_filled"] = [b in zero_bc for b in kept_bc]
        obs.index = [f"{lib}_{b}" for b in kept_bc]
        var = m["features"].set_index("gene_id", drop=False)
        adata = ad.AnnData(X=X, obs=obs, var=var.rename_axis(None))
        adata.var_names = pd.Index(var["gene_symbol"].values)
        adata.var_names_make_unique()
        adata.obsm["spatial"] = obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float)
        sf = json.loads(gzip.open(str(pre) + "scalefactors_json.json.gz", "rt").read())
        entry = {"scalefactors": sf}
        if load_images:
            entry["images"] = {"hires": _read_image(Path(str(pre) + "tissue_hires_image.png.gz")),
                               "lowres": _read_image(Path(str(pre) + "tissue_lowres_image.png.gz"))}
        adata.uns["spatial"] = {lib: entry}
        out[lib] = adata
    return out


def load_all(raw_dir: str | Path, **kwargs) -> dict[str, ad.AnnData]:
    """Load every library in the deposit, keyed by library id."""
    raw_dir = Path(raw_dir)
    out = {}
    for gsm in sorted(list_libraries(raw_dir)["gsm"].unique()):
        out.update(load_slide(raw_dir, gsm, **kwargs))
    return out
