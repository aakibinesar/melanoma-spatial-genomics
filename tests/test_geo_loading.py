"""Tests for src/geo_loading.py on tiny synthetic GEO-style folders that reproduce
the real deposit's quirks (shared slide matrix, missing barcodes, per-library
matrices). They do not need the real data."""
import gzip
import io
import json
import sys
from pathlib import Path

import numpy as np
import pytest
import scipy.io
import scipy.sparse as sp
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import geo_loading as gl  # noqa: E402

GRID = ["AAAA-1", "CCCC-1", "GGGG-1", "TTTT-1", "ACGT-1", "TGCA-1"]
GENES = [("ENSG1", "G1"), ("ENSG2", "G2"), ("ENSG3", "G3")]


def _gz(path, text_or_bytes):
    data = text_or_bytes.encode() if isinstance(text_or_bytes, str) else text_or_bytes
    with gzip.open(path, "wb") as fh:
        fh.write(data)


def _mtx_bytes(m):
    buf = io.BytesIO()
    scipy.io.mmwrite(buf, sp.coo_matrix(m), field="integer")
    return buf.getvalue()


def _png_bytes():
    buf = io.BytesIO()
    Image.fromarray(np.zeros((4, 4, 3), np.uint8)).save(buf, format="PNG")
    return buf.getvalue()


def _write_library(d, gsm, lib, matrix, barcodes, in_tissue, write_barcodes=True, write_features=True):
    pre = d / f"{gsm}_{lib}_"
    _gz(str(pre) + "matrix.mtx.gz", _mtx_bytes(matrix))
    if write_barcodes:
        _gz(str(pre) + "barcodes.tsv.gz", "\n".join(barcodes) + "\n")
    if write_features:
        _gz(str(pre) + "features.tsv.gz", "".join(f"{i}\t{s}\tGene Expression\n" for i, s in GENES))
    rows = ["barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres"]
    for k, b in enumerate(GRID):
        rows.append(f"{b},{int(b in in_tissue)},{k},{k},{100 + k},{200 + k}")
    _gz(str(pre) + "tissue_positions.csv.gz", "\n".join(rows) + "\n")
    _gz(str(pre) + "scalefactors_json.json.gz", json.dumps({"tissue_hires_scalef": 0.1, "tissue_lowres_scalef": 0.05, "spot_diameter_fullres": 10.0, "fiducial_diameter_fullres": 20.0}))
    _gz(str(pre) + "tissue_hires_image.png.gz", _png_bytes())
    _gz(str(pre) + "tissue_lowres_image.png.gz", _png_bytes())


@pytest.fixture
def geo_dir(tmp_path):
    counts = np.arange(1, 19).reshape(3, 6)          # genes x 6 spots, all non-zero
    # slide S1: two libraries share ONE slide-level matrix, split by in_tissue
    _write_library(tmp_path, "GSM1", "L1", counts, GRID, {"AAAA-1", "CCCC-1"})
    # L2 has the same data plus an extra all-zero column and NO barcodes/features files
    with_zero = np.hstack([counts, np.zeros((3, 1), int)])
    _write_library(tmp_path, "GSM1", "L2", with_zero, GRID, {"GGGG-1", "TTTT-1", "ACGT-1"},
                   write_barcodes=False, write_features=False)
    _gz(tmp_path / "GSE1_A2_features.tsv.gz", "".join(f"{i}\t{s}\tGene Expression\n" for i, s in GENES))
    # slide S2: a per-library matrix that already only holds that library's spots
    sub = counts[:, :2]
    _write_library(tmp_path, "GSM2", "L3", sub, GRID[:2], {"AAAA-1", "CCCC-1"})
    return tmp_path


def test_list_libraries(geo_dir):
    df = gl.list_libraries(geo_dir)
    assert list(df["library"]) == ["L1", "L2", "L3"]
    assert df.set_index("library").loc["L2", "has_barcodes"] == False  # noqa: E712
    assert df.set_index("library").loc["L1", "has_barcodes"] == True   # noqa: E712


def test_shared_matrix_is_split_by_in_tissue_without_double_counting(geo_dir):
    slide = gl.load_slide(geo_dir, "GSM1", load_images=False)
    a, b = slide["L1"], slide["L2"]
    assert set(a.obs["barcode"]) == {"AAAA-1", "CCCC-1"}
    assert set(b.obs["barcode"]) == {"GGGG-1", "TTTT-1", "ACGT-1"}
    assert not set(a.obs["barcode"]) & set(b.obs["barcode"])
    assert a.n_obs + b.n_obs == 5                      # not 2 * 6
    assert a.obs_names.is_unique and b.obs_names.is_unique


def test_matrix_orientation_and_values(geo_dir):
    a = gl.load_slide(geo_dir, "GSM1", load_images=False)["L1"]
    assert a.shape == (2, 3)                           # spots x genes
    # spot AAAA-1 is column 0 of the genes x spots matrix -> counts 1, 7, 13
    assert a[a.obs["barcode"] == "AAAA-1"].X.toarray().ravel().tolist() == [1, 7, 13]
    assert a.obsm["spatial"].shape == (2, 2)
    assert a.obsm["spatial"][0].tolist() == [200.0, 100.0]   # (x=col, y=row)


def test_missing_barcodes_recovered_from_sibling(geo_dir):
    b = gl.load_slide(geo_dir, "GSM1", load_images=False)["L2"]
    assert b.obs["barcodes_inferred"].all()
    assert b[b.obs["barcode"] == "GGGG-1"].X.toarray().ravel().tolist() == [3, 9, 15]
    assert list(b.var_names) == ["G1", "G2", "G3"]     # features came from the series-level file


def test_declared_in_tissue_barcode_with_no_nonzero_match_is_zero_filled(tmp_path):
    # TGCA-1's column (index 5) is all-zero in BOTH L1 and L2 - exactly the
    # real A2/A17 situation, where an all-zero column carries no fingerprint
    # to recover an identity from. L2 (no barcodes/features file) declares
    # TGCA-1 in_tissue anyway: since it can never be matched to a non-zero
    # column, it must appear as an explicit all-zero row, not be dropped.
    counts = np.arange(1, 19).reshape(3, 6)
    counts[:, 5] = 0
    _write_library(tmp_path, "GSM1", "L1", counts, GRID, {"AAAA-1", "CCCC-1"})
    _write_library(tmp_path, "GSM1", "L2", counts, GRID,
                   {"GGGG-1", "TTTT-1", "ACGT-1", "TGCA-1"},
                   write_barcodes=False, write_features=False)
    _gz(tmp_path / "GSE1_A2_features.tsv.gz", "".join(f"{i}\t{s}\tGene Expression\n" for i, s in GENES))

    b = gl.load_slide(tmp_path, "GSM1", load_images=False)["L2"]
    assert set(b.obs["barcode"]) == {"GGGG-1", "TTTT-1", "ACGT-1", "TGCA-1"}
    row = b[b.obs["barcode"] == "TGCA-1"]
    assert row.obs["expression_zero_filled"].item()
    assert row.X.toarray().ravel().tolist() == [0, 0, 0]
    assert not b[b.obs["barcode"] == "GGGG-1"].obs["expression_zero_filled"].item()


def test_inference_refuses_a_different_matrix():
    a = sp.csc_matrix(np.arange(1, 7).reshape(2, 3))
    b = sp.csc_matrix(np.arange(2, 8).reshape(2, 3))
    with pytest.raises(ValueError):
        gl.infer_missing_barcodes(b, a, ["x", "y", "z"])


def test_per_library_matrix_keeps_only_its_spots(geo_dir):
    c = gl.load_slide(geo_dir, "GSM2", load_images=False)["L3"]
    assert c.n_obs == 2 and set(c.obs["barcode"]) == {"AAAA-1", "CCCC-1"}


def test_images_and_scalefactors(geo_dir):
    a = gl.load_slide(geo_dir, "GSM1", load_images=True)["L1"]
    sp_entry = a.uns["spatial"]["L1"]
    assert sp_entry["scalefactors"]["tissue_lowres_scalef"] == 0.05
    assert sp_entry["images"]["hires"].shape == (4, 4, 3)
