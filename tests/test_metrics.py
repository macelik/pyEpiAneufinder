import pandas as pd

import pyEpiAneufinder as pea


def sample_result_table():
    return pd.DataFrame(
        {
            "seq": ["chr1", "chr1", "chr2", "chr2"],
            "start": [0, 100, 0, 100],
            "end": [100, 200, 100, 200],
            "cell1": [1, 2, 1, 0],
            "cell2": [1, 1, 2, 0],
        }
    )


def test_compute_cnv_burden_cell_returns_expected_columns():
    result = pea.compute_cnv_burden_cell(sample_result_table())
    assert list(result.columns) == ["barcodes", "cnv_burden"]
    assert len(result) == 2


def test_compute_aneuploidy_outputs_numeric_values():
    df = sample_result_table()
    assert isinstance(pea.compute_aneuploidy_across_sample(df), float)
    by_chr = pea.compute_aneuploidy_by_chr(df)
    assert list(by_chr.columns) == ["chr1", "chr2"]
    assert by_chr.shape == (1, 2)
    assert isinstance(pea.compute_heterogeneity_across_sample(df), float)
    heterogeneity_by_chr = pea.compute_heterogeneity_by_chr(df)
    assert list(heterogeneity_by_chr.columns) == ["chr1", "chr2"]
    assert heterogeneity_by_chr.shape == (1, 2)
