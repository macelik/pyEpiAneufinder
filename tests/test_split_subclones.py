import pandas as pd

import pyEpiAneufinder as pea


def test_split_subclones_returns_barcode_and_subclone_columns():
    df = pd.DataFrame(
        {
            "seq": ["chr1", "chr1", "chr2", "chr2"],
            "start": [0, 100, 0, 100],
            "end": [100, 200, 100, 200],
            "cell1": [1, 2, 1, 0],
            "cell2": [1, 1, 2, 0],
            "cell3": [2, 2, 0, 0],
        }
    )

    clones = pea.split_subclones(df, split_val=2)
    assert list(clones.columns) == ["barcode", "subclone"]
    assert set(clones["barcode"]) == {"cell1", "cell2", "cell3"}
    assert len(clones) == 3
