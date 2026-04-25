import pandas as pd
import pytest
from unittest.mock import patch

from pyEpiAneufinder.plotting import karyo_gainloss


@pytest.fixture
def sample_dataframe():
    return pd.DataFrame(
        {
            "seq": ["chr1", "chr1", "chr2", "chr2"],
            "start": [0, 1000000, 0, 1000000],
            "end": [1000000, 2000000, 1000000, 2000000],
            "cell1": [2, 1, 0, 1],
            "cell2": [2, 1, 1, 0],
            "cell3": [1, 1, 1, 1],
            "cell4": [0, 2, 2, 1],
        }
    )


@patch("matplotlib.pyplot.savefig")
@patch("pyEpiAneufinder.plotting.dendrogram")
def test_karyo_gainloss_saves_expected_output(mock_dendrogram, mock_savefig, sample_dataframe, tmp_path):
    mock_dendrogram.return_value = {"leaves": [0, 1, 2]}
    outpath = tmp_path / "karyogram.png"

    karyo_gainloss(sample_dataframe, str(outpath), "Title")

    mock_savefig.assert_called_once()
    assert mock_savefig.call_args[0][0] == str(outpath)


@patch("matplotlib.pyplot.savefig")
@patch("pyEpiAneufinder.plotting.dendrogram")
def test_karyo_gainloss_accepts_annotation_dataframe(mock_dendrogram, mock_savefig, sample_dataframe, tmp_path):
    mock_dendrogram.return_value = {"leaves": [0, 1, 2, 3]}
    outpath = tmp_path / "karyogram_annotated.png"
    annot_dt = pd.DataFrame(
        {"annot": ["a", "a", "b", "b"]},
        index=["cell1", "cell2", "cell3", "cell4"],
    )

    karyo_gainloss(sample_dataframe, str(outpath), "Annotated", annot_dt=annot_dt)

    mock_savefig.assert_called_once()
    assert mock_savefig.call_args[0][0] == str(outpath)


@patch("matplotlib.pyplot.savefig")
@patch("pyEpiAneufinder.plotting.dendrogram")
def test_karyo_gainloss_accepts_non_positional_coordinate_columns(mock_dendrogram, mock_savefig, tmp_path):
    mock_dendrogram.return_value = {"leaves": [0, 1]}
    df = pd.DataFrame(
        {
            "cell1": [1, 1],
            "seq": ["chr1", "chr1"],
            "cell2": [0, 2],
            "start": [0, 100],
            "cell3": [1, 0],
            "end": [100, 200],
        }
    )
    outpath = tmp_path / "karyogram_ordered.png"

    karyo_gainloss(df, str(outpath), "Ordered")

    mock_savefig.assert_called_once()
    assert mock_savefig.call_args[0][0] == str(outpath)
