import inspect

import pyEpiAneufinder as pea


def test_public_api_exports_documented_functions():
    assert callable(pea.epiAneufinder)
    assert callable(pea.split_subclones)
    assert callable(pea.karyo_gainloss)
    assert callable(pea.plot_single_cell_profile)
    assert callable(pea.compute_aneuploidy_across_sample)
    assert callable(pea.compute_aneuploidy_by_chr)
    assert callable(pea.compute_heterogeneity_across_sample)
    assert callable(pea.compute_heterogeneity_by_chr)
    assert callable(pea.compute_cnv_burden_cell)


def test_public_api_functions_have_nonempty_docstrings():
    public_functions = [
        pea.epiAneufinder,
        pea.split_subclones,
        pea.karyo_gainloss,
        pea.plot_single_cell_profile,
        pea.compute_aneuploidy_across_sample,
        pea.compute_aneuploidy_by_chr,
        pea.compute_heterogeneity_across_sample,
        pea.compute_heterogeneity_by_chr,
        pea.compute_cnv_burden_cell,
    ]
    for func in public_functions:
        assert inspect.getdoc(func)
