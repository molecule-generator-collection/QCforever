from pathlib import Path

import pytest

import compare_top3_spectra as MODULE


def test_load_supplied_workbook():
    workbook = Path(
        "/Users/wood/Library/CloudStorage/Dropbox/ResearchDoc/"
        "MitsuishiPaper/outputs/top3_examples_with_spectra.xlsx"
    )
    if not workbook.exists():
        pytest.skip("Reference workbook is not available")

    candidates, spectra = MODULE.load_workbook_data(workbook)

    assert len(candidates) == 18
    assert len(spectra) == 6
    assert {item.rank for item in candidates} == {1, 2, 3}
    assert all(spectra[item.sample_id] for item in candidates)


def test_result_record_exposes_difference_metrics():
    candidate = MODULE.Candidate("Top-3 correct", "sample:1", "CC", 2, "CO")
    qc_output = {
        "Spectrum_similarity": {
            "IR": {
                "similarity": 0.8,
                "dissimilarity": 2.0,
                "wasserstein_distance": 10.0,
            },
            "Raman": {
                "similarity": 0.6,
                "dissimilarity": 3.0,
                "wasserstein_distance": 20.0,
            },
        }
    }

    result = MODULE.result_record(candidate, "normal", qc_output)

    assert result["IR_difference"] == pytest.approx(0.2)
    assert result["Raman_difference"] == pytest.approx(0.4)
    assert result["mean_similarity"] == pytest.approx(0.7)
    assert result["mean_difference"] == pytest.approx(0.3)


def test_reference_files_are_qcforever_two_column_format(tmp_path):
    ir_path, raman_path = MODULE.write_reference_spectra(
        tmp_path, [(200.0, 1.5, 2.5), (100.0, 3.5, 4.5)]
    )

    assert ir_path.read_text().splitlines() == [
        "200.0000000000 1.5000000000",
        "100.0000000000 3.5000000000",
    ]
    assert raman_path.read_text().splitlines() == [
        "200.0000000000 2.5000000000",
        "100.0000000000 4.5000000000",
    ]


def test_candidate_seed_does_not_depend_on_processing_order():
    candidate = MODULE.Candidate("Top-3 correct", "sample:1", "CC", 2, "CO")

    assert MODULE.candidate_seed(candidate, 42) == MODULE.candidate_seed(candidate, 42)
    assert MODULE.candidate_seed(candidate, 42) != MODULE.candidate_seed(candidate, 43)


def test_qcforever_option_requires_xtb_conformer_search_and_optimization(tmp_path):
    option = MODULE.build_qcforever_option(
        tmp_path / "reference_ir.dat", tmp_path / "reference_raman.dat"
    )

    assert option.startswith("opt optconf=xtb freq=")
    assert option.endswith("reference_ir.dat," + str(tmp_path / "reference_raman.dat"))


def test_default_workbook_is_under_home_sumita_qcforever():
    args = MODULE.build_parser().parse_args([])

    assert args.workbook == Path(
        "/home/sumita/QCforever/top3_examples_with_spectra.xlsx"
    )
