import pytest

from qcforever.gaussian_run import GaussianRunPack
from qcforever.util import Spectrum_similarity


def test_merge_degenerate_peaks_sums_intensity_and_averages_position():
    positions, intensities = Spectrum_similarity.merge_degenerate_peaks(
        [200.8, 100.4, 100.0, 200.0], [4.0, 2.0, 1.0, 3.0], tolerance=1.0
    )

    assert positions == pytest.approx([100.2, 200.4])
    assert intensities == pytest.approx([3.0, 7.0])


@pytest.mark.parametrize("spectrum_type", ["IR", "Raman", "NMR"])
def test_identical_spectra_have_unit_similarity(spectrum_type):
    result = Spectrum_similarity.spectrum_similarity(
        [1.0, 2.0], [2.0, 1.0], [1.0, 2.0], [2.0, 1.0], spectrum_type
    )

    assert result["similarity"] == pytest.approx(1.0)
    assert result["dissimilarity"] == pytest.approx(0.0)
    assert result["wasserstein_distance"] == pytest.approx(0.0)


def test_nmr_nearby_peaks_are_degenerate():
    result = Spectrum_similarity.spectrum_similarity(
        [1.005], [2.0], [1.0, 1.01], [1.0, 1.0], "NMR"
    )

    assert result["similarity"] == pytest.approx(1.0)
    assert result["wasserstein_distance"] == pytest.approx(0.0)


def test_read_data_supports_comments_and_blank_lines(tmp_path):
    spectrum = tmp_path / "ir.dat"
    spectrum.write_text("# position intensity\n\n100.0 2.5 peak-a\n", encoding="utf-8")

    assert Spectrum_similarity.read_data(spectrum) == [[100.0], [2.5]]


def test_zero_intensity_spectrum_is_rejected():
    with pytest.raises(ValueError, match="positive-intensity"):
        Spectrum_similarity.spectrum_similarity(
            [100.0], [1.0], [100.0], [0.0], "IR"
        )


def test_gaussian_output_contains_similarity_dictionary(tmp_path, monkeypatch):
    class FakeLogParser:
        def Check_task(self):
            return {"freq_line": 0, "nmr_line": 1}, [[], [" 1 H Isotropic =  31.0"]]

    ir_file = tmp_path / "ir.dat"
    raman_file = tmp_path / "raman.dat"
    nmr_file = tmp_path / "nmr.dat"
    ir_file.write_text("100.0 3.0\n", encoding="utf-8")
    raman_file.write_text("100.0 7.0\n", encoding="utf-8")
    nmr_file.write_text("1.0 1.0\n", encoding="utf-8")

    monkeypatch.setattr(
        GaussianRunPack.gaussian_run.parse_log,
        "parse_log",
        lambda _filename: FakeLogParser(),
    )
    monkeypatch.setattr(
        GaussianRunPack.gaussian_run.Get_FreqPro,
        "Extract_Freq",
        lambda _lines: ([99.6, 100.4], [1.0, 2.0], [3.0, 4.0]) + (0.0,) * 7,
    )
    monkeypatch.setattr(
        GaussianRunPack.gaussian_run.Get_FreqPro,
        "Extract_vibvec",
        lambda _lines: ([], [], []),
    )
    monkeypatch.setattr(
        GaussianRunPack.gaussian_run.AtomInfo,
        "One_TMS_refer",
        lambda _element, _functional, _basis: 32.0,
    )

    calculation = GaussianRunPack.GaussianDFTRun.__new__(
        GaussianRunPack.GaussianDFTRun
    )
    calculation.functional = "b3lyp"
    calculation.basis = "sto-3g"
    calculation.ref_spectrum_paths = {
        "IR": str(ir_file),
        "Raman": str(raman_file),
        "NMR": str(nmr_file),
    }
    output = calculation.Extract_values(
        "job", {"freq": True, "nmr": True}, [], []
    )

    assert set(output["Spectrum_similarity"]) == {"IR", "Raman", "NMR"}
    assert output["Spectrum_similarity"]["IR"]["similarity"] == pytest.approx(1.0)
    assert output["Spectrum_similarity"]["Raman"]["similarity"] == pytest.approx(1.0)
    assert output["Spectrum_similarity"]["NMR"]["similarity"] == pytest.approx(1.0)
