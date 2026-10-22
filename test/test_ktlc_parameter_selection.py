"""KTLC selection tests; replace DFT evaluations, not the selection code."""

import subprocess
import sys
from types import SimpleNamespace

import pytest

from qcforever import gamess_run, gaussian_run
from qcforever.gamess_run.GamessRunPack import GamessDFTRun
from qcforever.gaussian_run.GaussianRunPack import GaussianDFTRun
from qcforever.util import check_resource, job_timeout


class FixedTrialOptimizer:
    """Exercise the production search loop using inexpensive fixed proposals."""

    def __init__(self, f, pbounds, verbose, random_state):
        assert pbounds == {"mu": (0, 1)}
        assert random_state == 1
        self.evaluate = f
        self.bounds = pbounds
        self.res = []
        self.extra_mu = 0.8

    def set_bounds(self, new_bounds):
        self.bounds = new_bounds

    def maximize(self, init_points, n_iter, **kwargs):
        assert init_points + n_iter == 1
        lower, upper = self.bounds["mu"]
        mu = lower if lower == upper else self.extra_mu
        target = self.evaluate(mu)
        assert isinstance(target, (int, float))
        assert target == target and abs(target) != float("inf")
        self.res.append({"params": {"mu": mu}, "target": target})
        if lower != upper:
            self.extra_mu += 0.001

    @property
    def max(self):
        return max(self.res, key=lambda trial: trial["target"])


@pytest.fixture
def optimizer(monkeypatch):
    instances = []

    def create_optimizer(**kwargs):
        result = FixedTrialOptimizer(**kwargs)
        instances.append(result)
        return result

    monkeypatch.setitem(sys.modules, "bayes_opt", SimpleNamespace(
        BayesianOptimization=create_optimizer,
        UtilityFunction=lambda **kwargs: kwargs,
    ))
    return instances


@pytest.fixture(params=["gaussian", "gamess"])
def calculation(request, monkeypatch, optimizer):
    """Allow each test to supply Koopmans residuals in eV or failed trials."""
    backend = request.param
    monkeypatch.setattr(check_resource, "respec_cores", lambda count: count)
    run_class = GaussianDFTRun if backend == "gaussian" else GamessDFTRun
    functional = "lc-blyp" if backend == "gaussian" else "LC-BLYP"
    runner = run_class(functional, "STO-3G", 1, "energy", "dummy.xyz")
    result = SimpleNamespace(
        runner=runner, backend=backend, residual=lambda mu: 0.1,
        writer=runner.make_input, optimizer=optimizer,
    )

    def residual():
        mu = runner.para_functional[0] if runner.para_functional else 0.0
        value = result.residual(mu)
        if isinstance(value, BaseException):
            raise value
        return value

    def execute(*args, **kwargs):
        mu = runner.para_functional[0] if runner.para_functional else 0.0
        return "error" if result.residual(mu) is None else "normal"

    if backend == "gaussian":
        monkeypatch.setattr(runner, "chain_job", lambda *args, **kwargs: None)
        monkeypatch.setattr(gaussian_run.Exe_Gaussian, "exe_Gaussian", execute)
        monkeypatch.setattr(job_timeout, "run", lambda *args, **kwargs: None)
        monkeypatch.setattr(runner, "Extract_values", lambda *args, **kwargs: {
            "satkoopmans": [residual(), 0.0],
        })
        result.tune = lambda: runner.LC_para_BOopt(
            "test", "xyz", ["He"], [0.0], [0.0], [0.0], 0, 1)
    else:
        monkeypatch.setattr(runner, "make_input", lambda *args, **kwargs: None)
        monkeypatch.setattr(gamess_run.Exe_Gamess, "exe_Gamess", execute)
        monkeypatch.setattr(runner, "Extract_values", lambda *args, **kwargs: {
            "vip": [residual()], "vea": [0.0],
            "Alpha_MO": [[0.0], [0.0]], "Beta_MO": [[0.0], [0.0]],
        })
        result.tune = lambda: runner.LC_para_BOopt(
            "test", ["He"], [0.0], [0.0], [0.0], 0, 1)
    return result


@pytest.mark.parametrize("best_error", [0.005, 0.01])
def test_initial_convergence_uses_best_evaluated_mu(calculation, best_error):
    calculation.residual = lambda mu: -best_error if mu == 0.35 else 0.2
    assert calculation.tune() == [0.35]
    report = calculation.runner.functional_parameter_optimization
    assert report["state"] == "converged"
    assert report["selected_mu"] == 0.35
    assert report["error_ev"] == pytest.approx(best_error)
    assert report["tolerance_ev"] == 0.01
    assert report["valid_evaluations"] == 3
    assert report["selection_reason"]
    assert [item["params"]["mu"] for item in calculation.optimizer[0].res] == [
        0.15, 0.35, 0.65]


def test_unconverged_best_is_used_by_next_real_input(calculation, tmp_path, monkeypatch):
    calculation.residual = lambda mu: -0.04 if mu == 0.35 else 0.3
    runner = calculation.runner
    runner.para_functional = calculation.tune()
    assert runner.para_functional == [0.35]
    report = runner.functional_parameter_optimization
    assert report["state"] == "best_evaluated"
    assert report["error_ev"] == pytest.approx(0.04)
    assert report["valid_evaluations"] == 52
    assert len(calculation.optimizer[0].res) == 52
    assert calculation.optimizer[0].res[-1]["params"]["mu"] != 0.35

    monkeypatch.chdir(tmp_path)
    if calculation.backend == "gaussian":
        calculation.writer("selected", 0, 1, Newinput=True,
                           Mol_atom=["He"], X=[0.0], Y=[0.0], Z=[0.0])
        text = (tmp_path / "selected.com").read_text()
        assert "3/107=0350000000" in text
        assert "3/108=0350000000" in text
    else:
        calculation.writer("ENERGY", 0, 1, "selected.inp",
                           Mol_atom=["He"], X=[0.0], Y=[0.0], Z=[0.0])
        assert "MU=0.35" in (tmp_path / "selected.inp").read_text()


def test_failed_trial_cannot_beat_valid_large_error(calculation):
    calculation.residual = lambda mu: 4.0 if mu == 0.35 else None
    assert calculation.tune() == [0.35]
    assert calculation.optimizer[0].max["target"] == -10
    report = calculation.runner.functional_parameter_optimization
    assert report["state"] == "best_evaluated"
    assert report["error_ev"] == 4.0
    assert report["valid_evaluations"] == 1


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), ValueError("bad output")])
def test_no_valid_evaluations_use_explicit_default(calculation, value):
    calculation.residual = lambda mu: value
    assert calculation.tune() == []
    report = calculation.runner.functional_parameter_optimization
    assert report["state"] == "default_no_valid"
    assert report["selected_mu"] is None
    assert report["error_ev"] is None
    assert report["valid_evaluations"] == 0
    assert report["selection_reason"]


@pytest.mark.parametrize("exception", [
    job_timeout.QCforeverTimeoutError("deadline"), KeyboardInterrupt(),
])
def test_timeout_and_cancellation_propagate(calculation, exception):
    calculation.residual = lambda mu: exception
    with pytest.raises(type(exception)):
        calculation.tune()


@pytest.mark.parametrize("calculation", ["gaussian"], indirect=True)
@pytest.mark.parametrize("exception", [
    subprocess.CalledProcessError(1, ["formchk"]),
    job_timeout.QCforeverTimeoutError("conversion deadline"), KeyboardInterrupt(),
])
def test_failed_conversion_never_reads_stale_output(calculation, monkeypatch, tmp_path, exception):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "test.fchk").write_text("stale successful output")
    parsed = []

    def convert(*args, **kwargs):
        assert args[0] == ["formchk", "test.chk"]
        assert kwargs["check"] is True
        raise exception

    monkeypatch.setattr(job_timeout, "run", convert)
    monkeypatch.setattr(calculation.runner, "Extract_values", lambda *args, **kwargs:
                        parsed.append(True) or {"satkoopmans": [0.005, 0.0]})
    if isinstance(exception, subprocess.CalledProcessError):
        assert calculation.tune() == []
        assert calculation.runner.functional_parameter_optimization["valid_evaluations"] == 0
    else:
        with pytest.raises(type(exception)):
            calculation.tune()
    assert parsed == []


@pytest.mark.parametrize("calculation", ["gamess"], indirect=True)
def test_gamess_zero_mu_keeps_existing_default_semantics(calculation, monkeypatch):
    original_maximize = FixedTrialOptimizer.maximize

    def propose_zero(self, *args, **kwargs):
        self.extra_mu = 0.0
        return original_maximize(self, *args, **kwargs)

    monkeypatch.setattr(FixedTrialOptimizer, "maximize", propose_zero)
    calculation.residual = lambda mu: 0.005 if mu == 0.0 else 0.1
    assert calculation.tune() == []
    report = calculation.runner.functional_parameter_optimization
    assert report["state"] == "converged"
    # GAMESS evaluates its built-in default here, not a literal MU=0 input.
    assert report["selected_mu"] is None
    assert report["error_ev"] == pytest.approx(0.005)


@pytest.mark.parametrize("calculation", ["gamess"], indirect=True)
def test_gamess_invalid_second_homo_cannot_hide_in_max(calculation, monkeypatch):
    monkeypatch.setattr(calculation.runner, "Extract_values", lambda *args: {
        "vip": [0.005], "vea": [0.0],
        "Alpha_MO": [[0.0], [0.0]], "Beta_MO": [[float("nan")], [0.0]],
    })
    assert calculation.tune() == []
    assert calculation.runner.functional_parameter_optimization["valid_evaluations"] == 0


def test_real_bayesian_optimizer_selects_observed_best(monkeypatch, tmp_path):
    """Use installed BO and the full search budget, without executing DFT."""
    monkeypatch.setattr(check_resource, "respec_cores", lambda count: count)
    runner = GaussianDFTRun("lc-blyp", "STO-3G", 1, "energy", "dummy.xyz")
    observations = []

    def extract_values(*args, **kwargs):
        mu = runner.para_functional[0]
        error_ev = 0.03 + (mu - 0.42)**2
        observations.append((mu, error_ev))
        return {"satkoopmans": [error_ev, 0.0]}

    monkeypatch.setattr(runner, "chain_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(gaussian_run.Exe_Gaussian, "exe_Gaussian",
                        lambda *args: "normal")
    monkeypatch.setattr(job_timeout, "run", lambda *args, **kwargs: None)
    monkeypatch.setattr(runner, "Extract_values", extract_values)
    runner.para_functional = runner.LC_para_BOopt(
        "test", "xyz", ["He"], [0.0], [0.0], [0.0], 0, 1)

    best_mu, best_error = min(observations, key=lambda trial: trial[1])
    assert len(observations) == 52
    assert runner.para_functional == [best_mu]
    report = runner.functional_parameter_optimization
    assert report["state"] == "best_evaluated"
    assert report["valid_evaluations"] == 52
    assert report["error_ev"] == best_error
    monkeypatch.chdir(tmp_path)
    runner.make_input("selected", 0, 1, Newinput=True,
                      Mol_atom=["He"], X=[0.0], Y=[0.0], Z=[0.0])
    expected_iop = gaussian_run.Gen_IOPline_functional.functional_para(
        "lc-blyp", [best_mu])
    assert expected_iop in (tmp_path / "selected.com").read_text()
