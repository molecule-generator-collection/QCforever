#!/usr/bin/env python3
"""Compare Top-3 predicted molecules with reference IR/Raman spectra.

The input workbook is expected to contain the ``Molecules`` and ``Spectra``
sheets used by ``top3_examples_with_spectra.xlsx``.  For every molecule and
Top-k predicted SMILES, this script:

1. creates a reproducible 3D SDF structure with RDKit;
2. writes the sample's reference IR and Raman peak files;
3. performs an xTB conformer search and DFT geometry optimization with
   ``opt optconf=xtb freq=IR.dat,Raman.dat``; and
4. incrementally saves similarity/difference metrics to CSV and full QC output
   to JSON, so a long calculation can be resumed safely.

Gaussian must be installed and configured in the execution environment.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sys
import zipfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any
from xml.etree import ElementTree


@dataclass(frozen=True)
class Candidate:
    selection_group: str
    sample_id: str
    reference_smiles: str
    rank: int
    predicted_smiles: str

    @property
    def key(self) -> str:
        return f"{self.sample_id}::top{self.rank}"


DEFAULT_WORKBOOK = Path("/home/sumita/QCforever/top3_examples_with_spectra.xlsx")


def _require_columns(headers: dict[str, int], required: set[str], sheet: str) -> None:
    missing = sorted(required - set(headers))
    if missing:
        raise ValueError(f"Sheet {sheet!r} is missing columns: {', '.join(missing)}")


def _excel_column_index(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference)
    if letters is None:
        raise ValueError(f"Invalid Excel cell reference: {reference}")
    index = 0
    for character in letters.group(0):
        index = index * 26 + ord(character) - ord("A") + 1
    return index - 1


def _parse_excel_scalar(value: str | None, cell_type: str | None) -> Any:
    if value is None:
        return None
    if cell_type in ("str", "inlineStr"):
        return value
    if cell_type == "b":
        return value == "1"
    try:
        number = float(value)
    except ValueError:
        return value
    return int(number) if number.is_integer() else number


def read_xlsx_sheet(workbook_path: Path, sheet_name: str) -> list[list[Any]]:
    """Read cell values from one XLSX sheet using only the standard library."""
    spreadsheet_ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    relationship_ns = (
        "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    )
    package_relationship_ns = (
        "http://schemas.openxmlformats.org/package/2006/relationships"
    )
    with zipfile.ZipFile(workbook_path) as archive:
        workbook_root = ElementTree.fromstring(archive.read("xl/workbook.xml"))
        relationship_root = ElementTree.fromstring(
            archive.read("xl/_rels/workbook.xml.rels")
        )
        relationships = {
            item.attrib["Id"]: item.attrib["Target"]
            for item in relationship_root.findall(
                f"{{{package_relationship_ns}}}Relationship"
            )
        }
        sheet_target = None
        for sheet in workbook_root.findall(
            f".//{{{spreadsheet_ns}}}sheet"
        ):
            if sheet.attrib.get("name") == sheet_name:
                rel_id = sheet.attrib[f"{{{relationship_ns}}}id"]
                sheet_target = relationships[rel_id]
                break
        if sheet_target is None:
            raise ValueError(f"Workbook has no sheet named {sheet_name!r}")
        sheet_path = sheet_target.lstrip("/")
        if not sheet_path.startswith("xl/"):
            sheet_path = f"xl/{sheet_path}"

        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared_root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in shared_root.findall(f"{{{spreadsheet_ns}}}si"):
                shared_strings.append(
                    "".join(
                        node.text or ""
                        for node in item.iter(f"{{{spreadsheet_ns}}}t")
                    )
                )

        sheet_root = ElementTree.fromstring(archive.read(sheet_path))
        rows: list[list[Any]] = []
        for row_node in sheet_root.findall(f".//{{{spreadsheet_ns}}}row"):
            values: dict[int, Any] = {}
            for cell in row_node.findall(f"{{{spreadsheet_ns}}}c"):
                column = _excel_column_index(cell.attrib["r"])
                cell_type = cell.attrib.get("t")
                if cell_type == "inlineStr":
                    text = "".join(
                        node.text or ""
                        for node in cell.iter(f"{{{spreadsheet_ns}}}t")
                    )
                    values[column] = text
                    continue
                value_node = cell.find(f"{{{spreadsheet_ns}}}v")
                raw_value = value_node.text if value_node is not None else None
                if cell_type == "s" and raw_value is not None:
                    values[column] = shared_strings[int(raw_value)]
                else:
                    values[column] = _parse_excel_scalar(raw_value, cell_type)
            if values:
                row = [None] * (max(values) + 1)
                for column, value in values.items():
                    row[column] = value
                rows.append(row)
        return rows


def load_workbook_data(
    workbook_path: Path, top_k: int = 3
) -> tuple[list[Candidate], dict[str, list[tuple[float, float, float]]]]:
    """Load Top-k candidates and reference peaks from the supplied workbook."""
    molecule_rows = iter(read_xlsx_sheet(workbook_path, "Molecules"))
    try:
        molecule_headers = {
            str(value): index for index, value in enumerate(next(molecule_rows)) if value
        }
        required_molecule_columns = {
            "selection_group",
            "sample_id",
            "canonical_smiles",
            *(f"predicted_smiles_top{rank}" for rank in range(1, top_k + 1)),
        }
        _require_columns(molecule_headers, required_molecule_columns, "Molecules")

        candidates: list[Candidate] = []
        for row in molecule_rows:
            sample_id = str(row[molecule_headers["sample_id"]] or "").strip()
            if not sample_id:
                continue
            for rank in range(1, top_k + 1):
                smiles = str(
                    row[molecule_headers[f"predicted_smiles_top{rank}"]] or ""
                ).strip()
                if not smiles:
                    continue
                candidates.append(
                    Candidate(
                        selection_group=str(
                            row[molecule_headers["selection_group"]] or ""
                        ),
                        sample_id=sample_id,
                        reference_smiles=str(
                            row[molecule_headers["canonical_smiles"]] or ""
                        ),
                        rank=rank,
                        predicted_smiles=smiles,
                    )
                )

        spectrum_rows = iter(read_xlsx_sheet(workbook_path, "Spectra"))
        spectrum_headers = {
            str(value): index for index, value in enumerate(next(spectrum_rows)) if value
        }
        _require_columns(
            spectrum_headers,
            {"sample_id", "frequency", "IR", "Raman"},
            "Spectra",
        )
        spectra: dict[str, list[tuple[float, float, float]]] = {}
        for row in spectrum_rows:
            sample_id = str(row[spectrum_headers["sample_id"]] or "").strip()
            if not sample_id:
                continue
            values = (
                row[spectrum_headers["frequency"]],
                row[spectrum_headers["IR"]],
                row[spectrum_headers["Raman"]],
            )
            if any(value is None for value in values):
                raise ValueError(f"Incomplete spectrum row for sample {sample_id}")
            peak = tuple(float(value) for value in values)
            if not all(math.isfinite(value) for value in peak):
                raise ValueError(f"Non-finite spectrum value for sample {sample_id}")
            if peak[1] < 0 or peak[2] < 0:
                raise ValueError(f"Negative spectrum intensity for sample {sample_id}")
            spectra.setdefault(sample_id, []).append(peak)

        missing_spectra = sorted({item.sample_id for item in candidates} - set(spectra))
        if missing_spectra:
            raise ValueError(
                "No reference spectrum found for: " + ", ".join(missing_spectra)
            )
        for peaks in spectra.values():
            peaks.sort(key=lambda peak: peak[0])
        return candidates, spectra
    except StopIteration as exc:
        raise ValueError("Molecules or Spectra sheet is empty") from exc


def safe_name(value: str) -> str:
    """Convert identifiers to portable directory names."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    return cleaned or "sample"


def write_reference_spectra(
    reference_dir: Path, peaks: list[tuple[float, float, float]]
) -> tuple[Path, Path]:
    """Write QCforever-compatible two-column IR and Raman peak files."""
    reference_dir.mkdir(parents=True, exist_ok=True)
    ir_path = (reference_dir / "reference_ir.dat").resolve()
    raman_path = (reference_dir / "reference_raman.dat").resolve()
    ir_path.write_text(
        "".join(f"{frequency:.10f} {ir:.10f}\n" for frequency, ir, _ in peaks),
        encoding="utf-8",
    )
    raman_path.write_text(
        "".join(
            f"{frequency:.10f} {raman:.10f}\n"
            for frequency, _, raman in peaks
        ),
        encoding="utf-8",
    )
    return ir_path, raman_path


def smiles_to_sdf(smiles: str, output_path: Path, seed: int) -> None:
    """Generate an H-complete, force-field-relaxed 3D structure."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    molecule = Chem.AddHs(molecule)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed) % 2147483647
    if AllChem.EmbedMolecule(molecule, params) != 0:
        params.useRandomCoords = True
        params.randomSeed = (int(seed) + 1) % 2147483647
        if AllChem.EmbedMolecule(molecule, params) != 0:
            raise RuntimeError(f"3D embedding failed for SMILES: {smiles}")

    try:
        if AllChem.MMFFHasAllMoleculeParams(molecule):
            AllChem.MMFFOptimizeMolecule(molecule, maxIters=1000)
        else:
            AllChem.UFFOptimizeMolecule(molecule, maxIters=1000)
    except (RuntimeError, ValueError):
        # ETKDG geometry is still a valid starting point for Gaussian.
        pass

    output_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_path))
    try:
        writer.write(molecule)
    finally:
        writer.close()


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if hasattr(value, "tolist"):
        return _jsonable(value.tolist())
    if hasattr(value, "item"):
        return _jsonable(value.item())
    return value


def candidate_seed(candidate: Candidate, base_seed: int) -> int:
    """Return an order-independent reproducible RDKit seed."""
    digest = hashlib.sha256(candidate.key.encode("utf-8")).digest()
    return (int.from_bytes(digest[:4], "big") + int(base_seed)) % 2147483647


def run_qcforever(
    sdf_path: Path,
    ir_path: Path,
    raman_path: Path,
    functional: str,
    basis: str,
    nproc: int,
    memory: str,
    calculation_timeout: int,
    job_timeout: int,
    stable: bool,
) -> dict[str, Any]:
    """Run one Gaussian/QCforever spectrum comparison."""
    from qcforever.gaussian_run import GaussianRunPack

    option = build_qcforever_option(ir_path, raman_path, stable)
    calculation = GaussianRunPack.GaussianDFTRun(
        functional,
        basis,
        nproc,
        option,
        sdf_path.name,
        restart=False,
        pklsave=True,
    )
    calculation.mem = memory
    calculation.timexe = calculation_timeout
    calculation.timejob = job_timeout

    original_directory = Path.cwd()
    try:
        os.chdir(sdf_path.parent)
        return calculation.run_gaussian()
    finally:
        os.chdir(original_directory)


def build_qcforever_option(
    ir_path: Path, raman_path: Path, stable: bool = False
) -> str:
    """Build the required conformer-search/optimization/frequency options."""
    option = f"opt optconf=xtb freq={ir_path},{raman_path}"
    return f"{option} stable" if stable else option


def write_calculated_spectrum(path: Path, qc_output: dict[str, Any]) -> None:
    frequencies = qc_output.get("freq", [])
    ir_values = qc_output.get("IR", [])
    raman_values = qc_output.get("Raman", [])
    if not frequencies:
        return
    with path.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["frequency", "IR", "Raman"])
        writer.writerows(zip(frequencies, ir_values, raman_values))


def result_record(
    candidate: Candidate,
    status: str,
    qc_output: dict[str, Any] | None = None,
    error: str = "",
) -> dict[str, Any]:
    similarities = (qc_output or {}).get("Spectrum_similarity", {})
    ir = similarities.get("IR", {})
    raman = similarities.get("Raman", {})
    ir_similarity = ir.get("similarity")
    raman_similarity = raman.get("similarity")
    available = [value for value in (ir_similarity, raman_similarity) if value is not None]
    mean_similarity = sum(available) / len(available) if available else None
    return {
        "key": candidate.key,
        **asdict(candidate),
        "status": status,
        "error": error,
        "IR_similarity": ir_similarity,
        "IR_difference": None if ir_similarity is None else 1.0 - ir_similarity,
        "IR_dissimilarity": ir.get("dissimilarity"),
        "IR_wasserstein_cm-1": ir.get("wasserstein_distance"),
        "Raman_similarity": raman_similarity,
        "Raman_difference": (
            None if raman_similarity is None else 1.0 - raman_similarity
        ),
        "Raman_dissimilarity": raman.get("dissimilarity"),
        "Raman_wasserstein_cm-1": raman.get("wasserstein_distance"),
        "mean_similarity": mean_similarity,
        "mean_difference": None if mean_similarity is None else 1.0 - mean_similarity,
        "qc_output": _jsonable(qc_output) if qc_output is not None else None,
    }


CSV_FIELDS = [
    "selection_group",
    "sample_id",
    "reference_smiles",
    "rank",
    "predicted_smiles",
    "status",
    "error",
    "IR_similarity",
    "IR_difference",
    "IR_dissimilarity",
    "IR_wasserstein_cm-1",
    "Raman_similarity",
    "Raman_difference",
    "Raman_dissimilarity",
    "Raman_wasserstein_cm-1",
    "mean_similarity",
    "mean_difference",
]


def save_results(output_dir: Path, records: list[dict[str, Any]]) -> None:
    """Atomically save compact CSV and full JSON outputs after each candidate."""
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "spectrum_comparison_results.json"
    json_tmp = output_dir / ".spectrum_comparison_results.json.tmp"
    json_tmp.write_text(
        json.dumps(records, ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )
    json_tmp.replace(json_path)

    csv_path = output_dir / "spectrum_comparison_results.csv"
    csv_tmp = output_dir / ".spectrum_comparison_results.csv.tmp"
    with csv_tmp.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    csv_tmp.replace(csv_path)


def load_previous_results(output_dir: Path) -> list[dict[str, Any]]:
    path = output_dir / "spectrum_comparison_results.json"
    if not path.exists():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError(f"Expected a JSON list in {path}")
    return data


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "workbook",
        type=Path,
        nargs="?",
        default=DEFAULT_WORKBOOK,
        help=f"input workbook (default: {DEFAULT_WORKBOOK})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("qcforever_spectrum_comparison"),
        help="calculation and result directory",
    )
    parser.add_argument("--functional", default="B3LYP")
    parser.add_argument(
        "--basis",
        default="6-31G(2df,p)",
        help="use the same method/basis as the reference spectra",
    )
    parser.add_argument("--nproc", type=int, default=8)
    parser.add_argument("--memory", default="8GB")
    parser.add_argument("--calculation-timeout", type=int, default=24 * 60 * 60)
    parser.add_argument("--job-timeout", type=int, default=48 * 60 * 60)
    parser.add_argument("--top-k", type=int, choices=(1, 2, 3), default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--stable", action="store_true", help="repair imaginary modes")
    parser.add_argument("--resume", action="store_true", help="skip existing result keys")
    parser.add_argument(
        "--sample-id",
        action="append",
        help="run only selected sample ID; may be repeated",
    )
    parser.add_argument("--limit", type=int, help="run only the first N candidates")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="prepare reference/SDF files without running Gaussian",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.nproc < 1:
        raise ValueError("--nproc must be at least 1")
    if args.limit is not None and args.limit < 1:
        raise ValueError("--limit must be at least 1")

    workbook_path = args.workbook.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    candidates, spectra = load_workbook_data(workbook_path, args.top_k)
    if args.sample_id:
        requested = set(args.sample_id)
        candidates = [item for item in candidates if item.sample_id in requested]
        missing = sorted(requested - {item.sample_id for item in candidates})
        if missing:
            raise ValueError("Unknown --sample-id: " + ", ".join(missing))
    if args.limit is not None:
        candidates = candidates[: args.limit]

    records = load_previous_results(output_dir) if args.resume else []
    record_indexes = {
        str(record.get("key")): index for index, record in enumerate(records)
    }
    completed = {
        str(record.get("key"))
        for record in records
        if record.get("IR_similarity") is not None
        and record.get("Raman_similarity") is not None
    }
    print(f"Loaded {len(candidates)} candidates from {workbook_path}")

    for index, candidate in enumerate(candidates, start=1):
        if candidate.key in completed:
            print(f"[{index}/{len(candidates)}] skip {candidate.key} (already saved)")
            continue

        sample_dir = output_dir / safe_name(candidate.sample_id)
        reference_paths = write_reference_spectra(
            sample_dir / "reference", spectra[candidate.sample_id]
        )
        candidate_dir = sample_dir / f"top{candidate.rank}"
        sdf_path = candidate_dir / "candidate.sdf"
        print(
            f"[{index}/{len(candidates)}] {candidate.key}: "
            f"{candidate.predicted_smiles}"
        )
        try:
            smiles_to_sdf(
                candidate.predicted_smiles,
                sdf_path,
                candidate_seed(candidate, args.seed),
            )
            if args.dry_run:
                record = result_record(candidate, "dry-run")
            else:
                qc_output = run_qcforever(
                    sdf_path=sdf_path,
                    ir_path=reference_paths[0],
                    raman_path=reference_paths[1],
                    functional=args.functional,
                    basis=args.basis,
                    nproc=args.nproc,
                    memory=args.memory,
                    calculation_timeout=args.calculation_timeout,
                    job_timeout=args.job_timeout,
                    stable=args.stable,
                )
                write_calculated_spectrum(candidate_dir / "calculated_spectrum.csv", qc_output)
                if "Spectrum_similarity" not in qc_output:
                    record = result_record(
                        candidate,
                        "error",
                        qc_output,
                        "QCforever output does not contain Spectrum_similarity; "
                        "inspect the saved qc_output and Gaussian status.",
                    )
                else:
                    record = result_record(
                        candidate, str(qc_output.get("log", "unknown")), qc_output
                    )
        except Exception as exc:  # keep the remaining long-running jobs alive
            record = result_record(
                candidate, "error", error=f"{type(exc).__name__}: {exc}"
            )
            print(f"  ERROR: {record['error']}", file=sys.stderr)
        if candidate.key in record_indexes:
            records[record_indexes[candidate.key]] = record
        else:
            record_indexes[candidate.key] = len(records)
            records.append(record)
        completed.add(candidate.key)
        save_results(output_dir, records)

    print(f"Results: {output_dir / 'spectrum_comparison_results.csv'}")
    print(f"Full QC output: {output_dir / 'spectrum_comparison_results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
