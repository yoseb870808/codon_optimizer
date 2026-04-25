"""Opt-in regression test against a real GenBank file in ``input/``.

Skipped automatically unless a ``.gb``/``.gbk`` file is present under
``input/`` and a target FASTA is also there. This is intentionally generic:
drop in any host genome and any FASTA of targets, and the test will run.

Run explicitly:
    python -m pytest tests/test_real_genome_regression.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from optimizer.pipeline import run_pipeline
from optimizer.utils import translate


INPUT_DIR = Path(__file__).resolve().parent.parent / "input"


def _first(glob: str) -> Path | None:
    matches = sorted(INPUT_DIR.glob(glob))
    return matches[0] if matches else None


def _preferred_xlsx() -> Path | None:
    """Prefer *_normalized.xlsx over plain .xlsx when both are present."""
    normalized = sorted(INPUT_DIR.glob("*_normalized.xlsx"))
    if normalized:
        return normalized[0]
    return _first("*.xlsx")


GENOME = _first("*.gbk") or _first("*.gb")
TARGETS = _first("*.fasta") or _first("*.fa")
RNASEQ = _preferred_xlsx()


pytestmark = pytest.mark.skipif(
    GENOME is None or TARGETS is None,
    reason="No host genome or target FASTA present under input/ (skip on CI)",
)


def test_real_pipeline(tmp_path):
    """End-to-end pipeline against whatever real data the user dropped in."""
    out = tmp_path / "real_run"
    result = run_pipeline(
        genome_path=GENOME,
        target_path=TARGETS,
        rnaseq_path=RNASEQ,
        mode="weighted",
        n_variants=3,
        output_dir=out,
        seed=42,
    )

    info = result["reference_info"]
    assert info["n_genes"] >= 1
    assert 0.20 < info["host_gc"] < 0.80

    assert len(result["results"]) >= 1

    for res in result["results"]:
        rec = next(v for v in res["variants"] if v["is_recommended"])
        body = rec["dna_sequence"]
        if body.endswith(("TAA", "TAG", "TGA")):
            body = body[:-3]
        assert translate(body) == res["protein_seq"]
        # ViennaRNA must have produced an MFE
        assert isinstance(rec["mfe_5prime"], float)

    for name in (
        "combined_recommended.fasta",
        "optimization_report.tsv",
        "codon_usage_table.tsv",
        "summary.md",
    ):
        assert (out / name).exists(), f"missing {name}"
