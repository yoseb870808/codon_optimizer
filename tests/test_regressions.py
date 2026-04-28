"""Regression tests for targeted maintenance fixes.

Each test here pins a specific bug and prevents reintroduction:

1. ``test_report_generator_imports_under_quote_nesting``
   The strftime format string inside f-strings used the same quote character
   as the outer f-string, which is invalid until Python 3.12 (PEP 701). The
   module must import on 3.10/3.11 too.

2. ``test_derive_seed_is_process_stable``
   ``_derive_seed`` mixed in ``hash(salt)``, which is randomized per process
   when ``PYTHONHASHSEED`` is unset. Two subprocesses now must produce the
   same child seed for the same inputs.

3. ``test_simulated_anneal_handles_minimal_cds``
   ``simulated_anneal`` crashed on a CDS short enough that no codon was
   mutable (e.g. one-codon CDS, or window // 3 == 1).

4. ``test_repair_sequence_motif_without_name``
   ``repair_sequence`` raised ``KeyError`` when a motif dict omitted ``name``;
   the compiled-motif path already fell back to the pattern, so the filter
   path must do the same.
"""

from __future__ import annotations

import importlib
import subprocess
import sys

import numpy as np

from optimizer.reference_builder import build_codon_usage_table
from optimizer.utils import AA_TO_CODONS


# ---------------------------------------------------------------------------
# (1) report_generator must import & generate output without SyntaxError
# ---------------------------------------------------------------------------

def test_report_generator_imports_under_quote_nesting(tmp_path):
    """Catch the old f-string quote-nesting SyntaxError at import time."""
    module = importlib.import_module("optimizer.report_generator")
    importlib.reload(module)

    # Exercise both writers that previously contained the bad f-strings.
    import pandas as pd

    cut = pd.DataFrame(
        [
            {"codon": "ATG", "amino_acid": "M", "count": 5, "frequency": 1.0, "w_i": 1.0},
            {"codon": "TAA", "amino_acid": "*", "count": 5, "frequency": 1.0, "w_i": 1.0},
        ]
    )
    module.write_cut_tsv(
        cut,
        reference_info={"mode": "test", "n_genes": 1, "n_total_cds": 1},
        output_dir=tmp_path,
        host_gc=0.5,
    )
    module.write_summary_md(
        all_results=[],
        reference_info={"mode": "test", "n_genes": 1, "n_total_cds": 1},
        config={"optimization": {"mode": "weighted", "n_variants": 1}},
        output_dir=tmp_path,
    )

    # Both files exist and contain the timestamp marker.
    assert (tmp_path / "codon_usage_table.tsv").exists()
    summary = (tmp_path / "summary.md").read_text()
    assert "_Generated " in summary and "Z_" in summary


# ---------------------------------------------------------------------------
# (2) Deterministic seed derivation across processes
# ---------------------------------------------------------------------------

def test_derive_seed_is_deterministic_in_process():
    from optimizer.sequence_optimizer import _derive_seed

    a = _derive_seed(42, 0, "weighted")
    b = _derive_seed(42, 0, "weighted")
    assert a == b
    assert _derive_seed(42, 0, "weighted") != _derive_seed(42, 0, "harmonized")
    assert _derive_seed(42, 0, "weighted") != _derive_seed(42, 1, "weighted")
    assert _derive_seed(None, 0, "weighted") is None


def test_derive_seed_is_process_stable():
    """A separate Python process must yield the same child seed.

    Run with ``PYTHONHASHSEED=random`` (the modern default) so that any
    accidental reintroduction of ``hash(salt)`` would diverge.
    """
    code = (
        "from optimizer.sequence_optimizer import _derive_seed;"
        "print(_derive_seed(42, 0, 'weighted'));"
        "print(_derive_seed(42, 5, 'harmonized'));"
        "print(_derive_seed(7, 3, 'repair'))"
    )

    def run_subprocess() -> str:
        env_extra = {"PYTHONHASHSEED": "random"}
        result = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            env={**dict(__import__("os").environ), **env_extra},
            check=True,
        )
        return result.stdout

    out_a = run_subprocess()
    out_b = run_subprocess()
    assert out_a == out_b, f"Seeds drifted between processes:\nA={out_a!r}\nB={out_b!r}"


def test_optimize_outputs_stable_across_runs():
    """End-to-end check: same seed → same sequence in repeated optimize() calls."""
    from optimizer.sequence_optimizer import optimize

    cut = build_codon_usage_table(
        [
            "ATG"
            + "".join(
                syn for codons in AA_TO_CODONS.values()
                for syn in codons if syn not in ("TAA", "TAG", "TGA")
            )
            + "TAA"
        ]
    )
    protein = "MAKLMAKLMAKL"
    v1 = optimize(protein, cut, mode="weighted", n_variants=2, host_gc=0.5, seed=123)
    v2 = optimize(protein, cut, mode="weighted", n_variants=2, host_gc=0.5, seed=123)
    assert [x["dna_sequence"] for x in v1] == [x["dna_sequence"] for x in v2]


# ---------------------------------------------------------------------------
# (3) simulated_anneal on very short CDS
# ---------------------------------------------------------------------------

def test_simulated_anneal_handles_minimal_cds():
    """A CDS with no mutable codon must not crash; returns unchanged."""
    from optimizer.mrna_repair import simulated_anneal

    cut = build_codon_usage_table(["ATGGCGAAACTGTAA" * 5])

    # One-codon CDS → mutable_positions is empty (codon 0 is Met, never mutable)
    one_codon = "ATG"
    new_dna, info = simulated_anneal(
        one_codon, cut, window=150, mfe_threshold=-30.0, iterations=50,
        rng=np.random.default_rng(0),
    )
    assert new_dna == one_codon
    assert info.attempts == 0
    assert info.method == "annealing"
    assert info.dna == one_codon

    # window // 3 == 1 forces mutation_window_codons=1 → empty range
    multi_codon = "ATG" + "GCG" * 20  # 21 codons
    new_dna, info = simulated_anneal(
        multi_codon, cut, window=2, mfe_threshold=-30.0, iterations=50,
        rng=np.random.default_rng(0),
    )
    assert info.attempts == 0


# ---------------------------------------------------------------------------
# (4) repair_sequence tolerates motifs without "name"
# ---------------------------------------------------------------------------

def test_repair_sequence_motif_without_name():
    """A motif dict missing 'name' must not raise KeyError in the filter path."""
    from optimizer.sequence_repair import find_motifs, repair_sequence

    cut = build_codon_usage_table(
        [
            "ATG"
            + "".join(
                syn for codons in AA_TO_CODONS.values()
                for syn in codons if syn not in ("TAA", "TAG", "TGA")
            )
            + "TAA"
        ]
    )

    # Embed an EcoRI site so repair has actual work to do.
    dna = "ATG" + "GAA" + "TTC" + "GCGGCGGCGGCGGCGGCG" + "TAA"
    assert len(dna) % 3 == 0

    motifs_no_name = [{"pattern": "GAATTC"}]  # NO "name" key
    # find_motifs already worked via _compile_motifs's fallback; smoke-check it.
    hits = find_motifs(dna, motifs_no_name)
    assert hits, "fixture should contain a GAATTC hit"
    assert hits[0].motif_name == "GAATTC"  # falls back to pattern

    new_dna, unresolved = repair_sequence(
        dna, cut, motifs_no_name, rng=np.random.default_rng(0)
    )
    # Translation preserved
    from optimizer.utils import translate
    assert translate(new_dna) == translate(dna)
    # The motif label flowed correctly through the repair filter path
    assert all(u["motif"] == "GAATTC" for u in unresolved)
