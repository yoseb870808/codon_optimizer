"""Belt-and-suspenders tests for the protein-identity invariant.

These tests exist to make the safety net explicit: every code path that
modifies a candidate DNA sequence must preserve the translated protein.
If any of these tests start failing, treat it as a critical regression.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from optimizer.pipeline import _assert_protein_identity, _expected_protein, run_pipeline
from optimizer.utils import translate


def test_assert_protein_identity_passes_on_match():
    dna = "ATGAAAACCTAA"  # M-K-T-stop
    expected = _expected_protein(dna)
    # Should not raise
    _assert_protein_identity(dna, expected, context="test")


def test_assert_protein_identity_raises_on_mismatch():
    dna = "ATGAAAACCTAA"        # M-K-T
    bogus_protein = "MKR"        # threonine→arginine swap (not synonymous)
    with pytest.raises(RuntimeError, match="PROTEIN INTEGRITY VIOLATION"):
        _assert_protein_identity(dna, bogus_protein, context="test")


def test_assert_protein_identity_raises_on_length_mismatch():
    dna = "ATGAAATAA"            # M-K
    bogus_protein = "MKT"        # one residue too long
    with pytest.raises(RuntimeError, match="PROTEIN INTEGRITY VIOLATION"):
        _assert_protein_identity(dna, bogus_protein, context="test")


def test_expected_protein_strips_trailing_stop():
    assert _expected_protein("ATGAAATAA") == "MK"
    assert _expected_protein("ATGAAA") == "MK"   # no trailing stop


@pytest.fixture
def small_inputs(tmp_path):
    """Synthetic mock inputs for the integration sweep."""
    from tests.conftest import (
        MOCK_PROTEIN_FASTA,
        _build_mock_expression,
        _build_mock_genome,
    )
    inp = tmp_path / "in"
    inp.mkdir()
    gbk = inp / "mock_genome.gbk"
    meta = _build_mock_genome(gbk)
    expr = inp / "mock_expression.xlsx"
    _build_mock_expression(expr, meta)
    fasta = inp / "mock_protein.fasta"
    fasta.write_text(MOCK_PROTEIN_FASTA)
    return {"genome": gbk, "rnaseq": expr, "target": fasta}


def test_full_pipeline_protein_invariant_holds(small_inputs, tmp_path):
    """Run the pipeline end-to-end and confirm no integrity violation fires."""
    cfg = tmp_path / "test_cfg.yaml"
    cfg.write_text(
        "reference:\n"
        "  min_cds_length: 500\n"
        "  top_expressed_pct: 30\n"
        "  min_reference_genes: 3\n"
        "  fallback_min_genes: 3\n"
    )
    out = tmp_path / "out"
    result = run_pipeline(
        genome_path=small_inputs["genome"],
        target_path=small_inputs["target"],
        rnaseq_path=small_inputs["rnaseq"],
        mode="weighted",
        n_variants=3,
        output_dir=out,
        config_path=cfg,
        seed=12345,
        mrna_repair_mode="both",
    )
    # Independent re-check: every recommended variant translates correctly
    for res in result["results"]:
        for v in res["variants"]:
            body = v["dna_sequence"]
            if body.endswith(("TAA", "TAG", "TGA")):
                body = body[:-3]
            assert translate(body) == res["protein_seq"]


def test_pipeline_aborts_on_corruption(small_inputs, tmp_path, monkeypatch):
    """Simulate a buggy repair that returns a non-synonymous mutation;
    confirm the pipeline raises RuntimeError before writing any output."""
    from optimizer import pipeline as pipeline_mod
    from optimizer.mrna_repair import RepairResult

    def _broken_repair(dna, codon_table, **kwargs):
        # Mutate the second codon to a different amino acid (M-X-... → M-A-... if K)
        # Just swap the codon at position 1 to "GCT" (Ala), which is unlikely to be synonymous
        broken = dna[:3] + "GCT" + dna[6:]
        info = RepairResult(
            dna=broken, initial_mfe=-50.0, final_mfe=-10.0,
            initial_cai=0.5, final_cai=0.5,
            attempts=1, succeeded=True, method="targeted",
        )
        return broken, info

    monkeypatch.setattr(pipeline_mod, "repair_5prime_structure", _broken_repair)

    cfg = tmp_path / "test_cfg.yaml"
    cfg.write_text(
        "reference:\n"
        "  min_cds_length: 500\n"
        "  top_expressed_pct: 30\n"
        "  min_reference_genes: 3\n"
        "  fallback_min_genes: 3\n"
        "qc:\n"
        "  mrna_mfe_threshold: -1.0\n"  # force every variant into repair
    )
    out = tmp_path / "out"
    with pytest.raises(RuntimeError, match="PROTEIN INTEGRITY VIOLATION"):
        run_pipeline(
            genome_path=small_inputs["genome"],
            target_path=small_inputs["target"],
            rnaseq_path=small_inputs["rnaseq"],
            mode="weighted",
            n_variants=2,
            output_dir=out,
            config_path=cfg,
            seed=42,
            mrna_repair_mode="targeted",
        )
