"""Tests for optimizer.report_generator."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from optimizer.metrics import compare_codons
from optimizer.reference_builder import build_codon_usage_table
from optimizer.report_generator import (
    write_combined_fasta,
    write_comparison_tsv,
    write_cut_tsv,
    write_optimized_fasta,
    write_report_tsv,
    write_summary_md,
)
from optimizer.sequence_optimizer import optimize


@pytest.fixture
def sample_variants():
    cut = build_codon_usage_table(["ATGGCGAAACTGTAA" * 40])
    protein = "MAKLMAKL"
    v = optimize(protein, cut, mode="weighted", n_variants=3, host_gc=0.5, seed=0)
    return protein, cut, v


@pytest.fixture
def sample_result_bundle(sample_variants):
    protein, cut, variants = sample_variants
    return {
        "header": "test_seq A test protein",
        "input_type": "protein",
        "length_aa": len(protein),
        "variants": variants,
        "qc": {
            "gc": {"within_range": True},
            "homopolymers": [],
            "repeats": [],
            "mrna_structure": {"mfe": -15.0, "has_strong_structure": False},
            "pass": True,
        },
        "original_cai": None,
        "original_gc": None,
        "host_gc": 0.5,
    }


class TestFASTAWriters:
    def test_per_seq_fasta(self, sample_variants, tmp_path):
        _, _, variants = sample_variants
        path = write_optimized_fasta("my_seq some description", variants, tmp_path)
        assert Path(path).exists()
        content = Path(path).read_text()
        # Three records
        assert content.count(">") == 3
        # Recommended first and contains metadata
        first_line = content.splitlines()[0]
        assert "optimized_recommended" in first_line
        assert "CAI=" in first_line and "GC=" in first_line and "score=" in first_line

    def test_combined_fasta(self, sample_result_bundle, tmp_path):
        bundles = [sample_result_bundle, {**sample_result_bundle, "header": "seq_two"}]
        path = write_combined_fasta(bundles, tmp_path)
        content = Path(path).read_text()
        assert content.count(">") == 2

    def test_output_dir_autocreate(self, sample_variants, tmp_path):
        target = tmp_path / "deep" / "nested" / "out"
        _, _, variants = sample_variants
        write_optimized_fasta("seq", variants, target)
        assert target.exists()


class TestTSVWriters:
    def test_report_tsv(self, sample_result_bundle, tmp_path):
        path = write_report_tsv([sample_result_bundle], tmp_path)
        df = pd.read_csv(path, sep="\t")
        assert len(df) == 1
        for col in (
            "sequence_name", "input_type", "length_aa",
            "recommended_cai", "recommended_gc", "composite_score",
            "optimization_mode", "qc_pass",
        ):
            assert col in df.columns

    def test_cut_tsv(self, sample_variants, tmp_path):
        _, cut, _ = sample_variants
        info = {"mode": "housekeeping", "n_genes": 10, "n_total_cds": 100}
        path = write_cut_tsv(cut, info, tmp_path, host_gc=0.52)
        text = Path(path).read_text()
        assert text.startswith("#")
        assert "Reference mode: housekeeping" in text
        assert "Host GC: 52.0%" in text
        # 64 codon rows + 1 header + comment lines
        data_lines = [l for l in text.splitlines() if not l.startswith("#") and l.strip()]
        assert len(data_lines) == 65  # 1 header + 64 codons

    def test_comparison_tsv(self, sample_variants, tmp_path):
        protein, cut, variants = sample_variants
        opt = variants[0]["dna_sequence"]
        # Build a fake "original" of the same length
        original = opt  # identical → all unchanged
        df = compare_codons(original, opt, cut)
        path = write_comparison_tsv("my_seq", df, tmp_path)
        out = pd.read_csv(path, sep="\t")
        for col in ("position", "amino_acid", "original_codon", "optimized_codon", "changed"):
            assert col in out.columns


class TestSummaryMD:
    def test_summary_contains_expected(self, sample_result_bundle, tmp_path):
        ref = {
            "mode": "rnaseq",
            "n_genes": 30,
            "n_total_cds": 500,
            "match_rate": 0.95,
            "host_gc": 0.48,
        }
        cfg = {"optimization": {"mode": "weighted", "n_variants": 3}}
        path = write_summary_md(
            [sample_result_bundle], ref, cfg, tmp_path,
            inputs={"genome": "mock.gbk", "target": "mock.fasta"},
        )
        text = Path(path).read_text(encoding="utf-8")
        assert "Run parameters" in text
        assert "Reference translatome" in text
        assert "rnaseq" in text
        assert "test_seq" in text
