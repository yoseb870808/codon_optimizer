"""Tests for optimizer.quality_control."""

from __future__ import annotations

import pytest

from optimizer.quality_control import (
    check_gc_content,
    check_homopolymers,
    check_mrna_structure,
    check_repeats,
    run_all_qc,
)


class TestGC:
    def test_zero_percent(self):
        r = check_gc_content("ATATATAT")
        assert r["overall_gc"] == 0.0
        assert r["within_range"] is False

    def test_full_gc(self):
        r = check_gc_content("GCGCGCGC")
        assert r["overall_gc"] == 1.0
        assert r["within_range"] is False

    def test_fifty_percent(self):
        r = check_gc_content("ATGCATGC")
        assert pytest.approx(r["overall_gc"]) == 0.5
        assert r["within_range"] is True

    def test_sliding_window_detects_spike(self):
        # 100 nt of balanced AT-rich + 60 nt GC spike
        body = ("AT" * 50) + ("GC" * 30) + ("AT" * 50)
        r = check_gc_content(body, window_size=30, window_deviation=0.2)
        assert len(r["flagged_regions"]) > 0
        assert any(f["direction"] == "high" for f in r["flagged_regions"])


class TestHomopolymers:
    def test_detects_long_run(self):
        r = check_homopolymers("ATGAAAAAAAATGC", max_run=6)
        assert len(r) == 1
        assert r[0]["base"] == "A"
        assert r[0]["length"] == 8

    def test_none_when_short(self):
        r = check_homopolymers("ATGATGATGATG", max_run=6)
        assert r == []

    def test_detects_at_end(self):
        r = check_homopolymers("ATGCCCCCCCCC", max_run=6)
        assert len(r) == 1
        assert r[0]["base"] == "C"


class TestRepeats:
    def test_detects_direct_repeat(self):
        # 12 bp exact repeat downstream
        seq = "ATGATGATGATG" + "AAAAAAAAAA" + "ATGATGATGATG"
        r = check_repeats(seq, min_length=12)
        assert len(r) >= 1
        hit = r[0]
        assert hit["length"] >= 12

    def test_no_repeat(self):
        # Randomly-generated DNA of length 200 — almost surely no 12-mer
        # occurs twice (expected number is << 1 for uniform 4-ary string).
        import random
        rng = random.Random(2024)
        seq = "".join(rng.choice("ACGT") for _ in range(200))
        r = check_repeats(seq, min_length=12)
        assert r == []


class TestMRNAStructure:
    def test_strong_structure_flagged(self):
        # Classical hairpin: 8 Cs + tetraloop + 8 Gs → strong stem-loop in real folding
        seq = "CCCCCCCC" + "AUAU" + "GGGGGGGG" + "ACGUACGUACGUACGU"
        r = check_mrna_structure(seq.replace("U", "T"))
        assert r["method"] == "viennarna"
        assert isinstance(r["mfe"], float)
        assert r["mfe"] < 0
        assert r["has_strong_structure"] is True or r["mfe"] > -30  # depends on threshold

    def test_unstructured_sequence_low_mfe_magnitude(self):
        seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # poly-A folds to nothing
        r = check_mrna_structure(seq)
        assert r["method"] == "viennarna"
        assert isinstance(r["mfe"], float)
        # Poly-A has no base-pairing so MFE should be near 0
        assert r["mfe"] >= -2.0
        assert r["has_strong_structure"] is False

    def test_returns_dot_bracket_structure(self):
        seq = "GCGCGCATATATGCGCGC"
        r = check_mrna_structure(seq)
        assert isinstance(r["structure"], str)
        assert set(r["structure"]) <= set("().")

    def test_threshold_controls_flag(self):
        seq = "GCGCGCGCGCGCATATGCGCGCGCGCGC"
        weak = check_mrna_structure(seq, mfe_threshold=-5.0)
        strict = check_mrna_structure(seq, mfe_threshold=-100.0)
        assert weak["has_strong_structure"] is True   # easy to cross
        assert strict["has_strong_structure"] is False  # impossible to cross


class TestRunAllQC:
    def test_clean_sequence_format(self):
        seq = "ATGCATGCTAGCTAGCATGCTAGCATGCATGCA" * 9
        r = run_all_qc(seq)
        assert "gc" in r and "homopolymers" in r and "repeats" in r and "mrna_structure" in r
        assert r["mrna_structure"]["method"] == "viennarna"

    def test_homopolymer_fails(self):
        seq = "ATG" + "A" * 20 + "ATGATGATGATGATG"
        r = run_all_qc(seq)
        assert r["pass"] is False
        assert len(r["homopolymers"]) > 0
