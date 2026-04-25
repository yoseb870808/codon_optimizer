"""Tests for optimizer.metrics."""

from __future__ import annotations

import math
import random

import pytest

from optimizer.metrics import (
    calculate_cai,
    calculate_nc,
    compare_codons,
    composite_score,
)
from optimizer.reference_builder import build_codon_usage_table
from optimizer.utils import AA_TO_CODONS


@pytest.fixture
def biased_cut():
    """CUT strongly biased: uses the preferred codon for each AA."""
    # Build one gene with dominant codons
    preferred = {
        "A": "GCG", "C": "TGC", "D": "GAC", "E": "GAA", "F": "TTC",
        "G": "GGC", "H": "CAC", "I": "ATC", "K": "AAA", "L": "CTG",
        "M": "ATG", "N": "AAC", "P": "CCG", "Q": "CAG", "R": "CGC",
        "S": "AGC", "T": "ACC", "V": "GTG", "W": "TGG", "Y": "TAC",
    }
    aa_sequence = "MKFLILVASAGTTLMITGNPKLRDEWYQHCFS"
    body = "".join(preferred[aa] for aa in aa_sequence)
    # Add a minor representation of other codons so w_i != 0 elsewhere
    minor = "".join(
        min(syns) for syns in AA_TO_CODONS.values() for _ in range(1)
    )
    return build_codon_usage_table([body + "TAA", minor + "TAA"])


class TestCAI:
    def test_cai_max_sequence_is_one(self, biased_cut):
        # Build a sequence using only the most-preferred codons → CAI = 1.0
        preferred = {
            aa: biased_cut[biased_cut["amino_acid"] == aa]
            .sort_values("w_i", ascending=False)
            .iloc[0]["codon"]
            for aa in AA_TO_CODONS
            if aa != "*"
        }
        protein = "MKFLILVASAGTTLMIT"
        seq = "".join(preferred[aa] for aa in protein) + "TAA"
        cai = calculate_cai(seq, biased_cut)
        assert pytest.approx(cai, abs=1e-9) == 1.0

    def test_cai_random_less_than_one(self, biased_cut):
        rng = random.Random(42)
        protein = "MKFLILVASAGTTLMITGN"
        # Randomly pick synonyms uniformly
        seq = "".join(rng.choice(AA_TO_CODONS[aa]) for aa in protein) + "TAA"
        cai = calculate_cai(seq, biased_cut)
        assert 0.0 < cai < 1.0

    def test_cai_all_met_returns_one(self, biased_cut):
        seq = "ATGATGATG"  # all Met, no synonymous contribution
        assert calculate_cai(seq, biased_cut) == 1.0

    def test_cai_bad_length_raises(self, biased_cut):
        with pytest.raises(ValueError):
            calculate_cai("ATGAA", biased_cut)


class TestNC:
    def test_nc_biased_near_20(self):
        # All codons the same preferred one → Nc near 20
        preferred = {
            "A": "GCG", "C": "TGC", "D": "GAC", "E": "GAA", "F": "TTC",
            "G": "GGC", "H": "CAC", "I": "ATC", "K": "AAA", "L": "CTG",
            "N": "AAC", "P": "CCG", "Q": "CAG", "R": "CGC", "S": "AGC",
            "T": "ACC", "V": "GTG", "Y": "TAC",
        }
        protein = "".join(preferred.keys()) * 20
        seq = "ATG" + "".join(preferred[aa] for aa in protein) + "TGG" + "TAA"
        nc = calculate_nc(seq)
        assert nc <= 25

    def test_nc_uniform_near_61(self):
        # Use every synonym roughly equally across many AAs
        codons = []
        codons.append("ATG")
        for aa, syns in AA_TO_CODONS.items():
            if aa in ("M", "W", "*"):
                continue
            # Use each synonym 5 times
            for c in syns:
                codons.extend([c] * 5)
        codons.append("TGG")
        codons.append("TAA")
        seq = "".join(codons)
        nc = calculate_nc(seq)
        # Uniform usage → near maximum. Allow generous tolerance.
        assert nc > 50


class TestCompareCodons:
    def test_compare_detects_changes(self, biased_cut):
        # Two Leu-only sequences, one CTG and one CTT
        original = "ATG" + "CTG" * 4 + "TAA"
        optimized = "ATG" + "CTT" * 4 + "TAA"
        df = compare_codons(original, optimized, biased_cut)
        # 4 positions of Leu should be flagged changed; Met and stop not changed
        changed = df["changed"].tolist()
        assert changed[0] is False or not changed[0]  # Met same
        assert all(changed[1:5])
        assert not changed[5]  # stop same

    def test_length_mismatch_raises(self, biased_cut):
        with pytest.raises(ValueError):
            compare_codons("ATGATG", "ATG", biased_cut)


class TestCompositeScore:
    def test_weighted_sum(self):
        # CAI=1, GC matches host, no MFE supplied → score = 0.5 + 0.3 + 0.2 = 1.0
        s = composite_score(cai=1.0, gc_content_value=0.5, host_gc=0.5)
        assert pytest.approx(s) == 1.0

    def test_gc_penalty(self):
        # GC way off host → gc_score low
        s_match = composite_score(cai=1.0, gc_content_value=0.5, host_gc=0.5)
        s_off = composite_score(cai=1.0, gc_content_value=0.9, host_gc=0.5)
        assert s_off < s_match

    def test_mfe_penalty(self):
        s_low = composite_score(
            cai=1.0, gc_content_value=0.5, host_gc=0.5, mfe=-5.0, mfe_threshold=-30.0
        )
        s_high = composite_score(
            cai=1.0, gc_content_value=0.5, host_gc=0.5, mfe=-40.0, mfe_threshold=-30.0
        )
        assert s_low > s_high

    def test_custom_weights(self):
        s = composite_score(
            cai=1.0,
            gc_content_value=0.5,
            host_gc=0.5,
            weights={"cai_weight": 0.8, "gc_weight": 0.1, "mrna_weight": 0.1},
        )
        assert pytest.approx(s) == 1.0
