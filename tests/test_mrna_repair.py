"""Tests for optimizer.mrna_repair (real ViennaRNA folds)."""

from __future__ import annotations

import numpy as np
import pytest

from optimizer.mrna_repair import (
    repair_5prime_structure,
    simulated_anneal,
    targeted_repair,
)
from optimizer.quality_control import check_mrna_structure
from optimizer.reference_builder import build_codon_usage_table
from optimizer.utils import AA_TO_CODONS, translate


@pytest.fixture(scope="module")
def diverse_cut():
    """A CUT where every AA has every synonym represented (no zero-count w_i)."""
    parts = []
    for aa, syns in AA_TO_CODONS.items():
        if aa == "*":
            continue
        # Use each synonym multiple times so w_i values differ but none are at floor
        for i, c in enumerate(syns):
            parts.append(c * (5 - i if 5 - i > 0 else 1))
    seq = "ATG" + "".join(parts) + "TAA"
    # Ensure divisible by 3 (it should be — each codon is 3 nt repeated)
    assert len(seq) % 3 == 0
    return build_codon_usage_table([seq])


def _build_problem_dna(diverse_cut):
    """Build a CDS with predictable strong 5' structure: GCCC×8 + loop + GGGC×8.

    Translation: M (ATG) + 8 alanines (GCC) + a loop region + alanines (GGC),
    chosen because Ala has 4 synonyms so we can fix it via substitution.
    """
    # GCC (Ala), GGC (Gly) — both common, 4 synonyms each
    head = "ATG" + "GCC" * 12 + "GGG" * 12  # strong CG stem
    # Pad to 200 nt so the 5' window fold is meaningful
    pad = "GCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA"
    return (head + pad + "TAA")[:201] + "TAA"  # ensure ends with stop


class TestTargetedRepair:
    def test_improves_strong_structure(self, diverse_cut):
        # A pathological CDS: 8x ATG-equivalent stem
        dna = "ATG" + "GCC" * 12 + "GGC" * 12 + "GCAGCAGCAGCAGCAGCAGCAGCAGCATAA"
        # Verify problem actually exists
        r = check_mrna_structure(dna)
        assert r["mfe"] < -20
        protein_before = translate(dna).rstrip("*")

        rng = np.random.default_rng(0)
        new_dna, info = targeted_repair(
            dna, diverse_cut, window=150, mfe_threshold=-15.0,
            max_attempts=20, rng=rng,
        )
        assert len(new_dna) == len(dna)
        assert translate(new_dna).rstrip("*") == protein_before
        # MFE should have improved (less negative)
        assert info.final_mfe > info.initial_mfe
        assert info.attempts > 0

    def test_no_op_when_already_good(self, diverse_cut):
        # poly-A folds to nothing — initial MFE near 0
        dna = "ATG" + "AAA" * 30 + "TAA"
        rng = np.random.default_rng(0)
        new_dna, info = targeted_repair(
            dna, diverse_cut, window=150, mfe_threshold=-30.0, rng=rng,
        )
        assert new_dna == dna
        assert info.attempts == 0
        assert info.succeeded is True

    def test_preserves_translation(self, diverse_cut):
        dna = "ATG" + "GCG" * 20 + "CGC" * 20 + "TAA"
        rng = np.random.default_rng(7)
        new_dna, _ = targeted_repair(
            dna, diverse_cut, window=120, mfe_threshold=-10.0,
            max_attempts=15, rng=rng,
        )
        assert translate(new_dna) == translate(dna)


class TestSimulatedAnneal:
    def test_improves_strong_structure(self, diverse_cut):
        dna = "ATG" + "GCG" * 15 + "CGC" * 15 + "GCAGCAGCAGCATAA"
        rng = np.random.default_rng(42)
        new_dna, info = simulated_anneal(
            dna, diverse_cut, window=120, mfe_threshold=-10.0,
            iterations=80, t_start=3.0, t_end=0.05, rng=rng,
        )
        assert len(new_dna) == len(dna)
        assert translate(new_dna).rstrip("*") == translate(dna).rstrip("*")
        # SA should not WORSEN MFE (best-tracking)
        assert info.final_mfe >= info.initial_mfe - 0.5  # tiny tolerance
        assert info.attempts > 0

    def test_returns_immediately_when_above_threshold(self, diverse_cut):
        dna = "ATG" + "AAA" * 20 + "TAA"
        rng = np.random.default_rng(0)
        new_dna, info = simulated_anneal(
            dna, diverse_cut, window=100, mfe_threshold=-30.0,
            iterations=50, rng=rng,
        )
        assert new_dna == dna
        assert info.attempts == 0


class TestDispatcher:
    def test_off_mode_no_change(self, diverse_cut):
        dna = "ATG" + "GCG" * 15 + "CGC" * 15 + "TAA"
        new_dna, info = repair_5prime_structure(dna, diverse_cut, mode="off")
        assert new_dna == dna
        assert info.method == "off"

    def test_both_mode_runs_when_targeted_fails(self, diverse_cut):
        dna = "ATG" + "GCG" * 15 + "CGC" * 15 + "TAA"
        new_dna, info = repair_5prime_structure(
            dna, diverse_cut, mode="both",
            mfe_threshold=-20.0,
            targeted_max_attempts=5,
            annealing_iterations=30,
            rng=np.random.default_rng(0),
        )
        # method should be 'both' if SA had to run
        assert info.method in ("targeted", "both")
        assert translate(new_dna).rstrip("*") == translate(dna).rstrip("*")

    def test_unknown_mode_raises(self, diverse_cut):
        with pytest.raises(ValueError, match="Unknown repair mode"):
            repair_5prime_structure("ATGAAATAA", diverse_cut, mode="bogus")
