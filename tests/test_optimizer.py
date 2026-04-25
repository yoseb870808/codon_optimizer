"""Tests for optimizer.sequence_optimizer."""

from __future__ import annotations

from collections import Counter

import pytest

from optimizer.metrics import calculate_cai
from optimizer.reference_builder import build_codon_usage_table
from optimizer.sequence_optimizer import optimize
from optimizer.utils import AA_TO_CODONS, CODON_TABLE_STANDARD, translate, chunk_codons


@pytest.fixture(scope="module")
def host_cut():
    # Strongly GC-rich biased host (prefers GC-ending codons where possible)
    preferred = {
        "A": "GCG", "C": "TGC", "D": "GAC", "E": "GAG", "F": "TTC",
        "G": "GGC", "H": "CAC", "I": "ATC", "K": "AAG", "L": "CTG",
        "M": "ATG", "N": "AAC", "P": "CCG", "Q": "CAG", "R": "CGC",
        "S": "AGC", "T": "ACC", "V": "GTG", "W": "TGG", "Y": "TAC",
    }
    protein = "MKFLILVASAGTTLMITGNPKLRDEWYQHCFSAGTTLMITGNPKLRDEWY"
    body = "".join(preferred[aa] for aa in protein)
    # Add some variation so alternatives have non-zero w_i
    alt = "".join(min(syns) for syns in AA_TO_CODONS.values() if len(syns) > 0)
    return build_codon_usage_table([body + "TAA", alt + "TAA"])


@pytest.fixture(scope="module")
def at_rich_source_cut():
    # A source organism biased toward AT-rich codons
    protein = "MKFLILVASAGTTLMITGNPKLRDEWYQHCFS"
    preferred_at = {
        "A": "GCA", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
        "G": "GGA", "H": "CAT", "I": "ATA", "K": "AAA", "L": "TTA",
        "M": "ATG", "N": "AAT", "P": "CCA", "Q": "CAA", "R": "AGA",
        "S": "TCA", "T": "ACA", "V": "GTA", "W": "TGG", "Y": "TAT",
    }
    body = "".join(preferred_at[aa] for aa in protein)
    return build_codon_usage_table([body + "TAA"])


@pytest.fixture
def short_protein():
    return "MKFLILVAS"


class TestMaxCai:
    def test_deterministic(self, host_cut, short_protein):
        v1 = optimize(short_protein, host_cut, mode="max_cai", n_variants=1, host_gc=0.5)
        v2 = optimize(short_protein, host_cut, mode="max_cai", n_variants=1, host_gc=0.5)
        assert v1[0]["dna_sequence"] == v2[0]["dna_sequence"]

    def test_cai_equals_one(self, host_cut, short_protein):
        v = optimize(short_protein, host_cut, mode="max_cai", host_gc=0.5)[0]
        assert pytest.approx(v["cai_score"], abs=1e-9) == 1.0

    def test_protein_preserved(self, host_cut, short_protein):
        v = optimize(short_protein, host_cut, mode="max_cai", host_gc=0.5)[0]
        # Strip trailing stop
        body = v["dna_sequence"][:-3] if CODON_TABLE_STANDARD.get(v["dna_sequence"][-3:]) == "*" else v["dna_sequence"]
        assert translate(body) == short_protein

    def test_n_variants_ignored(self, host_cut, short_protein):
        v = optimize(short_protein, host_cut, mode="max_cai", n_variants=5, host_gc=0.5)
        assert len(v) == 1


class TestWeighted:
    def test_three_variants(self, host_cut, short_protein):
        v = optimize(short_protein, host_cut, mode="weighted", n_variants=3, host_gc=0.5, seed=0)
        assert len(v) == 3

    def test_different_seeds_differ(self, host_cut, short_protein):
        v1 = optimize(short_protein, host_cut, mode="weighted", n_variants=3, host_gc=0.5, seed=1)
        v2 = optimize(short_protein, host_cut, mode="weighted", n_variants=3, host_gc=0.5, seed=2)
        seqs1 = {x["dna_sequence"] for x in v1}
        seqs2 = {x["dna_sequence"] for x in v2}
        assert seqs1 != seqs2

    def test_reproducibility(self, host_cut, short_protein):
        v1 = optimize(short_protein, host_cut, mode="weighted", n_variants=3, host_gc=0.5, seed=42)
        v2 = optimize(short_protein, host_cut, mode="weighted", n_variants=3, host_gc=0.5, seed=42)
        assert [x["dna_sequence"] for x in v1] == [x["dna_sequence"] for x in v2]

    def test_protein_preserved(self, host_cut):
        protein = "MKFLILVASAGTTLMITGNPKLRDEWYQHCFS"
        results = optimize(protein, host_cut, mode="weighted", n_variants=5, host_gc=0.5, seed=7)
        for v in results:
            body = v["dna_sequence"][:-3]
            assert translate(body) == protein

    def test_high_wi_preference(self, host_cut):
        # Over many Leu-only rounds, CTG should dominate since host_cut is biased toward CTG
        protein = "L" * 200
        v = optimize(protein, host_cut, mode="weighted", n_variants=1, host_gc=0.5, seed=3)[0]
        body = v["dna_sequence"][:-3]
        counts = Counter(chunk_codons(body))
        assert counts["CTG"] > counts["TTA"]
        assert counts["CTG"] > counts["CTC"]


class TestHarmonized:
    def test_requires_dna(self, host_cut, short_protein):
        with pytest.raises(ValueError, match="harmonized"):
            optimize(short_protein, host_cut, mode="harmonized", host_gc=0.5)

    def test_no_source_cut_estimates(self, host_cut, caplog):
        protein = "MKFLIL"
        # DNA with specific codons
        dna = "ATGAAATTCCTGATCCTG"
        import logging
        with caplog.at_level(logging.WARNING):
            v = optimize(
                protein, host_cut, mode="harmonized",
                original_dna=dna, host_gc=0.5, seed=7, n_variants=1,
            )
        assert len(v) == 1
        # Warning about source CUT
        assert any("source CUT" in m for m in caplog.messages)

    def test_rarity_preserved(self, host_cut, at_rich_source_cut):
        # Source uses rare codons (e.g. TTA for Leu). Harmonized host selection
        # should pick the codon that occupies the same rank in host CUT.
        protein = "MKFLILVASAGTTLMITGNPKLRDEWYQHCFS"
        # Build DNA using AT-rich preferred codons
        preferred_at = {
            "A": "GCA", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
            "G": "GGA", "H": "CAT", "I": "ATA", "K": "AAA", "L": "TTA",
            "M": "ATG", "N": "AAT", "P": "CCA", "Q": "CAA", "R": "AGA",
            "S": "TCA", "T": "ACA", "V": "GTA", "W": "TGG", "Y": "TAT",
        }
        dna = "".join(preferred_at[aa] for aa in protein)
        v = optimize(
            protein, host_cut, mode="harmonized",
            original_dna=dna, source_cut=at_rich_source_cut,
            host_gc=0.5, seed=0, n_variants=1,
        )[0]
        body = v["dna_sequence"][:-3]
        # Protein preserved
        assert translate(body) == protein

    def test_protein_preserved(self, host_cut, at_rich_source_cut):
        protein = "MKFLILVASAGTTL"
        dna = "ATG" + "AAATTCCTGATCCTGGTGGCCTCCGCCGGCACCACCCTG"
        assert len(dna) == len(protein) * 3
        v = optimize(
            protein, host_cut, mode="harmonized",
            original_dna=dna, source_cut=at_rich_source_cut,
            host_gc=0.5, seed=0, n_variants=3,
        )
        for var in v:
            body = var["dna_sequence"][:-3]
            assert translate(body) == protein


class TestVariantScoring:
    def test_sorted_by_composite_score(self, host_cut, short_protein):
        v = optimize(short_protein, host_cut, mode="weighted", n_variants=5, host_gc=0.5, seed=11)
        scores = [x["composite_score"] for x in v]
        assert scores == sorted(scores, reverse=True)

    def test_exactly_one_recommended(self, host_cut, short_protein):
        v = optimize(short_protein, host_cut, mode="weighted", n_variants=5, host_gc=0.5, seed=13)
        flagged = [x for x in v if x["is_recommended"]]
        assert len(flagged) == 1
        assert flagged[0] is v[0]


class TestInvalidInputs:
    def test_unknown_mode(self, host_cut, short_protein):
        with pytest.raises(ValueError, match="Unknown optimization mode"):
            optimize(short_protein, host_cut, mode="wtf")

    def test_empty_protein(self, host_cut):
        with pytest.raises(ValueError, match="Empty protein"):
            optimize("", host_cut, mode="max_cai")
