"""Tests for optimizer.sequence_repair."""

from __future__ import annotations

import numpy as np
import pytest

from optimizer.reference_builder import build_codon_usage_table
from optimizer.sequence_optimizer import optimize
from optimizer.sequence_repair import (
    DEFAULT_FORBIDDEN_MOTIFS,
    find_motifs,
    repair_sequence,
)
from optimizer.utils import translate


@pytest.fixture(scope="module")
def biased_cut():
    """A CUT biased toward one synonym per AA but with alternatives available."""
    protein_body = "MAKLMAKLMAKLMAKL" * 10 + "GSTCDERFHIPNWYQV"
    # Build two reference genes that together use every synonym at least once
    from optimizer.utils import AA_TO_CODONS
    alt = "".join(
        syn for codons in AA_TO_CODONS.values() for syn in codons if syn not in ("TAA", "TAG", "TGA")
    )
    return build_codon_usage_table([
        "ATGGCGAAACTGATGGCGAAACTGATGGCGAAACTG" * 5 + "TAA",
        "ATG" + alt + "TAA",
    ])


class TestFindMotifs:
    def test_literal_match(self):
        seq = "ATGGAATTCAAAATG"  # contains EcoRI
        hits = find_motifs(seq, [{"name": "EcoRI", "pattern": "GAATTC"}])
        assert len(hits) >= 1
        assert any(h.motif_name == "EcoRI" for h in hits)

    def test_reverse_complement(self):
        # RC of "GGTCTC" (BsaI) is "GAGACC"
        seq = "ATG" + "GAGACC" + "AAAATG"
        hits = find_motifs(seq, [{"name": "BsaI", "pattern": "GGTCTC"}])
        assert any(h.strand == "-" for h in hits)

    def test_no_match(self):
        seq = "ATGATGATGATGATGATG"
        hits = find_motifs(seq, [{"name": "EcoRI", "pattern": "GAATTC"}])
        assert hits == []

    def test_regex(self):
        seq = "ATGAAAATTTTTGCAT"
        hits = find_motifs(seq, [
            {"name": "longT", "pattern": "T{5,}", "regex": True,
             "also_scan_reverse_complement": False}
        ])
        assert len(hits) == 1


class TestRepairSequence:
    def test_removes_ecori_site(self, biased_cut):
        # A protein that, if naively encoded, includes GAATTC
        # Amino acids E-F: GAA-TTC → spans boundary produces "GAATTC"
        protein = "MEFLAK"
        # Start by building a DNA with an explicit EcoRI site
        dna = "ATG" + "GAA" + "TTC" + "CTT" + "GCT" + "AAA" + "TAA"
        assert "GAATTC" in dna
        rng = np.random.default_rng(0)
        new_dna, unresolved = repair_sequence(
            dna,
            biased_cut,
            [{"name": "EcoRI", "pattern": "GAATTC"}],
            rng=rng,
        )
        assert "GAATTC" not in new_dna
        assert unresolved == []
        # Translation preserved
        assert translate(new_dna).rstrip("*") == protein

    def test_preserves_translation(self, biased_cut):
        dna = "ATGGAATTCAAAAAAAAAAAAAAAATAA"  # two problems: EcoRI and polyA
        # Pad to codon-aligned length
        # We have: ATG-GAA-TTC-AAA-AAA-AAA-AAA-AAA-ATA-A... not divisible by 3. Let's use 27 nt.
        dna = "ATGGAATTCAAAAAAAAAAAAAAATAA"
        assert len(dna) % 3 == 0
        original_protein = translate(dna)
        rng = np.random.default_rng(1)
        new_dna, _ = repair_sequence(
            dna,
            biased_cut,
            [
                {"name": "EcoRI", "pattern": "GAATTC"},
                {"name": "polyA8", "pattern": "AAAAAAAA",
                 "also_scan_reverse_complement": False},
            ],
            rng=rng,
        )
        assert translate(new_dna) == original_protein

    def test_default_motifs_callable(self, biased_cut):
        dna = "ATG" + "GGTCTC" + "AAACTG" + "TAA"
        rng = np.random.default_rng(2)
        new_dna, _ = repair_sequence(dna, biased_cut, DEFAULT_FORBIDDEN_MOTIFS, rng=rng)
        assert "GGTCTC" not in new_dna


class TestOptimizeIntegration:
    def test_optimize_removes_default_motifs(self, biased_cut):
        # Weighted mode + oversampling + repair should rarely leave forbidden motifs
        protein = "MKFLILVASAGTTLMITGNPKLRDEWYQHCFSAGTTL"
        results = optimize(
            protein,
            biased_cut,
            mode="weighted",
            n_variants=3,
            host_gc=0.5,
            seed=42,
            config={
                "optimization": {
                    "mode": "weighted",
                    "n_variants": 3,
                    "oversample_factor": 3,
                    "repair_forbidden_motifs": True,
                    "use_default_forbidden": True,
                },
            },
        )
        assert len(results) == 3
        # At least the top-ranked variant should be clean
        top = results[0]
        assert len(top["forbidden_hits"]) == 0

    def test_forbidden_hits_field_present(self, biased_cut):
        protein = "MKFLIL"
        results = optimize(
            protein, biased_cut, mode="weighted", n_variants=1, host_gc=0.5, seed=0,
        )
        assert "forbidden_hits" in results[0]

    def test_per_variant_seed_reproducible(self, biased_cut):
        protein = "MKFLILVASAGTTLMITGNPKLRDEWYQHCFS"
        r1 = optimize(
            protein, biased_cut, mode="weighted",
            n_variants=3, host_gc=0.5, seed=1234,
        )
        r2 = optimize(
            protein, biased_cut, mode="weighted",
            n_variants=3, host_gc=0.5, seed=1234,
        )
        assert [v["dna_sequence"] for v in r1] == [v["dna_sequence"] for v in r2]

    def test_different_seed_changes_output(self, biased_cut):
        protein = "MKFLILVASAGTTLMITGNPKLRDEWYQHCFS"
        r1 = optimize(protein, biased_cut, mode="weighted",
                      n_variants=3, host_gc=0.5, seed=1)
        r2 = optimize(protein, biased_cut, mode="weighted",
                      n_variants=3, host_gc=0.5, seed=2)
        s1 = [v["dna_sequence"] for v in r1]
        s2 = [v["dna_sequence"] for v in r2]
        assert s1 != s2
