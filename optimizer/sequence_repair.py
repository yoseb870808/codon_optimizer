"""Remove forbidden motifs from an optimized DNA sequence by local recoding.

Given a variant sequence and a list of forbidden motifs (literal or regex),
this module walks the hits and rewrites the codons that overlap each motif
with synonymous alternatives drawn from the host CUT (weighted by ``w_i``).
Rewrite attempts are capped per motif; unresolved hits are reported so the
caller can surface a warning. Each repair keeps the protein translation
identical by construction.

Motifs are specified as dicts with keys:
  - ``name`` (str): human-readable label
  - ``pattern`` (str): either a literal DNA subsequence or a regex (if the
    ``regex`` flag is set)
  - ``regex`` (bool, optional): treat ``pattern`` as a regex
  - ``also_scan_reverse_complement`` (bool, optional): also search the reverse
    complement (important for restriction sites, which are typically symmetric
    but safer to confirm)
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .utils import CODON_TABLE_STANDARD, chunk_codons, reverse_complement

logger = logging.getLogger(__name__)


@dataclass
class MotifHit:
    motif_name: str
    start: int
    end: int
    matched: str
    strand: str  # "+" or "-"


def _compile_motifs(motifs: Iterable[dict]) -> list[dict]:
    """Return motifs as ``(name, compiled_regex, scan_rc)`` tuples."""
    compiled: list[dict] = []
    for m in motifs:
        name = m.get("name", m.get("pattern", "?"))
        pattern = m["pattern"]
        if not m.get("regex", False):
            pattern = re.escape(pattern.upper())
        compiled.append({
            "name": name,
            "regex": re.compile(pattern),
            "scan_rc": bool(m.get("also_scan_reverse_complement", True)),
        })
    return compiled


def find_motifs(sequence: str, motifs: Iterable[dict]) -> list[MotifHit]:
    """Return every motif hit (both strands if requested) in the sequence.

    Hits are returned sorted by start position on the forward strand.
    """
    seq = sequence.upper()
    compiled = _compile_motifs(motifs)
    hits: list[MotifHit] = []
    for m in compiled:
        for match in m["regex"].finditer(seq):
            hits.append(MotifHit(m["name"], match.start(), match.end(), match.group(0), "+"))
        if m["scan_rc"]:
            rc = reverse_complement(seq)
            for match in m["regex"].finditer(rc):
                # Map RC coordinates back to forward strand
                start = len(seq) - match.end()
                end = len(seq) - match.start()
                hits.append(MotifHit(m["name"], start, end, match.group(0), "-"))
    hits.sort(key=lambda h: h.start)
    return hits


def _codon_span(nt_start: int, nt_end: int) -> tuple[int, int]:
    """Return the (inclusive) codon index range that overlaps [nt_start, nt_end)."""
    first = nt_start // 3
    last = (nt_end - 1) // 3
    return first, last


def _sample_codon(
    aa: str,
    codon_table_by_aa: dict[str, tuple[list[str], np.ndarray]],
    rng: np.random.Generator,
    exclude: set[str] | None = None,
) -> str | None:
    """Sample a codon for ``aa`` weighted by w_i, optionally excluding some.

    Returns ``None`` if all candidates are excluded.
    """
    codons, probs = codon_table_by_aa[aa]
    if exclude:
        keep = np.array([c not in exclude for c in codons], dtype=bool)
        if not keep.any():
            return None
        probs_filtered = probs.copy()
        probs_filtered[~keep] = 0.0
        s = probs_filtered.sum()
        if s <= 0:
            candidates = [c for c, k in zip(codons, keep) if k]
            return rng.choice(candidates)
        probs_filtered /= s
        return rng.choice(codons, p=probs_filtered)
    return rng.choice(codons, p=probs)


def _by_aa(codon_table: pd.DataFrame) -> dict[str, tuple[list[str], np.ndarray]]:
    by_aa: dict[str, tuple[list[str], np.ndarray]] = {}
    for aa, group in codon_table.groupby("amino_acid"):
        codons = list(group["codon"])
        weights = group["w_i"].to_numpy(dtype=float)
        if weights.sum() <= 0:
            probs = np.ones_like(weights) / len(weights)
        else:
            probs = weights / weights.sum()
        by_aa[aa] = (codons, probs)
    return by_aa


def repair_sequence(
    dna: str,
    codon_table: pd.DataFrame,
    motifs: list[dict],
    rng: np.random.Generator | None = None,
    max_attempts_per_hit: int = 30,
) -> tuple[str, list[dict]]:
    """Rewrite codons to eliminate forbidden motifs, preserving translation.

    For each motif hit we locate the overlapping codons, then try up to
    ``max_attempts_per_hit`` synonymous substitutions (one codon at a time,
    starting from the middle of the overlap and working outward) until the
    motif is no longer detected in the window. Hits that cannot be removed
    after all attempts are recorded and returned.

    The protein translation is guaranteed unchanged. Stop codons and the
    single-synonym amino acids M and W are not touched.

    Args:
        dna: Input DNA sequence (must be divisible by 3).
        codon_table: Host CUT.
        motifs: List of motif specs (see module docstring).
        rng: Optional numpy Generator for deterministic retries.
        max_attempts_per_hit: Attempts per failing motif hit.

    Returns:
        Tuple of (possibly-modified DNA, unresolved hit descriptions).
    """
    if len(dna) % 3 != 0:
        raise ValueError(f"Sequence length {len(dna)} not divisible by 3")
    if not motifs:
        return dna, []

    rng = rng or np.random.default_rng()
    by_aa = _by_aa(codon_table)
    codons = chunk_codons(dna)
    protein = [CODON_TABLE_STANDARD.get(c, "X") for c in codons]
    unresolved: list[dict] = []

    # Iterate until no motif hits remain or no progress can be made
    previous_signature = None
    for _ in range(10):  # outer safety cap
        dna_now = "".join(codons)
        hits = find_motifs(dna_now, motifs)
        if not hits:
            break
        signature = tuple((h.motif_name, h.start, h.end) for h in hits)
        if signature == previous_signature:
            break
        previous_signature = signature

        for hit in hits:
            first_codon, last_codon = _codon_span(hit.start, hit.end)
            # Try centre-out ordering for codon candidates
            order = list(range(first_codon, last_codon + 1))
            mid = len(order) // 2
            ordered = sorted(order, key=lambda x: abs(x - order[mid]))

            solved = False
            for _ in range(max_attempts_per_hit):
                for idx in ordered:
                    aa = protein[idx]
                    if aa in ("M", "W", "*"):
                        continue
                    current = codons[idx]
                    exclude = {current}
                    new = _sample_codon(aa, by_aa, rng, exclude=exclude)
                    if new is None or new == current:
                        continue
                    # Tentatively swap
                    codons[idx] = new
                    # Re-scan just the local window
                    window_start = max(0, first_codon * 3 - 5)
                    window_end = min(len(dna_now), (last_codon + 1) * 3 + 5)
                    local = "".join(codons)[window_start:window_end]
                    still_hits = find_motifs(local, [
                        m for m in motifs if m["name"] == hit.motif_name
                    ])
                    if not still_hits:
                        solved = True
                        break
                    # Revert and try a different codon next round
                    codons[idx] = current
                if solved:
                    break
            if not solved:
                unresolved.append({
                    "motif": hit.motif_name,
                    "position": hit.start,
                    "length": hit.end - hit.start,
                    "matched": hit.matched,
                    "strand": hit.strand,
                })

    final_dna = "".join(codons)
    # One final pass — after all the swaps some earlier unresolved hits may
    # have been fixed by later edits
    remaining = find_motifs(final_dna, motifs)
    remaining_set = {(h.motif_name, h.start) for h in remaining}
    filtered_unresolved = [
        u for u in unresolved if (u["motif"], u["position"]) in remaining_set
    ]
    return final_dna, filtered_unresolved


DEFAULT_FORBIDDEN_MOTIFS: list[dict[str, Any]] = [
    # Type IIS sites commonly used in Golden Gate cloning
    {"name": "BsaI", "pattern": "GGTCTC"},
    {"name": "BsmBI", "pattern": "CGTCTC"},
    {"name": "SapI", "pattern": "GCTCTTC"},
    # Common Type II sites
    {"name": "EcoRI", "pattern": "GAATTC"},
    {"name": "BamHI", "pattern": "GGATCC"},
    {"name": "HindIII", "pattern": "AAGCTT"},
    {"name": "XhoI", "pattern": "CTCGAG"},
    {"name": "NdeI", "pattern": "CATATG"},
    {"name": "NcoI", "pattern": "CCATGG"},
    # Homopolymer runs (7+) that confound synthesis and can act as terminators
    # or replication-slippage sites. Matches the default QC threshold (max_run=6).
    {"name": "polyA7", "pattern": "A{7,}", "regex": True, "also_scan_reverse_complement": False},
    {"name": "polyT7", "pattern": "T{7,}", "regex": True, "also_scan_reverse_complement": False},
    {"name": "polyG7", "pattern": "G{7,}", "regex": True, "also_scan_reverse_complement": False},
    {"name": "polyC7", "pattern": "C{7,}", "regex": True, "also_scan_reverse_complement": False},
]
