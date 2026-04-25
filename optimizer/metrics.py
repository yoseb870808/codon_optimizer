"""Codon usage metrics and per-position sequence comparison."""

from __future__ import annotations

import logging
import math
from collections import Counter, defaultdict

import pandas as pd

from .utils import AA_TO_CODONS, CODON_TABLE_STANDARD, chunk_codons, gc_content

logger = logging.getLogger(__name__)


def _codon_to_wi(codon_table: pd.DataFrame) -> dict[str, float]:
    return dict(zip(codon_table["codon"], codon_table["w_i"]))


def _codon_to_aa(codon_table: pd.DataFrame) -> dict[str, str]:
    return dict(zip(codon_table["codon"], codon_table["amino_acid"]))


def calculate_cai(sequence: str, codon_table: pd.DataFrame) -> float:
    """Codon Adaptation Index (Sharp & Li, 1987).

    ``CAI = exp( (1/L) * Σ ln(w_i) )`` where the sum is over codons that
    belong to amino acid families with more than one synonym (Met, Trp, and
    stop codons are excluded). Codons with ``w_i == 0`` use a floor of 0.01
    (handled at table build time).

    Args:
        sequence: Nucleotide sequence, length divisible by 3.
        codon_table: CUT from :func:`build_codon_usage_table`.

    Returns:
        CAI as a float in (0, 1].
    """
    if len(sequence) % 3 != 0:
        raise ValueError(
            f"Sequence length {len(sequence)} is not divisible by 3"
        )

    wi_map = _codon_to_wi(codon_table)
    aa_map = _codon_to_aa(codon_table)

    log_sum = 0.0
    n = 0
    for codon in chunk_codons(sequence):
        if codon not in wi_map:
            continue  # unknown codon (e.g. containing N) — skip
        aa = aa_map[codon]
        if aa == "*":
            continue
        if aa in ("M", "W"):
            continue  # single-synonym families excluded per Sharp & Li
        w = wi_map[codon]
        if w <= 0:
            w = 0.01
        log_sum += math.log(w)
        n += 1

    if n == 0:
        return 1.0  # all-Met / all-Trp edge case
    return math.exp(log_sum / n)


_FAMILY_SIZES = {aa: len(codons) for aa, codons in AA_TO_CODONS.items() if aa != "*"}


def calculate_nc(sequence: str) -> float:
    """Effective Number of Codons (Wright, 1990).

    Uses Wright's F_cf homozygosity estimator aggregated by family size:

    ``Nc = 2 + 9/F̄_2 + 1/F̄_3 + 5/F̄_4 + 3/F̄_6``

    where ``F̄_k`` is the average ``F`` across all amino-acid families of
    synonym count ``k``. Families with zero usage in the sequence are
    excluded from their group average (if a whole group is missing, the
    corresponding term is dropped). When a family has only one observed
    codon, Wright's unbiased estimator returns F = 1, producing Nc = 20 in
    the extreme bias case.

    Range: roughly 20 (max bias) to 61 (no bias).
    """
    if len(sequence) % 3 != 0:
        raise ValueError(
            f"Sequence length {len(sequence)} is not divisible by 3"
        )

    codons = [c for c in chunk_codons(sequence) if c in CODON_TABLE_STANDARD]
    # Group codons by amino acid
    aa_codon_counts: dict[str, Counter] = defaultdict(Counter)
    for c in codons:
        aa = CODON_TABLE_STANDARD[c]
        if aa == "*":
            continue
        aa_codon_counts[aa][c] += 1

    # Compute F for each amino acid with > 1 synonym
    f_by_size: dict[int, list[float]] = defaultdict(list)
    for aa, fam_counts in aa_codon_counts.items():
        syns = AA_TO_CODONS[aa]
        n = len(syns)
        if n == 1:
            continue
        total = sum(fam_counts.values())
        if total <= 1:
            # Too few observations — treat as maximum homozygosity (F=1.0)
            f_by_size[n].append(1.0)
            continue
        ss = sum((fam_counts[c] / total) ** 2 for c in syns)
        f = (n * ss - 1) / (n - 1)
        # Clamp: occasionally f can be slightly negative due to floating error
        if f < 1.0 / n:
            f = 1.0 / n
        f_by_size[n].append(f)

    # Wright formula: Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6
    contributions = {2: 9, 3: 1, 4: 5, 6: 3}
    # Met and Trp contribute 1 each (family size 1, 2 total) — already the "2 +"

    nc = 2.0  # Met + Trp
    for size, weight in contributions.items():
        f_list = f_by_size.get(size)
        if not f_list:
            # Interpolate missing group from neighbours per Wright's convention:
            # if F_3 missing, estimate as (F_2 + F_4) / 2
            if size == 3 and f_by_size.get(2) and f_by_size.get(4):
                f_bar = (
                    sum(f_by_size[2]) / len(f_by_size[2])
                    + sum(f_by_size[4]) / len(f_by_size[4])
                ) / 2
            else:
                # Skip this term when truly unobservable
                continue
        else:
            f_bar = sum(f_list) / len(f_list)
        if f_bar <= 0:
            continue
        nc += weight / f_bar

    # Upper-bound Nc at 61 per the classical definition
    if nc > 61.0:
        nc = 61.0
    return nc


def compare_codons(
    original_dna: str,
    optimized_dna: str,
    codon_table: pd.DataFrame,
) -> pd.DataFrame:
    """Per-position codon comparison between two same-length coding sequences.

    Both sequences must encode the same protein and have identical length.

    Returns:
        DataFrame with columns ``position`` (1-based), ``amino_acid``,
        ``original_codon``, ``original_w_i``, ``optimized_codon``,
        ``optimized_w_i``, ``changed``.
    """
    if len(original_dna) != len(optimized_dna):
        raise ValueError(
            f"Length mismatch: original {len(original_dna)} vs optimized {len(optimized_dna)}"
        )
    if len(original_dna) % 3 != 0:
        raise ValueError("Sequence length is not divisible by 3")

    wi_map = _codon_to_wi(codon_table)
    aa_map = _codon_to_aa(codon_table)

    orig_codons = chunk_codons(original_dna)
    opt_codons = chunk_codons(optimized_dna)

    rows = []
    for i, (o, p) in enumerate(zip(orig_codons, opt_codons), start=1):
        aa_o = aa_map.get(o, "X")
        aa_p = aa_map.get(p, "X")
        if aa_o != aa_p:
            logger.warning(
                "Position %d: AA mismatch (original %s → %s, optimized %s → %s)",
                i, o, aa_o, p, aa_p,
            )
        rows.append(
            {
                "position": i,
                "amino_acid": aa_o,
                "original_codon": o,
                "original_w_i": wi_map.get(o, float("nan")),
                "optimized_codon": p,
                "optimized_w_i": wi_map.get(p, float("nan")),
                "changed": o != p,
            }
        )
    return pd.DataFrame(rows)


def composite_score(
    cai: float,
    gc_content_value: float,
    host_gc: float,
    mfe: float | None = None,
    mfe_threshold: float = -30.0,
    weights: dict[str, float] | None = None,
) -> float:
    """Weighted composite score used to rank variants.

    - ``cai_score`` = ``cai`` (already in [0, 1])
    - ``gc_score`` = ``1 - min(1, |gc - host_gc| / 0.5)``
    - ``mrna_score`` = ``1`` if no MFE supplied (neutral), else
      ``1 - min(1, |mfe| / |mfe_threshold|)``

    Default weights: ``cai=0.5, gc=0.3, mrna=0.2``.
    """
    w = {"cai_weight": 0.5, "gc_weight": 0.3, "mrna_weight": 0.2}
    if weights:
        w.update(weights)

    cai_score = max(0.0, min(1.0, cai))
    gc_score = 1.0 - min(1.0, abs(gc_content_value - host_gc) / 0.5)
    if mfe is None:
        mrna_score = 1.0
    else:
        denom = abs(mfe_threshold) if mfe_threshold else 1.0
        mrna_score = 1.0 - min(1.0, abs(mfe) / denom)
    return (
        w["cai_weight"] * cai_score
        + w["gc_weight"] * gc_score
        + w["mrna_weight"] * mrna_score
    )
