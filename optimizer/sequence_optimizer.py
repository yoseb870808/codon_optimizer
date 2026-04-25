"""Codon optimization engine: max_cai / weighted / harmonized modes."""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd

from .metrics import calculate_cai, composite_score
from .quality_control import check_mrna_structure
from .reference_builder import build_codon_usage_table
from .sequence_repair import DEFAULT_FORBIDDEN_MOTIFS, repair_sequence
from .utils import AA_TO_CODONS, CODON_TABLE_STANDARD, chunk_codons, gc_content

logger = logging.getLogger(__name__)

VALID_MODES = ("max_cai", "weighted", "harmonized")

# Stop codon to append after optimization. TAA tends to be the strongest
# terminator in most bacteria — a defensible default.
DEFAULT_STOP_CODON = "TAA"


def _codons_by_aa(codon_table: pd.DataFrame) -> dict[str, list[tuple[str, float]]]:
    """Return ``{aa: [(codon, w_i), ...]}`` sorted by w_i descending then codon."""
    by_aa: dict[str, list[tuple[str, float]]] = {}
    for aa, group in codon_table.groupby("amino_acid"):
        pairs = list(zip(group["codon"], group["w_i"]))
        pairs.sort(key=lambda x: (-x[1], x[0]))
        by_aa[aa] = pairs
    return by_aa


def optimize(
    protein_seq: str,
    codon_table: pd.DataFrame,
    mode: str = "weighted",
    n_variants: int = 3,
    host_gc: float | None = None,
    original_dna: str | None = None,
    source_cut: pd.DataFrame | None = None,
    seed: int | None = None,
    config: dict[str, Any] | None = None,
    append_stop: bool = True,
) -> list[dict]:
    """Generate one or more optimized DNA sequences for ``protein_seq``.

    Dispatches to the mode-specific helper, scores each variant with
    :func:`composite_score`, sorts variants by score descending, and marks the
    top variant as recommended.

    Args:
        protein_seq: Target protein (no stop codon).
        codon_table: Host CUT from :func:`build_codon_usage_table`.
        mode: ``"max_cai"`` | ``"weighted"`` | ``"harmonized"``.
        n_variants: Number of variants to generate (``max_cai`` always yields 1).
        host_gc: Host reference GC fraction (used for scoring).
        original_dna: Required for harmonized mode; provides source codons.
        source_cut: Source organism CUT for harmonized mode. If ``None`` and
            ``original_dna`` is given, a single-gene CUT is estimated from the
            input and a warning is logged.
        seed: Random seed for reproducibility in stochastic modes.
        config: Optional full config dict. Only the ``variant_scoring`` and
            ``optimization`` sections are read.
        append_stop: If True, append a stop codon (TAA) to every variant.

    Returns:
        List of dicts sorted by ``composite_score`` descending. Each dict
        carries ``dna_sequence``, ``cai_score``, ``gc_content``, ``mode``,
        ``variant_number`` (1-based input order), ``composite_score``, and
        ``is_recommended``.

    Raises:
        ValueError: On unknown mode, empty protein, or harmonized mode without DNA.
    """
    if mode not in VALID_MODES:
        raise ValueError(f"Unknown optimization mode: {mode}. Valid: {VALID_MODES}")
    if not protein_seq:
        raise ValueError("Empty protein sequence")
    protein_seq = protein_seq.upper().rstrip("*")

    cfg = config or {}
    opt_cfg = cfg.get("optimization", {})
    scoring_cfg = cfg.get("variant_scoring", {})
    weights = {
        "cai_weight": scoring_cfg.get("cai_weight", 0.5),
        "gc_weight": scoring_cfg.get("gc_weight", 0.3),
        "mrna_weight": scoring_cfg.get("mrna_weight", 0.2),
    }
    harm_noise = opt_cfg.get("harmonized_noise", 0.1)
    oversample = max(1, int(opt_cfg.get("oversample_factor", 3)))
    forbidden_motifs = cfg.get("forbidden_motifs") or []
    if not forbidden_motifs and opt_cfg.get("use_default_forbidden", True):
        forbidden_motifs = DEFAULT_FORBIDDEN_MOTIFS
    repair_enabled = bool(opt_cfg.get("repair_forbidden_motifs", True)) and bool(
        forbidden_motifs
    )

    by_aa = _codons_by_aa(codon_table)
    host_gc_value = host_gc if host_gc is not None else 0.5

    # For deterministic mode, oversampling is pointless
    if mode == "max_cai":
        candidates = [_optimize_max_cai(protein_seq, by_aa)]
    else:
        n_to_generate = n_variants * oversample
        candidates = _generate_candidates(
            mode=mode,
            protein_seq=protein_seq,
            by_aa=by_aa,
            codon_table=codon_table,
            source_cut=source_cut,
            original_dna=original_dna,
            n_variants=n_to_generate,
            base_seed=seed,
            harm_noise=harm_noise,
        )

    if append_stop:
        candidates = [s + DEFAULT_STOP_CODON for s in candidates]

    # Optional forbidden-motif repair
    unresolved_per_variant: list[list[dict]] = []
    if repair_enabled:
        repaired: list[str] = []
        for i, dna in enumerate(candidates):
            repair_rng = np.random.default_rng(
                _derive_seed(seed, i, salt="repair") if seed is not None else None
            )
            new_dna, unresolved = repair_sequence(
                dna, codon_table, forbidden_motifs, rng=repair_rng
            )
            repaired.append(new_dna)
            unresolved_per_variant.append(unresolved)
        candidates = repaired
    else:
        unresolved_per_variant = [[] for _ in candidates]

    # Per-variant 5' mRNA structure (used to score and to surface MFE)
    qc_cfg = cfg.get("qc", {})
    mfe_threshold = qc_cfg.get("mrna_mfe_threshold", -30.0)
    mrna_window = qc_cfg.get("mrna_5prime_window", 150)

    # Score each candidate
    scored: list[dict] = []
    for i, dna in enumerate(candidates):
        cai = calculate_cai(dna, codon_table)
        gc_v = gc_content(dna)
        mrna = check_mrna_structure(dna, window_5prime=mrna_window, mfe_threshold=mfe_threshold)
        mfe_value = mrna["mfe"]
        score = composite_score(
            cai=cai,
            gc_content_value=gc_v,
            host_gc=host_gc_value,
            mfe=mfe_value,
            mfe_threshold=mfe_threshold,
            weights=weights,
        )
        # Penalty for unresolved forbidden motifs (small, subtractive)
        n_unresolved = len(unresolved_per_variant[i])
        forbidden_penalty = 0.05 * n_unresolved
        final_score = score - forbidden_penalty
        scored.append(
            {
                "dna_sequence": dna,
                "cai_score": cai,
                "gc_content": gc_v,
                "mfe_5prime": mfe_value,
                "mrna_structure": mrna["structure"],
                "has_strong_structure": mrna["has_strong_structure"],
                "mode": mode,
                "variant_number": i + 1,
                "composite_score": final_score,
                "forbidden_hits": unresolved_per_variant[i],
                "is_recommended": False,
            }
        )

    scored.sort(
        key=lambda v: (
            # First: no unresolved forbidden motifs
            len(v["forbidden_hits"]) == 0,
            # Second: no strong 5' mRNA structure
            not v["has_strong_structure"],
            # Third: composite score
            v["composite_score"],
        ),
        reverse=True,
    )

    variants = scored[:max(1, n_variants)]
    # Renumber to reflect final output order
    for i, v in enumerate(variants, start=1):
        v["variant_number"] = i
    variants[0]["is_recommended"] = True
    return variants


def _derive_seed(base_seed: int | None, index: int, salt: str = "") -> int | None:
    """Derive a deterministic child seed from ``(base_seed, index, salt)``."""
    if base_seed is None:
        return None
    # Cheap mixing; stable across runs
    mix = (base_seed * 2654435761) ^ (index * 40503) ^ (hash(salt) & 0xFFFFFFFF)
    return mix & 0x7FFFFFFF


def _generate_candidates(
    mode: str,
    protein_seq: str,
    by_aa: dict[str, list[tuple[str, float]]],
    codon_table: pd.DataFrame,
    source_cut: pd.DataFrame | None,
    original_dna: str | None,
    n_variants: int,
    base_seed: int | None,
    harm_noise: float,
) -> list[str]:
    """Generate ``n_variants`` candidate DNA sequences in ``mode``."""
    candidates: list[str] = []
    if mode == "weighted":
        for i in range(n_variants):
            rng = np.random.default_rng(_derive_seed(base_seed, i, "weighted"))
            candidates.extend(_optimize_weighted(protein_seq, by_aa, 1, rng))
    elif mode == "harmonized":
        if not original_dna:
            raise ValueError("harmonized mode requires DNA input (original_dna)")
        if source_cut is None:
            logger.warning(
                "No source CUT provided; estimating from input gene itself (single-gene CUT)."
            )
            source_cut = build_codon_usage_table([original_dna])
        source_by_aa = _codons_by_aa(source_cut)
        for i in range(n_variants):
            rng = np.random.default_rng(_derive_seed(base_seed, i, "harmonized"))
            candidates.extend(
                _optimize_harmonized(
                    protein_seq, original_dna, by_aa, source_by_aa, 1, rng, harm_noise,
                )
            )
    else:
        raise ValueError(f"Unknown optimization mode: {mode}")
    return candidates


# ---------------------------------------------------------------------------
# Mode implementations
# ---------------------------------------------------------------------------

def _optimize_max_cai(
    protein_seq: str,
    by_aa: dict[str, list[tuple[str, float]]],
) -> str:
    out = []
    for aa in protein_seq:
        pairs = by_aa.get(aa)
        if not pairs:
            raise ValueError(f"Unknown amino acid '{aa}' in protein sequence")
        # Highest w_i, ties broken alphabetically — already sorted that way
        out.append(pairs[0][0])
    return "".join(out)


def _optimize_weighted(
    protein_seq: str,
    by_aa: dict[str, list[tuple[str, float]]],
    n_variants: int,
    rng: np.random.Generator,
) -> list[str]:
    # Pre-compute probability vectors per amino acid
    aa_probs: dict[str, tuple[list[str], np.ndarray]] = {}
    for aa, pairs in by_aa.items():
        if aa == "*":
            continue
        codons = [c for c, _ in pairs]
        weights = np.array([w for _, w in pairs], dtype=float)
        if weights.sum() <= 0:
            probs = np.ones_like(weights) / len(weights)
        else:
            probs = weights / weights.sum()
        aa_probs[aa] = (codons, probs)

    variants: list[str] = []
    for _ in range(max(1, n_variants)):
        out: list[str] = []
        for aa in protein_seq:
            if aa not in aa_probs:
                raise ValueError(f"Unknown amino acid '{aa}' in protein sequence")
            codons, probs = aa_probs[aa]
            out.append(rng.choice(codons, p=probs))
        variants.append("".join(out))
    return variants


def _optimize_harmonized(
    protein_seq: str,
    original_dna: str,
    host_by_aa: dict[str, list[tuple[str, float]]],
    source_by_aa: dict[str, list[tuple[str, float]]],
    n_variants: int,
    rng: np.random.Generator,
    noise_prob: float,
) -> list[str]:
    """Harmonized optimization: preserve per-position rarity rank from source.

    For each position we find the source codon's rank among its synonyms in
    the SOURCE CUT and pick the host codon at the same rank. With probability
    ``noise_prob`` we instead sample from the host weighted distribution
    (prevents over-rigidity in repetitive regions). Ties in source ranking
    are broken alphabetically to stay deterministic.
    """
    source_codons = chunk_codons(original_dna)
    # Strip trailing stop if present
    if source_codons and CODON_TABLE_STANDARD.get(source_codons[-1]) == "*":
        source_codons = source_codons[:-1]
    if len(source_codons) != len(protein_seq):
        raise ValueError(
            f"original_dna length {len(original_dna)} does not match protein "
            f"length {len(protein_seq)} (codons: {len(source_codons)})"
        )

    # Pre-compute ranks
    source_rank: dict[str, int] = {}
    for aa, pairs in source_by_aa.items():
        for rank, (c, _) in enumerate(pairs):
            source_rank[c] = rank

    host_weighted_probs: dict[str, tuple[list[str], np.ndarray]] = {}
    for aa, pairs in host_by_aa.items():
        if aa == "*":
            continue
        codons = [c for c, _ in pairs]
        weights = np.array([w for _, w in pairs], dtype=float)
        if weights.sum() <= 0:
            probs = np.ones_like(weights) / len(weights)
        else:
            probs = weights / weights.sum()
        host_weighted_probs[aa] = (codons, probs)

    variants: list[str] = []
    for _ in range(max(1, n_variants)):
        out: list[str] = []
        for i, aa in enumerate(protein_seq):
            host_pairs = host_by_aa.get(aa, [])
            if not host_pairs:
                raise ValueError(f"Unknown amino acid '{aa}' in protein sequence")
            host_codons = [c for c, _ in host_pairs]  # sorted by w_i desc

            if rng.random() < noise_prob:
                codons, probs = host_weighted_probs[aa]
                out.append(rng.choice(codons, p=probs))
                continue

            src = source_codons[i]
            rank = source_rank.get(src, 0)
            # Clip to host family size
            rank = min(rank, len(host_codons) - 1)
            out.append(host_codons[rank])
        variants.append("".join(out))
    return variants
