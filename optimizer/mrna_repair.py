"""Reduce strong 5' mRNA secondary structure via synonymous codon substitution.

Two algorithms are provided, both preserve protein translation by construction:

- :func:`targeted_repair` (Approach A) — folds the 5' window with ViennaRNA,
  parses the dot-bracket structure, and synonymously rewrites codons that
  participate in stems until the MFE crosses ``mfe_threshold`` (or attempts
  exhausted). Greedy: prefers substitutions that improve MFE the most while
  keeping CAI as high as possible.

- :func:`simulated_anneal` (Approach E) — stochastic Metropolis-Hastings
  search over the synonymous-codon neighbourhood of the input sequence,
  using a combined energy ``E = max(0, threshold - MFE) + λ·(1 - CAI)``.
  Slower but escapes local minima that block the greedy approach.

Both functions return ``(new_dna, info)`` where ``info`` carries
``initial_mfe``, ``final_mfe``, ``initial_cai``, ``final_cai``,
``attempts``, and ``succeeded`` (bool).
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .metrics import calculate_cai
from .quality_control import check_mrna_structure
from .utils import AA_TO_CODONS, CODON_TABLE_STANDARD, chunk_codons

logger = logging.getLogger(__name__)


@dataclass
class RepairResult:
    dna: str
    initial_mfe: float
    final_mfe: float
    initial_cai: float
    final_cai: float
    attempts: int
    succeeded: bool
    method: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "initial_mfe": self.initial_mfe,
            "final_mfe": self.final_mfe,
            "initial_cai": self.initial_cai,
            "final_cai": self.final_cai,
            "attempts": self.attempts,
            "succeeded": self.succeeded,
            "method": self.method,
            "delta_mfe": self.final_mfe - self.initial_mfe,
            "delta_cai": self.final_cai - self.initial_cai,
        }


def _wi_map(codon_table: pd.DataFrame) -> dict[str, float]:
    return dict(zip(codon_table["codon"], codon_table["w_i"]))


def _stem_codons_in_window(
    structure: str, window: int, n_codons_total: int
) -> list[int]:
    """Return codon indices that include any paired nt in the 5' window.

    Codons are 0-based. Only codons whose nt range falls within the window
    are considered. Codons containing the start codon (index 0) are
    excluded since ATG cannot be substituted.
    """
    paired_nt = {i for i, c in enumerate(structure[:window]) if c in "()"}
    codon_idx_set: set[int] = set()
    for nt in paired_nt:
        codon = nt // 3
        if codon == 0:
            continue
        if codon >= n_codons_total:
            continue
        codon_idx_set.add(codon)
    return sorted(codon_idx_set)


def _synonyms_excluding(codon: str) -> list[str]:
    aa = CODON_TABLE_STANDARD.get(codon)
    if aa is None or aa == "*":
        return []
    return [c for c in AA_TO_CODONS[aa] if c != codon]


def _replace_codon(dna: str, codon_idx: int, new_codon: str) -> str:
    start = codon_idx * 3
    return dna[:start] + new_codon + dna[start + 3:]


def _fold_5prime(dna: str, window: int, threshold: float) -> tuple[str, float]:
    """Fold the 5' window and return (structure, mfe)."""
    r = check_mrna_structure(dna, window_5prime=window, mfe_threshold=threshold)
    return r["structure"] or "", float(r["mfe"]) if r["mfe"] is not None else 0.0


# ---------------------------------------------------------------------------
# Approach A — targeted greedy repair
# ---------------------------------------------------------------------------

def targeted_repair(
    dna: str,
    codon_table: pd.DataFrame,
    window: int = 150,
    mfe_threshold: float = -30.0,
    max_attempts: int = 30,
    rng: np.random.Generator | None = None,
) -> tuple[str, RepairResult]:
    """Greedy 5' structure repair.

    Algorithm
    ---------
    1. Fold the first ``window`` nt; if MFE > threshold, return unchanged.
    2. Identify codons that overlap any paired position in the structure.
    3. For each such codon (random order), try every synonym in order of
       decreasing ``w_i`` (CAI-friendly first). Accept the first synonym
       that strictly improves MFE.
    4. Repeat the fold-and-substitute loop up to ``max_attempts`` times or
       until MFE > threshold.

    Args:
        dna: Input DNA sequence (length divisible by 3).
        codon_table: Host CUT.
        window: 5' window length to fold (nt).
        mfe_threshold: MFE above which structure is considered acceptable.
        max_attempts: Maximum substitution attempts.
        rng: Optional numpy Generator for deterministic order.

    Returns:
        ``(new_dna, RepairResult)``.
    """
    if len(dna) % 3 != 0:
        raise ValueError(f"DNA length {len(dna)} not divisible by 3")
    rng = rng or np.random.default_rng()
    wi = _wi_map(codon_table)

    structure, mfe = _fold_5prime(dna, window, mfe_threshold)
    initial_mfe = mfe
    initial_cai = calculate_cai(dna, codon_table)

    n_codons = len(dna) // 3
    current = dna
    attempts = 0

    if mfe > mfe_threshold:
        return current, RepairResult(
            dna=current, initial_mfe=initial_mfe, final_mfe=mfe,
            initial_cai=initial_cai, final_cai=initial_cai,
            attempts=0, succeeded=True, method="targeted",
        )

    for outer in range(max_attempts):
        candidates = _stem_codons_in_window(structure, window, n_codons)
        if not candidates:
            break
        rng.shuffle(candidates)

        improved = False
        for idx in candidates:
            attempts += 1
            current_codon = current[idx * 3:idx * 3 + 3]
            synonyms = _synonyms_excluding(current_codon)
            if not synonyms:
                continue
            # Try CAI-preferred synonyms first
            synonyms.sort(key=lambda c: -wi.get(c, 0.0))
            best_new_dna = None
            best_new_mfe = mfe
            best_new_struct = structure
            for syn in synonyms:
                tentative = _replace_codon(current, idx, syn)
                new_struct, new_mfe = _fold_5prime(tentative, window, mfe_threshold)
                if new_mfe > best_new_mfe + 1e-6:
                    best_new_mfe = new_mfe
                    best_new_dna = tentative
                    best_new_struct = new_struct
            if best_new_dna is not None:
                current = best_new_dna
                mfe = best_new_mfe
                structure = best_new_struct
                improved = True
                if mfe > mfe_threshold:
                    break
        if mfe > mfe_threshold or not improved:
            break

    final_cai = calculate_cai(current, codon_table)
    return current, RepairResult(
        dna=current,
        initial_mfe=initial_mfe,
        final_mfe=mfe,
        initial_cai=initial_cai,
        final_cai=final_cai,
        attempts=attempts,
        succeeded=mfe > mfe_threshold,
        method="targeted",
    )


# ---------------------------------------------------------------------------
# Approach E — simulated annealing
# ---------------------------------------------------------------------------

def simulated_anneal(
    dna: str,
    codon_table: pd.DataFrame,
    window: int = 150,
    mfe_threshold: float = -30.0,
    iterations: int = 200,
    t_start: float = 5.0,
    t_end: float = 0.05,
    cai_lambda: float = 3.0,
    mutation_window_codons: int | None = None,
    rng: np.random.Generator | None = None,
) -> tuple[str, RepairResult]:
    """Metropolis–Hastings search for higher 5' MFE while keeping CAI high.

    Energy
    ------
    ``E(seq) = max(0, threshold - MFE) + λ · (1 - CAI)``

    The first term is zero when the sequence already crosses the MFE
    threshold and grows linearly otherwise. The second term penalises CAI
    loss with weight ``cai_lambda`` (default 3.0). At each iteration we
    propose a synonymous substitution at a random codon position
    (restricted to the first ``mutation_window_codons`` codons by default,
    which is set to ``window // 3`` if not given), compute the new energy,
    and accept it with probability ``exp(-ΔE / T)``. Temperature decays
    geometrically from ``t_start`` to ``t_end``.

    Args:
        dna: Input DNA (length divisible by 3).
        codon_table: Host CUT.
        window: 5' window for MFE evaluation (nt).
        mfe_threshold: Target MFE.
        iterations: Number of MH steps.
        t_start: Starting temperature.
        t_end: Final temperature.
        cai_lambda: CAI penalty weight in the energy.
        mutation_window_codons: Restrict mutations to the first N codons.
            Defaults to ``window // 3`` (the same region used for MFE).
        rng: Optional numpy Generator.

    Returns:
        ``(best_dna, RepairResult)`` where best is the lowest-energy
        sequence encountered, not necessarily the final accepted state.
    """
    if len(dna) % 3 != 0:
        raise ValueError(f"DNA length {len(dna)} not divisible by 3")
    rng = rng or np.random.default_rng()
    n_codons = len(dna) // 3
    if mutation_window_codons is None:
        mutation_window_codons = max(1, window // 3)
    mutation_window_codons = min(mutation_window_codons, n_codons)
    # Skip codon 0 (Met) — never mutable
    mutable_positions = list(range(1, mutation_window_codons))

    structure, mfe = _fold_5prime(dna, window, mfe_threshold)
    initial_mfe = mfe
    cai = calculate_cai(dna, codon_table)
    initial_cai = cai

    # No mutable codon (e.g. CDS is one codon long, or window // 3 == 1):
    # nothing to optimize. Return unchanged with attempts=0.
    if not mutable_positions:
        return dna, RepairResult(
            dna=dna, initial_mfe=initial_mfe, final_mfe=mfe,
            initial_cai=initial_cai, final_cai=cai,
            attempts=0, succeeded=mfe > mfe_threshold, method="annealing",
        )

    def energy(mfe_val: float, cai_val: float) -> float:
        struct_pen = max(0.0, mfe_threshold - mfe_val)
        cai_pen = cai_lambda * max(0.0, 1.0 - cai_val)
        return struct_pen + cai_pen

    current = dna
    current_e = energy(mfe, cai)
    best = current
    best_e = current_e
    best_mfe = mfe
    best_cai = cai
    attempts = 0

    if mfe > mfe_threshold:
        # Already good — return as-is
        return current, RepairResult(
            dna=current, initial_mfe=initial_mfe, final_mfe=mfe,
            initial_cai=initial_cai, final_cai=cai,
            attempts=0, succeeded=True, method="annealing",
        )

    for step in range(iterations):
        attempts += 1
        t = t_start * (t_end / t_start) ** (step / max(1, iterations - 1))

        # Propose: pick a mutable codon, pick a random synonym
        idx = int(rng.choice(mutable_positions))
        cur_codon = current[idx * 3:idx * 3 + 3]
        synonyms = _synonyms_excluding(cur_codon)
        if not synonyms:
            continue
        new_codon = synonyms[int(rng.integers(0, len(synonyms)))]
        proposal = _replace_codon(current, idx, new_codon)
        prop_struct, prop_mfe = _fold_5prime(proposal, window, mfe_threshold)
        prop_cai = calculate_cai(proposal, codon_table)
        prop_e = energy(prop_mfe, prop_cai)

        delta = prop_e - current_e
        if delta <= 0 or rng.random() < math.exp(-delta / t):
            current = proposal
            current_e = prop_e
            cai = prop_cai
            mfe = prop_mfe
            if current_e < best_e:
                best = current
                best_e = current_e
                best_mfe = mfe
                best_cai = cai
                if best_mfe > mfe_threshold:
                    # Allow the schedule to continue — may still find smaller
                    # CAI penalty — but cheap early-out if CAI also healthy.
                    if best_cai >= initial_cai * 0.95:
                        break

    return best, RepairResult(
        dna=best,
        initial_mfe=initial_mfe,
        final_mfe=best_mfe,
        initial_cai=initial_cai,
        final_cai=best_cai,
        attempts=attempts,
        succeeded=best_mfe > mfe_threshold,
        method="annealing",
    )


# ---------------------------------------------------------------------------
# Combined dispatcher
# ---------------------------------------------------------------------------

REPAIR_MODES = ("off", "targeted", "annealing", "both")


def repair_5prime_structure(
    dna: str,
    codon_table: pd.DataFrame,
    mode: str = "targeted",
    window: int = 150,
    mfe_threshold: float = -30.0,
    targeted_max_attempts: int = 30,
    annealing_iterations: int = 200,
    annealing_t_start: float = 5.0,
    annealing_t_end: float = 0.05,
    cai_lambda: float = 3.0,
    rng: np.random.Generator | None = None,
) -> tuple[str, RepairResult]:
    """Dispatch by ``mode``: 'off' / 'targeted' / 'annealing' / 'both'.

    ``'both'`` runs targeted first; if it doesn't cross the threshold,
    follows up with simulated annealing starting from the targeted output.
    Returns the better of the two by final energy (lower MFE penalty +
    higher CAI = better).
    """
    if mode == "off":
        cai = calculate_cai(dna, codon_table)
        _, mfe = _fold_5prime(dna, window, mfe_threshold)
        return dna, RepairResult(
            dna=dna, initial_mfe=mfe, final_mfe=mfe,
            initial_cai=cai, final_cai=cai,
            attempts=0, succeeded=mfe > mfe_threshold, method="off",
        )
    if mode == "targeted":
        return targeted_repair(
            dna, codon_table, window=window, mfe_threshold=mfe_threshold,
            max_attempts=targeted_max_attempts, rng=rng,
        )
    if mode == "annealing":
        return simulated_anneal(
            dna, codon_table, window=window, mfe_threshold=mfe_threshold,
            iterations=annealing_iterations,
            t_start=annealing_t_start, t_end=annealing_t_end,
            cai_lambda=cai_lambda, rng=rng,
        )
    if mode == "both":
        new_dna, info_t = targeted_repair(
            dna, codon_table, window=window, mfe_threshold=mfe_threshold,
            max_attempts=targeted_max_attempts, rng=rng,
        )
        if info_t.succeeded:
            return new_dna, info_t
        new_dna_2, info_a = simulated_anneal(
            new_dna, codon_table, window=window, mfe_threshold=mfe_threshold,
            iterations=annealing_iterations,
            t_start=annealing_t_start, t_end=annealing_t_end,
            cai_lambda=cai_lambda, rng=rng,
        )
        # Combine attempt counts and report the SA-improved sequence
        info_a.initial_mfe = info_t.initial_mfe
        info_a.initial_cai = info_t.initial_cai
        info_a.attempts += info_t.attempts
        info_a.method = "both"
        return new_dna_2, info_a
    raise ValueError(f"Unknown repair mode: {mode}. Valid: {REPAIR_MODES}")
