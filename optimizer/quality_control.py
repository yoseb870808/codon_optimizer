"""Biophysical quality-control checks on optimized sequences."""

from __future__ import annotations

import logging
from typing import Any

from .utils import gc_content, reverse_complement

logger = logging.getLogger(__name__)

try:  # pragma: no cover — import side effect
    import RNA as _vrna  # type: ignore
    VIENNARNA_AVAILABLE = True
except ImportError:  # pragma: no cover
    _vrna = None
    VIENNARNA_AVAILABLE = False
    logger.info(
        "ViennaRNA not found; mRNA structure QC will use simple hairpin-scan fallback."
    )


def check_gc_content(
    sequence: str,
    target_gc: float | None = None,
    window_size: int = 50,
    acceptable_range: tuple[float, float] = (0.25, 0.65),
    window_deviation: float = 0.15,
) -> dict[str, Any]:
    """Compute overall GC and flag sliding-window hotspots.

    Args:
        sequence: Nucleotide sequence.
        target_gc: Reference GC fraction (host). If given, ``gc_deviation`` is
            |overall_gc - target_gc|.
        window_size: Sliding window size in nt.
        acceptable_range: Acceptable overall GC fraction, inclusive.
        window_deviation: Flag any window whose GC differs from the overall
            sequence GC by more than this.

    Returns:
        Dict with keys ``overall_gc``, ``target_gc``, ``gc_deviation``,
        ``within_range``, ``window_min``, ``window_max``, ``flagged_regions``.
    """
    overall = gc_content(sequence)
    window_min = window_max = overall
    flagged: list[dict[str, Any]] = []

    if len(sequence) >= window_size:
        values: list[float] = []
        for i in range(0, len(sequence) - window_size + 1):
            w = sequence[i:i + window_size]
            gc_w = gc_content(w)
            values.append(gc_w)
            if abs(gc_w - overall) > window_deviation:
                direction = "high" if gc_w > overall else "low"
                # Merge adjacent flagged windows
                if flagged and flagged[-1]["end"] >= i and flagged[-1]["direction"] == direction:
                    flagged[-1]["end"] = i + window_size
                else:
                    flagged.append({
                        "start": i,
                        "end": i + window_size,
                        "gc": gc_w,
                        "direction": direction,
                    })
        if values:
            window_min = min(values)
            window_max = max(values)

    return {
        "overall_gc": overall,
        "target_gc": target_gc,
        "gc_deviation": abs(overall - target_gc) if target_gc is not None else None,
        "within_range": acceptable_range[0] <= overall <= acceptable_range[1],
        "window_min": window_min,
        "window_max": window_max,
        "flagged_regions": flagged,
    }


def check_homopolymers(sequence: str, max_run: int = 6) -> list[dict[str, Any]]:
    """Detect homopolymer runs longer than ``max_run``.

    Returns one entry per run: ``{"position": 0-based start, "base": X, "length": N}``.
    """
    if not sequence:
        return []
    seq = sequence.upper()
    hits: list[dict[str, Any]] = []
    run_base = seq[0]
    run_start = 0
    for i in range(1, len(seq)):
        if seq[i] == run_base:
            continue
        run_len = i - run_start
        if run_len > max_run:
            hits.append({"position": run_start, "base": run_base, "length": run_len})
        run_base = seq[i]
        run_start = i
    # Trailing run
    run_len = len(seq) - run_start
    if run_len > max_run:
        hits.append({"position": run_start, "base": run_base, "length": run_len})
    return hits


def check_repeats(sequence: str, min_length: int = 12) -> list[dict[str, Any]]:
    """Detect exact direct repeats of length ≥ ``min_length``.

    Uses a suffix-index approach: build the set of all k-mers of size
    ``min_length`` and record the positions where each k-mer occurs. Any
    k-mer occurring ≥ 2 times is reported as a repeat pair (first two
    positions). This is O(n * k) which is fast enough for sequences of a
    few thousand nt typical in codon optimization.

    Returns one entry per repeat: ``{"pos1", "pos2", "sequence", "length"}``.
    """
    if min_length <= 0 or len(sequence) < 2 * min_length:
        return []
    seq = sequence.upper()
    positions: dict[str, list[int]] = {}
    for i in range(len(seq) - min_length + 1):
        kmer = seq[i:i + min_length]
        positions.setdefault(kmer, []).append(i)

    hits: list[dict[str, Any]] = []
    seen_pairs: set[tuple[int, int]] = set()
    for kmer, pos_list in positions.items():
        if len(pos_list) < 2:
            continue
        # Report the first two occurrences as a single hit per unique kmer
        pos1, pos2 = pos_list[0], pos_list[1]
        # Extend the hit as far right as it stays equal (maximal run)
        extend = min_length
        while (
            pos2 + extend < len(seq)
            and pos1 + extend < pos2
            and seq[pos1 + extend] == seq[pos2 + extend]
        ):
            extend += 1
        key = (pos1, pos2)
        if key in seen_pairs:
            continue
        seen_pairs.add(key)
        hits.append({
            "pos1": pos1,
            "pos2": pos2,
            "sequence": seq[pos1:pos1 + extend],
            "length": extend,
        })
    hits.sort(key=lambda h: (-h["length"], h["pos1"]))
    return hits


def _palindrome_scan(
    seq: str,
    min_stem: int = 5,
    max_stem: int = 20,
    max_loop: int = 20,
) -> list[dict[str, Any]]:
    """Find inverted-repeat pairs (stem-loop candidates).

    Search for pairs ``seq[i:i+L]`` and ``seq[j:j+L]`` where the second block
    is the reverse complement of the first and the two blocks are separated
    by a loop of 0..``max_loop`` nt. For each starting position we report the
    match with the shortest loop and greatest stem length.
    """
    seq = seq.upper()
    n = len(seq)
    hits: list[dict[str, Any]] = []
    reported: set[int] = set()

    for i in range(n - 2 * min_stem):
        if i in reported:
            continue
        for loop in range(0, max_loop + 1):
            j = i + min_stem + loop
            if j + min_stem > n:
                break
            # Does a min_stem match start here?
            if seq[i:i + min_stem] != reverse_complement(seq[j:j + min_stem]):
                continue
            # Extend as long as pair still reverse-complements
            best_len = min_stem
            for L in range(min_stem + 1, max_stem + 1):
                if i + L > j or j + L > n:
                    break
                if seq[i:i + L] == reverse_complement(seq[j:j + L]):
                    best_len = L
                else:
                    break
            hits.append({
                "position": i,
                "length": best_len,
                "sequence": seq[i:i + best_len],
                "loop": loop,
                "partner_position": j,
            })
            reported.add(i)
            break
    return hits


def check_mrna_structure(
    sequence: str,
    window_5prime: int = 150,
    mfe_threshold: float = -30.0,
) -> dict[str, Any]:
    """Predict secondary-structure strength of the 5' end.

    Primary path (ViennaRNA installed) runs RNA.fold over the first
    ``window_5prime`` nucleotides. Fallback does a palindrome scan for
    candidate hairpin stems. Flags ``has_strong_structure = True`` when
    MFE ≤ threshold (primary) or when any palindrome ≥ 8 bp is found
    (fallback).
    """
    window = sequence[:window_5prime]
    warning: str | None = None
    if not window:
        return {
            "method": "viennarna" if VIENNARNA_AVAILABLE else "fallback",
            "mfe": None,
            "structure": None,
            "has_strong_structure": False,
            "hairpins": [],
            "warning": "Empty sequence",
        }

    if VIENNARNA_AVAILABLE:
        try:
            structure, mfe = _vrna.fold(window.replace("T", "U"))
        except Exception as exc:  # pragma: no cover — robust to rare VRNA errors
            logger.warning("ViennaRNA fold failed (%s); falling back to palindrome scan", exc)
            hairpins = _palindrome_scan(window)
            return {
                "method": "fallback",
                "mfe": None,
                "structure": None,
                "has_strong_structure": any(h["length"] >= 8 for h in hairpins),
                "hairpins": hairpins,
                "warning": f"ViennaRNA error: {exc}",
            }
        return {
            "method": "viennarna",
            "mfe": mfe,
            "structure": structure,
            "has_strong_structure": mfe <= mfe_threshold,
            "hairpins": [],
            "warning": warning,
        }

    hairpins = _palindrome_scan(window)
    return {
        "method": "fallback",
        "mfe": None,
        "structure": None,
        "has_strong_structure": any(h["length"] >= 8 for h in hairpins),
        "hairpins": hairpins,
        "warning": "ViennaRNA not installed — using palindrome-based heuristic",
    }


def run_all_qc(
    sequence: str,
    target_gc: float | None = None,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run every QC check with thresholds from ``config["qc"]`` (or defaults).

    Returns a single dict with nested sub-reports keyed ``gc``, ``homopolymers``,
    ``repeats``, ``mrna_structure``, plus a top-level boolean ``pass`` that is
    True only when no critical flags fire (no homopolymer, no repeat, GC within
    range, no strong mRNA structure).
    """
    qc_cfg = (config or {}).get("qc", {})
    gc_range = tuple(qc_cfg.get("gc_range", (0.25, 0.65)))
    window_size = qc_cfg.get("gc_window_size", 50)
    window_dev = qc_cfg.get("gc_window_deviation", 0.15)
    max_homo = qc_cfg.get("max_homopolymer", 6)
    min_rep = qc_cfg.get("min_repeat_length", 12)
    mfe_thresh = qc_cfg.get("mrna_mfe_threshold", -30.0)
    five_prime = qc_cfg.get("mrna_5prime_window", 150)

    gc_result = check_gc_content(
        sequence,
        target_gc=target_gc,
        window_size=window_size,
        acceptable_range=gc_range,
        window_deviation=window_dev,
    )
    homo_result = check_homopolymers(sequence, max_run=max_homo)
    rep_result = check_repeats(sequence, min_length=min_rep)
    mrna_result = check_mrna_structure(
        sequence, window_5prime=five_prime, mfe_threshold=mfe_thresh
    )

    passed = (
        gc_result["within_range"]
        and not homo_result
        and not rep_result
        and not mrna_result["has_strong_structure"]
    )
    return {
        "gc": gc_result,
        "homopolymers": homo_result,
        "repeats": rep_result,
        "mrna_structure": mrna_result,
        "pass": passed,
    }
