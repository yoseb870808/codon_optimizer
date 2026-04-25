"""Biophysical quality-control checks on optimized sequences."""

from __future__ import annotations

import logging
from typing import Any

from .utils import gc_content

logger = logging.getLogger(__name__)

try:
    import RNA as _vrna  # type: ignore
except ImportError as exc:  # pragma: no cover — explicit hard requirement
    raise ImportError(
        "ViennaRNA is required for codon_optimizer. Install with one of:\n"
        "    pip install ViennaRNA           (Windows / Linux / macOS wheels available)\n"
        "    conda install -c bioconda viennarna\n"
        "Then verify with:  python -c \"import RNA; print(RNA.fold('GCGCGC'))\""
    ) from exc


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


def check_mrna_structure(
    sequence: str,
    window_5prime: int = 150,
    mfe_threshold: float = -30.0,
) -> dict[str, Any]:
    """Predict secondary-structure strength of the 5' end via ViennaRNA.

    Folds the first ``window_5prime`` nucleotides with ``RNA.fold`` and
    reports MFE (kcal/mol), the dot-bracket structure string, and whether
    the MFE crosses ``mfe_threshold`` (i.e. structure is strong enough to
    impede ribosome loading).

    Args:
        sequence: Full DNA sequence; only the first ``window_5prime`` nt
            are folded. T is transparently converted to U for folding.
        window_5prime: Length of the 5' window to fold.
        mfe_threshold: MFE (kcal/mol) at or below which
            ``has_strong_structure`` is True. Default -30.0.

    Returns:
        Dict with ``method`` (always ``"viennarna"``), ``mfe``,
        ``structure``, ``has_strong_structure``, and ``warning``.

    Raises:
        RuntimeError: If ViennaRNA fails to fold the window (rare; typically
            only on malformed input).
    """
    window = sequence[:window_5prime]
    if not window:
        return {
            "method": "viennarna",
            "mfe": None,
            "structure": None,
            "has_strong_structure": False,
            "warning": "Empty sequence",
        }

    try:
        structure, mfe = _vrna.fold(window.replace("T", "U"))
    except Exception as exc:
        raise RuntimeError(
            f"ViennaRNA fold failed on window of length {len(window)}: {exc}"
        ) from exc

    return {
        "method": "viennarna",
        "mfe": float(mfe),
        "structure": structure,
        "has_strong_structure": float(mfe) <= mfe_threshold,
        "warning": None,
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
