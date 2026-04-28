"""Output file generation (FASTA, TSV, markdown summary)."""

from __future__ import annotations

import logging
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

logger = logging.getLogger(__name__)


def _ensure_dir(output_dir: str | os.PathLike) -> Path:
    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    return path


def _safe_stem(header: str) -> str:
    """Turn a FASTA header into a filesystem-safe filename stem."""
    stem = header.split()[0] if header else "sequence"
    keep = []
    for ch in stem:
        if ch.isalnum() or ch in ("-", "_", "."):
            keep.append(ch)
        else:
            keep.append("_")
    return "".join(keep) or "sequence"


def _fmt_fasta_header(original_header: str, suffix: str, variant: dict) -> str:
    return (
        f">{original_header}_{suffix} "
        f"CAI={variant['cai_score']:.3f} "
        f"GC={variant['gc_content'] * 100:.1f}% "
        f"score={variant['composite_score']:.3f} "
        f"mode={variant['mode']}"
    )


def _wrap(seq: str, width: int = 70) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_optimized_fasta(
    header: str,
    variants: list[dict],
    output_dir: str | os.PathLike,
) -> str:
    """Write a per-sequence FASTA file containing every variant.

    Variants are assumed already sorted with the recommended one first. The
    recommended variant gets the ``_optimized_recommended`` suffix, others
    ``_variant_N``.

    Returns:
        Path (string) to the written file.
    """
    out = _ensure_dir(output_dir)
    stem = _safe_stem(header)
    path = out / f"{stem}_optimized.fasta"

    lines: list[str] = []
    for i, v in enumerate(variants):
        suffix = "optimized_recommended" if v["is_recommended"] else f"variant_{i + 1}"
        lines.append(_fmt_fasta_header(header, suffix, v))
        lines.append(_wrap(v["dna_sequence"]))

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def write_combined_fasta(
    all_results: list[dict],
    output_dir: str | os.PathLike,
) -> str:
    """Write a single FASTA containing only the recommended variant per input."""
    out = _ensure_dir(output_dir)
    path = out / "combined_recommended.fasta"

    lines: list[str] = []
    for res in all_results:
        header = res["header"]
        rec = next((v for v in res["variants"] if v["is_recommended"]), None)
        if rec is None:
            continue
        lines.append(_fmt_fasta_header(header, "optimized_recommended", rec))
        lines.append(_wrap(rec["dna_sequence"]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def write_report_tsv(
    all_results: list[dict],
    output_dir: str | os.PathLike,
) -> str:
    """Write a one-row-per-sequence TSV summarising the recommended variant."""
    out = _ensure_dir(output_dir)
    path = out / "optimization_report.tsv"

    rows: list[dict] = []
    for res in all_results:
        rec = next((v for v in res["variants"] if v["is_recommended"]), None)
        if rec is None:
            continue
        qc = res.get("qc", {})
        mrna = qc.get("mrna_structure", {})
        forbidden_hits = rec.get("forbidden_hits", [])
        rows.append(
            {
                "sequence_name": res["header"].split()[0],
                "input_type": res["input_type"],
                "length_aa": res["length_aa"],
                "n_variants": len(res["variants"]),
                "original_cai": res.get("original_cai"),
                "recommended_cai": rec["cai_score"],
                "cai_improvement": (
                    rec["cai_score"] - res["original_cai"]
                    if res.get("original_cai") is not None
                    else None
                ),
                "original_gc": res.get("original_gc"),
                "recommended_gc": rec["gc_content"],
                "host_gc": res.get("host_gc"),
                "homopolymer_count": len(qc.get("homopolymers", [])),
                "repeat_count": len(qc.get("repeats", [])),
                "mfe_5prime": rec.get("mfe_5prime", mrna.get("mfe")),
                "has_strong_structure": rec.get("has_strong_structure", mrna.get("has_strong_structure")),
                "forbidden_motif_count": len(forbidden_hits),
                "forbidden_motifs": ", ".join(sorted({h["motif"] for h in forbidden_hits})) or "",
                "composite_score": rec["composite_score"],
                "optimization_mode": rec["mode"],
                "qc_pass": qc.get("pass"),
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)
    return str(path)


def write_cut_tsv(
    codon_table: pd.DataFrame,
    reference_info: dict[str, Any],
    output_dir: str | os.PathLike,
    host_gc: float | None = None,
) -> str:
    """Write the custom Codon Usage Table as TSV with header comment lines."""
    out = _ensure_dir(output_dir)
    path = out / "codon_usage_table.tsv"

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
    header_lines = [
        f"# Codon Usage Table generated {timestamp}Z",
        f"# Reference mode: {reference_info.get('mode', 'unknown')}",
        f"# Reference genes: {reference_info.get('n_genes', 'unknown')}",
        f"# Total CDS in genome: {reference_info.get('n_total_cds', 'unknown')}",
    ]
    if host_gc is not None:
        header_lines.append(f"# Host GC: {host_gc * 100:.1f}%")
    if "match_rate" in reference_info:
        header_lines.append(f"# RNA-seq match rate: {reference_info['match_rate'] * 100:.1f}%")

    with open(path, "w", encoding="utf-8", newline="") as f:
        for line in header_lines:
            f.write(line + "\n")
        codon_table.to_csv(f, sep="\t", index=False)
    return str(path)


def write_comparison_tsv(
    header: str,
    comparison_df: pd.DataFrame,
    output_dir: str | os.PathLike,
) -> str:
    """Write a per-position comparison table (only meaningful for DNA input)."""
    out = _ensure_dir(output_dir)
    stem = _safe_stem(header)
    path = out / f"{stem}_codon_comparison.tsv"
    comparison_df.to_csv(path, sep="\t", index=False)
    return str(path)


def write_summary_md(
    all_results: list[dict],
    reference_info: dict,
    config: dict,
    output_dir: str | os.PathLike,
    inputs: dict[str, Any] | None = None,
) -> str:
    """Write a human-readable markdown summary of the run."""
    out = _ensure_dir(output_dir)
    path = out / "summary.md"

    inputs = inputs or {}
    lines: list[str] = []
    lines.append("# Codon Optimizer — Run Summary")
    lines.append("")
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
    lines.append(f"_Generated {timestamp}Z_")
    lines.append("")

    # Run parameters
    lines.append("## Run parameters")
    for key, value in inputs.items():
        lines.append(f"- **{key}**: `{value}`")
    lines.append(f"- **mode**: `{config.get('optimization', {}).get('mode', 'n/a')}`")
    lines.append(f"- **n_variants**: `{config.get('optimization', {}).get('n_variants', 'n/a')}`")
    lines.append("")

    # Reference translatome
    lines.append("## Reference translatome")
    lines.append(f"- Mode: **{reference_info.get('mode')}**")
    lines.append(f"- Genes used: **{reference_info.get('n_genes')}**")
    lines.append(f"- Total CDS in genome: **{reference_info.get('n_total_cds')}**")
    if "match_rate" in reference_info:
        lines.append(f"- RNA-seq match rate: **{reference_info['match_rate'] * 100:.1f}%**")
    if "host_gc" in reference_info:
        lines.append(f"- Host GC: **{reference_info['host_gc'] * 100:.1f}%**")
    lines.append("")

    # Per-sequence results
    lines.append("## Results")
    lines.append("")
    lines.append(
        "| Sequence | Input | Length (AA) | Original CAI | Recommended CAI | "
        "GC (%) | Score | QC pass |"
    )
    lines.append("|---|---|---|---|---|---|---|---|")
    for res in all_results:
        rec = next((v for v in res["variants"] if v["is_recommended"]), None)
        if rec is None:
            continue
        qc = res.get("qc", {})
        orig_cai_str = (
            f"{res['original_cai']:.3f}" if res.get("original_cai") is not None else "—"
        )
        lines.append(
            f"| {res['header'].split()[0]} | {res['input_type']} | "
            f"{res['length_aa']} | {orig_cai_str} | "
            f"{rec['cai_score']:.3f} | {rec['gc_content'] * 100:.1f} | "
            f"{rec['composite_score']:.3f} | "
            f"{'✓' if qc.get('pass') else '⚠'} |"
        )
    lines.append("")

    # Warnings aggregated
    warnings: list[str] = []
    for res in all_results:
        qc = res.get("qc", {})
        if qc.get("homopolymers"):
            warnings.append(
                f"{res['header'].split()[0]}: {len(qc['homopolymers'])} homopolymer run(s)"
            )
        if qc.get("repeats"):
            warnings.append(
                f"{res['header'].split()[0]}: {len(qc['repeats'])} repeat region(s)"
            )
        if qc.get("mrna_structure", {}).get("has_strong_structure"):
            warnings.append(
                f"{res['header'].split()[0]}: strong 5' mRNA structure predicted"
            )
    if warnings:
        lines.append("## QC flags")
        for w in warnings:
            lines.append(f"- {w}")
        lines.append("")
    else:
        lines.append("## QC flags")
        lines.append("_No flags raised._")
        lines.append("")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)
