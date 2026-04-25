"""End-to-end pipeline orchestration: config → CUT → optimize → QC → output."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any

import pandas as pd

from .input_handler import (
    load_config,
    parse_expression_data,
    parse_genbank,
    parse_target_sequences,
)
from .metrics import calculate_cai, compare_codons
from .quality_control import run_all_qc
from .reference_builder import (
    build_codon_usage_table,
    build_reference_set,
    compute_host_gc,
)
from .report_generator import (
    write_combined_fasta,
    write_comparison_tsv,
    write_cut_tsv,
    write_optimized_fasta,
    write_report_tsv,
    write_summary_md,
)
from .sequence_optimizer import optimize
from .utils import gc_content

logger = logging.getLogger(__name__)


def run_pipeline(
    genome_path: str | os.PathLike,
    target_path: str | os.PathLike,
    rnaseq_path: str | os.PathLike | None = None,
    mode: str | None = None,
    n_variants: int | None = None,
    output_dir: str | os.PathLike | None = None,
    config_path: str | os.PathLike | None = None,
    seed: int | None = None,
    source_genome_path: str | os.PathLike | None = None,
) -> dict[str, Any]:
    """Run the full codon-optimization pipeline and write outputs.

    Args:
        genome_path: Host GenBank file.
        target_path: FASTA of target protein or DNA sequences.
        rnaseq_path: Optional Excel file with RPKM expression data.
        mode: Optimization mode (``max_cai`` / ``weighted`` / ``harmonized``).
        n_variants: Number of variants per input sequence.
        output_dir: Directory to write outputs into. Created if missing.
        config_path: Optional user config YAML that merges on top of defaults.
        seed: Random seed for reproducibility.

    Returns:
        Dict with keys: ``reference_info``, ``codon_usage_table``, ``results``,
        ``output_paths``.
    """
    config = load_config(config_path)
    if mode is not None:
        config.setdefault("optimization", {})["mode"] = mode
    if n_variants is not None:
        config.setdefault("optimization", {})["n_variants"] = n_variants
    if output_dir is None:
        output_dir = config.get("output", {}).get("dir", "output")

    mode = config["optimization"]["mode"]
    n_variants = int(config["optimization"]["n_variants"])
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    logger.info("Parsing genome: %s", genome_path)
    all_cds = parse_genbank(genome_path)

    expression = None
    if rnaseq_path is not None:
        logger.info("Parsing RNA-seq: %s", rnaseq_path)
        expression = parse_expression_data(rnaseq_path)

    seqs, ref_info = build_reference_set(all_cds, expression, config)
    cut = build_codon_usage_table(seqs)
    host_gc = compute_host_gc(seqs)
    ref_info["host_gc"] = host_gc
    logger.info(
        "Reference: mode=%s, %d genes, host GC=%.1f%%",
        ref_info["mode"], ref_info["n_genes"], host_gc * 100,
    )

    # Build source CUT for harmonized mode if a source genome is supplied
    source_cut = None
    if source_genome_path is not None:
        logger.info("Parsing source genome: %s", source_genome_path)
        source_cds = parse_genbank(source_genome_path)
        source_sequences = [
            c["nucleotide_seq"] for c in source_cds.values()
            if c["length_nt"] >= 300  # keep anything reasonable
        ]
        if not source_sequences:
            raise ValueError(
                f"Source genome {source_genome_path} yielded no usable CDS"
            )
        source_cut = build_codon_usage_table(source_sequences)
        logger.info(
            "Source CUT built from %d CDS (source GC=%.1f%%)",
            len(source_sequences),
            compute_host_gc(source_sequences) * 100,
        )

    logger.info("Parsing targets: %s", target_path)
    targets = parse_target_sequences(target_path)

    results: list[dict[str, Any]] = []
    for idx, t in enumerate(targets):
        logger.info(
            "Optimizing [%d/%d]: %s (%s, %d AA)",
            idx + 1, len(targets), t["header"].split()[0], t["input_type"], len(t["protein_seq"]),
        )
        variants = optimize(
            protein_seq=t["protein_seq"],
            codon_table=cut,
            mode=mode,
            n_variants=n_variants,
            host_gc=host_gc,
            original_dna=t["original_dna"],
            source_cut=source_cut,
            seed=seed,
            config=config,
        )
        recommended = next(v for v in variants if v["is_recommended"])
        qc_report = run_all_qc(
            recommended["dna_sequence"],
            target_gc=host_gc,
            config=config,
        )
        original_cai = None
        original_gc = None
        if t["input_type"] == "dna" and t["original_dna"]:
            original_cai = calculate_cai(t["original_dna"], cut)
            original_gc = gc_content(t["original_dna"])
        results.append(
            {
                "header": t["header"],
                "input_type": t["input_type"],
                "protein_seq": t["protein_seq"],
                "original_dna": t["original_dna"],
                "length_aa": len(t["protein_seq"]),
                "variants": variants,
                "qc": qc_report,
                "original_cai": original_cai,
                "original_gc": original_gc,
                "host_gc": host_gc,
            }
        )

    # Write outputs
    output_paths: dict[str, str] = {}
    per_seq_paths: list[str] = []
    for res in results:
        per_seq_paths.append(
            write_optimized_fasta(res["header"], res["variants"], out_path)
        )
        if res["input_type"] == "dna" and res["original_dna"]:
            rec = next(v for v in res["variants"] if v["is_recommended"])
            # Trim trailing stop if present so lengths match
            opt_body = rec["dna_sequence"]
            orig = res["original_dna"]
            # Align by trimming both to common multiple of 3 and equal length
            min_len = min(len(opt_body), len(orig))
            min_len -= min_len % 3
            cmp_df = compare_codons(orig[:min_len], opt_body[:min_len], cut)
            write_comparison_tsv(res["header"], cmp_df, out_path)

    output_paths["per_sequence_fastas"] = per_seq_paths
    if config.get("output", {}).get("combined_fasta", True):
        output_paths["combined_fasta"] = write_combined_fasta(results, out_path)
    output_paths["report_tsv"] = write_report_tsv(results, out_path)
    output_paths["cut_tsv"] = write_cut_tsv(cut, ref_info, out_path, host_gc=host_gc)
    output_paths["summary_md"] = write_summary_md(
        results, ref_info, config, out_path,
        inputs={
            "genome": str(genome_path),
            "target": str(target_path),
            "rnaseq": str(rnaseq_path) if rnaseq_path else None,
            "seed": seed,
        },
    )

    return {
        "reference_info": ref_info,
        "codon_usage_table": cut,
        "results": results,
        "output_paths": output_paths,
        "config": config,
    }
