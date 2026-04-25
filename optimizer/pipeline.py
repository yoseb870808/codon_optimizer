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
import copy
import numpy as np

from .metrics import calculate_cai, compare_codons, composite_score
from .mrna_repair import REPAIR_MODES, repair_5prime_structure
from .quality_control import check_mrna_structure, run_all_qc
from .sequence_repair import find_motifs
from .utils import translate as _translate
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


def _repair_seed(base: int | None, target_idx: int, variant_idx: int) -> int | None:
    """Deterministic per-(target, variant) seed for the repair RNG."""
    if base is None:
        return None
    return ((base * 2654435761) ^ (target_idx * 786433) ^ (variant_idx * 40503)) & 0x7FFFFFFF


def _expected_protein(dna: str) -> str:
    """Translate ``dna``, dropping any single trailing stop codon."""
    aa = _translate(dna)
    return aa[:-1] if aa.endswith("*") else aa


def _assert_protein_identity(dna: str, expected_protein: str, context: str) -> None:
    """Hard runtime invariant: a translated DNA must equal the expected protein.

    Raises ``RuntimeError`` if the translation differs. This is a paranoid
    check that should never fire in normal operation — every code path that
    mutates a candidate sequence does so with synonymous substitutions only.
    If this fires, treat it as a critical bug and stop the pipeline.
    """
    actual = _expected_protein(dna)
    if actual != expected_protein:
        # Find the first divergent residue for the error message
        diff_at = -1
        for i, (a, b) in enumerate(zip(actual, expected_protein)):
            if a != b:
                diff_at = i
                break
        if diff_at < 0:
            diff_at = min(len(actual), len(expected_protein))
        raise RuntimeError(
            f"PROTEIN INTEGRITY VIOLATION at {context}: "
            f"expected {len(expected_protein)} AA but translation differs at position "
            f"{diff_at + 1} (got '{actual[diff_at:diff_at+10]}…' vs "
            f"expected '{expected_protein[diff_at:diff_at+10]}…'). "
            f"This is a critical bug — the optimizer must never produce a non-synonymous mutation."
        )


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
    mrna_repair_mode: str | None = None,
    write_pre_repair_outputs: bool | None = None,
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

    mrna_cfg = config.setdefault("mrna_repair", {})
    if mrna_repair_mode is not None:
        mrna_cfg["mode"] = mrna_repair_mode
    repair_mode = mrna_cfg.get("mode", "targeted")
    if repair_mode not in REPAIR_MODES:
        raise ValueError(
            f"mrna_repair.mode '{repair_mode}' invalid. Valid: {REPAIR_MODES}"
        )
    if write_pre_repair_outputs is None:
        write_pre_repair_outputs = mrna_cfg.get("write_pre_repair_outputs", True)
    if repair_mode == "off":
        write_pre_repair_outputs = False  # nothing to compare

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

    qc_cfg = config.get("qc", {})
    mfe_threshold = qc_cfg.get("mrna_mfe_threshold", -30.0)
    mrna_window = qc_cfg.get("mrna_5prime_window", 150)
    scoring_cfg = config.get("variant_scoring", {})
    scoring_weights = {
        "cai_weight": scoring_cfg.get("cai_weight", 0.5),
        "gc_weight": scoring_cfg.get("gc_weight", 0.3),
        "mrna_weight": scoring_cfg.get("mrna_weight", 0.2),
    }
    forbidden_motifs_cfg = config.get("forbidden_motifs") or []
    if not forbidden_motifs_cfg and config.get("optimization", {}).get("use_default_forbidden", True):
        from .sequence_repair import DEFAULT_FORBIDDEN_MOTIFS
        forbidden_motifs_cfg = DEFAULT_FORBIDDEN_MOTIFS

    pre_repair_results: list[dict[str, Any]] = []
    repaired_results: list[dict[str, Any]] = []

    for idx, t in enumerate(targets):
        logger.info(
            "Optimizing [%d/%d]: %s (%s, %d AA)",
            idx + 1, len(targets), t["header"].split()[0], t["input_type"], len(t["protein_seq"]),
        )
        pre_variants = optimize(
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

        # Invariant: every pre-repair variant must translate to the input protein
        for v in pre_variants:
            _assert_protein_identity(
                v["dna_sequence"], t["protein_seq"],
                context=f"post-optimize variant {v['variant_number']} of '{t['header'].split()[0]}'",
            )

        original_cai = None
        original_gc = None
        if t["input_type"] == "dna" and t["original_dna"]:
            original_cai = calculate_cai(t["original_dna"], cut)
            original_gc = gc_content(t["original_dna"])

        # Pre-repair record (always built; written only if requested)
        pre_recommended = next(v for v in pre_variants if v["is_recommended"])
        pre_qc = run_all_qc(
            pre_recommended["dna_sequence"], target_gc=host_gc, config=config,
        )
        pre_repair_results.append(
            {
                "header": t["header"],
                "input_type": t["input_type"],
                "protein_seq": t["protein_seq"],
                "original_dna": t["original_dna"],
                "length_aa": len(t["protein_seq"]),
                "variants": pre_variants,
                "qc": pre_qc,
                "original_cai": original_cai,
                "original_gc": original_gc,
                "host_gc": host_gc,
            }
        )

        # Repair pass — only on variants flagged for strong structure
        if repair_mode == "off":
            repaired_variants = pre_variants
        else:
            repaired_variants = []
            for v_idx, v in enumerate(pre_variants):
                v_copy = copy.deepcopy(v)
                if v_copy["has_strong_structure"]:
                    rng = np.random.default_rng(_repair_seed(seed, idx, v_idx))
                    new_dna, info = repair_5prime_structure(
                        v_copy["dna_sequence"],
                        cut,
                        mode=repair_mode,
                        window=mrna_window,
                        mfe_threshold=mfe_threshold,
                        targeted_max_attempts=int(mrna_cfg.get("targeted_max_attempts", 30)),
                        annealing_iterations=int(mrna_cfg.get("annealing_iterations", 200)),
                        annealing_t_start=float(mrna_cfg.get("annealing_t_start", 5.0)),
                        annealing_t_end=float(mrna_cfg.get("annealing_t_end", 0.05)),
                        cai_lambda=float(mrna_cfg.get("cai_lambda", 3.0)),
                        rng=rng,
                    )
                    # === Invariant 1: synonymous-only mutation ===
                    _assert_protein_identity(
                        new_dna, t["protein_seq"],
                        context=f"post-mRNA-repair variant {v['variant_number']} of '{t['header'].split()[0]}'",
                    )

                    v_copy["dna_sequence"] = new_dna
                    v_copy["cai_score"] = info.final_cai
                    v_copy["gc_content"] = gc_content(new_dna)
                    v_copy["mfe_5prime"] = info.final_mfe
                    v_copy["has_strong_structure"] = info.final_mfe <= mfe_threshold
                    v_copy["repair_info"] = info.as_dict()

                    # Refresh dot-bracket structure (was stale)
                    new_struct = check_mrna_structure(
                        new_dna, window_5prime=mrna_window, mfe_threshold=mfe_threshold,
                    )
                    v_copy["mrna_structure"] = new_struct["structure"]

                    # Re-check forbidden motifs — repair could have introduced new ones
                    if forbidden_motifs_cfg:
                        new_hits = find_motifs(new_dna, forbidden_motifs_cfg)
                        v_copy["forbidden_hits"] = [
                            {
                                "motif": h.motif_name,
                                "position": h.start,
                                "length": h.end - h.start,
                                "matched": h.matched,
                                "strand": h.strand,
                            }
                            for h in new_hits
                        ]
                    # Recompute composite_score with the post-repair sequence
                    forbidden_penalty = 0.05 * len(v_copy["forbidden_hits"])
                    v_copy["composite_score"] = composite_score(
                        cai=v_copy["cai_score"],
                        gc_content_value=v_copy["gc_content"],
                        host_gc=host_gc,
                        mfe=v_copy["mfe_5prime"],
                        mfe_threshold=mfe_threshold,
                        weights=scoring_weights,
                    ) - forbidden_penalty
                else:
                    v_copy["repair_info"] = {
                        "method": "skipped_already_good",
                        "initial_mfe": v_copy["mfe_5prime"],
                        "final_mfe": v_copy["mfe_5prime"],
                        "initial_cai": v_copy["cai_score"],
                        "final_cai": v_copy["cai_score"],
                        "attempts": 0,
                        "succeeded": True,
                        "delta_mfe": 0.0,
                        "delta_cai": 0.0,
                    }
                repaired_variants.append(v_copy)
            # Re-rank after repair (a previously-rejected variant might now win)
            repaired_variants.sort(
                key=lambda v: (
                    len(v["forbidden_hits"]) == 0,
                    not v["has_strong_structure"],
                    v["composite_score"],
                ),
                reverse=True,
            )
            # Renumber so variant_number reflects post-repair rank (1 = best)
            for new_i, v in enumerate(repaired_variants, start=1):
                v["variant_number"] = new_i
                v["is_recommended"] = False
            repaired_variants[0]["is_recommended"] = True

        post_recommended = next(v for v in repaired_variants if v["is_recommended"])
        post_qc = run_all_qc(
            post_recommended["dna_sequence"], target_gc=host_gc, config=config,
        )
        repaired_results.append(
            {
                "header": t["header"],
                "input_type": t["input_type"],
                "protein_seq": t["protein_seq"],
                "original_dna": t["original_dna"],
                "length_aa": len(t["protein_seq"]),
                "variants": repaired_variants,
                "qc": post_qc,
                "original_cai": original_cai,
                "original_gc": original_gc,
                "host_gc": host_gc,
            }
        )

    # Final results = repaired (or pre, if mode == off)
    results = repaired_results

    # === Final integrity sweep ===
    # Independently re-translate every variant in every output (pre + post)
    # against the input protein. This is the last guardrail before writing
    # any output file. If anything mutated non-synonymously anywhere in the
    # pipeline, we abort here rather than produce a misleading FASTA.
    for bucket_name, bucket in (("post-repair", repaired_results), ("pre-repair", pre_repair_results)):
        for res in bucket:
            for v in res["variants"]:
                _assert_protein_identity(
                    v["dna_sequence"], res["protein_seq"],
                    context=f"final {bucket_name} sweep — variant {v['variant_number']} of '{res['header'].split()[0]}'",
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
            opt_body = rec["dna_sequence"]
            orig = res["original_dna"]
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
            "mrna_repair_mode": repair_mode,
        },
    )

    # Pre-repair outputs in a sibling directory + side-by-side comparison
    if write_pre_repair_outputs and repair_mode != "off":
        pre_dir = out_path / "pre_repair"
        pre_dir.mkdir(exist_ok=True)
        pre_per_seq: list[str] = []
        for res in pre_repair_results:
            pre_per_seq.append(
                write_optimized_fasta(res["header"], res["variants"], pre_dir)
            )
        output_paths["pre_repair_fastas"] = pre_per_seq
        output_paths["pre_repair_combined_fasta"] = write_combined_fasta(
            pre_repair_results, pre_dir
        )
        output_paths["pre_repair_report_tsv"] = write_report_tsv(
            pre_repair_results, pre_dir
        )

        # Side-by-side comparison
        comp_rows = []
        for pre_res, post_res in zip(pre_repair_results, repaired_results):
            pre_rec = next(v for v in pre_res["variants"] if v["is_recommended"])
            post_rec = next(v for v in post_res["variants"] if v["is_recommended"])
            info = post_rec.get("repair_info", {})
            comp_rows.append({
                "sequence_name": pre_res["header"].split()[0],
                "repair_mode": repair_mode,
                "repair_method_used": info.get("method", "n/a"),
                "repair_attempts": info.get("attempts", 0),
                "repair_succeeded": info.get("succeeded"),
                "pre_cai": pre_rec["cai_score"],
                "post_cai": post_rec["cai_score"],
                "delta_cai": post_rec["cai_score"] - pre_rec["cai_score"],
                "pre_gc": pre_rec["gc_content"],
                "post_gc": post_rec["gc_content"],
                "pre_mfe_5prime": pre_rec["mfe_5prime"],
                "post_mfe_5prime": post_rec["mfe_5prime"],
                "delta_mfe": (post_rec["mfe_5prime"] or 0.0) - (pre_rec["mfe_5prime"] or 0.0),
                "pre_strong_structure": pre_rec["has_strong_structure"],
                "post_strong_structure": post_rec["has_strong_structure"],
            })
        comp_path = out_path / "mrna_repair_comparison.tsv"
        pd.DataFrame(comp_rows).to_csv(comp_path, sep="\t", index=False)
        output_paths["mrna_repair_comparison_tsv"] = str(comp_path)

    return {
        "reference_info": ref_info,
        "codon_usage_table": cut,
        "results": results,
        "pre_repair_results": pre_repair_results,
        "output_paths": output_paths,
        "config": config,
    }
