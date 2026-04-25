"""Build reference translatome and custom Codon Usage Table."""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Any

import pandas as pd

from .utils import AA_TO_CODONS, CODON_TABLE_STANDARD, chunk_codons, gc_content

logger = logging.getLogger(__name__)

DEFAULT_HOUSEKEEPING_KEYWORDS: dict[str, list[str]] = {
    "ribosomal": ["ribosomal protein", "30S ribosomal", "50S ribosomal"],
    "translation": ["elongation factor", "translation initiation factor"],
    "chaperones": ["chaperone", "chaperonin", "GroEL", "GroES", "DnaK"],
    "core_metabolism": [
        "glyceraldehyde-3-phosphate dehydrogenase",
        "enolase",
        "phosphoglycerate kinase",
        "pyruvate kinase",
        "fructose-bisphosphate aldolase",
    ],
    "rna_polymerase": ["DNA-directed RNA polymerase"],
}

MIN_WI_FLOOR = 0.01


def _flatten_keywords(kw_map: dict[str, list[str]]) -> list[str]:
    out: list[str] = []
    for vals in kw_map.values():
        out.extend(vals)
    return out


def build_reference_set(
    all_cds: dict[str, dict],
    expression_data: pd.DataFrame | None = None,
    config: dict[str, Any] | None = None,
) -> tuple[list[str], dict]:
    """Select the reference gene set that feeds the CUT.

    Two modes, selected automatically by whether ``expression_data`` is given:

    - **Mode A (RNA-seq)**: match locus_tags between CDS and expression, keep
      CDS ≥ ``min_cds_length``, take the top ``top_expressed_pct``%. If that
      yields fewer than ``min_reference_genes`` entries, expand to top 20% and
      warn. Match rate below 80 % produces a warning; below 10 % is fatal.
    - **Mode B (housekeeping)**: filter CDS by product keywords and length. If
      fewer than ``fallback_min_genes`` are found, fall back to all CDS
      satisfying the length filter and warn.

    Args:
        all_cds: Output of :func:`parse_genbank`.
        expression_data: Output of :func:`parse_expression_data` (or ``None``).
        config: Full config dict. Missing keys use sensible defaults.

    Returns:
        A pair ``(sequences, info)``. ``sequences`` is the list of CDS
        nucleotide sequences selected for CUT construction. ``info`` is a
        metadata dict with keys ``mode``, ``n_genes``, ``n_total_cds``,
        ``match_rate`` (Mode A only), and ``gene_list``.

    Raises:
        ValueError: If Mode A match rate < 10 %, or if no genes are selected.
    """
    cfg = (config or {}).get("reference", {})
    min_len = cfg.get("min_cds_length", 500)
    top_pct = cfg.get("top_expressed_pct", 10)
    min_ref = cfg.get("min_reference_genes", 30)
    fallback_min = cfg.get("fallback_min_genes", 20)

    hk_keywords_cfg = (config or {}).get("housekeeping_keywords", DEFAULT_HOUSEKEEPING_KEYWORDS)
    hk_flat = _flatten_keywords(hk_keywords_cfg)

    n_total = len(all_cds)
    if expression_data is not None and not expression_data.empty:
        info = _reference_mode_rnaseq(
            all_cds, expression_data, min_len, top_pct, min_ref
        )
    else:
        info = _reference_mode_housekeeping(all_cds, min_len, hk_flat, fallback_min)

    info["n_total_cds"] = n_total

    sequences = [all_cds[tag]["nucleotide_seq"] for tag in info["gene_list"]]
    if not sequences:
        raise ValueError("Reference gene set is empty — cannot build CUT")
    info["n_genes"] = len(sequences)
    return sequences, info


def _reference_mode_rnaseq(
    all_cds: dict[str, dict],
    expression: pd.DataFrame,
    min_len: int,
    top_pct: float,
    min_ref: int,
) -> dict:
    gb_tags = set(all_cds.keys())
    expr_tags = set(expression["locus_tag"].astype(str))
    matched = gb_tags & expr_tags
    n_expr = len(expr_tags)

    match_rate = len(matched) / n_expr if n_expr else 0.0
    if match_rate < 0.10:
        raise ValueError(
            f"Too few locus_tag matches between GenBank and expression file "
            f"({match_rate:.1%}). Check that both files are from the same organism."
        )
    if match_rate < 0.80:
        logger.warning(
            "Only %.1f%% of expression locus_tags matched the GenBank file",
            match_rate * 100,
        )

    # Filter expression to matched + length filter
    expr_matched = expression[expression["locus_tag"].astype(str).isin(matched)].copy()
    expr_matched["length_nt"] = expr_matched["locus_tag"].map(
        lambda t: all_cds[t]["length_nt"]
    )
    expr_matched = expr_matched[expr_matched["length_nt"] >= min_len]
    expr_matched = expr_matched.sort_values("RPKM_average", ascending=False)

    if expr_matched.empty:
        raise ValueError(
            f"No matched CDS is ≥ {min_len} nt after length filtering."
        )

    n_filtered = len(expr_matched)
    n_top = max(1, int(round(n_filtered * top_pct / 100.0)))
    top_selection = expr_matched.head(n_top)

    if len(top_selection) < min_ref:
        expanded_pct = max(top_pct * 2, 20)
        n_expanded = max(1, int(round(n_filtered * expanded_pct / 100.0)))
        top_selection = expr_matched.head(n_expanded)
        logger.warning(
            "Top %d%% yielded only %d genes (< %d). Expanded to top %d%% → %d genes.",
            top_pct,
            n_top,
            min_ref,
            expanded_pct,
            len(top_selection),
        )

    gene_list = top_selection["locus_tag"].astype(str).tolist()
    return {
        "mode": "rnaseq",
        "match_rate": match_rate,
        "gene_list": gene_list,
    }


def _reference_mode_housekeeping(
    all_cds: dict[str, dict],
    min_len: int,
    keywords: list[str],
    fallback_min: int,
) -> dict:
    lc_keywords = [k.lower() for k in keywords]
    matched: list[str] = []
    for tag, entry in all_cds.items():
        if entry["length_nt"] < min_len:
            continue
        product = (entry.get("product") or "").lower()
        if any(kw in product for kw in lc_keywords):
            matched.append(tag)

    if len(matched) >= fallback_min:
        return {"mode": "housekeeping", "gene_list": matched}

    logger.warning(
        "Only %d housekeeping genes found (< %d). Falling back to all CDS ≥ %d nt.",
        len(matched),
        fallback_min,
        min_len,
    )
    fallback = [
        tag for tag, entry in all_cds.items() if entry["length_nt"] >= min_len
    ]
    if not fallback:
        fallback = list(all_cds.keys())
    return {"mode": "all_cds_fallback", "gene_list": fallback}


def build_codon_usage_table(reference_sequences: list[str]) -> pd.DataFrame:
    """Build a 64-row Codon Usage Table from the reference translatome.

    For each codon we compute ``count`` (absolute occurrences), ``frequency``
    (fraction of the encoding amino acid family), ``w_i`` (relative
    adaptiveness = count / max_count_within_family, Sharp & Li 1987), and
    ``rscu`` (relative synonymous codon usage). Codons with zero observed
    count receive ``w_i = MIN_WI_FLOOR`` so CAI does not collapse to zero.
    ATG and TGG have ``w_i = 1.0`` unconditionally.

    Args:
        reference_sequences: A list of coding nucleotide sequences.

    Returns:
        DataFrame with columns ``codon``, ``amino_acid``, ``count``,
        ``frequency``, ``w_i``, ``rscu``; sorted by amino acid then w_i
        descending.
    """
    counts: dict[str, int] = defaultdict(int)
    for seq in reference_sequences:
        for codon in chunk_codons(seq):
            if codon in CODON_TABLE_STANDARD:
                counts[codon] += 1

    rows = []
    for aa, codons in AA_TO_CODONS.items():
        family_total = sum(counts[c] for c in codons)
        family_max = max((counts[c] for c in codons), default=0)
        n_syn = len(codons)

        for codon in codons:
            c = counts[codon]
            if family_total == 0:
                freq = 0.0
                rscu = 0.0
            else:
                freq = c / family_total
                rscu = (c * n_syn) / family_total

            if n_syn == 1:
                w = 1.0
            elif family_max == 0:
                w = MIN_WI_FLOOR
            else:
                w = c / family_max
                if w == 0.0:
                    w = MIN_WI_FLOOR

            rows.append(
                {
                    "codon": codon,
                    "amino_acid": aa,
                    "count": c,
                    "frequency": freq,
                    "w_i": w,
                    "rscu": rscu,
                }
            )

    df = pd.DataFrame(rows)
    df = df.sort_values(by=["amino_acid", "w_i"], ascending=[True, False])
    df = df.reset_index(drop=True)
    return df


def compute_host_gc(reference_sequences: list[str]) -> float:
    """Return the GC fraction of all reference sequences concatenated."""
    if not reference_sequences:
        return 0.0
    joined = "".join(reference_sequences)
    return gc_content(joined)
