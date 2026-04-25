"""Input parsing and validation for GenBank, Excel, FASTA, and config YAML."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from Bio import SeqIO

from .utils import (
    AMBIGUOUS_AMINO_ACIDS,
    DNA_BASES,
    STANDARD_AMINO_ACIDS,
    has_internal_stop,
    is_dna,
    is_protein,
    translate,
)

logger = logging.getLogger(__name__)

REQUIRED_EXPRESSION_COLUMNS = [
    "locus_tag",
    "RPKM_rep1",
    "RPKM_rep2",
    "RPKM_rep3",
    "RPKM_average",
]

# Common column-name variants the parser auto-renames to canonical form.
# Match is case-insensitive against the lowercased column header.
EXPRESSION_COLUMN_ALIASES: dict[str, list[str]] = {
    "locus_tag": ["locus_tag", "locus", "gene_id", "geneid", "gene", "id", "tag"],
    "RPKM_rep1": ["rpkm_rep1", "rpkm_wt1", "rpkm_1", "rpkm1",
                  "rep1", "wt1", "replicate_1", "replicate1"],
    "RPKM_rep2": ["rpkm_rep2", "rpkm_wt2", "rpkm_2", "rpkm2",
                  "rep2", "wt2", "replicate_2", "replicate2"],
    "RPKM_rep3": ["rpkm_rep3", "rpkm_wt3", "rpkm_3", "rpkm3",
                  "rep3", "wt3", "replicate_3", "replicate3"],
    "RPKM_average": ["rpkm_average", "rpkm_avg", "average_rpkm", "avg_rpkm",
                     "mean_rpkm", "rpkm_mean", "average", "mean"],
}


def _read_expression_table(path: Path) -> pd.DataFrame:
    """Dispatch to the right pandas reader based on file extension."""
    ext = path.suffix.lower()
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(path)
    if ext == ".csv":
        # Quoted commas in numbers (e.g. "1,503.54") survive read_csv when the
        # field is properly quoted; we coerce them later.
        return pd.read_csv(path)
    if ext in (".tsv", ".txt"):
        return pd.read_csv(path, sep="\t")
    # Last-ditch: sniff the delimiter
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception as exc:
        raise ValueError(
            f"Cannot determine format of {path}. "
            f"Supported extensions: .xlsx, .xls, .csv, .tsv, .txt"
        ) from exc


def _canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename column variants to the canonical schema (in-place safe copy)."""
    lower_to_actual = {c.lower().strip(): c for c in df.columns}
    rename_map: dict[str, str] = {}
    for canonical, aliases in EXPRESSION_COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in lower_to_actual:
                actual = lower_to_actual[alias]
                if actual != canonical:
                    rename_map[actual] = canonical
                break
    if rename_map:
        df = df.rename(columns=rename_map)
        logger.info("Renamed expression columns: %s", rename_map)
    return df


def _coerce_numeric(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Strip thousands separators and coerce to float; drop unparseable rows."""
    df = df.copy()
    for col in columns:
        # Strip commas/whitespace for any non-numeric column (handles object,
        # pandas StringDtype, mixed types). Pure numeric columns pass through.
        if not pd.api.types.is_numeric_dtype(df[col]):
            df[col] = (
                df[col]
                .astype(str)
                .str.replace(",", "", regex=False)
                .str.strip()
            )
        df[col] = pd.to_numeric(df[col], errors="coerce")
    n_dropped = int(df[columns].isna().any(axis=1).sum())
    if n_dropped:
        logger.warning(
            "Dropped %d row(s) with non-numeric expression values", n_dropped
        )
        df = df.dropna(subset=columns).reset_index(drop=True)
    return df


def parse_genbank(filepath: str | os.PathLike) -> dict[str, dict]:
    """Parse a GenBank file and extract valid CDS features.

    Walks every record in the file (supporting multi-contig GenBank) and
    emits one entry per CDS keyed by ``locus_tag``. Features missing a
    ``locus_tag``, not divisible by 3, or containing an internal stop codon
    are skipped with a logged warning.

    Args:
        filepath: Path to a ``.gb`` or ``.gbk`` file.

    Returns:
        Mapping of ``locus_tag`` to a dict with keys: ``locus_tag``, ``gene``,
        ``product``, ``nucleotide_seq``, ``amino_acid_seq``, ``length_nt``,
        ``length_aa``, ``location``, ``contig``.

    Raises:
        FileNotFoundError: If ``filepath`` does not exist.
        ValueError: If the file yields zero valid CDS features.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"GenBank file not found: {filepath}")

    result: dict[str, dict] = {}
    n_skipped_no_tag = 0
    n_skipped_length = 0
    n_skipped_pseudogene = 0
    n_skipped_dup = 0

    for record in SeqIO.parse(str(path), "genbank"):
        record_seq = record.seq
        contig_id = record.id or record.name or "unknown"
        for feature in record.features:
            if feature.type != "CDS":
                continue

            locus_tag_list = feature.qualifiers.get("locus_tag")
            if not locus_tag_list:
                n_skipped_no_tag += 1
                logger.warning(
                    "Skipping CDS without locus_tag at %s on %s",
                    feature.location,
                    contig_id,
                )
                continue
            locus_tag = locus_tag_list[0]

            try:
                nt_seq = str(feature.extract(record_seq)).upper()
            except Exception as exc:
                logger.warning(
                    "Skipping CDS %s: could not extract sequence (%s)",
                    locus_tag,
                    exc,
                )
                continue

            if len(nt_seq) == 0 or len(nt_seq) % 3 != 0:
                n_skipped_length += 1
                logger.warning(
                    "Skipping CDS %s: length %d not divisible by 3",
                    locus_tag,
                    len(nt_seq),
                )
                continue

            aa_seq = translate(nt_seq)
            # Strip a single trailing stop codon for the internal-stop check
            aa_body = aa_seq[:-1] if aa_seq.endswith("*") else aa_seq
            if has_internal_stop(aa_body):
                n_skipped_pseudogene += 1
                logger.warning(
                    "Skipping CDS %s: internal stop codon (likely pseudogene)",
                    locus_tag,
                )
                continue
            if "X" in aa_body:
                logger.warning(
                    "CDS %s contains ambiguous codons (X); keeping but flagged",
                    locus_tag,
                )

            gene = feature.qualifiers.get("gene", [None])[0]
            product = feature.qualifiers.get("product", [None])[0]

            if locus_tag in result:
                n_skipped_dup += 1
                logger.warning("Duplicate locus_tag %s; keeping first occurrence", locus_tag)
                continue

            result[locus_tag] = {
                "locus_tag": locus_tag,
                "gene": gene,
                "product": product,
                "nucleotide_seq": nt_seq,
                "amino_acid_seq": aa_seq,
                "length_nt": len(nt_seq),
                "length_aa": len(aa_seq.rstrip("*")),
                "location": str(feature.location),
                "contig": contig_id,
            }

    if not result:
        raise ValueError(
            f"No valid CDS features found in {filepath} "
            f"(skipped: {n_skipped_no_tag} without locus_tag, "
            f"{n_skipped_length} with invalid length, "
            f"{n_skipped_pseudogene} pseudogenes)"
        )

    logger.info(
        "Parsed %d CDS from %s (skipped %d no-tag, %d length, %d pseudogene, %d dup)",
        len(result),
        path.name,
        n_skipped_no_tag,
        n_skipped_length,
        n_skipped_pseudogene,
        n_skipped_dup,
    )
    return result


def parse_expression_data(filepath: str | os.PathLike) -> pd.DataFrame:
    """Parse an RNA-seq expression file (xlsx / xls / csv / tsv / txt).

    Required columns (canonical names): ``locus_tag``, ``RPKM_rep1``,
    ``RPKM_rep2``, ``RPKM_rep3``, ``RPKM_average``.

    Common variants are auto-renamed before validation:
      * Replicate columns: ``RPKM_WT1``/``WT1``/``rep1``/``replicate_1`` →
        ``RPKM_rep1`` (and similarly for 2, 3).
      * Average column: ``RPKM_avg``/``mean_rpkm``/``average`` →
        ``RPKM_average``.
      * Locus column: ``locus``/``gene_id``/``geneid``/``gene`` →
        ``locus_tag``.

    Numeric columns may contain thousands separators (``"1,503.54"``);
    they are stripped automatically. Rows with unparseable expression
    values are dropped with a logged warning.

    Args:
        filepath: Path to ``.xlsx``, ``.xls``, ``.csv``, ``.tsv``, or
            ``.txt`` (tab-separated).

    Returns:
        DataFrame sorted by ``RPKM_average`` descending.

    Raises:
        FileNotFoundError: If ``filepath`` does not exist.
        ValueError: If required columns cannot be resolved after aliasing,
            or any expression value is negative.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Expression file not found: {filepath}")

    df = _read_expression_table(path)
    df = _canonicalize_columns(df)

    missing = [c for c in REQUIRED_EXPRESSION_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Expression file {filepath} missing required columns: {missing}.\n"
            f"Found columns: {list(df.columns)}.\n"
            f"Recognized aliases for each canonical column:\n"
            + "\n".join(f"  {k}: {v}" for k, v in EXPRESSION_COLUMN_ALIASES.items())
        )

    df = _coerce_numeric(df, ["RPKM_rep1", "RPKM_rep2", "RPKM_rep3", "RPKM_average"])

    for col in ["RPKM_rep1", "RPKM_rep2", "RPKM_rep3", "RPKM_average"]:
        if (df[col] < 0).any():
            raise ValueError(f"Column {col} contains negative values in {filepath}")

    df["locus_tag"] = df["locus_tag"].astype(str).str.strip()
    df = df[REQUIRED_EXPRESSION_COLUMNS]  # canonical column order
    df = df.sort_values("RPKM_average", ascending=False).reset_index(drop=True)
    return df


def parse_target_sequences(filepath: str | os.PathLike) -> list[dict]:
    """Parse a FASTA file and auto-detect protein vs DNA per record.

    Args:
        filepath: Path to a ``.fasta`` or ``.fa`` file.

    Returns:
        List of dicts with keys: ``header``, ``input_type`` (``"protein"`` or
        ``"dna"``), ``protein_seq``, ``original_dna`` (None for protein input).

    Raises:
        FileNotFoundError: If ``filepath`` does not exist.
        ValueError: If the file contains no valid sequences, or if a sequence
            contains ambiguous amino acids, or if a DNA sequence has a length
            not divisible by 3, or if a DNA sequence translates to X.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Target FASTA not found: {filepath}")

    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"FASTA file {filepath} contains no sequences")

    result: list[dict] = []
    for record in records:
        header = record.description or record.id
        raw_seq = str(record.seq).strip()
        if not raw_seq:
            raise ValueError(f"Empty sequence for FASTA record '{header}'")

        if is_dna(raw_seq):
            dna = raw_seq.upper()
            if "N" in dna:
                logger.warning("DNA sequence '%s' contains N bases", header)
            if len(dna) % 3 != 0:
                raise ValueError(
                    f"DNA sequence '{header}' length {len(dna)} is not divisible by 3"
                )
            if not dna.startswith("ATG"):
                logger.warning("DNA sequence '%s' does not start with ATG", header)
            protein = translate(dna)
            # Strip one trailing stop
            protein_body = protein[:-1] if protein.endswith("*") else protein
            if "X" in protein_body:
                raise ValueError(
                    f"DNA sequence '{header}' translates to ambiguous residues (X) — "
                    f"cannot optimize"
                )
            if has_internal_stop(protein_body):
                raise ValueError(
                    f"DNA sequence '{header}' has an internal stop codon"
                )
            result.append(
                {
                    "header": header,
                    "input_type": "dna",
                    "protein_seq": protein_body,
                    "original_dna": dna,
                }
            )
        else:
            # Not pure DNA — treat as protein. Reject ambiguous residues first.
            protein = raw_seq.upper().rstrip("*")
            ambig = set(protein) & AMBIGUOUS_AMINO_ACIDS
            if ambig:
                raise ValueError(
                    f"Protein sequence '{header}' contains ambiguous amino acids: "
                    f"{sorted(ambig)}"
                )
            extra = set(protein) - STANDARD_AMINO_ACIDS
            if extra:
                raise ValueError(
                    f"Protein sequence '{header}' contains non-standard residues: "
                    f"{sorted(extra)}"
                )
            if not protein:
                raise ValueError(f"Empty protein body for FASTA record '{header}'")
            result.append(
                {
                    "header": header,
                    "input_type": "protein",
                    "protein_seq": protein,
                    "original_dna": None,
                }
            )

    return result


_DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent.parent / "config_default.yaml"


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge ``override`` into ``base`` (base is not mutated)."""
    out = dict(base)
    for key, value in override.items():
        if (
            key in out
            and isinstance(out[key], dict)
            and isinstance(value, dict)
        ):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def load_config(filepath: str | os.PathLike | None = None) -> dict[str, Any]:
    """Load the default YAML config and merge a user config on top of it.

    Args:
        filepath: Path to a user YAML file to merge on top of defaults. If
            ``None``, only defaults are returned.

    Returns:
        Merged config dict.

    Raises:
        FileNotFoundError: If ``filepath`` is given but does not exist.
    """
    if not _DEFAULT_CONFIG_PATH.exists():
        raise FileNotFoundError(
            f"Default config missing at {_DEFAULT_CONFIG_PATH}"
        )
    with open(_DEFAULT_CONFIG_PATH, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f) or {}

    if filepath is None:
        return config

    user_path = Path(filepath)
    if not user_path.exists():
        raise FileNotFoundError(f"Config file not found: {filepath}")
    with open(user_path, "r", encoding="utf-8") as f:
        user_config = yaml.safe_load(f) or {}
    return _deep_merge(config, user_config)
