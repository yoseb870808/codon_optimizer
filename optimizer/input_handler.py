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
    """Parse an RNA-seq Excel file with fixed column structure.

    Required columns (exact names): ``locus_tag``, ``RPKM_rep1``, ``RPKM_rep2``,
    ``RPKM_rep3``, ``RPKM_average``. Additional columns are kept.

    Args:
        filepath: Path to a ``.xlsx`` file.

    Returns:
        DataFrame sorted by ``RPKM_average`` descending.

    Raises:
        FileNotFoundError: If ``filepath`` does not exist.
        ValueError: If required columns are missing or contain non-numeric /
            negative expression values.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Expression file not found: {filepath}")

    df = pd.read_excel(path)
    missing = [c for c in REQUIRED_EXPRESSION_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Expression file {filepath} missing required columns: {missing}"
        )

    if df.columns[0] != "locus_tag":
        raise ValueError(
            f"Expression file first column must be 'locus_tag', got '{df.columns[0]}'"
        )
    if df.columns[-1] != "RPKM_average":
        # Allow extra cols after average? Spec says last must be RPKM_average.
        # We enforce strictly per spec.
        raise ValueError(
            f"Expression file last column must be 'RPKM_average', got '{df.columns[-1]}'"
        )

    for col in ["RPKM_rep1", "RPKM_rep2", "RPKM_rep3", "RPKM_average"]:
        if not pd.api.types.is_numeric_dtype(df[col]):
            raise ValueError(f"Column {col} must be numeric in {filepath}")
        if (df[col] < 0).any():
            raise ValueError(f"Column {col} contains negative values in {filepath}")

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
