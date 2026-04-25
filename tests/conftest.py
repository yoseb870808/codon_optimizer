"""Shared pytest fixtures — build mock GenBank, Excel, and FASTA files on the fly."""

from __future__ import annotations

import random
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation
from Bio.SeqRecord import SeqRecord

FIXTURE_DIR = Path(__file__).parent / "fixtures"
FIXTURE_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Helpers to synthesize realistic-ish CDS sequences
# ---------------------------------------------------------------------------

# A small pool of codons per amino acid, weighted toward one "preferred" codon
# per AA to produce meaningful w_i spread (not uniform / random).
PREFERRED_CODONS = {
    "A": ("GCG", ["GCT", "GCC", "GCA", "GCG"]),
    "C": ("TGC", ["TGT", "TGC"]),
    "D": ("GAC", ["GAT", "GAC"]),
    "E": ("GAA", ["GAA", "GAG"]),
    "F": ("TTC", ["TTT", "TTC"]),
    "G": ("GGC", ["GGT", "GGC", "GGA", "GGG"]),
    "H": ("CAC", ["CAT", "CAC"]),
    "I": ("ATC", ["ATT", "ATC", "ATA"]),
    "K": ("AAA", ["AAA", "AAG"]),
    "L": ("CTG", ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]),
    "M": ("ATG", ["ATG"]),
    "N": ("AAC", ["AAT", "AAC"]),
    "P": ("CCG", ["CCT", "CCC", "CCA", "CCG"]),
    "Q": ("CAG", ["CAA", "CAG"]),
    "R": ("CGC", ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]),
    "S": ("AGC", ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]),
    "T": ("ACC", ["ACT", "ACC", "ACA", "ACG"]),
    "V": ("GTG", ["GTT", "GTC", "GTA", "GTG"]),
    "W": ("TGG", ["TGG"]),
    "Y": ("TAC", ["TAT", "TAC"]),
}


def _make_cds(length_aa: int, bias_preferred: float = 0.6, seed: int = 0) -> str:
    """Build a synthetic CDS: ATG + random body + TAA.

    Each body codon is the preferred codon with probability ``bias_preferred``
    and a uniform choice among synonyms otherwise. This produces a reference
    translatome with non-trivial codon usage bias.
    """
    rng = random.Random(seed)
    aa_pool = list(PREFERRED_CODONS.keys())
    aa_pool.remove("M")
    aa_pool.remove("W")  # keep rare to avoid distorting

    body_aa = [rng.choice(aa_pool) for _ in range(length_aa - 2)]
    codons = ["ATG"]
    for aa in body_aa:
        preferred, synonyms = PREFERRED_CODONS[aa]
        if rng.random() < bias_preferred:
            codons.append(preferred)
        else:
            codons.append(rng.choice(synonyms))
    codons.append("TAA")
    return "".join(codons)


def _add_cds_feature(
    record: SeqRecord,
    start: int,
    nt_seq: str,
    locus_tag: str | None,
    product: str | None = None,
    gene: str | None = None,
    complement: bool = False,
) -> int:
    """Insert ``nt_seq`` into ``record.seq`` at ``start`` and add a CDS feature.

    If ``complement`` is True, the reverse-complement of ``nt_seq`` is written
    into the record sequence so that extracting the feature on the minus strand
    recovers ``nt_seq``. Returns the next free position.
    """
    if complement:
        placed = str(Seq(nt_seq).reverse_complement())
        strand = -1
    else:
        placed = nt_seq
        strand = +1

    end = start + len(placed)
    new_seq = record.seq[:start] + Seq(placed) + record.seq[end:]
    record.seq = new_seq

    qualifiers: dict[str, list[str]] = {}
    if locus_tag is not None:
        qualifiers["locus_tag"] = [locus_tag]
    if product is not None:
        qualifiers["product"] = [product]
    if gene is not None:
        qualifiers["gene"] = [gene]

    feature = SeqFeature(
        FeatureLocation(start, end, strand=strand),
        type="CDS",
        qualifiers=qualifiers,
    )
    record.features.append(feature)
    return end


def _build_mock_genome(path: Path) -> dict:
    """Build a mock GenBank file and return metadata about its contents."""
    # Size generously so every placement fits inside the record seq.
    total_len = 60_000
    record = SeqRecord(
        Seq("A" * total_len),
        id="MOCK001",
        name="MOCK001",
        description="Mock genome for codon optimizer tests",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )

    cursor = 100
    gap = 30

    metadata: dict = {
        "locus_tags": [],
        "housekeeping_tags": [],
        "skipped_expected": [],
        "complement_tag": None,
        "pseudogene_tag": None,
        "bad_length_tag": None,
        "no_tag_position": None,
    }

    # 10 regular CDS GENE_001..GENE_010, lengths varying
    aa_lengths = [200, 700, 450, 350, 180, 250, 600, 320, 280, 190]
    for i, aa_len in enumerate(aa_lengths, start=1):
        tag = f"GENE_{i:03d}"
        seq = _make_cds(aa_len, seed=i)
        cursor = _add_cds_feature(record, cursor, seq, locus_tag=tag)
        cursor += gap
        metadata["locus_tags"].append(tag)

    # 3 housekeeping CDS (length >= 500 nt — use aa=200 → 603 nt)
    hk_specs = [
        ("GENE_HK_RPS1", "30S ribosomal protein S1"),
        ("GENE_HK_EFTU", "elongation factor Tu"),
        ("GENE_HK_GROEL", "chaperonin GroEL"),
    ]
    for tag, product in hk_specs:
        seq = _make_cds(220, seed=hash(tag) % 10_000)
        cursor = _add_cds_feature(record, cursor, seq, locus_tag=tag, product=product)
        cursor += gap
        metadata["housekeeping_tags"].append(tag)

    # 1 CDS on the complement strand
    comp_tag = "GENE_COMP"
    comp_seq = _make_cds(150, seed=9999)
    cursor = _add_cds_feature(
        record, cursor, comp_seq, locus_tag=comp_tag, complement=True
    )
    cursor += gap
    metadata["complement_tag"] = comp_tag
    metadata["complement_seq"] = comp_seq

    # 1 CDS with length not divisible by 3 — craft manually
    bad_tag = "GENE_BAD_LEN"
    bad_seq = "ATGAAAGTT" + "AAAA"  # 13 nt, not divisible by 3
    end = cursor + len(bad_seq)
    record.seq = record.seq[:cursor] + Seq(bad_seq) + record.seq[end:]
    record.features.append(
        SeqFeature(
            FeatureLocation(cursor, end, strand=+1),
            type="CDS",
            qualifiers={"locus_tag": [bad_tag]},
        )
    )
    cursor = end + gap
    metadata["bad_length_tag"] = bad_tag

    # 1 pseudogene with internal stop
    pseudo_tag = "GENE_PSEUDO"
    # Build: ATG + TAA + a bunch of codons + TAA → internal stop
    pseudo_seq = "ATG" + "TAA" + "GCGGCGGCGGCG" + "TAA"
    assert len(pseudo_seq) % 3 == 0
    end = cursor + len(pseudo_seq)
    record.seq = record.seq[:cursor] + Seq(pseudo_seq) + record.seq[end:]
    record.features.append(
        SeqFeature(
            FeatureLocation(cursor, end, strand=+1),
            type="CDS",
            qualifiers={"locus_tag": [pseudo_tag]},
        )
    )
    cursor = end + gap
    metadata["pseudogene_tag"] = pseudo_tag

    # 1 CDS without locus_tag
    notag_seq = _make_cds(60, seed=4242)
    end = cursor + len(notag_seq)
    record.seq = record.seq[:cursor] + Seq(notag_seq) + record.seq[end:]
    record.features.append(
        SeqFeature(
            FeatureLocation(cursor, end, strand=+1),
            type="CDS",
            qualifiers={},  # intentionally no locus_tag
        )
    )
    metadata["no_tag_position"] = (cursor, end)
    cursor = end + gap

    # Extend the record to ensure seq is long enough (pad with As)
    if len(record.seq) < total_len:
        record.seq = record.seq + Seq("A" * (total_len - len(record.seq)))

    SeqIO.write([record], str(path), "genbank")
    return metadata


def _build_mock_expression(path: Path, genome_meta: dict) -> None:
    """Build a mock Excel expression file aligned with the mock genome."""
    rows = []
    # GENE_001..GENE_010 — GENE_001-003 highest, decreasing from there
    expressions = [5000, 4500, 4000, 200, 150, 80, 50, 20, 5, 0.5]
    for tag, rpkm_avg in zip(genome_meta["locus_tags"], expressions):
        r1 = rpkm_avg * 0.95
        r2 = rpkm_avg * 1.05
        r3 = rpkm_avg * 1.00
        rows.append(
            {
                "locus_tag": tag,
                "RPKM_rep1": r1,
                "RPKM_rep2": r2,
                "RPKM_rep3": r3,
                "RPKM_average": rpkm_avg,
            }
        )
    # 2 mismatched locus_tags not present in the GenBank file
    for fake_tag in ("MISMATCH_001", "MISMATCH_002"):
        rows.append(
            {
                "locus_tag": fake_tag,
                "RPKM_rep1": 10,
                "RPKM_rep2": 11,
                "RPKM_rep3": 9,
                "RPKM_average": 10,
            }
        )
    df = pd.DataFrame(rows)
    df.to_excel(path, index=False)


MOCK_PROTEIN_FASTA = """>test_protein_short A small test protein
MKFLILVAS
>test_protein_medium A medium length protein
MKFLILVASAGTTLMITGNPKLRDEWYQHCFS
"""

MOCK_DNA_FASTA = """>test_dna_ecoli_style An E.coli-like CDS (GC-rich codons)
ATGAAATTCCTGATCCTGGTGGCCTCCGCCGGCACCACCCTGATGATCACCGGCAACCCGAAACTGCGCGATGAATGGTACCAGCACTGCTTCTCC
>test_dna_at_rich An AT-rich CDS
ATGAAATTTTTAATTTTAATTGCATCAGCAGGTACAACATTAATGATTACAGGTAATCCAAAATTAAGAGATGAATGGTATAAA
"""


def _write_mock_fastas(proto_path: Path, dna_path: Path) -> None:
    proto_path.write_text(MOCK_PROTEIN_FASTA)
    dna_path.write_text(MOCK_DNA_FASTA)


# ---------------------------------------------------------------------------
# Session-scoped fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def mock_genome_path() -> Path:
    path = FIXTURE_DIR / "mock_genome.gbk"
    meta_path = FIXTURE_DIR / "mock_genome.meta.json"
    if not path.exists() or not meta_path.exists():
        meta = _build_mock_genome(path)
        import json
        meta_path.write_text(json.dumps({k: v for k, v in meta.items() if k != "complement_seq"}))
        # save complement seq separately as FASTA
        (FIXTURE_DIR / "complement_seq.txt").write_text(meta["complement_seq"])
    return path


@pytest.fixture(scope="session")
def mock_genome_meta(mock_genome_path) -> dict:
    import json
    meta_path = FIXTURE_DIR / "mock_genome.meta.json"
    data = json.loads(meta_path.read_text())
    data["complement_seq"] = (FIXTURE_DIR / "complement_seq.txt").read_text()
    return data


@pytest.fixture(scope="session")
def mock_expression_path(mock_genome_meta) -> Path:
    path = FIXTURE_DIR / "mock_expression.xlsx"
    if not path.exists():
        _build_mock_expression(path, mock_genome_meta)
    return path


@pytest.fixture(scope="session")
def mock_protein_fasta_path() -> Path:
    path = FIXTURE_DIR / "mock_protein.fasta"
    if not path.exists():
        path.write_text(MOCK_PROTEIN_FASTA)
    return path


@pytest.fixture(scope="session")
def mock_dna_fasta_path() -> Path:
    path = FIXTURE_DIR / "mock_dna.fasta"
    if not path.exists():
        path.write_text(MOCK_DNA_FASTA)
    return path


@pytest.fixture(scope="session")
def project_root() -> Path:
    return Path(__file__).resolve().parent.parent
