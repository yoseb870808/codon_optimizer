"""Shared utilities for the codon optimizer."""

from __future__ import annotations

CODON_TABLE_STANDARD: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA_TO_CODONS: dict[str, list[str]] = {}
for _codon, _aa in CODON_TABLE_STANDARD.items():
    AA_TO_CODONS.setdefault(_aa, []).append(_codon)
for _aa in AA_TO_CODONS:
    AA_TO_CODONS[_aa].sort()

STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AMINO_ACIDS = set("BJOUXZ")
DNA_BASES = set("ATCGN")

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def translate(dna_seq: str) -> str:
    """Translate DNA sequence to protein using the standard genetic code.

    Args:
        dna_seq: DNA nucleotide sequence (case-insensitive). Must be divisible
            by 3 (trailing partial codon is ignored with no warning here).

    Returns:
        Amino acid sequence. Unknown codons (containing N or otherwise absent
        from the standard table) translate to ``X``. Stop codons translate to
        ``*``.
    """
    seq = dna_seq.upper()
    out: list[str] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        out.append(CODON_TABLE_STANDARD.get(codon, "X"))
    return "".join(out)


def is_protein(seq: str) -> bool:
    """Return True if the sequence contains only standard amino acids or stop.

    A sequence of only ``ATCGN`` is ambiguous (could be DNA); this function
    returns False in that case so :func:`is_dna` wins. Empty strings return
    False.
    """
    if not seq:
        return False
    s = seq.upper().rstrip("*")
    if not s:
        return False
    if set(s) <= DNA_BASES:
        return False
    return set(s) <= STANDARD_AMINO_ACIDS


def is_dna(seq: str) -> bool:
    """Return True if the sequence contains only DNA bases (ATCGN)."""
    if not seq:
        return False
    return set(seq.upper()) <= DNA_BASES


def gc_content(seq: str) -> float:
    """Return the GC fraction of a nucleotide sequence.

    Returns 0.0 for an empty sequence. Non-ACGT bases (including N) are
    counted in the denominator but contribute 0 to the numerator.
    """
    if not seq:
        return 0.0
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    return gc / len(s)


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence (case preserved)."""
    return seq.translate(_COMPLEMENT)[::-1]


def has_internal_stop(protein_seq: str) -> bool:
    """Return True if the protein has a stop codon before the final position."""
    if "*" not in protein_seq:
        return False
    idx = protein_seq.find("*")
    return idx < len(protein_seq) - 1


def chunk_codons(dna_seq: str) -> list[str]:
    """Split a DNA sequence into a list of 3-letter codons (trailing partials dropped)."""
    s = dna_seq.upper()
    n = len(s) - (len(s) % 3)
    return [s[i:i + 3] for i in range(0, n, 3)]
