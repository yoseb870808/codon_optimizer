# Architecture — Codon Optimizer for Non-Model Organisms

## 1. System Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                        CLI (run_optimizer.py)                       │
│  --genome  --rnaseq(opt)  --target  --mode  --variants  --config   │
└─────────────┬───────────────────────────────────────────────────────┘
              │
              ▼
┌─────────────────────────┐
│  Module 1: input_handler │ → Parse GenBank, Excel, FASTA
│  - parse_genbank()       │   Validate all inputs
│  - parse_expression()    │   Auto-detect protein vs DNA
│  - parse_targets()       │
│  - load_config()         │
└────────────┬────────────┘
             │
             ▼
┌──────────────────────────────┐
│  Module 2: reference_builder  │ → Build reference translatome
│  - build_reference_set()      │   Mode A: RNA-seq top 10%
│  - build_codon_usage_table()  │   Mode B: Housekeeping fallback
│  - compute_host_gc()          │   Calculate w_i (Sharp & Li 1987)
└────────────┬─────────────────┘
             │
             ▼
┌──────────────────────────────────┐
│  Module 3: sequence_optimizer     │ → Generate optimized sequences
│  - optimize_max_cai()             │   3 modes: max_cai, weighted, harmonized
│  - optimize_weighted()            │   Weighted: stochastic sampling ∝ w_i
│  - optimize_harmonized()          │   Harmonized: preserve source rarity
│  - select_best_variants()         │   Pick top 3 by composite score
└────────────┬─────────────────────┘
             │
             ▼
┌────────────────────────────────┐
│  Module 4: quality_control      │ → Biophysical QC
│  - check_gc_content()           │   GC overall + sliding window
│  - check_homopolymers()         │   Runs of identical bases
│  - check_repeats()              │   Direct repeats
│  - check_mrna_structure()       │   5' MFE via ViennaRNA (or fallback)
│  - run_all_qc()                 │
└────────────┬───────────────────┘
             │
             ▼
┌────────────────────────────────┐
│  Module 5: metrics              │ → CAI, Nc, comparison
│  - calculate_cai()              │   Sharp & Li 1987
│  - calculate_nc()               │   Wright 1990
│  - compare_codons()             │   Per-position original vs optimized
│  - composite_score()            │   0.5*CAI + 0.3*GC + 0.2*MFE
└────────────┬───────────────────┘
             │
             ▼
┌────────────────────────────────┐
│  Module 6: report_generator     │ → Write outputs
│  - write_fasta()                │   FASTA with 3 variants
│  - write_report_tsv()           │   Metrics table (Excel-openable)
│  - write_cut_tsv()              │   Codon usage table
│  - write_comparison_tsv()       │   Per-position codon comparison
│  - write_summary_md()           │   Human-readable summary
└─────────────────────────────────┘
```

---

## 2. Input Specifications

### 2.1 GenBank File (required)

Standard NCBI-annotated GenBank format (`.gb` or `.gbk`). The parser extracts CDS features.

**Required qualifiers per CDS:**
- `locus_tag` — unique identifier (primary key for matching with RNA-seq)
- Nucleotide sequence derivable from feature location (start, end, strand, join)

**Optional qualifiers (used when available):**
- `gene` — gene name
- `product` — protein product description (used for housekeeping gene identification)
- `translation` — amino acid sequence (used for validation)

**Validation rules:**
- Reject CDS without `locus_tag` (log warning, skip)
- Reject CDS with nucleotide length not divisible by 3 (log warning, skip)
- Reject CDS with internal stop codons when translated (log warning, skip; these are pseudogenes)
- Accept CDS on either strand (handle complement)
- Accept join features (concatenate exon sequences)
- Accept multi-record GenBank files (multiple contigs/chromosomes/plasmids)

### 2.2 RNA-seq Expression File (optional, Excel .xlsx)

Fixed column structure:

| Column | Type | Description |
|--------|------|-------------|
| `locus_tag` | string | Must match GenBank CDS locus_tag values |
| `RPKM_rep1` | float | Replicate 1 normalized expression |
| `RPKM_rep2` | float | Replicate 2 normalized expression |
| `RPKM_rep3` | float | Replicate 3 normalized expression |
| `RPKM_average` | float | Average across replicates (used for ranking) |

**Validation rules:**
- First column must be named exactly `locus_tag`
- Last column must be named exactly `RPKM_average`
- All RPKM values must be numeric and ≥ 0
- Report match rate between Excel locus_tags and GenBank locus_tags
- Warn if match rate < 80%
- Fail if match rate < 10% (likely wrong file or wrong organism)

### 2.3 Target Sequences (required, FASTA)

Standard FASTA format. One file, one or multiple sequences.

**Auto-detection per sequence:**
- **Protein:** sequence contains only `ACDEFGHIKLMNPQRSTVWY*` (20 standard amino acids + stop)
- **DNA:** sequence contains only `ATCGN` (and optionally lowercase)
- If DNA: translate to protein using standard genetic code, store both original DNA and translation
- Strip trailing `*` from protein sequences before optimization (re-add appropriate stop codon after)

**Validation rules:**
- Reject empty sequences
- Reject sequences with ambiguous amino acids (B, J, O, U, X, Z) — these cannot be unambiguously reverse-translated
- Warn if DNA sequence has `N` bases (translate as X, then fail on X)
- Warn if DNA sequence does not start with ATG

### 2.4 Configuration (optional, YAML)

See `config_default.yaml` for all parameters. User config merges with and overrides defaults.

---

## 3. Module Specifications

### 3.1 Module 1: input_handler.py

```python
"""Input parsing and validation for all file types."""

def parse_genbank(filepath: str) -> dict[str, dict]:
    """
    Parse GenBank file and extract all valid CDS features.

    Args:
        filepath: Path to .gb or .gbk file.

    Returns:
        Dict keyed by locus_tag. Each value is a dict:
        {
            "locus_tag": str,
            "gene": str or None,
            "product": str or None,
            "nucleotide_seq": str,     # coding sequence (sense strand, ATG...stop)
            "amino_acid_seq": str,     # translated sequence
            "length_nt": int,
            "length_aa": int,
            "location": str,           # human-readable location string
            "contig": str,             # source record ID
        }

    Raises:
        FileNotFoundError: If filepath doesn't exist.
        ValueError: If file contains no valid CDS features.
    """

def parse_expression_data(filepath: str) -> pd.DataFrame:
    """
    Parse RNA-seq Excel file with fixed column structure.

    Expected columns: locus_tag, RPKM_rep1, RPKM_rep2, RPKM_rep3, RPKM_average

    Args:
        filepath: Path to .xlsx file.

    Returns:
        DataFrame sorted by RPKM_average descending.

    Raises:
        FileNotFoundError: If filepath doesn't exist.
        ValueError: If required columns are missing.
    """

def parse_target_sequences(filepath: str) -> list[dict]:
    """
    Parse multi-FASTA file. Auto-detect protein vs DNA per sequence.

    Args:
        filepath: Path to .fasta or .fa file.

    Returns:
        List of dicts:
        {
            "header": str,             # FASTA header (without >)
            "input_type": "protein" | "dna",
            "protein_seq": str,        # amino acid sequence (always present)
            "original_dna": str | None # original DNA if input was DNA
        }

    Raises:
        FileNotFoundError: If filepath doesn't exist.
        ValueError: If file contains no valid sequences.
    """

def load_config(filepath: str = None) -> dict:
    """
    Load YAML config and merge with defaults.

    Args:
        filepath: Path to user config YAML. None uses defaults only.

    Returns:
        Merged config dict.
    """
```

### 3.2 Module 2: reference_builder.py

```python
"""Build reference translatome and custom Codon Usage Table."""

# Housekeeping gene keywords for Mode B (no RNA-seq)
HOUSEKEEPING_KEYWORDS = {
    "ribosomal": [
        "ribosomal protein", "30S ribosomal", "50S ribosomal"
    ],
    "translation": [
        "elongation factor", "translation initiation factor"
    ],
    "chaperones": [
        "chaperone", "chaperonin", "GroEL", "GroES", "DnaK"
    ],
    "core_metabolism": [
        "glyceraldehyde-3-phosphate dehydrogenase", "enolase",
        "phosphoglycerate kinase", "pyruvate kinase",
        "fructose-bisphosphate aldolase"
    ],
    "rna_polymerase": [
        "DNA-directed RNA polymerase"
    ]
}

def build_reference_set(
    all_cds: dict[str, dict],
    expression_data: pd.DataFrame | None = None,
    config: dict = None
) -> tuple[list[str], dict]:
    """
    Select reference gene set for CUT construction.

    Mode A (RNA-seq provided):
        1. Match locus_tags between CDS and expression data
        2. Filter: CDS length >= 500 nt
        3. Sort by RPKM_average descending
        4. Take top 10%
        5. If < 30 genes, expand to top 20% with warning

    Mode B (no RNA-seq):
        1. Filter: CDS length >= 500 nt
        2. Match product descriptions against HOUSEKEEPING_KEYWORDS
        3. If < 20 genes found, use ALL CDS >= 500 nt with warning

    Args:
        all_cds: Dict from parse_genbank()
        expression_data: DataFrame from parse_expression_data() or None
        config: Config dict with thresholds

    Returns:
        Tuple of:
        - list of nucleotide sequences (reference translatome)
        - dict with metadata: {
            "mode": "rnaseq" | "housekeeping" | "all_cds_fallback",
            "n_genes": int,
            "n_total_cds": int,
            "match_rate": float (for RNA-seq mode),
            "gene_list": list of locus_tags used
          }
    """

def build_codon_usage_table(reference_sequences: list[str]) -> pd.DataFrame:
    """
    Calculate custom Codon Usage Table from reference translatome.

    For each of the 64 codons, calculate:
        - count: raw occurrence count across all reference sequences
        - amino_acid: which AA this codon encodes
        - frequency: count / total synonymous codons for that AA
        - w_i: relative adaptiveness = count / max_count_for_that_AA (Sharp & Li 1987)
        - rscu: relative synonymous codon usage = (count * n_synonyms) / total_for_AA

    Special cases:
        - ATG (Met) and TGG (Trp): w_i = 1.0 (single codon, no synonyms)
        - Stop codons (TAA, TAG, TGA): include in table, calculate w_i from reference usage
        - Codons with zero count: w_i = 0.01 (floor to avoid log(0) in CAI)

    Args:
        reference_sequences: List of nucleotide CDS sequences.

    Returns:
        DataFrame with columns: [codon, amino_acid, count, frequency, w_i, rscu]
        64 rows (one per codon), sorted by amino_acid then w_i descending.
    """

def compute_host_gc(reference_sequences: list[str]) -> float:
    """
    Calculate overall GC content of the reference translatome.
    Used as target GC% for variant scoring.

    Returns:
        GC fraction (0.0 to 1.0)
    """
```

### 3.3 Module 3: sequence_optimizer.py

```python
"""Codon optimization engine with three modes."""

def optimize(
    protein_seq: str,
    codon_table: pd.DataFrame,
    mode: str = "weighted",
    n_variants: int = 3,
    host_gc: float = None,
    original_dna: str = None,
    source_cut: pd.DataFrame = None,
    seed: int = None
) -> list[dict]:
    """
    Main optimization dispatch. Routes to mode-specific function.

    Args:
        protein_seq: Target protein sequence (no stop codon).
        codon_table: Host CUT from build_codon_usage_table().
        mode: "max_cai" | "weighted" | "harmonized"
        n_variants: Number of variant sequences to generate.
        host_gc: Host GC fraction (for variant scoring).
        original_dna: Original DNA sequence (required for harmonized mode).
        source_cut: Source organism CUT (required for harmonized mode).
        seed: Random seed for reproducibility.

    Returns:
        List of variant dicts, sorted by composite score descending:
        {
            "dna_sequence": str,
            "cai_score": float,
            "gc_content": float,
            "mode": str,
            "variant_number": int,
            "composite_score": float,
            "is_recommended": bool   # True for top variant
        }
    """

def _optimize_max_cai(protein_seq: str, codon_table: pd.DataFrame) -> str:
    """
    Deterministic: always select codon with highest w_i for each AA.
    Single output (n_variants ignored).

    Algorithm:
        For each amino acid in sequence:
            Select codon where w_i == max(w_i for all synonyms of this AA)
            If tie: select first alphabetically (deterministic)
    """

def _optimize_weighted(
    protein_seq: str,
    codon_table: pd.DataFrame,
    n_variants: int,
    rng: np.random.Generator
) -> list[str]:
    """
    Stochastic: sample codons proportional to w_i.

    Algorithm:
        For each variant:
            For each amino acid in sequence:
                Get all synonymous codons {c1, c2, ..., cn} with w_i {w1, w2, ..., wn}
                Normalize: p_j = w_j / sum(w_all)
                Sample one codon from distribution with probabilities p
                Append to growing sequence
            Store complete sequence

    Returns n_variants DNA sequences.
    """

def _optimize_harmonized(
    protein_seq: str,
    original_dna: str,
    codon_table: pd.DataFrame,
    source_cut: pd.DataFrame,
    n_variants: int,
    rng: np.random.Generator
) -> list[str]:
    """
    Preserve codon rarity profile from source organism.

    Algorithm:
        For each position i in protein sequence:
            1. Identify the source codon at position i in original_dna
            2. Calculate the percentile rank of that codon among its synonyms
               in the SOURCE CUT (e.g., this codon is the 3rd most common for Leu)
            3. Select the codon at the same percentile rank in the HOST CUT
               (e.g., pick the 3rd most common Leu codon in the host)
            4. Add stochastic noise: with probability 0.1, sample from weighted
               distribution instead (prevents over-rigidity)

    Requires original_dna (same length as protein_seq * 3).
    Requires source_cut (built from source organism genome, or estimated from
    the input gene's own codon usage if genome unavailable).

    If source_cut is None: estimate source codon usage from the original_dna
    sequence itself (single-gene CUT). Log warning that this is approximate.
    """

def _score_variant(
    dna_seq: str,
    codon_table: pd.DataFrame,
    host_gc: float,
    config: dict
) -> float:
    """
    Composite score for variant ranking.

    score = cai_weight * CAI
          + gc_weight * (1.0 - abs(gc_variant - gc_host) / 0.5)
          + mrna_weight * mrna_stability_score

    Where mrna_stability_score:
        If ViennaRNA available: 1.0 - min(1.0, abs(MFE) / abs(threshold))
        If not: 1.0 (neutral, no penalty)

    Weights from config (default: 0.5, 0.3, 0.2)

    Returns: float (higher = better)
    """
```

### 3.4 Module 4: quality_control.py

```python
"""Biophysical quality control checks on optimized sequences."""

def check_gc_content(
    sequence: str,
    target_gc: float = None,
    window_size: int = 50,
    acceptable_range: tuple = (0.25, 0.65),
    window_deviation: float = 0.15
) -> dict:
    """
    GC content analysis.

    Returns:
    {
        "overall_gc": float,
        "target_gc": float or None,
        "gc_deviation": float (absolute difference from target),
        "within_range": bool,
        "window_min": float,
        "window_max": float,
        "flagged_regions": list of {start, end, gc, direction ("high"/"low")}
    }
    """

def check_homopolymers(sequence: str, max_run: int = 6) -> list[dict]:
    """
    Detect runs of identical nucleotides exceeding max_run.

    Returns: list of {"position": int, "base": str, "length": int}
    Empty list if no violations.
    """

def check_repeats(sequence: str, min_length: int = 12) -> list[dict]:
    """
    Detect direct repeats exceeding min_length.

    Algorithm: sliding window comparison. For each position, search downstream
    for matching subsequences. Report pairs.

    Returns: list of {
        "pos1": int, "pos2": int,
        "sequence": str, "length": int
    }
    """

def check_mrna_structure(
    sequence: str,
    window_5prime: int = 150,
    mfe_threshold: float = -30.0
) -> dict:
    """
    Assess mRNA secondary structure in the 5' region.

    Primary method (ViennaRNA installed):
        Import RNA module from ViennaRNA
        Calculate MFE of first `window_5prime` nucleotides using RNA.fold()
        Report MFE, structure string, and flag if MFE < threshold

    Fallback (ViennaRNA not installed):
        Scan for palindromic sequences (potential hairpin stems) in first 150 nt
        Report count and locations of potential hairpins
        Set has_strong_structure = True if any palindrome > 8 bp found

    Returns:
    {
        "method": "viennarna" | "fallback",
        "mfe": float or None,
        "structure": str or None,
        "has_strong_structure": bool,
        "hairpins": list (fallback only),
        "warning": str or None
    }
    """

def run_all_qc(
    sequence: str,
    target_gc: float = None,
    config: dict = None
) -> dict:
    """
    Run all QC checks and return combined report.

    Returns: {
        "gc": check_gc_content result,
        "homopolymers": check_homopolymers result,
        "repeats": check_repeats result,
        "mrna_structure": check_mrna_structure result,
        "pass": bool (True if no critical flags)
    }
    """
```

### 3.5 Module 5: metrics.py

```python
"""Codon usage metrics and sequence comparison."""

import math

def calculate_cai(sequence: str, codon_table: pd.DataFrame) -> float:
    """
    Codon Adaptation Index (Sharp & Li, 1987).

    CAI = exp( (1/L) * Σ ln(w_i) )  for i = 1 to L codons

    Where:
        L = number of codons (excluding Met, Trp, stop — these have no synonyms)
        w_i = relative adaptiveness from CUT

    If any w_i = 0: use floor value (0.01) to avoid log(0).

    Args:
        sequence: DNA sequence (must be divisible by 3).
        codon_table: CUT DataFrame with 'codon' and 'w_i' columns.

    Returns:
        CAI value (float, 0.0 to 1.0).
    """

def calculate_nc(sequence: str) -> float:
    """
    Effective Number of Codons (Wright, 1990).

    Measures codon usage bias independent of any reference set.
    Range: 20 (maximum bias, one codon per AA) to 61 (no bias, uniform usage).

    Uses the F_cf (homozygosity) method:
        For each amino acid family with n synonymous codons:
            F = (n * Σ(p_i²) - 1) / (n - 1)
        Where p_i = frequency of codon i among synonymous codons in this sequence.

    Nc = 2 + 9/F_2 + 1/F_3 + 5/F_4 + 3/F_6
    (for families of size 2, 3, 4, 6)

    If a family has zero usage in the sequence, use F = 0 (maximum diversity assumption).
    """

def compare_codons(
    original_dna: str,
    optimized_dna: str,
    codon_table: pd.DataFrame
) -> pd.DataFrame:
    """
    Per-position codon comparison between original and optimized sequence.

    Both sequences must encode the same protein (same length, same AA at each position).

    Returns DataFrame:
        [position, amino_acid, original_codon, original_w_i,
         optimized_codon, optimized_w_i, changed]

    'changed' is boolean: True if codon differs between original and optimized.
    """

def composite_score(
    cai: float,
    gc_content: float,
    host_gc: float,
    mfe: float = None,
    mfe_threshold: float = -30.0,
    weights: dict = None
) -> float:
    """
    Combined quality score for variant ranking.

    Default weights: cai=0.5, gc=0.3, mrna=0.2

    Components:
        cai_score = cai (already 0-1)
        gc_score = 1.0 - min(1.0, abs(gc_content - host_gc) / 0.5)
        mrna_score = 1.0 if no ViennaRNA, else 1.0 - min(1.0, abs(mfe) / abs(mfe_threshold))

    Returns: weighted sum (float, higher = better)
    """
```

### 3.6 Module 6: report_generator.py

```python
"""Output file generation."""

def write_optimized_fasta(
    header: str,
    variants: list[dict],
    output_dir: str
) -> str:
    """
    Write FASTA file with all variants for one input sequence.

    Header format:
        >[original_header]_optimized_recommended CAI=0.87 GC=42.1% score=0.82
        ATGATG...
        >[original_header]_variant_2 CAI=0.85 GC=43.2% score=0.79
        ATGATG...
        >[original_header]_variant_3 CAI=0.84 GC=41.8% score=0.78
        ATGATG...

    Recommended variant (highest composite score) listed first.

    Returns: path to written file.
    """

def write_combined_fasta(
    all_results: list[dict],
    output_dir: str
) -> str:
    """
    Write single FASTA with only the recommended variant per input sequence.

    Returns: path to written file.
    """

def write_report_tsv(
    all_results: list[dict],
    output_dir: str
) -> str:
    """
    Write metrics report as TSV (Excel-openable).

    Columns:
        sequence_name, input_type, length_aa, n_variants,
        original_cai, recommended_cai, cai_improvement,
        original_gc, recommended_gc, host_gc,
        homopolymer_count, repeat_count,
        mfe_5prime, has_strong_structure,
        composite_score, optimization_mode

    One row per input sequence (metrics for the recommended variant).

    Returns: path to written file.
    """

def write_cut_tsv(
    codon_table: pd.DataFrame,
    reference_info: dict,
    output_dir: str
) -> str:
    """
    Write the computed Codon Usage Table as TSV.

    Includes header comment lines (starting with #) documenting:
        # Reference mode: rnaseq | housekeeping | all_cds_fallback
        # Reference genes: N
        # Host GC: XX.X%

    Returns: path to written file.
    """

def write_comparison_tsv(
    header: str,
    comparison_df: pd.DataFrame,
    output_dir: str
) -> str:
    """
    Write per-position codon comparison for one sequence.

    Only generated when input was DNA (original codons available).

    Returns: path to written file.
    """

def write_summary_md(
    all_results: list[dict],
    reference_info: dict,
    config: dict,
    output_dir: str
) -> str:
    """
    Human-readable markdown summary of the entire run.

    Sections:
        - Run parameters (genome file, RNA-seq file, mode, n_variants)
        - Reference translatome summary (mode, gene count, GC%)
        - Per-sequence results table
        - QC flag summary
        - Recommendations and warnings

    Returns: path to written file.
    """
```

### 3.7 utils.py

```python
"""Shared utilities."""

CODON_TABLE_STANDARD = {
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

# Reverse mapping: amino acid → list of synonymous codons
AA_TO_CODONS = {}  # Build from CODON_TABLE_STANDARD at module load

def translate(dna_seq: str) -> str:
    """Translate DNA to protein using standard genetic code."""

def is_protein(seq: str) -> bool:
    """Check if sequence is protein (only standard AA characters + *)."""

def is_dna(seq: str) -> bool:
    """Check if sequence is DNA (only ATCGN, case-insensitive)."""

def gc_content(seq: str) -> float:
    """Calculate GC fraction of a nucleotide sequence."""

def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
```

---

## 4. CLI Specification (run_optimizer.py)

```
usage: run_optimizer.py [-h] --genome GENOME --target TARGET
                        [--rnaseq RNASEQ]
                        [--mode {max_cai,weighted,harmonized}]
                        [--variants N] [--output DIR]
                        [--config YAML]
                        [--locus-col COL] [--expr-col COL]
                        [--seed INT]

Codon Optimizer for Non-Model Organisms

required arguments:
  --genome GENOME       Host genome GenBank file (.gb/.gbk)
  --target TARGET       Sequences to optimize (FASTA, protein or DNA)

optional arguments:
  --rnaseq RNASEQ       RNA-seq expression data (Excel .xlsx)
  --mode MODE           Optimization mode (default: weighted)
  --variants N          Number of variants to generate (default: 3)
  --output DIR          Output directory (default: ./output/)
  --config YAML         Custom config file (default: config_default.yaml)
  --seed INT            Random seed for reproducibility
  -h, --help            Show this help message
```

**Pipeline execution order in CLI:**
1. Load config
2. Parse genome → extract CDS
3. If RNA-seq: parse Excel, match locus_tags, filter top 10%
4. Else: filter housekeeping genes from CDS annotations
5. Build CUT from reference set
6. Compute host GC from reference set
7. Parse target FASTA → auto-detect per sequence
8. For each target: optimize → QC → score → select top 3
9. For each target: calculate CAI of original (if DNA input) and optimized
10. Write all output files
11. Print summary to console

---

## 5. Edge Cases and Error Handling

| Situation | Response |
|-----------|----------|
| GenBank file with zero CDS | Raise ValueError with message |
| Excel has locus_tags not found in GenBank | Log each mismatch as warning, continue with matched set |
| Match rate < 80% | Log warning: "Only X% of locus_tags matched" |
| Match rate < 10% | Raise ValueError: "Too few matches, check if files are from the same organism" |
| Top 10% yields < 30 genes (RNA-seq mode) | Expand to top 20%, log warning |
| < 20 housekeeping genes found (no RNA-seq) | Fall back to all CDS ≥ 500 nt, log warning |
| Target protein has unusual length (< 10 AA) | Log warning, optimize anyway |
| Target DNA not divisible by 3 | Raise ValueError with position info |
| Harmonized mode but no original DNA | Raise ValueError: "harmonized mode requires DNA input" |
| ViennaRNA not installed | Log info: "ViennaRNA not found, using simple hairpin detection" |
| Codon with zero count in reference | Set w_i = 0.01 (floor) |
| Output directory doesn't exist | Create it |
| Output file already exists | Overwrite with warning |
