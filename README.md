# Codon Optimizer for Non-Model Organisms

A generic, organism-agnostic codon optimization tool. Give it a host
genome (GenBank) and — optionally — RNA-seq expression data; it builds a
custom Codon Usage Table (CUT), reverse-translates one or more target
proteins or genes into host-optimized DNA, scrubs forbidden motifs,
runs biophysical QC, and writes a structured report.

**No hardcoded organisms. No fake gene names.** Everything is driven
by the input files you provide.

---

## Highlights

- **GenBank → CUT pipeline.** Parses CDS features, handles the complement
  strand and multi-contig records, skips pseudogenes and malformed
  features, and builds the reference translatome either from the
  top-expressed genes (RNA-seq mode) or annotated housekeeping genes
  (no-RNA-seq mode).
- **Three optimization modes.**
  - `max_cai` — deterministic best-codon-per-AA (CAI = 1.0).
  - `weighted` *(default)* — stochastic sampling proportional to the
    host relative adaptiveness `w_i` (Sharp & Li 1987), producing
    diverse variants.
  - `harmonized` — preserves the source organism's per-codon rarity
    rank (pair with `--source-genome` for a proper source CUT).
- **Forbidden-motif scrubbing.** Restriction sites (BsaI, BsmBI, SapI,
  EcoRI, BamHI, HindIII, XhoI, NdeI, NcoI) and synthesis-unfriendly
  homopolymer runs (≥7 of A/T/G/C) are rewritten synonymously after
  optimization. Protein translation is preserved by construction.
- **QC-aware oversampling.** The engine generates
  `n_variants × oversample_factor` candidates and returns the cleanest
  `n_variants` ranked by `(no-forbidden-motifs, composite_score)`.
- **Biophysical QC.** Overall + sliding-window GC, homopolymer runs,
  direct repeats, and 5′ mRNA secondary structure (ViennaRNA primary,
  palindrome fallback).
- **Reproducible.** `--seed N` makes every variant bit-reproducible via
  per-variant child seeds.
- **Structured outputs** suitable for wet-lab hand-off: per-sequence
  multi-variant FASTA, combined recommended-only FASTA, metrics report
  TSV, custom CUT TSV, per-position codon comparison TSV, and a
  markdown run summary.

---

## Installation

```bash
git clone <repo-url> codon_optimizer
cd codon_optimizer
python -m venv .venv && source .venv/bin/activate   # or .venv\Scripts\activate on Windows
pip install -r requirements.txt

# Optional — ViennaRNA enables MFE-based mRNA structure QC:
conda install -c bioconda viennarna
# or: pip install ViennaRNA
```

The tool runs without ViennaRNA and falls back to a palindrome-based
hairpin heuristic.

Verify the install:

```bash
python -m pytest tests/ -v
python run_optimizer.py --help
```

---

## Quick start

```bash
python run_optimizer.py \
    --genome  host_genome.gbk \
    --target  my_proteins.fasta \
    --rnaseq  rnaseq_expression.xlsx \
    --mode    weighted \
    --variants 3 \
    --seed    42 \
    --output  output/
```

**Minimum (no expression data):**

```bash
python run_optimizer.py --genome host.gbk --target protein.fasta
```

**Harmonized with a source genome:**

```bash
python run_optimizer.py \
    --genome        host.gbk \
    --source-genome ecoli_K12.gbk \
    --target        my_gene_from_ecoli.fasta \
    --mode          harmonized
```

---

## CLI reference

| Flag | Required | Description |
|---|---|---|
| `--genome` | ✓ | Host genome GenBank file (`.gb`/`.gbk`). Multi-record supported. |
| `--target` | ✓ | Multi-FASTA of target proteins or DNA sequences. Auto-detected per record. |
| `--rnaseq` |  | Excel `.xlsx` with columns `locus_tag, RPKM_rep1, RPKM_rep2, RPKM_rep3, RPKM_average`. |
| `--source-genome` |  | Source organism GenBank — builds a full source CUT for harmonized mode. |
| `--mode` |  | `max_cai`, `weighted`, or `harmonized` (default: from config, usually `weighted`). |
| `--variants` |  | Number of variants to return per input sequence. |
| `--output` |  | Output directory (created if missing). Default `./output/`. |
| `--config` |  | User YAML that merges on top of `config_default.yaml`. |
| `--seed` |  | Random seed for reproducibility. |
| `--verbose` / `-v` |  | Debug-level logging. |

---

## Inputs

### Host GenBank (`--genome`)

Standard NCBI annotation. Every CDS must carry a `locus_tag`. Features
without a `locus_tag`, or whose nucleotide length is not divisible by
three, or that contain an internal stop codon (pseudogenes) are logged
as warnings and skipped. Joins and complement-strand features are
handled.

### RNA-seq expression (`--rnaseq`)

Excel `.xlsx` with the exact columns, in this order:

```
locus_tag | RPKM_rep1 | RPKM_rep2 | RPKM_rep3 | RPKM_average
```

Locus tags must match the GenBank CDS tags. Match rate < 80 % triggers
a warning; < 10 % is fatal (almost certainly a file mismatch).

Need to reshape an RNA-seq file that doesn't match this schema? See
`scripts/prepare_sr7_expression.py` for an example that renames columns
and coerces comma-quoted numbers to floats.

### Target sequences (`--target`)

Multi-FASTA. Each record is auto-detected as protein or DNA:

- **Protein** — contains only standard amino acids (`ACDEFGHIKLMNPQRSTVWY`)
  plus optional trailing `*`. Ambiguous residues (`B`, `J`, `O`, `U`, `X`,
  `Z`) are rejected.
- **DNA** — contains only `ATCGN`. Must be divisible by three. Must not
  translate to any `X`. Internal stop codons are rejected.

DNA input enables per-position codon comparison and harmonized mode.

### Config YAML (`--config`)

Your YAML merges recursively over
[`config_default.yaml`](config_default.yaml). Override only the keys you
care about; defaults fill the rest.

**Example `my_config.yaml`:**

```yaml
optimization:
  mode: weighted
  n_variants: 5
  oversample_factor: 5
  repair_forbidden_motifs: true

forbidden_motifs:
  - { name: BsaI,  pattern: GGTCTC }
  - { name: BsmBI, pattern: CGTCTC }
  - { name: SapI,  pattern: GCTCTTC }
  - { name: my_stem_loop, pattern: "GC{5,}GC", regex: true }

qc:
  max_homopolymer: 6
  min_repeat_length: 12
```

---

## Outputs

Everything goes into `--output`. Per-run layout:

```
output/
├── combined_recommended.fasta       # One recommended variant per input
├── <gene>_optimized.fasta           # All variants for <gene>, recommended first
├── <gene>_codon_comparison.tsv      # Per-position diff vs. input (DNA input only)
├── codon_usage_table.tsv            # 64-row host CUT + comment header
├── optimization_report.tsv          # One row per input: CAI, GC, QC, motif flags
└── summary.md                       # Human-readable markdown run summary
```

`optimization_report.tsv` columns:

```
sequence_name, input_type, length_aa, n_variants,
original_cai, recommended_cai, cai_improvement,
original_gc, recommended_gc, host_gc,
homopolymer_count, repeat_count, mfe_5prime, has_strong_structure,
forbidden_motif_count, forbidden_motifs,
composite_score, optimization_mode, qc_pass
```

---

## Optimization modes

### `max_cai`

Always picks the codon with the highest `w_i` per amino acid (ties
broken alphabetically). CAI of the output = 1.0 against the host CUT.
Produces a single deterministic variant. Useful when you want the
theoretical ceiling, but the resulting sequences can be repetitive;
forbidden-motif repair may then swap individual codons.

### `weighted` *(default)*

At each position, sample a codon from the amino acid's synonyms with
probabilities proportional to `w_i`. Each requested variant uses an
independent RNG stream derived from `--seed`, so `--seed 42 --variants 3`
always gives the same three sequences across runs.

### `harmonized`

Requires DNA input. For each source codon, rank its position among
synonymous peers in the source CUT (high rank = commonly used in the
source); pick the host codon at the same rank. This preserves
rhythm-dependent effects such as co-translational folding pauses. With
probability `harmonized_noise` (default 0.1) the engine samples from
the weighted host distribution instead, which prevents over-rigidity
in repetitive regions.

When `--source-genome` is omitted, the source CUT is estimated from the
single input gene and a warning is logged — OK in a pinch but noisy for
small genes.

---

## Architecture

```
optimizer/
├── input_handler.py        Parse GenBank, Excel, FASTA, config
├── reference_builder.py    Reference translatome + CUT construction
├── sequence_optimizer.py   Three modes + oversampling + ranking
├── sequence_repair.py      Forbidden-motif detection + local repair
├── quality_control.py      GC, homopolymer, repeat, mRNA structure
├── metrics.py              CAI, Nc, codon comparison, composite score
├── report_generator.py     FASTA, TSV, markdown writers
├── pipeline.py             End-to-end orchestration
└── utils.py                Genetic code, translate, GC, RC

tests/                      105 tests across 8 suites
scripts/                    One-off helpers (e.g. RNA-seq reshape)
```

---

## Testing

```bash
python -m pytest tests/ -v
```

105 tests covering every module plus end-to-end and CLI integration
cases. Mock GenBank / Excel / FASTA fixtures are generated on the fly
into `tests/fixtures/`.

---

## Worked example: *Priestia megaterium* SR7

A real-world run is shipped under `input/`:

- `Pmegaterium_SR7_tss.gbk` — annotated *B. megaterium* SR7 genome (3 replicons, 5,531 valid CDS)
- `SR7_anaerobic_RNA.xlsx` — anaerobic RNA-seq triplicates

First reshape the expression file to the expected schema:

```bash
python scripts/prepare_sr7_expression.py
```

Then build the SR7 CUT standalone (optional, for inspection):

```bash
python scripts/build_sr7_cut.py
# → sr7_cut/codon_usage_table.tsv + reference_gene_list.tsv
```

Or just run the full pipeline on a FASTA of targets
(e.g. `input/MVA.fasta`, 7 MVA-pathway enzymes from *E. coli* and *S. cerevisiae*):

```bash
python run_optimizer.py \
    --genome   input/Pmegaterium_SR7_tss.gbk \
    --rnaseq   input/SR7_anaerobic_RNA_normalized.xlsx \
    --target   input/MVA.fasta \
    --variants 3 --seed 42 \
    --output   output_MVA/
```

Representative results on the 7 MVA enzymes:

| Gene        | Length (AA) | Original → Recommended CAI | GC (%) | QC |
|---|---|---|---|---|
| atoB (E. coli) | 394  | 0.58 → 0.72 | 44.9 | ✓ |
| ERG8         | 451  | 0.66 → 0.73 | 40.3 | ✓ |
| ERG12        | 443  | 0.60 → 0.71 | 38.5 | ✓ |
| ERG13        | 491  | 0.63 → 0.72 | 39.0 | ✓ |
| ERG19        | 396  | 0.60 → 0.73 | 40.1 | ✓ |
| HMG1         | 1054 | 0.63 → 0.72 | 39.7 | ✓ |
| IDI1         | 288  | 0.67 → 0.74 | 37.3 | ⚠ (one 12-bp repeat) |

---

## Citation

If you use this tool in published work, please cite:

- Sharp, P.M. & Li, W.-H. (1987). The codon adaptation index — a measure
  of directional synonymous codon usage bias. *Nucleic Acids Research*
  **15**(3), 1281–1295.
- Wright, F. (1990). The 'effective number of codons' used in a gene.
  *Gene* **87**(1), 23–29.
- Lorenz, R. *et al.* (2011). ViennaRNA Package 2.0.
  *Algorithms for Molecular Biology* **6**, 26.

---

## License

MIT — see [`LICENSE`](LICENSE).
