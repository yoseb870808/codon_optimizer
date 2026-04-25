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
- **5′ mRNA structure repair.** Strong stem-loops in the first 150 nt
  reduce ribosome loading; the optimizer detects them with ViennaRNA
  (real ΔG calculation) and removes them through three combinable
  strategies — *Approach A* (targeted, dot-bracket-driven codon repair),
  *Approach B* (5′-region temperature flattening to preserve the natural
  ribosome ramp), and *Approach E* (Metropolis–Hastings simulated
  annealing). Both pre- and post-repair outputs are written for
  side-by-side comparison.
- **QC-aware oversampling.** The engine generates
  `n_variants × oversample_factor` candidates and returns the cleanest
  `n_variants` ranked by `(no-forbidden-motifs, no-strong-structure, composite_score)`.
- **Biophysical QC.** Overall + sliding-window GC, homopolymer runs,
  direct repeats, and 5′ mRNA MFE via ViennaRNA `RNA.fold()` (mandatory
  dependency, real thermodynamic calculation — no fallback heuristic).
- **Reproducible.** `--seed N` makes every variant bit-reproducible via
  per-variant child seeds.
- **Protein-identity invariant** asserted at three independent layers
  (post-optimize, post-repair, end-of-pipeline sweep). Any non-synonymous
  mutation aborts the pipeline rather than silently emitting a bad FASTA.
- **Structured outputs** suitable for wet-lab hand-off: per-sequence
  multi-variant FASTA, combined recommended-only FASTA, metrics report
  TSV, custom CUT TSV, per-position codon comparison TSV, mRNA-repair
  side-by-side TSV, and a markdown run summary.

---

## Installation

```bash
git clone https://github.com/yoseb870808/codon_optimizer.git
cd codon_optimizer
python -m venv .venv && source .venv/bin/activate   # or .venv\Scripts\activate on Windows
pip install -r requirements.txt
```

`requirements.txt` already pins ViennaRNA. The pip wheel works on
Windows / Linux / macOS for Python 3.10–3.13. If the wheel is missing
for your platform, fall back to:

```bash
conda install -c bioconda viennarna
```

ViennaRNA is **mandatory** — the pipeline imports `RNA` at startup and
fails fast with install instructions if it's missing.

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
| `--mrna-repair` |  | `off`, `targeted`, `annealing`, or `both` (default: `both`). Controls how strong 5′ structure is removed. |
| `--no-pre-repair-outputs` |  | Skip writing the `pre_repair/` subdir + `mrna_repair_comparison.tsv`. Default is to write both. |
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

Any of `.xlsx`, `.xls`, `.csv`, `.tsv`, or `.txt` (tab-separated).
Required columns:

```
locus_tag | RPKM_rep1 | RPKM_rep2 | RPKM_rep3 | RPKM_average
```

The parser auto-renames common alternates so most real-world files
work without preprocessing:

| Canonical | Aliases recognised (case-insensitive) |
|---|---|
| `locus_tag` | `locus`, `gene_id`, `geneid`, `gene`, `id`, `tag` |
| `RPKM_rep1/2/3` | `RPKM_WT1/2/3`, `RPKM_1/2/3`, `WT1/2/3`, `replicate_1/2/3`, `rep1/2/3` |
| `RPKM_average` | `RPKM_avg`, `mean_rpkm`, `average_rpkm`, `mean`, `average` |

Numeric columns may contain thousands separators (`"1,503.54"`); they
are stripped automatically. Rows with unparseable values are dropped
with a logged warning.

Locus tags must match the GenBank CDS tags. Match rate < 80 % triggers
a warning; < 10 % is fatal (almost certainly a file mismatch).

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
├── combined_recommended.fasta       # Final recommended variant per input (post-repair)
├── <gene>_optimized.fasta           # All variants for <gene>, recommended first
├── <gene>_codon_comparison.tsv      # Per-position diff vs. input (DNA input only)
├── codon_usage_table.tsv            # 64-row host CUT + comment header
├── optimization_report.tsv          # One row per input: CAI, GC, MFE, QC, motif flags
├── mrna_repair_comparison.tsv       # Side-by-side: pre vs post repair MFE/CAI per input
├── summary.md                       # Human-readable markdown run summary
└── pre_repair/                      # Same files as above, computed BEFORE the repair pass
    ├── combined_recommended.fasta
    ├── <gene>_optimized.fasta
    └── optimization_report.tsv
```

The `pre_repair/` subdir lets you see what the optimizer would have
recommended without the 5′ structure repair step. Use
`--no-pre-repair-outputs` to skip writing it, and use
`--mrna-repair off` to disable the repair pass entirely.

`optimization_report.tsv` columns:

```
sequence_name, input_type, length_aa, n_variants,
original_cai, recommended_cai, cai_improvement,
original_gc, recommended_gc, host_gc,
homopolymer_count, repeat_count, mfe_5prime, has_strong_structure,
forbidden_motif_count, forbidden_motifs,
composite_score, optimization_mode, qc_pass
```

`mrna_repair_comparison.tsv` columns:

```
sequence_name, repair_mode, repair_method_used, repair_attempts,
repair_succeeded, pre_cai, post_cai, delta_cai,
pre_gc, post_gc, pre_mfe_5prime, post_mfe_5prime, delta_mfe,
pre_strong_structure, post_strong_structure
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

## 5′ mRNA structure repair

Strong secondary structure in the first ~150 nt of an mRNA blocks
ribosome loading and is one of the strongest predictors of low protein
expression in bacteria (Kudla et al. 2009, Goodman et al. 2013). The
optimizer measures the 5′ minimum free energy with ViennaRNA on every
candidate, flags variants where MFE ≤ `qc.mrna_mfe_threshold`
(default −30 kcal/mol), and tries to repair them through a combination
of three strategies — all of which preserve the protein translation by
construction.

| Strategy | What it does | When it kicks in |
|---|---|---|
| **A — targeted repair** | Parses the dot-bracket structure, identifies codons that participate in stems, synonymously rewrites them in CAI-preferred order, refolds, repeats. | After candidate generation, only on flagged variants. |
| **B — 5′-region temperature** | Flattens the per-AA codon distribution for the first N codons (default 30). Less probability mass on the most-biased codons preserves the natural slow-codon "ramp" (Tuller et al. 2010), which prevents most strong 5′ structure from forming in the first place. | During candidate generation in `weighted` mode. |
| **E — simulated annealing** | Metropolis–Hastings search over synonymous-codon space with energy `max(0, threshold − MFE) + λ·(1 − CAI)`. Slower but escapes local minima that block A. | Used by `--mrna-repair both` only when A fails to cross the threshold. |

Configure via `--mrna-repair {off,targeted,annealing,both}` or the
`mrna_repair` block in `config_default.yaml`. Default is `both`, which
is the right choice for most users: A handles 95 % of cases cheaply, E
catches the rest.

The `pre_repair/` subdir + `mrna_repair_comparison.tsv` show exactly
what each strategy did for every gene. Example:

```
sequence_name  repair_method_used  pre_mfe  post_mfe  delta_mfe  delta_cai
ERG8_Scer      targeted            -46.1    -29.4     +16.7      +0.001
atoB_Ecoli     targeted            -29.4    -28.9     +0.5        0.000
ERG12_Scer     skipped_already_good
…
```

### Wet-lab complements

The optimizer fixes **CDS-internal** structure. For maximum protein
yield, pair it with two upstream interventions:

- **RBS Calculator** (Salis Lab) — design the 5′ UTR / Shine-Dalgarno
  for your target translation initiation rate.
- **Self-cleaving ribozyme** (RiboJ; Lou et al. 2012) inserted between
  promoter and 5′ UTR — creates a defined, sequence-clean 5′ end and
  decouples the CDS from upstream context. Doesn't touch CDS-internal
  structure, but eliminates promoter-context-induced variability.

### Protein-identity invariant

Every code path that mutates a candidate sequence does so via
**synonymous substitution only**. This is enforced at three independent
layers:

1. After every per-variant repair call (`pipeline.py` line ~270)
2. After each `optimize()` call (line ~205)
3. End-of-pipeline sweep over every variant in every output bucket
   (line ~360) — runs *just before* any FASTA is written.

A non-synonymous mutation anywhere in the pipeline raises
`RuntimeError: PROTEIN INTEGRITY VIOLATION` and aborts. The
`tests/test_protein_integrity.py` suite includes a deliberate
corruption test that confirms the assertion fires correctly.

---

## Architecture

```
optimizer/
├── input_handler.py        Parse GenBank, Excel, FASTA, config
├── reference_builder.py    Reference translatome + CUT construction
├── sequence_optimizer.py   Three modes + oversampling + ranking + 5'-region temperature
├── sequence_repair.py      Forbidden-motif detection + local repair
├── mrna_repair.py          Targeted + simulated-annealing 5' MFE repair
├── quality_control.py      GC, homopolymer, repeat, mRNA structure (mandatory ViennaRNA)
├── metrics.py              CAI, Nc, codon comparison, composite score
├── report_generator.py     FASTA, TSV, markdown writers
├── pipeline.py             End-to-end orchestration + protein-identity invariants
└── utils.py                Genetic code, translate, GC, RC

tests/                      122 tests across 11 suites
```

---

## Testing

```bash
python -m pytest tests/ -v
```

122 tests covering every module plus end-to-end and CLI integration
cases. Mock GenBank / Excel / FASTA fixtures are generated on the fly
into `tests/fixtures/`. ViennaRNA is a hard dependency, so all tests
require a working `RNA.fold()`.

---

## Worked example

A two-step run on any host genome:

```bash
# 1. Verify your install
python -m pytest tests/ -v
python run_optimizer.py --help

# 2. Run the optimizer on your inputs
python run_optimizer.py \
    --genome   input/host.gbk \
    --rnaseq   input/expression.xlsx \
    --target   input/my_proteins.fasta \
    --variants 3 --seed 42 \
    --output   output/
```

After the run, inspect:

- `output/summary.md` — human-readable run summary
- `output/optimization_report.tsv` — per-input metrics (CAI, GC, MFE, QC flags)
- `output/codon_usage_table.tsv` — the host CUT used for optimization
- `output/combined_recommended.fasta` — recommended variant per input
- `output/<gene>_optimized.fasta` — all variants per input
- `output/<gene>_codon_comparison.tsv` — per-position diff vs. input (DNA input only)

Two example FASTAs are bundled under `examples/` for a quick smoke test:

```bash
python run_optimizer.py \
    --genome  <your_host.gbk> \
    --target  examples/example_protein.fasta \
    --output  output_smoke/
```

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
