# Changelog

All notable changes to this project will be documented in this file.
This project follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] — 2026-04-24

### Added
- **Forbidden-motif avoidance.** Optimized sequences are automatically
  scanned for restriction sites and long homopolymers; offending codons are
  re-sampled synonymously until the motif is gone. Default blacklist covers
  BsaI, BsmBI, SapI, EcoRI, BamHI, HindIII, XhoI, NdeI, NcoI, plus
  polyA/T/G/C runs of 7 or more. Users may override via
  `forbidden_motifs` in the config YAML.
- **QC-aware oversampling.** `optimization.oversample_factor` (default 3)
  generates `n_variants * oversample_factor` candidates; the pipeline
  keeps the cleanest `n_variants` ranked by (no-forbidden-motifs,
  composite score).
- **Per-variant reproducibility.** Each candidate uses a deterministic
  child seed derived from `(--seed, index, salt)`, so a single run is
  bit-reproducible and two variants never collide on the same RNG stream.
- **Source genome support.** `--source-genome <GenBank>` builds a
  full source CUT for harmonized mode instead of estimating from the
  single input gene.
- Forbidden-motif columns in the report TSV (`forbidden_motif_count`,
  `forbidden_motifs`).
- 11 new tests covering motif detection, repair, and reproducibility.

### Changed
- `optimization_report.tsv` schema gained `forbidden_motif_count` and
  `forbidden_motifs` columns (after `has_strong_structure`).

## [0.1.0] — 2026-04-24

### Added
- Initial release. Six modules (`input_handler`, `reference_builder`,
  `sequence_optimizer`, `quality_control`, `metrics`, `report_generator`)
  plus `pipeline.py` and the `run_optimizer.py` CLI.
- Three optimization modes: `max_cai`, `weighted`, `harmonized`.
- Biophysical QC: GC content, sliding-window GC, homopolymer detection,
  direct-repeat detection, 5′ mRNA structure (ViennaRNA primary,
  palindrome fallback).
- Output bundle: per-sequence FASTA with variants, combined
  recommended-only FASTA, metrics TSV, custom CUT TSV, codon comparison
  TSV (DNA input), and a markdown summary.
- 94 tests across 6 module-level suites plus end-to-end + CLI integration.
