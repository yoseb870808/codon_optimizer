# Changelog

All notable changes to this project will be documented in this file.
This project follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.2] — 2026-04-28

### Fixed
- **CI red on Python 3.10 / 3.11**:
  ``optimizer/report_generator.py`` had ``strftime("…")`` calls nested
  inside f-strings using the same quote character. PEP 701 lifted that
  restriction in 3.12, so 3.10 and 3.11 raised ``SyntaxError`` at module
  import — taking down 6 of 12 GitHub Actions jobs. Fixed by extracting
  the timestamp into a local before each f-string.
- **Non-deterministic seeds across processes.**
  ``sequence_optimizer._derive_seed`` mixed in ``hash(salt)``, which
  Python randomizes per process under the default
  ``PYTHONHASHSEED=random``. Replaced with ``zlib.crc32`` so child
  seeds are stable across runs, processes, and platforms.
- **``simulated_anneal`` crash on minimal CDS.**
  When ``mutable_positions`` was empty (one-codon CDS, or
  ``window // 3 == 1``), ``rng.choice([])`` raised ``ValueError``.
  Now returns the input unchanged with ``attempts=0``.
- **``KeyError`` in motif repair when motif dict omits ``"name"``.**
  ``_compile_motifs`` already fell back to the pattern; the filter path
  in ``repair_sequence`` did not. Extracted ``_motif_label`` helper and
  used it in both places.

### Docs
- ``CLAUDE.md``, ``CONTRIBUTING.md``, ``architecture.md``, and
  ``test_plan.md`` no longer mention a ViennaRNA fallback. The runtime
  has required ViennaRNA since 0.3.0; the docs now match.

### Tests
- 6 new regression tests in ``tests/test_regressions.py``, including a
  cross-process subprocess test that pins ``_derive_seed`` stability
  against any future reintroduction of ``hash()``.

## [0.4.1] — 2026-04-25

### Added
- **RNA-seq parser accepts more file formats and column names.**
  - File extensions: `.xlsx`, `.xls`, `.csv`, `.tsv`, `.txt` (tab-separated).
    Previously `.xlsx` only.
  - Column auto-rename for common alternates (`RPKM_WT1` → `RPKM_rep1`,
    `gene_id` → `locus_tag`, `mean_rpkm` → `RPKM_average`, etc.).
    Full alias map in `EXPRESSION_COLUMN_ALIASES`.
  - Thousands separators (`"1,503.54"`) are stripped automatically.
  - Helpful error message that lists all recognized aliases when a
    required column can't be resolved.
- 6 new tests covering txt/csv parsing, column aliasing, and comma
  stripping (`test_input_handler.py`).

### Fixed
- `_coerce_numeric` now handles pandas 3.x `StringDtype` columns, not
  just legacy `object` dtype.
- `test_different_seeds_differ` now uses a longer protein so the
  reproducibility check is statistically robust.

## [0.4.0] — 2026-04-25

### Safety
- **Protein-identity invariant enforced at three independent layers**:
  (1) immediately after every per-variant repair call, (2) immediately
  after `optimize()` returns, (3) a final end-of-pipeline sweep over
  every pre- and post-repair variant before any FASTA is written. Any
  non-synonymous mutation aborts the pipeline with
  ``RuntimeError: PROTEIN INTEGRITY VIOLATION`` instead of silently
  emitting a corrupted sequence.
- **Post-repair staleness fixed.** After mRNA repair the pipeline
  refreshes the dot-bracket structure, recomputes the composite score
  with the new MFE, re-checks forbidden motifs (repair could introduce
  new restriction sites), and renumbers `variant_number` so the report
  TSV and FASTA agree on which variant is which.
- 6 new tests in ``tests/test_protein_integrity.py`` including a
  deliberate-corruption test that confirms the assertion fires when
  fed a non-synonymous mutation.

### Added
- **Approach A — targeted 5' structure repair.** New
  ``optimizer/mrna_repair.py`` parses the ViennaRNA dot-bracket structure,
  identifies codons that participate in 5'-window stems, and synonymously
  rewrites them in CAI-preferred order until MFE crosses
  ``mrna_mfe_threshold`` or the attempt cap is reached.
- **Approach B — N-terminal codon-temperature flattening.** New config
  knobs ``optimization.n_terminal_temperature`` (default 1.5) and
  ``optimization.n_terminal_codons`` (default 30). The first N codons are
  sampled with a flatter probability distribution, preserving the natural
  ribosome ramp (Tuller et al. 2010) and dramatically reducing the rate
  at which strong 5' structure forms in the first place.
- **Approach E — simulated-annealing 5' structure repair.** Metropolis–
  Hastings search over synonymous-codon space with energy
  ``E = max(0, threshold - MFE) + λ · (1 - CAI)``. Slower than A but
  escapes local minima.
- **Dispatcher.** ``mrna_repair.mode = off | targeted | annealing | both``.
  ``both`` runs A first, falls back to E only if A fails.
- **Dual outputs.** Pipeline writes a ``pre_repair/`` subdirectory with
  the un-repaired variants alongside the final outputs, plus a
  side-by-side ``mrna_repair_comparison.tsv`` showing
  ``pre_mfe_5prime / post_mfe_5prime / delta_mfe / delta_cai`` per gene.
  Disable with ``--no-pre-repair-outputs``.
- New CLI flags ``--mrna-repair {off,targeted,annealing,both}`` and
  ``--no-pre-repair-outputs``.
- 8 new tests in ``tests/test_mrna_repair.py`` (real ``RNA.fold`` calls).

### Changed
- Variant ranking still uses ``(no-forbidden, no-strong-structure, score)``
  but is re-applied after the repair pass — a candidate previously
  rejected for its structure may now win.
- ``optimization_report.tsv`` columns unchanged; the new
  ``mrna_repair_comparison.tsv`` carries the repair-specific data.

## [0.3.0] — 2026-04-25

### Changed (BREAKING)
- **ViennaRNA is now mandatory.** The palindrome-scan fallback has been
  removed. Importing `optimizer.quality_control` raises `ImportError`
  with install instructions if `RNA` (the ViennaRNA Python binding) is
  missing. Install with `pip install ViennaRNA` (Windows / Linux / macOS
  wheels available for Python 3.10–3.13) or `conda install -c bioconda
  viennarna`.
- Variant ranking now uses the actual `RNA.fold()` MFE: the optimizer
  ranks candidates by `(no-forbidden-motifs, no-strong-structure,
  composite_score)`. The composite score now consumes the MFE term
  (previously a no-op when ViennaRNA was unavailable).
- `mfe_5prime` and `has_strong_structure` columns in
  `optimization_report.tsv` now always carry real values.

### Added
- Per-variant fields `mfe_5prime`, `mrna_structure` (dot-bracket), and
  `has_strong_structure` on every variant dict returned by `optimize()`.
- Tests now exercise real `RNA.fold` calls on poly-A, GC-rich
  hairpins, and threshold sensitivity.

### Removed
- `_palindrome_scan` and the `VIENNARNA_AVAILABLE` flag in
  `quality_control.py`.

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
