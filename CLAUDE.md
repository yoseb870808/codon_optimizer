# CLAUDE.md — Codon Optimizer Development Harness

## Project Identity
**Name:** Codon Optimizer for Non-Model Organisms
**Type:** Standalone Python CLI tool (pip-installable, GitHub-distributable)
**Location:** `D:\Dropbox\2020_MIT\B_megaterium_SR7\ALE\_SR7_ALE_computation_analysis\Sec3_literature_search\codon_optimizer\`

## What This Tool Does
Takes a host organism's GenBank genome file + optional RNA-seq expression data → builds a custom Codon Usage Table → reverse-translates target protein or DNA sequences into host-optimized DNA using weighted stochastic codon selection → runs biophysical quality control → outputs optimized FASTA + metrics report.

**No organism-specific logic. No hardcoded gene names. Fully generic.**

## Read Order
1. This file (`CLAUDE.md`) — development rules and workflow
2. `architecture.md` — full technical specification (modules, functions, algorithms, I/O)
3. `test_plan.md` — verification criteria per module

---

## Development Workflow: Plan → Work → Review → Gate

You build this project **one module at a time**. For each module:

### PLAN
- Read the module spec in `architecture.md`
- Identify inputs, outputs, dependencies on other modules
- Identify edge cases from `test_plan.md`

### WORK
- Implement the module
- Write tests for the module (in `tests/`)
- Write docstrings for every public function

### REVIEW
- Run the module's tests: `python -m pytest tests/test_<module>.py -v`
- Verify edge cases pass
- Verify no regressions: `python -m pytest tests/ -v`
- Check: does the implementation match `architecture.md`?

### GATE
- **All tests pass** → proceed to next module
- **Any test fails** → fix before moving on. Do NOT skip. Do NOT comment out failing tests.

---

## Build Sequence

Determine the optimal implementation order yourself based on the dependency graph in `architecture.md`. The logical dependency chain is:

```
input_handler (no dependencies)
    → reference_builder (depends on input_handler)
        → metrics (depends on reference_builder for CUT format)
        → sequence_optimizer (depends on reference_builder + metrics)
            → quality_control (depends on sequence_optimizer output format)
                → report_generator (depends on all above)
                    → run_optimizer.py CLI (integrates everything)
```

You may reorder within this chain if you find a better approach. You may refactor module boundaries if the architecture benefits. Document any deviations from `architecture.md` in comments.

After all modules are built and tested individually, run a **full integration test** that exercises the complete pipeline end-to-end using the example files in `examples/`.

---

## Guardrails — Hard Rules

### NEVER do these:
- **Never invent biological data.** No fake gene names, no made-up titers, no fabricated accession numbers. Test fixtures use clearly synthetic data (e.g., `ATGATGATGATG` repeats, `test_gene_001`).
- **Never hardcode organism-specific logic.** No "if Bacillus then X". The tool works for ANY organism given a GenBank file.
- **Never skip tests.** Every module has tests. Tests must pass before moving to the next module.
- **Never silently swallow errors.** Parse failures, missing files, mismatched locus_tags — all must produce clear error messages with actionable guidance.
- **Never assume input format.** Validate GenBank structure, Excel columns, FASTA format. Fail fast with helpful messages.

### ALWAYS do these:
- **Write tests first** for edge cases listed in `test_plan.md`, then implement to pass them.
- **Handle the ViennaRNA optional dependency** gracefully. The tool must work fully without it (fallback to simple hairpin detection). Use try/except import.
- **Use type hints** on all function signatures.
- **Write docstrings** for all public functions (Google style).
- **Log warnings** for recoverable issues (few locus_tag matches, short reference set, etc.) using Python `logging` module.
- **Validate inputs early.** Check file existence, format, and content before processing.

---

## Code Quality Standards

- **Python 3.10+** compatibility
- **Type hints** on all public function parameters and return types
- **Docstrings** on all public functions (Google style: Args, Returns, Raises)
- **Logging** via `logging` module (not print statements) for warnings and info
- **Error handling** with custom exceptions where appropriate
- **No global state** — all functions receive their dependencies as arguments
- `black` formatting compatible (88 char line length)
- Imports organized: stdlib → third-party → local

---

## Testing Standards

- Use `pytest` as test framework
- Test fixtures in `tests/fixtures/` (mock GenBank, mock Excel, mock FASTA)
- Each module has its own test file: `tests/test_<module>.py`
- Integration test: `tests/test_integration.py` (full pipeline, runs last)
- **Minimum coverage targets:**
  - input_handler: all input formats + all error cases
  - reference_builder: both modes (RNA-seq and housekeeping fallback)
  - sequence_optimizer: all 3 modes (max_cai, weighted, harmonized)
  - quality_control: all QC checks including ViennaRNA fallback
  - metrics: CAI calculation verified against known values
  - report_generator: output file existence and format validation
- Create test fixtures during the build (mock .gbk, .xlsx, .fasta files) — do NOT rely on real SR7 data for tests

---

## Dependencies

See `requirements.txt` for the full list. Key constraints:
- `biopython` for GenBank/FASTA parsing
- `pandas` + `openpyxl` for Excel I/O
- `numpy` for stochastic sampling
- `ViennaRNA` for mRNA structure (optional, graceful fallback)
- `pyyaml` for config
- `pytest` for testing

---

## Output Expectations

When the full build is complete, this directory should contain:

```
codon_optimizer/
├── CLAUDE.md                    # This file (read-only reference)
├── architecture.md              # Technical spec (read-only reference)
├── test_plan.md                 # Test spec (read-only reference)
├── README.md                    # User documentation (YOU write this)
├── requirements.txt             # Dependencies
├── config_default.yaml          # Default configuration
├── pyproject.toml               # Package metadata (pip-installable)
├── run_optimizer.py             # CLI entry point
├── optimizer/
│   ├── __init__.py
│   ├── input_handler.py
│   ├── reference_builder.py
│   ├── sequence_optimizer.py
│   ├── quality_control.py
│   ├── metrics.py
│   ├── report_generator.py
│   └── utils.py
├── tests/
│   ├── conftest.py              # Shared fixtures
│   ├── fixtures/                # Mock input files
│   │   ├── mock_genome.gbk
│   │   ├── mock_expression.xlsx
│   │   ├── mock_protein.fasta
│   │   └── mock_dna.fasta
│   ├── test_input_handler.py
│   ├── test_reference_builder.py
│   ├── test_optimizer.py
│   ├── test_qc.py
│   ├── test_metrics.py
│   ├── test_report.py
│   └── test_integration.py
├── examples/
│   ├── example_protein.fasta
│   └── example_dna.fasta
├── data/                        # Empty (populated at runtime)
├── output/                      # Empty (populated at runtime)
└── .gitignore
```

Every `.py` file must be importable without errors. Every test must pass. `python run_optimizer.py --help` must print usage information. The README must explain installation and usage for a scientist who has never used the tool before.

---

## Final Checklist (run before declaring done)

```
[ ] python -m pytest tests/ -v                    → ALL PASS
[ ] python run_optimizer.py --help                 → prints usage
[ ] python run_optimizer.py --genome examples/... --target examples/... → produces output/
[ ] output/ contains: FASTA, report TSV, CUT TSV, comparison TSV, summary MD
[ ] README.md exists with installation + usage instructions
[ ] No hardcoded organism names anywhere in optimizer/
[ ] ViennaRNA import failure handled gracefully (tool still works)
[ ] All public functions have type hints + docstrings
```
