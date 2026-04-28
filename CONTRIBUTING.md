# Contributing

Thanks for your interest! A few lightweight conventions keep the repo tidy.

## Setup

```bash
git clone <repo-url> codon_optimizer
cd codon_optimizer
python -m venv .venv
source .venv/bin/activate            # .venv\Scripts\activate on Windows
pip install -r requirements.txt
pip install -e .
```

## Workflow

1. Open an issue describing the bug or the feature.
2. Branch off `main` (`fix/<slug>` or `feat/<slug>`).
3. Write the test first (mandatory — see below).
4. Implement until tests pass.
5. Open a PR. CI runs the full `pytest` suite.

## Testing

```bash
python -m pytest tests/ -v
```

- Every public function in `optimizer/` has tests (see the module-level
  test files in `tests/`).
- Add your test to the appropriate `tests/test_<module>.py` or create a
  new one if the scope is new.
- End-to-end behavior belongs in `tests/test_integration.py`.
- Real-data regression tests belong in `tests/test_real_genome_regression.py`
  (auto-skipped when inputs are absent).
- Never delete a test to "make it pass." If a test is wrong, fix the
  test and document why in the PR description.

## Code style

- Python 3.10+ with type hints on every public signature.
- Google-style docstrings on public functions.
- `logging`, not `print`, for warnings and info.
- `black` compatible (88-char lines). Imports: stdlib → third-party → local.
- Custom exceptions go in the module that raises them.
- No global state — pass dependencies as arguments.

## Scope guardrails

- No hardcoded organism-specific logic. Everything must work for an
  arbitrary GenBank.
- Do not invent biological data in tests. Test fixtures are synthetic
  (e.g., `GENE_001`, `MOCK001`).
- ViennaRNA is a mandatory runtime dependency. The pipeline imports
  ``RNA`` at startup and raises a clear ``ImportError`` if missing — do
  not reintroduce a try/except fallback.
- Input validation fails fast with a clear, actionable message.

## Commit messages

Conventional-Commits-style tags help the changelog:

```
feat: add forbidden-motif repair
fix: handle CDS with zero length
docs: expand harmonized-mode section
test: cover oversampling reproducibility
```
