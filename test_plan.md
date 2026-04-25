# Test Plan — Codon Optimizer

## General Rules
- Use `pytest` as framework
- Test fixtures in `tests/fixtures/` (create mock files during test setup or store as static files)
- Each module has its own test file
- Run full suite after every module: `python -m pytest tests/ -v`
- Integration test (`test_integration.py`) runs last and exercises the full pipeline

---

## Test Fixtures to Create

### mock_genome.gbk
A minimal GenBank file with ~20 CDS features on a single 10,000 bp contig. Include:
- 10 CDS with locus_tags `GENE_001` through `GENE_010`, various lengths (200-2000 nt)
- 3 CDS with known product annotations matching housekeeping keywords ("30S ribosomal protein S1", "elongation factor Tu", "chaperonin GroEL")
- 1 CDS with length not divisible by 3 (should be skipped)
- 1 CDS with an internal stop codon (pseudogene, should be skipped)
- 1 CDS on the complement strand
- 1 CDS without a locus_tag (should be skipped)
- All sequences should use real-ish codons (not random) so CUT calculations produce meaningful w_i values
- Use synthetic gene names and sequences — do NOT use real organism data

### mock_expression.xlsx
Excel file with columns: locus_tag, RPKM_rep1, RPKM_rep2, RPKM_rep3, RPKM_average
- Include `GENE_001` through `GENE_010` (matching mock_genome.gbk)
- Include 2 locus_tags NOT in the GenBank file (to test mismatch handling)
- Expression values ranging from 0.5 to 5000 RPKM
- Make `GENE_001`, `GENE_002`, `GENE_003` the highest expressed

### mock_protein.fasta
```
>test_protein_short A small test protein
MKFLILVAS
>test_protein_medium A medium length protein
MKFLILVASAGTTLMITGNPKLRDEWYQHCFS
```

### mock_dna.fasta
```
>test_dna_ecoli_style An E.coli-like CDS (GC-rich codons)
ATGAAATTCCTGATCCTGGTGGCCTCCGCCGGCACCACCCTGATGATCACCGGCAACCCGAAACTGCGCGATGAATGGTACCAGCACTGCTTCTCC
>test_dna_at_rich An AT-rich CDS
ATGAAATTTTTAATTTTAATTGCATCAGCAGGTACAACATTAATGATTACAGGTAATCCAAAATTAAGAGATGAATGGTATAAA
```

---

## Module Test Specifications

### test_input_handler.py

| Test | Input | Expected |
|------|-------|----------|
| Parse valid GenBank | mock_genome.gbk | Dict with 10+ valid CDS entries |
| Skip CDS without locus_tag | mock_genome.gbk | Excluded from results, warning logged |
| Skip CDS not divisible by 3 | mock_genome.gbk | Excluded, warning logged |
| Skip pseudogene (internal stop) | mock_genome.gbk | Excluded, warning logged |
| Handle complement strand CDS | mock_genome.gbk | Correct reverse-complement sequence extracted |
| FileNotFoundError for missing GenBank | "nonexistent.gbk" | Raises FileNotFoundError |
| Parse valid Excel | mock_expression.xlsx | DataFrame with correct columns, sorted by RPKM_average desc |
| Missing required columns | Excel without RPKM_average | Raises ValueError |
| Auto-detect protein FASTA | mock_protein.fasta | input_type = "protein" |
| Auto-detect DNA FASTA | mock_dna.fasta | input_type = "dna", protein_seq populated from translation |
| Multi-sequence FASTA | mock_protein.fasta | Returns list with 2 entries |
| Empty FASTA | empty file | Raises ValueError |
| Config loading | config_default.yaml | Returns dict with all expected keys |
| Config merging | custom.yaml overriding one key | Merged value overrides default |

### test_reference_builder.py

| Test | Scenario | Expected |
|------|----------|----------|
| Mode A: RNA-seq reference set | mock CDS + mock expression | Top 10% by RPKM, all ≥ 500 nt |
| Mode A: expansion to 20% | If top 10% < 30 genes | Expands, warning logged |
| Mode A: match rate logging | 2 mismatched locus_tags | Warning logged, match rate < 100% |
| Mode B: housekeeping detection | mock CDS with annotated products | Finds ribosomal + EF-Tu + GroEL genes |
| Mode B: fallback to all CDS | < 20 housekeeping genes | Uses all CDS ≥ 500 nt, warning logged |
| CUT calculation: codon counts | Known sequence "ATGATGATG" (3x Met) | ATG count = 3 |
| CUT calculation: w_i for Leu | Reference with known Leu codon distribution | w_i correct per Sharp & Li |
| CUT calculation: w_i = 1.0 for Met | Any reference set | ATG w_i = 1.0 |
| CUT calculation: zero-count codon | Codon absent from reference | w_i = 0.01 (floor) |
| CUT shape | Any reference set | 64 rows, columns: codon, amino_acid, count, frequency, w_i, rscu |
| Host GC calculation | Known sequences | Correct GC fraction |

### test_optimizer.py

| Test | Scenario | Expected |
|------|----------|----------|
| max_cai: deterministic | Same input twice | Identical output |
| max_cai: always highest w_i | Short protein, known CUT | Every codon has w_i = 1.0 in CUT |
| max_cai: CAI = 1.0 | Any input | CAI of output = 1.0 against same CUT |
| weighted: stochastic | Same input, different seeds | Different sequences |
| weighted: 3 variants | n_variants=3 | Returns list of length 3 |
| weighted: high-w_i preference | Many runs with known CUT | High-w_i codons appear more frequently |
| weighted: all valid codons | Any input | Output translates to same protein as input |
| weighted: reproducibility | Same seed | Same output |
| harmonized: rarity preservation | DNA input with known rare codons | Corresponding rare codons in host CUT selected |
| harmonized: requires DNA | Protein input only | Raises ValueError |
| harmonized: no source CUT | DNA input, source_cut=None | Estimates from input gene, warning logged |
| variant scoring | 3 variants with known CAI/GC/MFE | Correct composite scores, correct ranking |
| recommended flag | 3 variants | Exactly one has is_recommended=True (highest score) |
| protein output integrity | Any optimization | All variants translate to the same protein as input |

### test_qc.py

| Test | Scenario | Expected |
|------|----------|----------|
| GC overall | "ATATATAT" (0% GC) | overall_gc = 0.0, within_range = False |
| GC overall | "GCGCGCGC" (100% GC) | overall_gc = 1.0, within_range = False |
| GC overall | "ATGCATGC" (50% GC) | overall_gc = 0.5, within_range = True |
| GC sliding window | Sequence with local GC spike | flagged_regions contains the spike region |
| Homopolymer detection | "ATGAAAAAAAATGC" (8x A) | One entry: position, base="A", length=8 |
| Homopolymer none | "ATGATGATGATG" | Empty list |
| Repeat detection | "ATGATGATGATG" with 12+ bp repeat downstream | Detected |
| Repeat none | Non-repetitive sequence | Empty list |
| mRNA structure (ViennaRNA) | GC-rich 5' region | MFE reported, possibly flagged |
| mRNA structure (fallback) | Palindromic 5' sequence | Hairpin detected |
| mRNA structure (fallback) | Non-structured 5' | No hairpin |
| run_all_qc | Clean sequence | pass = True |
| run_all_qc | Sequence with homopolymer | pass = False, homopolymer reported |

### test_metrics.py

| Test | Scenario | Expected |
|------|----------|----------|
| CAI of max_cai sequence | Optimized with max_cai mode | CAI = 1.0 |
| CAI of random sequence | Random codons | CAI < 1.0 (typically 0.1-0.3) |
| CAI excludes Met/Trp | Sequence of all Met | Returns appropriate value (edge case) |
| CAI with zero-count codon | Sequence using codon with w_i=0.01 | CAI reduced but not zero |
| Nc of biased sequence | All-max-CAI sequence | Nc near 20 |
| Nc of uniform sequence | Each synonym used equally | Nc near 61 |
| Codon comparison | Two known sequences | Correct per-position diff, changed flags accurate |
| Composite score calculation | Known CAI, GC, MFE values | Correct weighted sum |

### test_report.py

| Test | Scenario | Expected |
|------|----------|----------|
| FASTA output | 3 variants, one recommended | File exists, 3 sequences, recommended first |
| FASTA headers | Generated output | Headers contain CAI, GC%, score |
| Combined FASTA | 2 input sequences | File exists, 2 recommended sequences |
| Report TSV | Full pipeline output | File exists, correct columns, correct row count |
| CUT TSV | Generated CUT | File exists, 64 rows, comment headers present |
| Comparison TSV | DNA input | File exists, correct columns |
| Summary MD | Full pipeline | File exists, contains run parameters |
| Output directory creation | Non-existent output dir | Created automatically |

### test_integration.py

| Test | Scenario | Expected |
|------|----------|----------|
| Full pipeline with RNA-seq | mock_genome + mock_expression + mock_protein | All output files generated, no errors |
| Full pipeline without RNA-seq | mock_genome + mock_protein (no Excel) | Falls back to housekeeping, all outputs generated |
| Full pipeline DNA input | mock_genome + mock_dna | Comparison TSV generated (has original codons) |
| Full pipeline multi-FASTA | mock_genome + FASTA with 2 sequences | Per-sequence outputs + combined FASTA |
| Full pipeline max_cai mode | --mode max_cai | 1 variant per sequence, CAI = 1.0 |
| Full pipeline harmonized | mock_genome + mock_dna --mode harmonized | Runs without error, rarity profile preserved |
| CLI --help | run_optimizer.py --help | Prints usage, exit code 0 |
| CLI missing --genome | run_optimizer.py --target x.fasta | Error message, exit code != 0 |
| CLI missing --target | run_optimizer.py --genome x.gbk | Error message, exit code != 0 |
| Reproducibility | Same inputs + --seed 42 twice | Identical outputs |

---

## Quality Gate Summary

| Gate | When | Criterion |
|------|------|-----------|
| Gate 1 | After input_handler | All test_input_handler tests pass |
| Gate 2 | After reference_builder | Gates 1 + all test_reference_builder pass |
| Gate 3 | After metrics | Gates 1-2 + all test_metrics pass |
| Gate 4 | After sequence_optimizer | Gates 1-3 + all test_optimizer pass |
| Gate 5 | After quality_control | Gates 1-4 + all test_qc pass |
| Gate 6 | After report_generator | Gates 1-5 + all test_report pass |
| Gate 7 | After CLI + integration | ALL tests pass, including test_integration |
