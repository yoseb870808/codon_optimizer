"""End-to-end integration tests for the codon optimizer pipeline."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from optimizer.pipeline import run_pipeline
from optimizer.utils import translate


# Use relaxed config so the tiny mock genome can build a reference set.
def _write_test_config(tmp_path: Path) -> Path:
    cfg = tmp_path / "test_config.yaml"
    cfg.write_text(
        "reference:\n"
        "  min_cds_length: 500\n"
        "  top_expressed_pct: 30\n"
        "  min_reference_genes: 3\n"
        "  fallback_min_genes: 3\n"
    )
    return cfg


class TestFullPipeline:
    def test_rnaseq_mode(
        self, mock_genome_path, mock_expression_path, mock_protein_fasta_path, tmp_path
    ):
        cfg = _write_test_config(tmp_path)
        out = tmp_path / "out"
        result = run_pipeline(
            genome_path=mock_genome_path,
            target_path=mock_protein_fasta_path,
            rnaseq_path=mock_expression_path,
            mode="weighted",
            n_variants=3,
            output_dir=out,
            config_path=cfg,
            seed=7,
        )
        assert result["reference_info"]["mode"] == "rnaseq"
        assert (out / "combined_recommended.fasta").exists()
        assert (out / "optimization_report.tsv").exists()
        assert (out / "codon_usage_table.tsv").exists()
        assert (out / "summary.md").exists()

        df = pd.read_csv(out / "optimization_report.tsv", sep="\t")
        assert len(df) == 2  # two protein sequences

    def test_housekeeping_mode(
        self, mock_genome_path, mock_protein_fasta_path, tmp_path
    ):
        cfg = _write_test_config(tmp_path)
        out = tmp_path / "hk_out"
        result = run_pipeline(
            genome_path=mock_genome_path,
            target_path=mock_protein_fasta_path,
            rnaseq_path=None,
            mode="weighted",
            n_variants=3,
            output_dir=out,
            config_path=cfg,
            seed=7,
        )
        assert result["reference_info"]["mode"] in ("housekeeping", "all_cds_fallback")
        assert (out / "summary.md").exists()

    def test_dna_input_creates_comparison(
        self, mock_genome_path, mock_expression_path, mock_dna_fasta_path, tmp_path
    ):
        cfg = _write_test_config(tmp_path)
        out = tmp_path / "dna_out"
        run_pipeline(
            genome_path=mock_genome_path,
            target_path=mock_dna_fasta_path,
            rnaseq_path=mock_expression_path,
            mode="weighted",
            n_variants=3,
            output_dir=out,
            config_path=cfg,
            seed=7,
        )
        comparison_files = list(out.glob("*_codon_comparison.tsv"))
        assert len(comparison_files) == 2

    def test_max_cai_mode(
        self, mock_genome_path, mock_expression_path, mock_protein_fasta_path, tmp_path
    ):
        cfg = _write_test_config(tmp_path)
        out = tmp_path / "maxcai_out"
        result = run_pipeline(
            genome_path=mock_genome_path,
            target_path=mock_protein_fasta_path,
            rnaseq_path=mock_expression_path,
            mode="max_cai",
            n_variants=3,
            output_dir=out,
            config_path=cfg,
            seed=7,
        )
        for res in result["results"]:
            # max_cai always yields exactly 1 variant
            assert len(res["variants"]) == 1
            assert pytest.approx(res["variants"][0]["cai_score"], abs=1e-9) == 1.0

    def test_harmonized_mode(
        self, mock_genome_path, mock_expression_path, mock_dna_fasta_path, tmp_path
    ):
        cfg = _write_test_config(tmp_path)
        out = tmp_path / "harm_out"
        result = run_pipeline(
            genome_path=mock_genome_path,
            target_path=mock_dna_fasta_path,
            rnaseq_path=mock_expression_path,
            mode="harmonized",
            n_variants=3,
            output_dir=out,
            config_path=cfg,
            seed=7,
        )
        for res in result["results"]:
            rec = next(v for v in res["variants"] if v["is_recommended"])
            body = rec["dna_sequence"][:-3]
            assert translate(body) == res["protein_seq"]

    def test_reproducibility(
        self, mock_genome_path, mock_expression_path, mock_protein_fasta_path, tmp_path
    ):
        cfg = _write_test_config(tmp_path)
        out1 = tmp_path / "repro1"
        out2 = tmp_path / "repro2"
        r1 = run_pipeline(
            genome_path=mock_genome_path,
            target_path=mock_protein_fasta_path,
            rnaseq_path=mock_expression_path,
            mode="weighted",
            n_variants=3,
            output_dir=out1,
            config_path=cfg,
            seed=42,
        )
        r2 = run_pipeline(
            genome_path=mock_genome_path,
            target_path=mock_protein_fasta_path,
            rnaseq_path=mock_expression_path,
            mode="weighted",
            n_variants=3,
            output_dir=out2,
            config_path=cfg,
            seed=42,
        )
        s1 = [v["dna_sequence"] for res in r1["results"] for v in res["variants"]]
        s2 = [v["dna_sequence"] for res in r2["results"] for v in res["variants"]]
        assert s1 == s2


class TestCLI:
    def test_help(self, project_root):
        result = subprocess.run(
            [sys.executable, str(project_root / "run_optimizer.py"), "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "usage" in result.stdout.lower()
        assert "--genome" in result.stdout

    def test_missing_genome(self, project_root):
        result = subprocess.run(
            [sys.executable, str(project_root / "run_optimizer.py"), "--target", "x.fasta"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0

    def test_missing_target(self, project_root):
        result = subprocess.run(
            [sys.executable, str(project_root / "run_optimizer.py"), "--genome", "x.gbk"],
            capture_output=True, text=True,
        )
        assert result.returncode != 0

    def test_end_to_end_cli(
        self,
        project_root,
        mock_genome_path,
        mock_expression_path,
        mock_protein_fasta_path,
        tmp_path,
    ):
        cfg = _write_test_config(tmp_path)
        out = tmp_path / "cli_out"
        result = subprocess.run(
            [
                sys.executable, str(project_root / "run_optimizer.py"),
                "--genome", str(mock_genome_path),
                "--target", str(mock_protein_fasta_path),
                "--rnaseq", str(mock_expression_path),
                "--mode", "weighted",
                "--variants", "3",
                "--output", str(out),
                "--config", str(cfg),
                "--seed", "7",
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert (out / "summary.md").exists()
        assert (out / "optimization_report.tsv").exists()
