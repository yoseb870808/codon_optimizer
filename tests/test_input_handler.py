"""Tests for optimizer.input_handler."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import pytest

from optimizer import input_handler
from optimizer.input_handler import (
    load_config,
    parse_expression_data,
    parse_genbank,
    parse_target_sequences,
)


class TestParseGenBank:
    def test_parse_returns_valid_cds(self, mock_genome_path, mock_genome_meta):
        cds = parse_genbank(mock_genome_path)
        # 10 GENE_ + 3 housekeeping + 1 complement = 14 valid
        assert len(cds) >= 10
        for tag in mock_genome_meta["locus_tags"]:
            assert tag in cds

    def test_parse_skips_no_locus_tag(self, mock_genome_path, caplog):
        with caplog.at_level(logging.WARNING):
            cds = parse_genbank(mock_genome_path)
        # the CDS without locus_tag should not appear
        assert not any(
            v.get("location", "").startswith("[")
            and v.get("locus_tag") is None
            for v in cds.values()
        )
        assert any("without locus_tag" in msg for msg in caplog.messages)

    def test_parse_skips_bad_length(self, mock_genome_path, mock_genome_meta, caplog):
        with caplog.at_level(logging.WARNING):
            cds = parse_genbank(mock_genome_path)
        assert mock_genome_meta["bad_length_tag"] not in cds
        assert any(
            "not divisible by 3" in msg for msg in caplog.messages
        )

    def test_parse_skips_pseudogene(self, mock_genome_path, mock_genome_meta, caplog):
        with caplog.at_level(logging.WARNING):
            cds = parse_genbank(mock_genome_path)
        assert mock_genome_meta["pseudogene_tag"] not in cds
        assert any("internal stop" in msg for msg in caplog.messages)

    def test_parse_handles_complement(self, mock_genome_path, mock_genome_meta):
        cds = parse_genbank(mock_genome_path)
        tag = mock_genome_meta["complement_tag"]
        assert tag in cds
        # extracted sense strand should match what we originally encoded
        assert cds[tag]["nucleotide_seq"] == mock_genome_meta["complement_seq"]

    def test_parse_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            parse_genbank("nonexistent_file.gbk")

    def test_cds_entries_have_all_fields(self, mock_genome_path):
        cds = parse_genbank(mock_genome_path)
        entry = next(iter(cds.values()))
        for field in (
            "locus_tag", "gene", "product", "nucleotide_seq",
            "amino_acid_seq", "length_nt", "length_aa", "location", "contig",
        ):
            assert field in entry

    def test_parse_detects_housekeeping_products(self, mock_genome_path, mock_genome_meta):
        cds = parse_genbank(mock_genome_path)
        for tag in mock_genome_meta["housekeeping_tags"]:
            assert cds[tag]["product"] is not None


class TestParseExpression:
    def test_parse_valid_excel(self, mock_expression_path):
        df = parse_expression_data(mock_expression_path)
        assert list(df.columns)[0] == "locus_tag"
        assert list(df.columns)[-1] == "RPKM_average"
        assert "RPKM_rep1" in df.columns
        # sorted descending
        assert df["RPKM_average"].iloc[0] >= df["RPKM_average"].iloc[-1]

    def test_missing_required_columns(self, tmp_path):
        bad = tmp_path / "bad.xlsx"
        pd.DataFrame({"locus_tag": ["X"], "other": [1]}).to_excel(bad, index=False)
        with pytest.raises(ValueError, match="missing required columns"):
            parse_expression_data(bad)

    def test_missing_file(self):
        with pytest.raises(FileNotFoundError):
            parse_expression_data("nope.xlsx")

    def test_negative_values_rejected(self, tmp_path):
        bad = tmp_path / "neg.xlsx"
        df = pd.DataFrame({
            "locus_tag": ["X"],
            "RPKM_rep1": [-1.0],
            "RPKM_rep2": [1.0],
            "RPKM_rep3": [1.0],
            "RPKM_average": [0.33],
        })
        df.to_excel(bad, index=False)
        with pytest.raises(ValueError, match="negative"):
            parse_expression_data(bad)

    def test_tab_separated_txt(self, tmp_path):
        """Parser accepts .txt with tab-separated values."""
        path = tmp_path / "expr.txt"
        path.write_text(
            "locus_tag\tRPKM_rep1\tRPKM_rep2\tRPKM_rep3\tRPKM_average\n"
            "GENE_001\t100.5\t110.3\t90.7\t100.5\n"
            "GENE_002\t1500.0\t1600.0\t1400.0\t1500.0\n"
        )
        df = parse_expression_data(path)
        assert len(df) == 2
        assert df["RPKM_average"].iloc[0] == 1500.0   # sorted desc

    def test_csv_format(self, tmp_path):
        path = tmp_path / "expr.csv"
        path.write_text(
            "locus_tag,RPKM_rep1,RPKM_rep2,RPKM_rep3,RPKM_average\n"
            "GENE_A,5.0,5.0,5.0,5.0\n"
        )
        df = parse_expression_data(path)
        assert len(df) == 1

    def test_auto_renames_wt_columns(self, tmp_path):
        """WT-style replicate column names are auto-renamed to rep1/2/3."""
        path = tmp_path / "wt.txt"
        path.write_text(
            "locus_tag\tRPKM_WT1\tRPKM_WT2\tRPKM_WT3\tRPKM_average\n"
            "GENE_001\t10\t20\t30\t20\n"
        )
        df = parse_expression_data(path)
        assert "RPKM_rep1" in df.columns
        assert "RPKM_WT1" not in df.columns

    def test_auto_renames_other_aliases(self, tmp_path):
        """Other recognized aliases (gene_id, replicate_1, mean_rpkm) are renamed."""
        path = tmp_path / "alt.txt"
        path.write_text(
            "gene_id\treplicate_1\treplicate_2\treplicate_3\tmean_rpkm\n"
            "G1\t1\t2\t3\t2\n"
        )
        df = parse_expression_data(path)
        for col in ("locus_tag", "RPKM_rep1", "RPKM_rep2", "RPKM_rep3", "RPKM_average"):
            assert col in df.columns

    def test_strips_thousands_separators(self, tmp_path):
        """Comma-quoted numbers like '1,503.54' are coerced to floats."""
        path = tmp_path / "commas.csv"
        path.write_text(
            "locus_tag,RPKM_rep1,RPKM_rep2,RPKM_rep3,RPKM_average\n"
            'G1,"1,503.54","63,035.45","14,340.54","26,293.18"\n'
        )
        df = parse_expression_data(path)
        assert df["RPKM_average"].iloc[0] == pytest.approx(26293.18, abs=0.1)
        assert df["RPKM_rep2"].iloc[0] == pytest.approx(63035.45, abs=0.1)

    def test_helpful_error_lists_aliases(self, tmp_path):
        """Missing-column error message lists known aliases."""
        path = tmp_path / "bad.csv"
        path.write_text("foo,bar\n1,2\n")
        with pytest.raises(ValueError, match="Recognized aliases"):
            parse_expression_data(path)


class TestParseTargets:
    def test_parse_protein_fasta(self, mock_protein_fasta_path):
        targets = parse_target_sequences(mock_protein_fasta_path)
        assert len(targets) == 2
        for t in targets:
            assert t["input_type"] == "protein"
            assert t["original_dna"] is None
            assert "*" not in t["protein_seq"]

    def test_parse_dna_fasta(self, mock_dna_fasta_path):
        targets = parse_target_sequences(mock_dna_fasta_path)
        assert len(targets) == 2
        for t in targets:
            assert t["input_type"] == "dna"
            assert t["original_dna"] is not None
            assert len(t["original_dna"]) % 3 == 0

    def test_empty_fasta_raises(self, tmp_path):
        empty = tmp_path / "empty.fasta"
        empty.write_text("")
        with pytest.raises(ValueError):
            parse_target_sequences(empty)

    def test_ambiguous_aa_rejected(self, tmp_path):
        bad = tmp_path / "ambig.fasta"
        bad.write_text(">bad\nMKXFZBV\n")
        with pytest.raises(ValueError, match="ambiguous"):
            parse_target_sequences(bad)

    def test_missing_file(self):
        with pytest.raises(FileNotFoundError):
            parse_target_sequences("missing.fasta")

    def test_dna_not_div_by_3(self, tmp_path):
        bad = tmp_path / "short.fasta"
        bad.write_text(">bad\nATGAAAT\n")
        with pytest.raises(ValueError, match="not divisible by 3"):
            parse_target_sequences(bad)


class TestLoadConfig:
    def test_default_has_expected_keys(self):
        cfg = load_config()
        for section in ("reference", "optimization", "variant_scoring", "qc", "output"):
            assert section in cfg
        assert cfg["optimization"]["mode"] in ("max_cai", "weighted", "harmonized")

    def test_config_merging(self, tmp_path):
        user = tmp_path / "custom.yaml"
        user.write_text("optimization:\n  n_variants: 7\n")
        cfg = load_config(user)
        assert cfg["optimization"]["n_variants"] == 7
        # other defaults preserved
        assert "cai_weight" in cfg["variant_scoring"]

    def test_missing_user_config_raises(self):
        with pytest.raises(FileNotFoundError):
            load_config("no_such_config.yaml")
