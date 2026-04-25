#!/usr/bin/env python3
"""Command-line entry point for the codon optimizer."""

from __future__ import annotations

import argparse
import logging
import sys

from optimizer.pipeline import run_pipeline


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Codon Optimizer for Non-Model Organisms",
    )
    parser.add_argument("--genome", required=True, help="Host genome GenBank file (.gb/.gbk)")
    parser.add_argument("--target", required=True, help="Sequences to optimize (FASTA: protein or DNA)")
    parser.add_argument("--rnaseq", default=None, help="RNA-seq expression data (Excel .xlsx)")
    parser.add_argument(
        "--source-genome",
        default=None,
        help="Source organism GenBank (for harmonized mode — builds a full source CUT)",
    )
    parser.add_argument(
        "--mode",
        choices=["max_cai", "weighted", "harmonized"],
        default=None,
        help="Optimization mode (default: from config — 'weighted')",
    )
    parser.add_argument("--variants", type=int, default=None, help="Number of variants to generate")
    parser.add_argument("--output", default=None, help="Output directory (default: ./output/)")
    parser.add_argument("--config", default=None, help="Custom config YAML (merged over defaults)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose logging",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    try:
        result = run_pipeline(
            genome_path=args.genome,
            target_path=args.target,
            rnaseq_path=args.rnaseq,
            mode=args.mode,
            n_variants=args.variants,
            output_dir=args.output,
            config_path=args.config,
            seed=args.seed,
            source_genome_path=args.source_genome,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    paths = result["output_paths"]
    print("\nOptimization complete.")
    print(f"  Reference mode:  {result['reference_info']['mode']}")
    print(f"  Reference genes: {result['reference_info']['n_genes']}")
    print(f"  Host GC:         {result['reference_info']['host_gc'] * 100:.1f}%")
    print(f"  Sequences:       {len(result['results'])}")
    print()
    print("Output files:")
    for key, value in paths.items():
        if isinstance(value, list):
            for v in value:
                print(f"  {key}: {v}")
        else:
            print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
