"""
cli.py
AlleleSelect command-line interface.

Usage:
    alleleselect --variant c.575G>A --transcript ENST00000360228.10 --output R192Q_candidates/
    alleleselect --variant c.575G>A --output R192Q_candidates/ --no-blast --no-rnafold

Run 'alleleselect --help' for full options.
"""

import argparse
import os
import sys


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="alleleselect",
        description=(
            "AlleleSelect: Allele-Selective ASO Design Pipeline for CACNA1A GoF Mutations.\n"
            "Xiu Lab | thexiulab.org | github.com/axshoe/alleleselect"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--variant", "-v",
        required=True,
        help="HGVS coding notation for the variant (e.g. c.575G>A for R192Q).",
    )
    parser.add_argument(
        "--transcript", "-t",
        default="ENST00000360228",
        help="Ensembl transcript ID. Default: ENST00000360228 (CACNA1A canonical).",
    )
    parser.add_argument(
        "--output", "-o",
        default="alleleselect_output/",
        help="Output directory for CSV and HTML report. Default: alleleselect_output/",
    )
    parser.add_argument(
        "--aso-lengths",
        nargs="+",
        type=int,
        default=[18, 19, 20, 21, 22],
        help="ASO lengths to generate. Default: 18 19 20 21 22.",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=30,
        help="Nucleotide flank on each side of mutation site for window sliding. Default: 30.",
    )
    parser.add_argument(
        "--concentration-uM",
        type=float,
        default=1.0,
        help="Strand concentration in micromolar for Tm calculation. Default: 1.0.",
    )
    parser.add_argument(
        "--no-blast",
        action="store_true",
        help="Skip BLASTn off-target check (use if BLAST+ or GENCODE FASTA not available).",
    )
    parser.add_argument(
        "--no-rnafold",
        action="store_true",
        help="Skip RNAfold accessibility scoring (use if ViennaRNA not installed).",
    )
    parser.add_argument(
        "--top-n-blast",
        type=int,
        default=50,
        help="Number of top candidates to run BLASTn on. Default: 50.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed progress messages.",
    )
    return parser


def run(args) -> None:
    from alleleselect.sequence.hgvs_parser import parse_hgvs_coding, validate_against_cds
    from alleleselect.sequence.fetcher import fetch_cds_sequence, generate_mutant_cds
    from alleleselect.scoring.allele_selectivity import generate_candidate_windows
    from alleleselect.scoring.accessibility import run_rnafold, compute_window_accessibility, check_rnafold_available
    from alleleselect.scoring.offtarget import run_blast_offtarget
    from alleleselect.scoring.splice import flag_splice_risk, get_splice_positions_for_r192q
    from alleleselect.modification.annotator import annotate_all_candidates
    from alleleselect.output.report import save_csv, save_html_report

    os.makedirs(args.output, exist_ok=True)
    log = print if args.verbose else lambda *a, **k: None

    # 1. Parse HGVS variant
    print(f"[AlleleSelect] Parsing variant: {args.variant}")
    parsed = parse_hgvs_coding(args.variant)
    log(f"  Position: {parsed['position']}, Ref: {parsed['ref']}, Alt: {parsed['alt']}")

    # 2. Fetch CDS sequence from Ensembl
    print(f"[AlleleSelect] Fetching CDS from Ensembl ({args.transcript})...")
    wt_cds = fetch_cds_sequence(args.transcript)
    validate_against_cds(parsed, wt_cds)
    mut_cds = generate_mutant_cds(wt_cds, parsed)
    log(f"  CDS length: {len(wt_cds)} nt")
    log(f"  Context: ...{wt_cds[parsed['position']-4:parsed['position']+3]}... (WT)")
    log(f"  Context: ...{mut_cds[parsed['position']-4:parsed['position']+3]}... (Mut)")

    # 3. Generate candidate windows
    print(f"[AlleleSelect] Generating ASO candidates (lengths: {args.aso_lengths}, flank: {args.flank})...")
    candidates = generate_candidate_windows(
        wt_cds, mut_cds,
        mutation_pos=parsed["position"],
        aso_lengths=args.aso_lengths,
        flank=args.flank,
    )
    print(f"  {len(candidates)} candidates generated.")

    # 4. mRNA accessibility scoring
    if not args.no_rnafold:
        print("[AlleleSelect] Running RNAfold accessibility scoring...")
        import tempfile as _tempfile
        import os as _os
        from alleleselect.sequence.fetcher import extract_window
        window_seq, window_start = extract_window(wt_cds, parsed["position"], flank=200)

        # Explicit work_dir so we can diagnose what files RNAfold actually writes
        rnafold_work_dir = _tempfile.mkdtemp(prefix="alleleselect_rnafold_diag_")
        print(f"  [diag] RNAfold work_dir: {rnafold_work_dir}")
        print(f"  [diag] window_seq length: {len(window_seq)} nt")

        rnafold_result = run_rnafold(window_seq, work_dir=rnafold_work_dir)

        try:
            files_in_workdir = _os.listdir(rnafold_work_dir)
        except Exception:
            files_in_workdir = ["(could not list)"]
        print(f"  [diag] files in work_dir after RNAfold: {files_in_workdir}")
        print(f"  [diag] bp_matrix size: {len(rnafold_result['bp_matrix'])} pairs")
        if rnafold_result['per_base_unpaired']:
            sample = rnafold_result['per_base_unpaired'][:5]
            print(f"  [diag] first 5 accessibility values: {[round(v,3) for v in sample]}")

        for c in candidates:
            acc = compute_window_accessibility(
                rnafold_result["per_base_unpaired"],
                c["mRNA_start"],
                c["mRNA_end"],
                cds_start_in_sequence=1 - window_start,
            )
            c["accessibility_score"] = acc
        log(f"  RNAfold complete. Mean accessibility: {sum(c['accessibility_score'] for c in candidates)/len(candidates):.3f}")
    else:
        print("[AlleleSelect] Skipping RNAfold (--no-rnafold). Setting accessibility = 0.5.")
        for c in candidates:
            c["accessibility_score"] = 0.5

    # 5. Off-target assessment
    if not args.no_blast:
        print(f"[AlleleSelect] Running BLASTn off-target check (top {args.top_n_blast} candidates)...")
        candidates = run_blast_offtarget(candidates, top_n=args.top_n_blast)
    else:
        print("[AlleleSelect] Skipping BLASTn (--no-blast). Setting off_target_count = -1.")
        for c in candidates:
            c["off_target_count"] = -1
            c["off_target_genes"] = []

    # 6. Splice site flagging
    print("[AlleleSelect] Flagging splice site proximity...")
    splice_positions = get_splice_positions_for_r192q()
    candidates = flag_splice_risk(candidates, splice_positions)

    # 7. Gapmer modification annotation
    print("[AlleleSelect] Annotating gapmer modification patterns...")
    candidates = annotate_all_candidates(candidates)

    # 8. Re-rank after all scoring
    # Priority: top_candidate AND accessibility > 0.65 AND off_target = 0
    def priority_key(c):
        asr = c.get("allele_selectivity_ratio", 0)
        acc = c.get("accessibility_score", 0)
        ot_penalty = 100 if c.get("off_target_count", -1) > 0 else 0
        sr_penalty = 50 if c.get("splice_risk") == "Y" else 0
        return asr - (acc * 0.5) + ot_penalty + sr_penalty

    candidates.sort(key=priority_key)

    # 9. Save outputs
    variant_label = f"CACNA1A {args.variant}"
    csv_path = os.path.join(args.output, "candidates.csv")
    html_path = os.path.join(args.output, "report.html")

    save_csv(candidates, csv_path)
    save_html_report(candidates, variant_label, html_path)

    # Print top 5 summary
    print(f"\n[AlleleSelect] Complete. Outputs in: {args.output}/")
    print(f"  Top 5 candidates:")
    for i, c in enumerate(candidates[:5], 1):
        print(
            f"  {i}. {c.get('ASO_ID','?')} | {c.get('aso_seq','')} | "
            f"ASR={c.get('allele_selectivity_ratio',0):.3f} | "
            f"Access={c.get('accessibility_score',0):.2f} | "
            f"OT={c.get('off_target_count',-1)} | "
            f"Splice={c.get('splice_risk','?')}"
        )


def main():
    parser = build_parser()
    args = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(f"[AlleleSelect] Error: {e}", file=sys.stderr)
        if "--verbose" in sys.argv or "-v" in sys.argv:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()