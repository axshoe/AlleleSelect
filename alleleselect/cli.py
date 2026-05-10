"""
cli.py
AlleleSelect command-line interface — v4

New in v4:
  --extra-mismatch      Engineered mismatch mode (van Roon-Mom / Ostergaard 2013 Fig 7)
  --diff-accessibility  Differential WT vs. mutant RNAfold scoring (Aguti & Zhou 2024)
  --recommend-mods      Chemical modification recommendations at SNP-flanking positions
  --gene / --no-splice-check  (v3, unchanged)

All other flags unchanged from v3.

Usage:
    alleleselect --variant c.575G>A --transcript ENST00000360228.10 --output demo/R192Q_output/
    alleleselect --variant c.575G>A --output demo/R192Q_output/ --extra-mismatch --recommend-mods
    alleleselect --variant c.689C>T --transcript ENST00000436367.6 --gene ATXN1 --no-splice-check --output atxn1/
"""

import argparse
import os
import sys


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="alleleselect",
        description=(
            "AlleleSelect v4: Allele-Selective ASO Design Pipeline for Dominant Neurological Mutations.\n"
            "Xiu Lab | thexiulab.org | github.com/axshoe/alleleselect"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--variant", "-v", required=True,
        help="HGVS coding notation (e.g. c.575G>A).")
    parser.add_argument("--transcript", "-t", default="ENST00000360228",
        help="Ensembl transcript ID. Default: ENST00000360228 (CACNA1A canonical).")
    parser.add_argument("--gene", "-g", default=None,
        help="Gene name for output labeling (e.g. CACNA1A, ATXN1).")
    parser.add_argument("--output", "-o", default="alleleselect_output/",
        help="Output directory. Default: alleleselect_output/")
    parser.add_argument("--aso-lengths", nargs="+", type=int, default=[18, 19, 20, 21, 22],
        help="ASO lengths to generate. Default: 18 19 20 21 22.")
    parser.add_argument("--flank", type=int, default=30,
        help="Nucleotide flank around mutation for window sliding. Default: 30.")
    parser.add_argument("--concentration-uM", type=float, default=1.0,
        help="Strand concentration in uM for Tm calculation. Default: 1.0.")
    parser.add_argument("--no-blast", action="store_true",
        help="Skip BLASTn off-target check.")
    parser.add_argument("--no-rnafold", action="store_true",
        help="Skip RNAfold accessibility scoring.")
    parser.add_argument("--no-splice-check", action="store_true",
        help="Skip splice site proximity check. Recommended for non-CACNA1A transcripts.")
    parser.add_argument("--top-n-blast", type=int, default=50,
        help="Number of top candidates to BLAST. Default: 50.")
    # ── v4 new flags ──────────────────────────────────────────────────────────
    parser.add_argument("--diff-accessibility", action="store_true",
        help=(
            "Run RNAfold on BOTH wildtype and mutant sequences and compute "
            "differential accessibility (mutant - wildtype). Candidates where "
            "mutant mRNA is more accessible score higher. Motivated by Aguti & Zhou 2024 "
            "(PMID 38993932). Requires RNAfold. Adds ~60s runtime."
        ))
    parser.add_argument("--extra-mismatch", action="store_true",
        help=(
            "Generate engineered-mismatch variants for top 5 candidates. "
            "Introduces one deliberate additional mismatch against the wildtype allele "
            "2-3 positions from the SNP, giving WT two mismatches and mutant one. "
            "Achieves >100-fold discrimination in Ostergaard 2013 Fig 7. "
            "Van Roon-Mom (LUMC) proposed this for AlleleSelect."
        ))
    parser.add_argument("--mismatch-type", default="auto",
        help=(
            "Base to introduce as deliberate extra mismatch: auto (most destabilizing) "
            "or A/C/G/T. Default: auto."
        ))
    parser.add_argument("--recommend-mods", action="store_true",
        help=(
            "Recommend chemical modifications at gap positions flanking the SNP "
            "(2S-dT, FRNA, S-cEt) to suppress minor RNase H cleavage sites. "
            "Based on Ostergaard 2013 Figures 3, 5, 7. "
            "Khvorova (UMass) identified this as the key remaining gap in v2."
        ))
    parser.add_argument("--verbose", action="store_true",
        help="Print detailed progress messages.")
    return parser


def run(args) -> None:
    from alleleselect.sequence.hgvs_parser import parse_hgvs_coding, validate_against_cds
    from alleleselect.sequence.fetcher import fetch_cds_sequence, generate_mutant_cds, extract_window
    from alleleselect.scoring.allele_selectivity import generate_candidate_windows
    from alleleselect.scoring.accessibility import (
        run_rnafold, compute_window_accessibility, check_rnafold_available,
        compute_differential_accessibility,
    )
    from alleleselect.scoring.offtarget import run_blast_offtarget
    from alleleselect.scoring.splice import (
        flag_splice_risk, get_splice_positions_for_r192q,
        get_splice_positions_for_transcript,
    )
    from alleleselect.modification.annotator import annotate_all_candidates
    from alleleselect.output.report import save_csv, save_html_report
    from alleleselect.scoring.snp_position import (
        score_snp_position, screen_toxic, composite_score,
        recommend_gap_modifications,
    )

    os.makedirs(args.output, exist_ok=True)
    log = print if args.verbose else lambda *a, **k: None

    gene_label = args.gene if args.gene else args.transcript.split(".")[0]

    # 1. Parse HGVS variant
    print(f"[AlleleSelect] Parsing variant: {args.variant}")
    parsed = parse_hgvs_coding(args.variant)
    log(f"  Position: {parsed['position']}, Ref: {parsed['ref']}, Alt: {parsed['alt']}")

    # 2. Fetch CDS
    print(f"[AlleleSelect] Fetching CDS from Ensembl ({args.transcript})...")
    wt_cds  = fetch_cds_sequence(args.transcript)
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

    # 4. mRNA accessibility (wildtype)
    if not args.no_rnafold:
        print("[AlleleSelect] Running RNAfold accessibility scoring...")
        import tempfile as _tempfile
        window_seq, window_start = extract_window(wt_cds, parsed["position"], flank=200)
        rnafold_work_dir = _tempfile.mkdtemp(prefix="alleleselect_rnafold_")
        rnafold_result   = run_rnafold(window_seq, work_dir=rnafold_work_dir)
        for c in candidates:
            acc = compute_window_accessibility(
                rnafold_result["per_base_unpaired"],
                c["mRNA_start"], c["mRNA_end"],
                cds_start_in_sequence=1 - window_start,
            )
            c["accessibility_score"] = acc
        log(f"  Mean accessibility: {sum(c['accessibility_score'] for c in candidates)/len(candidates):.3f}")
    else:
        print("[AlleleSelect] Skipping RNAfold (--no-rnafold). accessibility = 0.5.")
        for c in candidates:
            c["accessibility_score"] = 0.5

    # 4b. Differential accessibility (v4)
    if args.diff_accessibility and not args.no_rnafold:
        import tempfile as _tempfile
        diff_work = _tempfile.mkdtemp(prefix="alleleselect_diff_")
        candidates = compute_differential_accessibility(
            wt_cds, mut_cds, parsed["position"], candidates,
            flank=200, work_dir=diff_work,
        )
    else:
        for c in candidates:
            c["diff_accessibility"] = 0.0
            c["wt_accessibility"]   = c.get("accessibility_score", 0.5)
            c["mut_accessibility"]  = c.get("accessibility_score", 0.5)

    # 5. Off-target
    if not args.no_blast:
        print(f"[AlleleSelect] Running BLASTn off-target check (top {args.top_n_blast})...")
        candidates = run_blast_offtarget(
            candidates, top_n=args.top_n_blast, gene_name=gene_label
        )
    else:
        print("[AlleleSelect] Skipping BLASTn (--no-blast).")
        from alleleselect.scoring.offtarget import _set_unscreened
        for c in candidates:
            _set_unscreened(c)

    # 6. Splice site flagging
    print("[AlleleSelect] Flagging splice site proximity...")
    if args.no_splice_check:
        splice_positions = []
        log("  Splice check disabled (--no-splice-check).")
    elif args.transcript.startswith("ENST00000360228") and args.gene in (None, "CACNA1A"):
        splice_positions = get_splice_positions_for_r192q()
    else:
        splice_positions = get_splice_positions_for_transcript(
            transcript_id=args.transcript,
            mutation_pos=parsed["position"],
            flank=args.flank,
        )
    candidates = flag_splice_risk(candidates, splice_positions)

    # 7. Gapmer modification annotation
    print("[AlleleSelect] Annotating gapmer modification patterns...")
    candidates = annotate_all_candidates(candidates)

    # 7b. SNP position scoring + toxicity + composite (v2/v4)
    print("[AlleleSelect] Scoring SNP position and screening toxic sequences (v4)...")
    SNP_CDS_POS = parsed["position"]
    WING_LEN    = 5

    for c in candidates:
        aso_seq   = c.get("aso_seq", "")
        win_start = c.get("mRNA_start", 0)
        aso_len   = len(aso_seq)

        pos_result = score_snp_position(
            aso_len=aso_len, window_start=win_start,
            snp_cds_pos=SNP_CDS_POS, wing_len=WING_LEN,
        )
        c["snp_pos_in_aso"] = pos_result["snp_pos_in_aso"]
        c["snp_pos_score"]  = pos_result["snp_pos_score"]
        c["snp_region"]     = pos_result["snp_region"]

        tox = screen_toxic(aso_seq)
        c["tox_summary"] = tox["summary"]
        c["tox_serious"] = tox["serious"]
        c["tox_warning"] = tox["warning"]
        c["tox_flags"]   = "; ".join(
            f"{f['motif']}: {f['reason']}" for f in tox["flags"]
        ) if tox["flags"] else ""

        c["composite_score"] = composite_score(
            asr=c.get("allele_selectivity_ratio", 0.0),
            accessibility=c.get("accessibility_score", 0.5),
            snp_pos_score=c["snp_pos_score"],
            tox_serious=c["tox_serious"],
            diff_accessibility=c.get("diff_accessibility", 0.0),
        )

        # Chemical modification recommendation (v4, if requested)
        if args.recommend_mods:
            mod_rec = recommend_gap_modifications(
                aso_seq=aso_seq,
                snp_pos_in_aso=c["snp_pos_in_aso"],
                wing_len=WING_LEN,
            )
            c["mod_recommendation"] = mod_rec["synthesis_note"]
            c["mod_primary"]        = (
                mod_rec["recommendations"][0]["modification"]
                if mod_rec["recommendations"] else "N/A"
            )
            c["mod_primary_pos"]    = (
                mod_rec["recommendations"][0]["position_1indexed"]
                if mod_rec["recommendations"] else -1
            )
        else:
            c["mod_recommendation"] = ""
            c["mod_primary"]        = ""
            c["mod_primary_pos"]    = -1

    # 8. Rank by composite score
    candidates.sort(key=lambda c: -c.get("composite_score", 0.0))

    # 8b. Engineered mismatch candidates (v4, if requested)
    em_candidates = []
    if args.extra_mismatch:
        print("[AlleleSelect] Generating engineered-mismatch variants (top 5 candidates)...")
        from alleleselect.scoring.engineered_mismatch import (
            generate_mismatched_candidates, annotate_em_candidates,
        )
        em_raw = generate_mismatched_candidates(
            top_candidates=candidates[:5],
            wt_cds=wt_cds,
            snp_cds_pos=SNP_CDS_POS,
            extra_mismatch_type=args.mismatch_type,
        )
        em_candidates = annotate_em_candidates(em_raw, SNP_CDS_POS, wing_len=WING_LEN)
        # Re-score composite for em candidates
        for c in em_candidates:
            c["composite_score"] = composite_score(
                asr=c.get("allele_selectivity_ratio", 0.0),
                accessibility=c.get("accessibility_score", 0.5),
                snp_pos_score=c.get("snp_pos_score", 0.0),
                tox_serious=c.get("tox_serious", False),
                diff_accessibility=c.get("diff_accessibility", 0.0),
            )
        em_candidates.sort(key=lambda c: -c.get("composite_score", 0.0))
        print(f"  {len(em_candidates)} engineered-mismatch variants generated.")

    # 9. Save outputs
    variant_label = f"{gene_label} {args.variant}"
    csv_path  = os.path.join(args.output, "candidates.csv")
    html_path = os.path.join(args.output, "report.html")

    save_csv(candidates, csv_path)

    if em_candidates:
        em_csv = os.path.join(args.output, "candidates_engineered_mismatch.csv")
        save_csv(em_candidates, em_csv)
        print(f"Engineered mismatch CSV saved: {em_csv}")

    save_html_report(candidates, variant_label, html_path)

    # 10. Print top 5 summary
    print(f"\n[AlleleSelect] Complete. Outputs in: {args.output}/")
    print(f"  Top 5 candidates (ranked by composite score):")
    for i, c in enumerate(candidates[:5], 1):
        diff_str = f" | Diff={c.get('diff_accessibility',0):+.3f}" if args.diff_accessibility else ""
        mod_str  = f" | Mod={c.get('mod_primary','')}" if args.recommend_mods else ""
        ot_str   = (f"OT={c.get('off_target_count',-1)}"
                    f"(min_mm={c.get('min_off_target_mismatches',-1)},"
                    f"risk={c.get('ot_risk_level','?')})")
        print(
            f"  {i}. {c.get('ASO_ID','?')} | {c.get('aso_seq','')} | "
            f"Composite={c.get('composite_score',0):.4f} | "
            f"ASR={c.get('allele_selectivity_ratio',0):.3f} | "
            f"SNPpos={c.get('snp_pos_in_aso','?')}({c.get('snp_region','?')}) | "
            f"PosScore={c.get('snp_pos_score',0):.3f} | "
            f"Access={c.get('accessibility_score',0):.3f}"
            f"{diff_str} | "
            f"Tox={c.get('tox_summary','?')} | "
            f"{ot_str} | "
            f"Splice={c.get('splice_risk','?')}"
            f"{mod_str}"
        )

    if em_candidates:
        print(f"\n  Top 3 engineered-mismatch variants:")
        for i, c in enumerate(em_candidates[:3], 1):
            print(
                f"  EM{i}. {c.get('ASO_ID','?')} | {c.get('aso_seq','')} | "
                f"Composite={c.get('composite_score',0):.4f} | "
                f"ASR={c.get('allele_selectivity_ratio',0):.3f} | "
                f"Extra_mm@pos{c.get('extra_mismatch_pos','?')}={c.get('extra_mismatch_base','?')} | "
                f"Tox={c.get('tox_summary','?')}"
            )


def main():
    parser = build_parser()
    args   = parser.parse_args()
    try:
        run(args)
    except Exception as e:
        print(f"[AlleleSelect] Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()