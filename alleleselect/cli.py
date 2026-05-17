"""
cli.py
AlleleSelect command-line interface — v5

New in v5:
  --junction-mode       Junction/exon-skipping ASO design mode (Aguti & Zhou 2024)
  --exon-skip INT       Exon number skipped in mutant (auto-fetches junction sequences)
  --mut-fasta FILE      Manual mutant junction sequence (alternative to --exon-skip)
  --wt-fasta FILE       Manual wildtype region sequence
  --junction-center INT Position of junction center in mut-fasta (0-indexed)
  Wing position fix: SNP in wing now scores 0.10 (wing_caution) not 0.0
                     Based on Aguti (UCL) + Elgersma (Erasmus MC) feedback

New in v4:
  --extra-mismatch      Engineered mismatch mode (van Roon-Mom / Ostergaard 2013 Fig 7)
  --diff-accessibility  Differential WT vs. mutant RNAfold scoring (Aguti & Zhou 2024)
  --recommend-mods      Chemical modification recommendations at SNP-flanking positions
  --gene / --no-splice-check  (v3, unchanged)

Usage:
    # Standard SNP mode (CACNA1A R192Q):
    alleleselect --variant c.575G>A --transcript ENST00000360228.10 --output demo/R192Q_output/

    # Junction mode (COL6A3 exon 16 skipping):
    alleleselect --junction-mode --gene COL6A3 --wt-transcript ENST00000295550.9 --exon-skip 16 --output demo/COL6A3_junction/

    # Junction mode with manual FASTA files:
    alleleselect --junction-mode --gene COL6A3 --mut-fasta mut_junction.fa --wt-fasta wt_region.fa --junction-center 50 --output demo/COL6A3_junction/
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
    parser.add_argument("--variant", "-v", default=None,
        help="HGVS coding notation (e.g. c.575G>A). Required unless --junction-mode is used.")
    parser.add_argument("--transcript", "-t", default="ENST00000360228",
        help="Ensembl transcript ID. Default: ENST00000360228 (CACNA1A canonical).")
    parser.add_argument("--gene", "-g", default=None,
        help="Gene name for output labeling (e.g. CACNA1A, COL6A3, ATXN1).")
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
    # ── v5: junction mode ─────────────────────────────────────────────────────
    parser.add_argument("--junction-mode", action="store_true",
        help=(
            "Junction/exon-skipping design mode. For mutations that cause a novel "
            "exon-exon junction in the mutant mRNA (e.g. splice site mutations). "
            "The ASO targets the novel junction which is absent from wildtype mRNA. "
            "Use with --exon-skip or --mut-fasta/--wt-fasta."
        ))
    parser.add_argument("--wt-transcript", default=None,
        help="Ensembl transcript ID for junction mode (auto-fetches exon sequences).")
    parser.add_argument("--exon-skip", type=int, default=None,
        help="Exon number (1-indexed) that is skipped in the mutant mRNA.")
    parser.add_argument("--mut-fasta", default=None,
        help="FASTA file with mutant junction sequence (alternative to --exon-skip).")
    parser.add_argument("--wt-fasta", default=None,
        help="FASTA file with wildtype region sequence (for cross-reactivity check).")
    parser.add_argument("--junction-center", type=int, default=None,
        help="0-indexed position of junction in --mut-fasta sequence.")
    parser.add_argument("--junction-flank", type=int, default=100,
        help="Nucleotides from each exon to include in the junction. Default: 100.")
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


def run_junction(args) -> None:
    """
    Junction mode: design ASOs for exon-skipping / novel junction mutations.
    Motivated by Aguti & Zhou 2024 (COL6A3 c.6210+1G>A exon 16 skipping).
    """
    from alleleselect.sequence.junction_fetcher import (
        build_junction_sequences, load_fasta_sequence,
        generate_junction_candidates, score_junction_composite,
    )
    from alleleselect.scoring.accessibility import run_rnafold, compute_window_accessibility
    from alleleselect.scoring.offtarget import run_blast_offtarget
    from alleleselect.modification.annotator import annotate_all_candidates
    from alleleselect.scoring.snp_position import screen_toxic
    from alleleselect.output.report import save_csv, save_html_report

    import tempfile as _tempfile

    os.makedirs(args.output, exist_ok=True)
    log = print if args.verbose else lambda *a, **k: None

    gene_label = args.gene if args.gene else "GENE"
    print(f"[AlleleSelect] Junction mode: designing ASOs for {gene_label} novel junction")

    # Step 1: Get sequences
    if args.exon_skip:
        transcript = args.wt_transcript or args.transcript
        print(f"[AlleleSelect] Fetching exon sequences from Ensembl ({transcript}, exon {args.exon_skip} skipped)...")
        mut_seq, wt_seq, junction_center = build_junction_sequences(
            transcript_id=transcript,
            skipped_exon=args.exon_skip,
            junction_flank=args.junction_flank,
        )
        log(f"  Mutant junction: {len(mut_seq)} nt, junction at position {junction_center}")
        log(f"  Wildtype region: {len(wt_seq)} nt")
    else:
        print(f"[AlleleSelect] Loading sequences from FASTA files...")
        mut_seq = load_fasta_sequence(args.mut_fasta)
        wt_seq  = load_fasta_sequence(args.wt_fasta) if args.wt_fasta else ""
        junction_center = args.junction_center

    # Step 2: Generate candidates
    print(f"[AlleleSelect] Generating junction ASO candidates (lengths: {args.aso_lengths}, flank: {args.flank})...")
    candidates = generate_junction_candidates(
        mut_junction_seq=mut_seq,
        wt_region_seq=wt_seq,
        junction_center=junction_center,
        aso_lengths=args.aso_lengths,
        flank=args.flank,
        gene_label=gene_label,
    )
    print(f"  {len(candidates)} candidates generated.")
    n_specific = sum(1 for c in candidates if c["junction_specificity"] == 1.0)
    print(f"  {n_specific} candidates with no wildtype cross-reactivity (junction-specific).")

    # Step 3: RNAfold accessibility on mutant junction
    if not args.no_rnafold:
        print("[AlleleSelect] Running RNAfold accessibility scoring (mutant junction)...")
        rnafold_dir = _tempfile.mkdtemp(prefix="alleleselect_junc_rnafold_")
        rnafold_result = run_rnafold(mut_seq, work_dir=rnafold_dir)
        for c in candidates:
            acc = compute_window_accessibility(
                rnafold_result["per_base_unpaired"],
                c["mRNA_start"], c["mRNA_end"],
                cds_start_in_sequence=0,
            )
            c["accessibility_score"] = acc
        log(f"  Mean accessibility: {sum(c['accessibility_score'] for c in candidates)/len(candidates):.3f}")

    # Step 4: Toxicity screening
    print("[AlleleSelect] Screening toxic sequences...")
    for c in candidates:
        tox = screen_toxic(c["aso_seq"])
        c["tox_summary"] = tox["summary"]
        c["tox_serious"] = tox["serious"]
        c["tox_warning"] = tox["warning"]
        c["tox_flags"]   = "; ".join(
            f"{f['motif']}: {f['reason']}" for f in tox["flags"]
        ) if tox["flags"] else ""

    # Step 5: Gapmer annotation
    print("[AlleleSelect] Annotating gapmer modification patterns...")
    candidates = annotate_all_candidates(candidates)

    # Step 6: Score and rank
    candidates = score_junction_composite(candidates)
    candidates.sort(key=lambda c: -c.get("composite_score", 0.0))

    # Step 7: Off-target BLAST (optional)
    if not args.no_blast:
        print(f"[AlleleSelect] Running BLASTn off-target check (top {args.top_n_blast})...")
        candidates = run_blast_offtarget(candidates, top_n=args.top_n_blast, gene_name=gene_label)

    # Step 8: Save outputs
    variant_label = f"{gene_label} junction (exon {args.exon_skip} skip)" if args.exon_skip else f"{gene_label} junction"
    csv_path  = os.path.join(args.output, "candidates_junction.csv")
    html_path = os.path.join(args.output, "report_junction.html")
    save_csv(candidates, csv_path)
    save_html_report(candidates, variant_label, html_path)

    print(f"\n[AlleleSelect] Junction mode complete. Outputs in: {args.output}/")
    print(f"  Top 5 junction candidates:")
    for i, c in enumerate(candidates[:5], 1):
        print(
            f"  {i}. {c.get('ASO_ID','?')} | {c.get('aso_seq','')} | "
            f"Composite={c.get('composite_score',0):.4f} | "
            f"JuncSpec={c.get('junction_specificity',0):.1f} | "
            f"JuncPos={c.get('junction_pos_score',0):.3f} | "
            f"Access={c.get('accessibility_score',0):.3f} | "
            f"Tox={c.get('tox_summary','?')} | "
            f"OT={c.get('off_target_count',-1)} | "
            f"WTcross={c.get('wt_cross_reactivity',0)}"
        )


def main():
    parser = build_parser()
    args   = parser.parse_args()

    # Validate arguments
    if args.junction_mode:
        # Junction mode: needs either exon-skip+wt-transcript or mut-fasta+wt-fasta
        if args.exon_skip and not args.wt_transcript:
            parser.error("--exon-skip requires --wt-transcript")
        if args.mut_fasta and args.junction_center is None:
            parser.error("--mut-fasta requires --junction-center")
        if not args.exon_skip and not args.mut_fasta:
            parser.error("--junction-mode requires either --exon-skip or --mut-fasta")
    else:
        if not args.variant:
            parser.error("--variant is required unless --junction-mode is used")

    try:
        if args.junction_mode:
            run_junction(args)
        else:
            run(args)
    except Exception as e:
        print(f"[AlleleSelect] Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()