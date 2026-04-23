import os
from cbtiseq_tools.config import load_paths
from cbtiseq_tools.modules.reformat_fastq import reformat_fastq
from cbtiseq_tools.modules.process_bam import process_bam
from cbtiseq_tools.modules.tag_bam import tag_bam
from cbtiseq_tools.modules.alignment_bam import alignment
from cbtiseq_tools.modules.summary_bam import generate_summary


def run_pipeline(args):

    paths = load_paths()

    os.makedirs(args.workdir, exist_ok=True)

    formatted_r1 = f"{args.workdir}/{args.sample}_formatted_R1.fastq"
    formatted_r2 = f"{args.workdir}/{args.sample}_formatted_R2.fastq"

    print("[STEP 1] Reformat FASTQ")

    reformat_fastq(
        args.r1,
        args.r2,
        formatted_r1,
        formatted_r2,
        args.bcA,
        args.umiA,
        args.bcB,
        args.umiB
    )

    print("[STEP 2] Process BAM")

    filtered_bam = process_bam(
        args.workdir,
        args.sample,
        formatted_r1,
        formatted_r2,
        paths["picard"],
        paths["dropseq"],
        args.bcA,
        args.umiA,
        args.bcB,
        args.umiB
    )
    
    print("[STEP 3] Tag BAM and PolyA trim")

    merged_bam = tag_bam(
        args.workdir,
        args.sample,
        filtered_bam,
        paths["dropseq"]
    )

    print("[STEP 4] Alignment and DGE")

    final_bam, dge_file, dge_summary = alignment(
        args.workdir,
        args.sample,
        merged_bam,
        paths,
        args.whitelist
    )

    print(f"[FINAL BAM] {final_bam}")
    print(f"[DGE MATRIX] {dge_file}")
    print(f"[SUMMARY] {dge_summary}")
    
    print("[STEP 5] Generate enhanced summary")

    enhanced_summary = f"{args.workdir}/{args.sample}_enhance_summary.txt"

    generate_summary(
        raw_bam=filtered_bam,
        trimmed_bam=merged_bam,
        dge_summary=dge_summary,
        output_file=enhanced_summary
    )

    print(f"[SUMMARY REPORT] {enhanced_summary}")

print("[PIPELINE COMPLETED SUCCESSFULLY]")
