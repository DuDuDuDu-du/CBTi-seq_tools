import os
import subprocess
from cbtiseq_tools.utils import run_cmd

def parse_length(pos):
    if pos is None:
        return 0
    start, end = map(int, pos.split("-"))
    return end - start + 1

def process_bam(
    workdir,
    sample,
    formatted_r1,
    formatted_r2,
    picard_path,
    dropseq_path,
    bcA=None,
    umiA=None,
    bcB=None,
    umiB=None
):
    os.makedirs(f"{workdir}/tmp", exist_ok=True)

    # 自动计算长度
    umi_len = parse_length(umiA) + parse_length(umiB)
    bc_len = parse_length(bcA) + parse_length(bcB)

    umi_range = f"1-{umi_len}"
    bc_range = f"{umi_len + 1}-{umi_len + bc_len}"

    raw_bam = f"{workdir}/unaligned_raw_{sample}.bam"
    cell_bam = f"{workdir}/unaligned_tagged_{sample}_Cell.bam"
    umi_bam = f"{workdir}/unaligned_tagged_{sample}_CellMolecular.bam"
    filtered_bam = f"{workdir}/unaligned_tagged_{sample}_filtered.bam"

    # Step 1 FastqToSam
    cmd1 = f"""
    java -Xmx4000m -jar {picard_path} FastqToSam \
      TMP_DIR={workdir}/tmp \
      F1={formatted_r1} \
      F2={formatted_r2} \
      O={raw_bam} \
      SM={sample} \
      RG={sample}
    """
    run_cmd(cmd1, "FastqToSam")

    # Step 2 Cell barcode tagging
    cmd2 = f"""
    {dropseq_path}/TagBamWithReadSequenceExtended \
      SUMMARY={workdir}/cell_barcode_summary.txt \
      BASE_RANGE={bc_range} \
      BASE_QUALITY=10 \
      BARCODED_READ=1 \
      DISCARD_READ=false \
      TAG_NAME=XC \
      NUM_BASES_BELOW_QUALITY=1 \
      INPUT={raw_bam} \
      OUTPUT={cell_bam} \
      TMP_DIR={workdir}/tmp
    """
    run_cmd(cmd2, "Cell barcode tagging")

    # Step 3 UMI tagging
    cmd3 = f"""
    {dropseq_path}/TagBamWithReadSequenceExtended \
      SUMMARY={workdir}/umi_summary.txt \
      BASE_RANGE={umi_range} \
      BASE_QUALITY=10 \
      BARCODED_READ=1 \
      DISCARD_READ=false \
      TAG_NAME=XM \
      NUM_BASES_BELOW_QUALITY=1 \
      INPUT={cell_bam} \
      OUTPUT={umi_bam} \
      TMP_DIR={workdir}/tmp
    """
    run_cmd(cmd3, "UMI tagging")

    # Step 4 Filter low quality
    cmd4 = f"""
    {dropseq_path}/FilterBam \
      TAG_REJECT=XQ \
      INPUT={umi_bam} \
      OUTPUT={filtered_bam} \
      TMP_DIR={workdir}/tmp
    """
    run_cmd(cmd4, "Filter BAM")

    return filtered_bam
