import os
import subprocess
from cbtiseq_tools.utils import run_cmd

def alignment(workdir, sample, input_bam, paths, whitelist=None):
    tmp_dir = f"{workdir}/tmp"
    os.makedirs(tmp_dir, exist_ok=True)

    read1_fastq = f"{workdir}/{sample}_read1.fastq"
    read2_fastq = f"{workdir}/{sample}_read2.fastq"

    aligned_bam = f"{workdir}/aligned.sorted.{sample}.bam"
    merged_bam = f"{workdir}/{sample}_merged.bam"
    final_bam = f"{workdir}/{sample}_final.bam"

    dge_file = f"{workdir}/{sample}_UMI1.txt.gz"
    dge_summary = f"{workdir}/{sample}_UMI1.summary.txt"

    # Step 1 BAM -> paired FASTQ
    cmd1 = f"""
    java -Xmx500m -jar {paths["picard"]} SamToFastq \
      INPUT={input_bam} \
      TMP_DIR={tmp_dir} \
      FASTQ={read1_fastq} \
      SECOND_END_FASTQ={read2_fastq}
    """
    run_cmd(cmd1, "SamToFastq")

    # Step 2 STAR alignment
    cmd2 = f"""
    {paths["star"]} \
      --readFilesIn {read1_fastq} {read2_fastq} \
      --genomeDir {paths["genome_dir"]} \
      --outFileNamePrefix {workdir}/star. \
      --runThreadN 8
    """
    run_cmd(cmd2, "STAR alignment")

    # Step 3 SortSam
    cmd3 = f"""
    java -Xmx4000m -jar {paths["picard"]} SortSam \
      INPUT={workdir}/star.Aligned.out.sam \
      OUTPUT={aligned_bam} \
      SORT_ORDER=queryname \
      TMP_DIR={tmp_dir}
    """
    run_cmd(cmd3, "SortSam")

    # Step 4 Merge alignment
    cmd4 = f"""
    java -Xmx4000m -jar {paths["picard"]} MergeBamAlignment \
      REFERENCE_SEQUENCE={paths["reference_fasta"]} \
      UNMAPPED_BAM={input_bam} \
      ALIGNED_BAM={aligned_bam} \
      INCLUDE_SECONDARY_ALIGNMENTS=false \
      PAIRED_RUN=true \
      CLIP_ADAPTERS=false \
      TMP_DIR={tmp_dir} \
      OUTPUT={merged_bam}
    """
    run_cmd(cmd4, "MergeBamAlignment")

    # Step 5 Gene annotation
    cmd5 = f"""
    {paths["dropseq"]}/TagReadWithGeneFunction \
      O={final_bam} \
      ANNOTATIONS_FILE={paths["refflat"]} \
      TMP_DIR={tmp_dir} \
      INPUT={merged_bam}
    """
    run_cmd(cmd5, "Gene annotation")

    # Step 6 DGE
    dge_cmd = f"""
    {paths["dropseq"]}/DigitalExpression \
      INPUT={final_bam} \
      O={dge_file} \
      TMP_DIR={tmp_dir} \
      SUMMARY={dge_summary} \
      MIN_NUM_TRANSCRIPTS_PER_CELL=1
    """

    if whitelist:
        dge_cmd += f" CELL_BC_FILE={whitelist}"

    run_cmd(dge_cmd, "Digital expression matrix")

    return final_bam, dge_file, dge_summary
