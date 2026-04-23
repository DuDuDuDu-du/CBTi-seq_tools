import os
import subprocess
import pysam
from cbtiseq_tools.utils import run_cmd

def split_bam_by_pair(input_bam, out_r1, out_r2):
    print("[RUNNING] Split BAM into read1/read2")

    with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as bam:
        with pysam.AlignmentFile(out_r1, "wb", header=bam.header) as f1, \
             pysam.AlignmentFile(out_r2, "wb", header=bam.header) as f2:

            count1 = 0
            count2 = 0

            for read in bam:
                if read.is_read1:
                    f1.write(read)
                    count1 += 1
                elif read.is_read2:
                    f2.write(read)
                    count2 += 1

    print(f"[DONE] Split complete: R1={count1}, R2={count2}\n")


def filter_paired_reads(read1_bam, read2_bam, out_r1, out_r2):
    print("[RUNNING] Filter paired reads")

    read2_dict = {}

    with pysam.AlignmentFile(read2_bam, "rb", check_sq=False) as r2:
        for read in r2:
            read2_dict[read.query_name] = read

    paired_count = 0

    with pysam.AlignmentFile(read1_bam, "rb", check_sq=False) as r1, \
         pysam.AlignmentFile(out_r1, "wb", header=r1.header) as out1, \
         pysam.AlignmentFile(out_r2, "wb", header=r1.header) as out2:

        for read1 in r1:
            qname = read1.query_name

            if qname in read2_dict:
                read2 = read2_dict[qname]

                # 复制 XC / XM
                if read2.has_tag("XC"):
                    read1.set_tag("XC", read2.get_tag("XC"))

                if read2.has_tag("XM"):
                    read1.set_tag("XM", read2.get_tag("XM"))

                out1.write(read1)
                out2.write(read2)

                paired_count += 1

    print(f"[DONE] Paired reads retained: {paired_count}\n")


def merge_paired_bam(read1_bam, read2_bam, output_bam):
    print("[RUNNING] Merge paired BAM")

    read2_dict = {}

    with pysam.AlignmentFile(read2_bam, "rb", check_sq=False) as r2:
        for read in r2:
            read2_dict[read.query_name] = read

    merged_count = 0

    with pysam.AlignmentFile(read1_bam, "rb", check_sq=False) as r1, \
         pysam.AlignmentFile(output_bam, "wb", header=r1.header) as out:

        for read1 in r1:
            qname = read1.query_name

            if qname in read2_dict:
                out.write(read1)
                out.write(read2_dict[qname])
                merged_count += 1

    print(f"[DONE] Merge complete: {merged_count} pairs\n")


def tag_bam(workdir, sample, filtered_bam, dropseq_path):
    tmp_dir = f"{workdir}/tmp"

    read1_bam = f"{workdir}/{sample}_read1.bam"
    read2_bam = f"{workdir}/{sample}_read2.bam"

    trimmed_r2 = f"{workdir}/{sample}_trimmed_read2.bam"

    paired_r1 = f"{workdir}/{sample}_paired_read1.bam"
    paired_r2 = f"{workdir}/{sample}_paired_read2.bam"

    merged_bam = f"{workdir}/unaligned_mc_tagged_polyA_{sample}_filtered.bam"

    # Step 1 split
    split_bam_by_pair(filtered_bam, read1_bam, read2_bam)

    # Step 2 PolyA trimming
    cmd = f"""
    {dropseq_path}/PolyATrimmer \
      OUTPUT={trimmed_r2} \
      OUTPUT_SUMMARY={workdir}/polyA_report.txt \
      TMP_DIR={tmp_dir} \
      MISMATCHES=0 \
      NUM_BASES=6 \
      NEW=true \
      INPUT={read2_bam}
    """
    run_cmd(cmd, "PolyA trimming")

    # Step 3 pair filter + tag sync
    filter_paired_reads(
        read1_bam,
        trimmed_r2,
        paired_r1,
        paired_r2
    )

    # Step 4 merge
    merge_paired_bam(
        paired_r1,
        paired_r2,
        merged_bam
    )

    return merged_bam
